"""Step 02: lipocalin genes in the RefSeq gene annotations (run inside WSL: run_wsl.sh step02_annotated.py).

Why: NCBI's RefSeq pipeline has already predicted genes in 17 of the genomes. Its proteins are the quickest place to
look for MUPs and their relatives. Pseudogenes have no protein, so they are listed from the annotation by name.

What it does, per species in config 'refseq':
1. Reads the GFF3 annotation: every gene (protein-coding or not, pseudogenes included) with its location, and for
   protein-coding genes the longest protein and NCBI's 'model_evidence' note (how much RNA-seq supports the model).
   NCBI's assembly report gives chromosome names and sequence lengths.
2. hmmsearch (HMMER 3) with the Pfam lipocalin profile PF00061 over the longest protein of every gene. A gene is a
   lipocalin candidate if its protein has a PF00061 domain with independent E-value < 1e-3. The domain's envelope
   (the stretch of protein that matches the profile) is kept for the tree, padded by 10 residues each side.
3. Every gene whose name or description says major urinary, odorant-binding, aphrodisin or lipocalin is listed too,
   so pseudogenes and genes without a protein are not missed.

Outputs (derived/annotated/)
- genes_<slug>.tsv.gz      every gene: gene_id, symbol, description, biotype, seqid, chrom, start, end, strand
- lipocalins.tsv           lipocalin candidates of all species (hmm hits + name matches), one row per gene
- lipocalin_domains.faa    PF00061 envelope (+-10 aa) of each candidate protein; full proteins in lipocalins_full.faa
- hmmsearch.domtbl         raw HMMER domain table
"""
from __future__ import annotations

import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

from vm_common import (fasta_records, gff_attrs, load_config, lower_priority, paths, read_tsv, refseq_files, run,
                       unquote, write_fasta, write_tsv)

NAME_RE = re.compile(r"major urinary|odorant[- ]binding|aphrodisin|lipocalin|alpha-2u|\bMup\d|\bObp\d", re.I)
EVID_RE = re.compile(r"(\d+)% coverage of the annotated genomic feature by RNAseq alignments(?:, including (\d+) "
                     r"samples? with support for all annotated introns)?")


def assembly_report(path: Path) -> dict:
    """RefSeq accession -> (chromosome/scaffold name, length)."""
    out = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            name, role, mol, refseq_acc, length, ucsc = f[0], f[1], f[2], f[6], f[8], f[9]
            if refseq_acc and refseq_acc != "na":
                chrom = mol if role == "assembled-molecule" and mol != "na" else name
                out[refseq_acc] = (chrom, int(length) if length.isdigit() else 0)
    return out


def parse_gff(path: Path):
    genes, mrna_gene, mrna_evid, prot_gene, product = {}, {}, {}, {}, {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqid, _src, typ, start, end, _s, strand, _p, attr = f
            if typ in ("gene", "pseudogene"):
                a = gff_attrs(attr)
                gid = ""
                for x in a.get("Dbxref", "").split(","):
                    if x.startswith("GeneID:"):
                        gid = x[7:]
                genes[a["ID"]] = {"gene_id": gid, "symbol": a.get("Name", ""), "description": unquote(a.get("description", "")),
                                  "biotype": a.get("gene_biotype", typ), "pseudo": a.get("pseudo", "") == "true",
                                  "seqid": seqid, "start": int(start), "end": int(end), "strand": strand}
            elif typ in ("mRNA", "transcript"):
                a = gff_attrs(attr)
                mrna_gene[a["ID"]] = a.get("Parent", "")
                if "product" in a and a.get("Parent", "") not in product:
                    product[a["Parent"]] = re.sub(r",? transcript variant \S+$", "", unquote(a["product"]))
                if "model_evidence" in a:
                    mrna_evid[a["ID"]] = unquote(a["model_evidence"])
            elif typ == "CDS":
                a = gff_attrs(attr)
                pid = a.get("protein_id")
                if pid and pid not in prot_gene:
                    par = a.get("Parent", "")
                    prot_gene[pid] = (mrna_gene.get(par, par), par)
    for gkey, g in genes.items():          # older RefSeq GFFs (e.g. MicOch1.0) have no gene 'description'
        if not g["description"] and gkey in product:
            g["description"] = product[gkey]
    return genes, prot_gene, mrna_evid


def evidence_summary(note: str) -> tuple[int | None, int | None]:
    if not note:
        return None, None
    m = EVID_RE.search(note)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2)) if m.group(2) else 0


def main() -> None:
    lower_priority()
    cfg = load_config()
    P = paths(cfg)
    out = P["derived"] / "annotated"
    out.mkdir(parents=True, exist_ok=True)
    threads = str(cfg["threads"])

    hmm = out / "PF00061.hmm"
    if not hmm.exists():
        with gzip.open(P["raw"] / f"pfam_{cfg['pfam_hmm']}" / f"{cfg['pfam_hmm']}.hmm.gz", "rt") as fh:
            hmm.write_text(fh.read())

    all_prot, info = [], {}            # longest protein per gene, all species
    name_hits = []
    for sp in cfg["refseq"]:
        f = refseq_files(cfg, sp["accession"])
        rep = assembly_report(f["report"])
        genes, prot_gene, evid = parse_gff(f["gff"])
        rows = []
        for gkey, g in genes.items():
            chrom, slen = rep.get(g["seqid"], (g["seqid"], 0))
            g["chrom"], g["seqlen"] = chrom, slen
            rows.append(g)
        write_tsv(out / f"genes_{sp['slug']}.tsv", rows,
                  ["gene_id", "symbol", "description", "biotype", "pseudo", "seqid", "chrom", "seqlen", "start", "end",
                   "strand"])
        # compress
        raw = (out / f"genes_{sp['slug']}.tsv").read_bytes()
        with gzip.open(out / f"genes_{sp['slug']}.tsv.gz", "wb") as fh:
            fh.write(raw)
        (out / f"genes_{sp['slug']}.tsv").unlink()

        best = {}
        for pid, desc, seq in fasta_records(f["protein"]):
            gkey, mrna = prot_gene.get(pid, (None, None))
            if gkey is None or gkey not in genes:
                continue
            if gkey not in best or len(seq) > len(best[gkey][1]):
                best[gkey] = (pid, seq.rstrip("*"), mrna, desc)
        for gkey, (pid, seq, mrna, desc) in best.items():
            g = genes[gkey]
            key = f"{sp['slug']}|{g['gene_id']}"
            cov, nsamp = evidence_summary(evid.get(mrna, ""))
            info[key] = dict(g, slug=sp["slug"], species=sp["name"], common=sp["common"], group=sp["group"],
                             protein_id=pid, length=len(seq), rnaseq_cov=cov, rnaseq_intron_samples=nsamp,
                             model_evidence=evid.get(mrna, ""), seq=seq)
            all_prot.append((key, seq))
        for gkey, g in genes.items():
            if NAME_RE.search(g["description"]) or NAME_RE.search(g["symbol"]):
                key = f"{sp['slug']}|{g['gene_id']}"
                name_hits.append(key)
                if key not in info:
                    info[key] = dict(g, slug=sp["slug"], species=sp["name"], common=sp["common"], group=sp["group"],
                                     protein_id="", length=0, rnaseq_cov=None, rnaseq_intron_samples=None,
                                     model_evidence="", seq="")
        print(f"{sp['slug']}: {len(genes)} genes, {len(best)} proteins", flush=True)

    allfa = out / "all_longest_proteins.faa"
    write_fasta(allfa, all_prot)
    domtbl = out / "hmmsearch.domtbl"
    if not domtbl.exists() or "--redo-hmm" in sys.argv:      # the proteins only change when RefSeq files do
        run(["hmmsearch", "--cpu", threads, "-E", "1e-3", "--domE", "1e-3", "--domtblout", domtbl, "-o", "/dev/null",
             hmm, allfa])
    dom = defaultdict(list)
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            x = line.split()
            dom[x[0]].append({"i_evalue": float(x[12]), "score": float(x[13]), "env_from": int(x[19]),
                              "env_to": int(x[20])})
    keep = set(dom) | set(name_hits)
    rows, dom_fa, full_fa = [], [], []
    for key in sorted(keep):
        r = info[key]
        d = sorted(dom.get(key, []), key=lambda z: z["i_evalue"])
        r = dict(r, hmm_evalue=d[0]["i_evalue"] if d else "", hmm_score=d[0]["score"] if d else "",
                 n_domains=len(d), env_from=d[0]["env_from"] if d else "", env_to=d[0]["env_to"] if d else "",
                 name_match=key in name_hits, key=key)
        rows.append(r)
        if d and r["seq"]:
            a = max(1, min(z["env_from"] for z in d[:1]) - 10)
            b = min(len(r["seq"]), max(z["env_to"] for z in d[:1]) + 10)
            dom_fa.append((key, r["seq"][a - 1:b]))
            full_fa.append((key, r["seq"]))
    fields = ["key", "slug", "species", "common", "group", "gene_id", "symbol", "description", "biotype", "pseudo",
              "seqid", "chrom", "seqlen", "start", "end", "strand", "protein_id", "length", "hmm_evalue", "hmm_score",
              "n_domains", "env_from", "env_to", "name_match", "rnaseq_cov", "rnaseq_intron_samples", "model_evidence"]
    write_tsv(out / "lipocalins.tsv", rows, fields)
    write_fasta(out / "lipocalin_domains.faa", dom_fa)
    write_fasta(out / "lipocalins_full.faa", full_fa)
    print(f"lipocalin candidates: {len(rows)} genes ({len(dom_fa)} with a PF00061 domain)", flush=True)


if __name__ == "__main__":
    main()
