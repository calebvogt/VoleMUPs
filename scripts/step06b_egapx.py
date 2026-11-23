"""Step 06b: sanity check against the Dennis lab's EGAPx gene annotation of the 2026 vole assemblies
(Abuelanin et al. 2026; github.com/mydennislab/2026-voles-assembly, annotation/). Used for analysis only, not
redistributed (that repository has no license).

1. Lipocalin-family genes (description: urinary / odorant / aphrodisin / lipocalin) in the prairie and meadow vole
   EGAPx annotations, by chromosome.
2. Prairie vole: each EGAPx gene matched by position to the DNA-search loci (step 04) on the same assembly (its chrN
   sequences are the GenBank hap1 chromosomes, same lengths and orientation; checked here).
3. Prairie vole liver RNA-seq (ERR17322872) counted against the EGAPx MUP and OBP transcripts (exons spliced from the
   hap1 assembly), with the same k-mer method as step 06, so the two gene sets can be compared.

Outputs: <data_root>/derived/egapx/ (lipocalin_genes.tsv, prairie_vs_dnasearch.tsv, prairie_liver_counts.tsv)
usage: python scripts/step06b_egapx.py
"""
from __future__ import annotations

import gzip
import re

import requests

import step06_expression as s6
from vm_common import fasta_records, genome_fasta, load_config, now, paths, read_tsv, revcomp, write_tsv

REPO = "https://raw.githubusercontent.com/mydennislab/2026-voles-assembly/main/annotation"
FILES = {"prairie": "prairie_hap1_egapx.gff.gz", "meadow": "meadow_hap1_egapx.gff.gz"}
PRAIRIE_ASM = "GCA_984569965.1"
PRAIRIE_GENOME = "microtus_ochrogaster_ucd1"
LIVER_RUN = "ERR17322872"
LIPO = re.compile(r"urinary|odorant|aphrodisin|lipocalin", re.I)


def attrs(s: str) -> dict:
    return dict(kv.split("=", 1) for kv in s.split(";") if "=" in kv)


def fetch(P, sess) -> dict:
    out = P["raw"] / "egapx_dennis2026"
    out.mkdir(parents=True, exist_ok=True)
    got = {}
    for sp, name in FILES.items():
        dest = out / name
        if not dest.exists():
            r = sess.get(f"{REPO}/{name}", timeout=300)
            r.raise_for_status()
            dest.write_bytes(r.content)
        got[sp] = dest
    readme = out / "README.txt"
    if not readme.exists():
        readme.write_text(f"EGAPx gene annotations of the 2026 vole hap1 assemblies (Abuelanin et al. 2026, Cell Genomics), "
                          f"from {REPO}/, downloaded {now()} by scripts/step06b_egapx.py. The repository has no "
                          "license: for analysis only, not redistributed.\n", encoding="utf-8")
    return got


def parse(path) -> tuple[dict, dict, list]:
    """genes (lipocalin-family only), mRNA -> gene, and seq-region lengths."""
    genes, mrna, exons, lens = {}, {}, {}, {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                _t, seq, _a, b = line.split()
                lens[seq] = int(b)
                continue
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9:
                continue
            a = attrs(x[8])
            if x[2] in ("gene", "pseudogene") and LIPO.search(a.get("description", "")):
                genes[a["ID"]] = {"gene": a.get("Name", a["ID"]), "type": x[2], "description": a.get("description", ""),
                                  "chrom": x[0], "start": int(x[3]), "end": int(x[4]), "strand": x[6]}
            elif x[2] in ("mRNA", "transcript") and a.get("Parent") in genes:
                mrna[a["ID"]] = a["Parent"]
            elif x[2] == "exon" and a.get("Parent") in mrna:
                exons.setdefault(a["Parent"], []).append((int(x[3]), int(x[4])))
    return genes, {m: (g, sorted(exons.get(m, []))) for m, g in mrna.items()}, lens


def chrom_map(sess, lens: dict) -> dict:
    """EGAPx chrN -> GenBank accession of the same hap1 chromosome; lengths must agree."""
    rep = sess.get(f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{PRAIRIE_ASM}/sequence_reports",
                   params={"page_size": 1000}, timeout=120).json()["reports"]
    out = {}
    for r in rep:
        name = "chr" + str(r.get("chr_name", ""))
        if name in lens and r.get("role") == "assembled-molecule":
            if int(r["length"]) != lens[name]:
                raise RuntimeError(f"{name}: EGAPx length {lens[name]} != GenBank {r['length']}")
            out[name] = r["genbank_accession"]
    return out


_CACHE = {}


def egapx_models(cfg, P) -> tuple[dict, dict]:
    """Prairie vole EGAPx MUP and OBP transcripts (exons spliced from the hap1 assembly): gene key -> (label, seq), and
    gene key -> MUP / OBP. Families come from the overlapping DNA-search locus (or the description where none)."""
    if "prairie" in _CACHE:
        return _CACHE["prairie"]
    sess = requests.Session()
    sess.headers["User-Agent"] = cfg["user_agent"]
    files = fetch(P, sess)
    genes, mrna, lens = parse(files["prairie"])
    cmap = chrom_map(sess, lens)
    loci = [l for l in read_tsv(P["derived"] / "genomes" / PRAIRIE_GENOME / "loci.tsv") if l["family"] in ("MUP", "OBP")]
    fam_of = {}
    for gid, g in genes.items():
        acc = cmap.get(g["chrom"], "")
        hit = [l for l in loci if l["contig"] == acc and int(l["start"]) <= g["end"] and int(l["end"]) >= g["start"]]
        fams = {l["family"] for l in hit}
        if len(fams) == 1:
            fam_of[gid] = fams.pop()
        elif not hit and re.search("major urinary", g["description"], re.I):
            fam_of[gid] = "MUP"
        elif not hit and g["chrom"] == "chrX" and re.search("odorant", g["description"], re.I):
            fam_of[gid] = "OBP"
    best = {}
    for m, (gid, ex) in mrna.items():
        if gid in fam_of and ex:
            ln = sum(b - a + 1 for a, b in ex)
            if gid not in best or ln > best[gid][1]:
                best[gid] = (ex, ln)
    contigs = {cmap[genes[g]["chrom"]] for g in best}
    g4 = next(x for x in cfg["genomes"] if x["slug"] == PRAIRIE_GENOME)
    seqs = {}
    for name, _d, seq in fasta_records(genome_fasta(cfg, g4)):
        if name in contigs:
            seqs[name] = seq.upper()
            if len(seqs) == len(contigs):
                break
    tx, fams = {}, {}
    for gid, (ex, _ln) in best.items():
        g = genes[gid]
        s = "".join(seqs[cmap[g["chrom"]]][a - 1:b] for a, b in ex)
        if g["strand"] == "-":
            s = revcomp(s)
        key = f"egapx:{g['gene']}"
        tx[key] = (f"{g['gene']}|{g['description']} (EGAPx)", s)
        fams[key] = fam_of[gid]
    _CACHE["prairie"] = (tx, fams)
    return tx, fams


def main() -> None:
    cfg = load_config()
    P = paths(cfg)
    sess = requests.Session()
    sess.headers["User-Agent"] = cfg["user_agent"]
    files = fetch(P, sess)
    out = P["derived"] / "egapx"
    out.mkdir(parents=True, exist_ok=True)

    # 1. lipocalin genes in both annotations
    parsed = {sp: parse(f) for sp, f in files.items()}
    rows = []
    for sp, (genes, _m, _l) in parsed.items():
        for gid, g in genes.items():
            rows.append({"species": sp, "id": gid, **g})
    write_tsv(out / "lipocalin_genes.tsv", rows)
    for sp in parsed:
        n = sum(1 for r in rows if r["species"] == sp)
        obpx = sum(1 for r in rows if r["species"] == sp and r["chrom"] == "chrX" and re.search("odorant", r["description"], re.I))
        mup = sum(1 for r in rows if r["species"] == sp and re.search("urinary", r["description"], re.I))
        print(f"{sp}: {n} lipocalin-family genes; MUP-named {mup}; X-linked OBP-named {obpx}", flush=True)

    # 2. prairie EGAPx genes vs DNA-search loci (same assembly)
    genes, mrna, lens = parsed["prairie"]
    cmap = chrom_map(sess, lens)
    loci = [l for l in read_tsv(P["derived"] / "genomes" / PRAIRIE_GENOME / "loci.tsv") if l["family"] in ("MUP", "OBP")]
    used = set()
    cmp_rows = []
    for gid, g in sorted(genes.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["start"])):
        acc = cmap.get(g["chrom"], "")
        hit = [l for l in loci if l["contig"] == acc and int(l["start"]) <= g["end"] and int(l["end"]) >= g["start"]]
        for l in hit:
            used.add(l["locus"])
        cmp_rows.append({"egapx_gene": g["gene"], "description": g["description"], "type": g["type"],
                         "chrom": g["chrom"], "genbank": acc, "start": g["start"], "end": g["end"],
                         "dnasearch_loci": ",".join(f"{l['locus']}:{l['family']}:{l['status']}" for l in hit)})
    for l in loci:
        if l["locus"] not in used:
            cmp_rows.append({"egapx_gene": "", "description": "no EGAPx gene", "type": "", "chrom": "",
                             "genbank": l["contig"], "start": l["start"], "end": l["end"],
                             "dnasearch_loci": f"{l['locus']}:{l['family']}:{l['status']}"})
    write_tsv(out / "prairie_vs_dnasearch.tsv", cmp_rows)
    for fam in ("MUP", "OBP"):
        fl = [l for l in loci if l["family"] == fam and l["status"] != "fragment"]
        matched = [l for l in fl if l["locus"] in used]
        print(f"prairie {fam}: {len(fl)} DNA-search loci (non-fragment), {len(matched)} overlap an EGAPx gene "
              f"({sum(l['status'] == 'intact' for l in matched)} of {sum(l['status'] == 'intact' for l in fl)} intact)",
              flush=True)

    # 3. prairie liver reads against the EGAPx MUP and OBP transcripts
    run = next(r for r in s6.RUNS if r["run"] == LIVER_RUN)
    tx = s6.transcripts(cfg, P, run["slug"])
    allg = read_tsv(P["derived"] / "summary" / "genes.tsv")
    clade = {g["id"].split("|")[-1]: g["clade"] for g in allg if g["slug"] == run["slug"]}
    tx = {k: v for k, v in tx.items() if clade.get(k) not in ("MUP", "OBP")}
    etx, efam = egapx_models(cfg, P)
    tx.update(etx)
    famkey = {k: {"MUP": "MUPFAM", "OBP": "OBPFAM"}[v] for k, v in efam.items()}
    fam = {k: famkey.get(k, k) for k in tx}
    fq = P["raw"] / f"ena_{LIVER_RUN}" / f"{LIVER_RUN}_1.head.fastq.gz"
    n, c = s6.count_reads(fq, tx, fam)
    crow = []
    for k, (lab, seq) in list(tx.items()) + [(f, (s6.FAMILY_LABEL[f], "")) for f in s6.FAMILY_LABEL]:
        fm = k if k in s6.FAMILY_LABEL else fam.get(k, k)
        crow.append({"gene_id": k, "symbol": lab.split("|")[0], "family": {"MUPFAM": "MUP", "OBPFAM": "OBP"}.get(fm, ""),
                     "reads": c[k], "reads_total": n, "rpm": round(c[k] / n * 1e6, 1)})
    write_tsv(out / "prairie_liver_counts.tsv", crow)
    for f1, name in (("MUPFAM", "MUP"), ("OBPFAM", "OBP")):
        tot = sum(r["reads"] for r in crow if r["family"] == name)
        ng = sum(1 for k in famkey.values() if k == f1)
        print(f"prairie liver vs EGAPx {name} ({ng} transcripts): {tot} reads, {tot / n * 1e6:.0f} RPM", flush=True)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
