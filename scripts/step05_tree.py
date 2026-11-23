"""Step 05: lipocalin gene trees, to sort every candidate into the MUP clade, the OBP/aphrodisin clade or another
lipocalin group (WSL/Mac: blastp, MAFFT, trimAl, IQ-TREE).

Why: gene names given by annotation pipelines ('major urinary protein 4-like') come from the single best BLAST
match and can mislead in a family this large. A gene tree groups each protein with its closest relatives across all
species, so a vole gene that is a true MUP should fall inside the clade that holds the mouse and rat MUPs and the
horse, cat and dog MUP allergens (Logan et al. 2008 showed most mammals keep one ancestral Mup).

What it does
1. Candidates: the PF00061 domain of every RefSeq lipocalin candidate (Step 02), the reference proteins from UniProt,
   and the translations of lipocalin loci found in unannotated assemblies (Step 04; intact or disrupted, with
   premature stops and frameshift-positions masked as X; partial loci shorter than 120 aa are left out).
2. Cytosolic fatty-acid-binding proteins (FABP/CRABP/CRBP, which share the Pfam profile but are not secreted
   lipocalins) are removed: blastp of every candidate against the mouse + human lipocalins; a candidate whose best
   hit is one of these is dropped.
3. Full tree: MAFFT (--auto) alignment, trimAl (-gappyout) trimming, IQ-TREE 3 maximum likelihood with ModelFinder
   (-mset LG,WAG,JTT) and 1000 ultrafast bootstrap replicates. Rooted on the ApoD lipocalins (an old group, outside
   the MUP/OBP/LCN groups in published lipocalin trees).
4. Clades: MUP clade = the smallest clade holding all reference MUPs (mouse, rat, horse, cat, dog); OBP clade = the
   smallest clade holding hamster aphrodisin, the three bank vole OBPs and mouse Obp1a/Obp1b (the X-linked OBP1 /
   aphrodisin family); OBP2 clade = the smallest clade holding mouse and human OBP2A/OBP2B. Every candidate is
   labelled MUP, OBP, OBP2 or other, with its nearest labelled reference (blastp).
5. Focused tree for figures: the MUP, OBP and OBP2 clades plus their nearest outgroup lipocalins, realigned with MAFFT
   L-INS-i and re-estimated as in 3.

Outputs (derived/tree/): candidates.tsv, full/ (alignment, iqtree files), focused/ (same), classification.tsv
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

from vm_common import fasta_records, load_config, lower_priority, paths, read_tsv, run, write_fasta, write_tsv

FABP = re.compile(r"^(Fabp\d*|FABP\d*|Crabp\d|CRABP\d|Rbp[1257]|RBP[1257]|Pmp2|PMP2)$")
MUP_REF = re.compile(r"(Mus_musculus\|Mup\d|Rattus_norvegicus\|.*|Equ_c_1|Fel_d_4|Can_f_6)")


def short_label(r: dict) -> str:
    sym = r["symbol"] or ""
    return sym


def candidates(cfg, P) -> tuple[list[tuple[str, str]], list[dict]]:
    lip = {r["key"]: r for r in read_tsv(P["derived"] / "annotated" / "lipocalins.tsv")}
    seqs, meta = [], []
    for key, _d, s in fasta_records(P["derived"] / "annotated" / "lipocalin_domains.faa"):
        r = lip[key]
        tid = f"{r['slug']}|{r['gene_id']}"
        seqs.append((tid, s))
        meta.append({"tip": tid, "slug": r["slug"], "source": "refseq", "symbol": r["symbol"],
                     "description": r["description"], "chrom": r["chrom"], "seqid": r["seqid"], "start": r["start"],
                     "end": r["end"], "strand": r["strand"], "length": r["length"], "status": "annotated",
                     "rnaseq_cov": r["rnaseq_cov"], "rnaseq_intron_samples": r["rnaseq_intron_samples"]})
    for name, _d, s in fasta_records(P["raw"] / "uniprot_reference" / "reference_lipocalins.fasta"):
        _up, lab, acc = name.split("|", 2)[0], name.split("|")[1], name.split("|")[3]
        tid = f"uniprot|{acc}"
        seqs.append((tid, s))
        meta.append({"tip": tid, "slug": "uniprot", "source": "uniprot", "symbol": name.split("|")[2],
                     "description": lab, "length": len(s), "status": "reference"})
    gdir = P["derived"] / "genomes"
    for g in cfg["genomes"]:
        f = gdir / g["slug"] / "loci.tsv"
        if not f.exists():
            continue
        for l in read_tsv(f):
            t = l["translation"].upper().replace("*", "X")
            if len(t.replace("X", "")) < 120:
                continue
            tid = f"{g['slug']}|locus{l['locus']}"
            seqs.append((tid, t))
            meta.append({"tip": tid, "slug": g["slug"], "source": "genome", "symbol": f"locus{l['locus']}",
                         "description": f"miniprot locus, best query {l['best_query']}", "chrom": l["contig"],
                         "seqid": l["contig"], "start": l["start"], "end": l["end"], "strand": l["strand"],
                         "length": len(t), "status": l["status"], "frameshifts": l["frameshifts"],
                         "stops": l["stops"]})
    return seqs, meta


def blast_label(P, seqs, meta, out: Path, threads) -> dict:
    """Best blastp hit of every candidate among the mouse + human lipocalins (by symbol)."""
    ref = [(t, s) for (t, s), m in zip(seqs, meta) if m["slug"] in ("mus_musculus", "homo_sapiens")]
    sym = {m["tip"]: m["symbol"] for m in meta}
    write_fasta(out / "ref.faa", ref)
    write_fasta(out / "cand.faa", seqs)
    run(["makeblastdb", "-in", out / "ref.faa", "-dbtype", "prot", "-out", out / "refdb"])
    tsv = out / "blast.tsv"
    run(["blastp", "-query", out / "cand.faa", "-db", out / "refdb", "-outfmt", "6 qseqid sseqid pident bitscore",
         "-max_target_seqs", "5", "-evalue", "1e-3", "-num_threads", threads, "-out", tsv])
    best = {}
    for line in open(tsv):
        q, s, pid, bits = line.split("\t")
        if q == s:
            continue
        if q not in best or float(bits) > best[q][2]:
            best[q] = (sym.get(s, s), float(pid), float(bits))
    return best


def tree(aln_in: Path, prefix: Path, threads, focused: bool = False) -> Path:
    """MAFFT + trimAl + IQ-TREE. The full family tree (~700 proteins) is only used to sort proteins into clades, so it
    uses LG+G4 with IQ-TREE's fast search and SH-aLRT branch support; the focused tree uses MAFFT L-INS-i,
    ModelFinder (LG/WAG/JTT families) and 1000 ultrafast bootstraps + 1000 SH-aLRT replicates."""
    aln = prefix.with_suffix(".aln.faa")
    if not aln.exists() or aln.stat().st_mtime < aln_in.stat().st_mtime:
        with open(str(aln) + ".part", "w") as fh:
            args = ["mafft", "--thread", threads, "--quiet"] + (["--localpair", "--maxiterate", "1000"] if focused
                                                                 else ["--auto"])
            run(args + [aln_in], stdout=fh)
        Path(str(aln) + ".part").replace(aln)
        for f in prefix.parent.glob(prefix.name + ".*"):
            if f.suffix not in (".faa",):
                f.unlink()
    trim = prefix.with_suffix(".trim.faa")
    run(["trimal", "-in", aln, "-out", trim, "-gappyout"])
    if not Path(str(prefix) + ".treefile").exists():
        model = ["-m", "MFP", "-mset", "LG,WAG,JTT", "-B", "1000", "-alrt", "1000"] if focused else             ["-m", "LG+G4", "--fast", "-alrt", "1000"]
        run(["iqtree3", "-s", trim, *model, "-T", threads, "--seed", "1", "--prefix", prefix, "-quiet", "-redo"])
    return Path(str(prefix) + ".treefile")


MUP_POSITION: set[str] = set()     # tips of the lipocalins inside any genome's MUP neighbourhood (filled in main)


def write_if_changed(path: Path, records: list[tuple[str, str]]) -> None:
    """Rewrite an input FASTA only when its sequences change, so that a rerun does not redo a tree (the tree steps
    redo the alignment and the tree when the input file is newer than the alignment)."""
    if path.exists() and [(n, s) for n, _d, s in fasta_records(path)] == list(records):
        return
    write_fasta(path, records)


def clade_sets(t, meta: list[dict]) -> dict:
    """MUP, OBP and OBP2 clades of a rooted tree (Bio.Phylo), from the reference proteins in it.

    MUP: the smallest clade holding the mouse MUPs, the horse, cat and dog MUP allergens, and every lipocalin that sits
    in the MUP neighbourhood (between the genes that flank the mouse cluster) of any genome. Anchoring on position as
    well as on the references keeps the clade stable when the deep branching order inside it is uncertain: in the
    focused tree, the hamster and vole MUPs branch just outside the mouse + horse + carnivore MUPs, which the root
    then splits apart; the clade so defined still holds nothing but MUP-named genes (checked in the report)."""
    tips = {c.name: c for c in t.get_terminals()}

    def pick(test):
        return [tips[m["tip"]] for m in meta if m["tip"] in tips and test(m)]
    refs = {
        "MUP": pick(lambda m: (m["slug"] == "mus_musculus" and re.match(r"^Mup\d+$", m["symbol"])) or
                    (m["slug"] == "uniprot" and m["symbol"] in ("Equ_c_1_MUP", "Fel_d_4_MUP", "Can_f_6_MUP")) or
                    m["tip"] in MUP_POSITION),
        # the X-linked OBP1 / aphrodisin family: hamster aphrodisin, the bank vole OBPs and mouse Obp1a/Obp1b
        "OBP": pick(lambda m: (m["slug"] == "uniprot" and (m["symbol"].startswith("Obp") or m["symbol"] == "aphrodisin"))
                    or (m["slug"] == "mus_musculus" and m["symbol"] in ("Obp1a", "Obp1b"))),
        "OBP2": pick(lambda m: m["slug"] in ("mus_musculus", "homo_sapiens") and m["symbol"].upper() in ("OBP2A", "OBP2B")),
    }
    out = {}
    for k, r in refs.items():
        if len(r) >= 2:
            ca = t.common_ancestor(r)
            out[k] = (ca, {c.name for c in ca.get_terminals()})
    return out


def label(tips, sets) -> dict:
    lab = {}
    for name in tips:
        lab[name] = "other"
        for k in ("OBP2", "OBP", "MUP"):          # most specific last wins
            if k in sets and name in sets[k][1]:
                lab[name] = k
    return lab


def root_on(t, names: list[str]) -> None:
    tips = [c for c in t.get_terminals() if c.name in set(names)]
    try:
        if len(tips) == 1:
            t.root_with_outgroup(tips[0])
        elif tips:
            t.root_with_outgroup(t.common_ancestor(tips))
        else:
            t.root_at_midpoint()
    except Exception:
        t.root_at_midpoint()


def info_of(sets) -> dict:
    d = {}
    for k, (ca, names) in sets.items():
        d[f"{k}_clade_size"] = len(names)
        # IQ-TREE writes 'SH-aLRT/UFBoot' (e.g. 98.5/100) as the node label, which Bio.Phylo keeps as the name
        d[f"{k}_support"] = ca.confidence if ca.confidence is not None else (ca.name or "")
    return d


def main() -> None:
    from Bio import Phylo
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-focused", action="store_true", help="full tree + classification only")
    args = ap.parse_args()
    lower_priority()
    cfg = load_config()
    P = paths(cfg)
    th = str(cfg["threads"])
    out = P["derived"] / "tree"
    (out / "full").mkdir(parents=True, exist_ok=True)
    (out / "focused").mkdir(parents=True, exist_ok=True)
    seqs, meta = candidates(cfg, P)
    seqd = dict(seqs)
    # lipocalins inside the MUP neighbourhood of any genome (RefSeq genes from step 03, DNA loci from step 04)
    for g in read_tsv(P["derived"] / "regions" / "region_genes.tsv"):
        if g["region"] == "MUP" and g["inside"] == "True" and g["lipocalin"] == "True":
            MUP_POSITION.add(f"{g['slug']}|{g['gene_id']}")
    for gen in cfg["genomes"]:
        f = P["derived"] / "genomes" / gen["slug"] / "loci.tsv"
        if f.exists():
            for l in read_tsv(f):
                if "MUP" in (l.get("in_region") or "") and l["status"] != "fragment":
                    MUP_POSITION.add(f"{gen['slug']}|locus{l['locus']}")
    best = blast_label(P, seqs, meta, out, th)
    keep = []
    for (tid, s), m in zip(seqs, meta):
        b = best.get(tid)
        m["best_ref"], m["best_ref_pid"] = (b[0], round(b[1], 1)) if b else ("", "")
        ref_sym = m["symbol"] if m["slug"] in ("mus_musculus", "homo_sapiens") else m["best_ref"]
        m["fabp_like"] = bool(FABP.match(ref_sym or "")) or bool(FABP.match(m["symbol"] or ""))
        if not m["fabp_like"] and len(s.replace("X", "")) >= 100 and m["source"] != "genome":
            keep.append((tid, s))
    keep_ids = {t for t, _ in keep}
    # 1. full family tree (RefSeq + UniProt), rooted on ApoD
    write_if_changed(out / "full" / "input.faa", keep)
    full = tree(out / "full" / "input.faa", out / "full" / "lipocalins", th)
    t = Phylo.read(str(full), "newick")
    root_on(t, [m["tip"] for m in meta if m["symbol"].upper() == "APOD" and m["tip"] in keep_ids])
    sets = clade_sets(t, meta)
    lab_full = label([c.name for c in t.get_terminals()], sets)
    print("full tree:", info_of(sets), flush=True)
    Phylo.write(t, str(out / "full" / "lipocalins.rooted.nwk"), "newick")
    # 2. focused tree: the MUP and OBP clades and everything between them, the OBP2 clade, 20 mouse/human lipocalins
    #    from the rest of the family as outgroup, and the DNA-search loci of the MUP/OBP families (all MUP-family loci;
    #    the OBP/OBP2-family loci of the UC Davis prairie vole only, to keep the tree small)
    both = t.common_ancestor([sets["MUP"][0], sets["OBP"][0]])
    focus = {c.name for c in both.get_terminals()} | (sets["OBP2"][1] if "OBP2" in sets else set())
    rest = [c.name for c in t.get_terminals() if c.name not in focus]
    sym_of = {m["tip"]: m["symbol"].upper() for m in meta}
    mh = [x for x in rest if x.split("|")[0] in ("mus_musculus", "homo_sapiens")]
    outg = ([x for x in mh if sym_of.get(x) == "APOD"] + [x for x in mh if sym_of.get(x) != "APOD"])[:20]
    q_clade = {m["symbol"]: lab_full.get(m["tip"]) for m in meta
               if m["slug"] in ("mus_musculus", "rattus_norvegicus", "uniprot") and m["tip"] in lab_full}
    loci = []
    for m in meta:
        if m["source"] != "genome" or m["fabp_like"] or m["status"] == "fragment":
            continue
        qc = q_clade.get(m["description"].rsplit("|", 1)[-1], "")
        if qc == "MUP" or (qc in ("OBP", "OBP2") and m["slug"].startswith("microtus_ochrogaster_ucd1")):
            loci.append(m["tip"])
    sel = focus | set(outg) | set(loci)
    fsets, lab_foc, foc = {}, {}, None
    if not args.skip_focused:
        write_if_changed(out / "focused" / "input.faa", [(x, seqd[x]) for x in seqd if x in sel])
        write_tsv(out / "focused" / "outgroup.tsv", [{"tip": x} for x in outg])
        foc = tree(out / "focused" / "input.faa", out / "focused" / "mup_obp", th, focused=True)
        tf = Phylo.read(str(foc), "newick")
        root_on(tf, [x for x in outg if sym_of.get(x) == "APOD"][:1] or outg[:1])   # root on mouse ApoD
        fsets = clade_sets(tf, meta)
        lab_foc = label([c.name for c in tf.get_terminals()], fsets)
        Phylo.write(tf, str(out / "focused" / "mup_obp.rooted.nwk"), "newick")
        print("focused tree:", info_of(fsets), flush=True)
    rows = []
    for m in meta:
        m["in_full_tree"] = m["tip"] in keep_ids
        m["clade_full"] = lab_full.get(m["tip"], "")
        m["clade_focused"] = lab_foc.get(m["tip"], "")
        m["clade"] = m["clade_focused"] or m["clade_full"]
        m["in_focused_tree"] = m["tip"] in sel
        m["focused_outgroup"] = m["tip"] in outg
        rows.append(m)
    write_tsv(out / "classification.tsv", rows,
              ["tip", "slug", "source", "symbol", "description", "chrom", "seqid", "start", "end", "strand", "length",
               "status", "frameshifts", "stops", "rnaseq_cov", "rnaseq_intron_samples", "best_ref", "best_ref_pid",
               "fabp_like", "in_full_tree", "clade_full", "clade_focused", "clade", "in_focused_tree",
               "focused_outgroup"])
    write_tsv(out / "clades.tsv", [dict({f"full_{k}": v for k, v in info_of(sets).items()},
                                        **{f"focused_{k}": v for k, v in info_of(fsets).items()},
                                        focused_tips=len(sel), outgroup_tips=len(outg), genome_loci=len(loci))])
    print(f"focused tree: {len(sel)} tips -> {foc}", flush=True)


if __name__ == "__main__":
    main()
