"""Step 03: the MUP and OBP gene neighbourhoods in every annotated genome (synteny), and a DNA-level search there.

Why: a gene family can be missed by a gene annotation (pseudogenes, unusual genes) or lie in an assembly gap. Genes
that sit next to each other in the mouse usually sit next to each other in other rodents too ('synteny'), so the
mouse's neighbours of the Mup cluster mark where a vole's Mup genes should be. Looking at that stretch of DNA directly
shows whether the voles have MUP genes there, broken copies (pseudogenes), or a hole in the assembly.

Modes (run in order; --define and --fetch need internet and run anywhere, --search needs miniprot: WSL/Mac):
  --define  1. Query proteins for the DNA search (derived/queries/queries.faa): seven mouse MUPs chosen to span the
               family (central Mup1, peripheral Mup3-6, Mup20 = darcin, Mup21), the rat MUPs, the horse, cat and dog
               MUP allergens, the mouse and bank vole OBPs, hamster aphrodisin, water vole arvicolin, and every other
               secreted mouse lipocalin as a decoy, so that each DNA hit can be assigned to the family it matches best.
            2. Three mouse lipocalin clusters (anchors): MUP (Mup genes, chr 4), OBP1 (Obp1a/Obp1b, chr X) and LCN2
               (the chr 2 lipocalin cluster that holds Obp2a/Obp2b). For each, the 8 nearest protein-coding genes on
               each side (named genes only; no lipocalins, no LOC/Gm/Rik placeholders) are the 'flanking genes'.
            3. In each other species, the flanking genes are found by gene symbol (RefSeq names orthologs after their
               mouse or human gene). The region is the stretch between the innermost left and right flanking genes
               found on the same chromosome or scaffold. If only one side is found there (the scaffold ends), the
               region runs from that side to the scaffold end (an 'open' region: the assembly breaks there).
  --fetch   The DNA of each region (+-50 kb) from NCBI (E-utilities efetch) into raw/ncbi_regions/.
  --search  Assembly gaps (runs of >= 10 N) in each region, and miniprot (protein-to-genome alignment, splice-aware)
            of the query proteins to each region. Each locus is assigned the family of its best-scoring query, and
            called intact (no frameshift, no premature stop, >= 80% of the query aligned), disrupted (frameshift or
            stop: a possible pseudogene, or an assembly error) or partial.

Outputs (derived/regions/): anchors.tsv, flanks.tsv, regions.tsv, region_genes.tsv, gaps.tsv, region_loci.tsv
"""
from __future__ import annotations

import argparse
import gzip
import re
import time
from pathlib import Path

from vm_common import (cluster_loci, fasta_records, load_config, lower_priority, n_runs, parse_miniprot_gff, paths,
                       read_tsv, run, write_fasta, write_tsv)

PLACEHOLDER = re.compile(r"^(LOC\d+|Gm\d+|\d+Rik|[A-Z]\d+Rik|.*Rik$|ENS|TRNA|Trna|Mir\d|Snor|Rn\d|Vmn|Olfr|Or\d)", re.I)
LIPO_SYM = re.compile(r"^(Mup|Obp|Lcn|Ptgds|Rbp4|Ambp|Apod|Apom|Orm|C8g|Paep|Lcnl|Fabp|Crabp|Rbp[1-7]|Pmp2)", re.I)
FABP_SYM = re.compile(r"^(Fabp|Crabp|Rbp[1-3567]$|Pmp2)", re.I)
K_FLANK = 8
MARGIN = 50_000
MAX_REGION = 4_000_000
ANCHORS = {
    "MUP": {"desc": "Mup cluster (mouse chr 4)", "match": re.compile(r"^Mup(\d+|-ps\d+)$")},
    "OBP1": {"desc": "Obp1a/Obp1b (mouse chr X)", "match": re.compile(r"^Obp1[ab]$")},
    "LCN2": {"desc": "chr 2 lipocalin cluster with Obp2a/Obp2b (mouse chr 2)",
             "match": re.compile(r"^(Obp2[ab]|Lcn(3|5|6|8|9|10|11|12|15))$")},
}
MOUSE_MUP_QUERIES = ["Mup1", "Mup3", "Mup4", "Mup5", "Mup6", "Mup20", "Mup21"]


def read_genes(P, slug: str) -> list[dict]:
    rows = []
    with gzip.open(P["derived"] / "annotated" / f"genes_{slug}.tsv.gz", "rt", encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            r["start"], r["end"] = int(r["start"]), int(r["end"])
            r["seqlen"] = int(r["seqlen"] or 0)
            rows.append(r)
    return rows


# ----------------------------------------------------------------------------------------------------------------------
# --define
# ----------------------------------------------------------------------------------------------------------------------
def build_queries(cfg, P) -> None:
    lip = read_tsv(P["derived"] / "annotated" / "lipocalins.tsv")
    seqs = {k: s for k, _d, s in fasta_records(P["derived"] / "annotated" / "lipocalins_full.faa")}
    recs = []
    for r in lip:
        if r["key"] not in seqs or not r["hmm_evalue"]:
            continue
        sym, desc = r["symbol"], r["description"]
        if r["slug"] == "mus_musculus":
            if FABP_SYM.match(sym) or int(r["length"]) > 400:
                continue
            if sym.startswith("Mup"):
                if sym not in MOUSE_MUP_QUERIES:
                    continue
                fam = "MUP"
            elif sym.startswith("Obp"):
                fam = "OBP"
            else:
                fam = "other"
            recs.append((f"Q|{fam}|Mus_musculus|{sym}", seqs[r["key"]]))
        elif r["slug"] == "rattus_norvegicus" and re.search(r"major urinary|alpha-2u", desc, re.I):
            recs.append((f"Q|MUP|Rattus_norvegicus|{sym}", seqs[r["key"]]))
    for name, _d, seq in fasta_records(P["raw"] / "uniprot_reference" / "reference_lipocalins.fasta"):
        lab = name.split("|")
        fam = "MUP" if "MUP" in lab[2] else "OBP"
        recs.append((f"Q|{fam}|{lab[1]}|{lab[2]}", seq))
    # rat: keep at most 4 distinct MUPs (they are near-identical)
    out, rat = [], 0
    for n, s in recs:
        if "|Rattus_norvegicus|" in n:
            rat += 1
            if rat > 4:
                continue
        out.append((n, s))
    write_fasta(P["derived"] / "queries" / "queries.faa", out)
    print(f"queries: {len(out)} ({sum('|MUP|' in n for n, _ in out)} MUP, {sum('|OBP|' in n for n, _ in out)} OBP)")


def flanking(genes: list[dict], anchor_genes: list[dict]) -> tuple[list[dict], list[dict], tuple[str, int, int]]:
    seqid = anchor_genes[0]["seqid"]
    lo, hi = min(g["start"] for g in anchor_genes), max(g["end"] for g in anchor_genes)
    pc = [g for g in genes if g["seqid"] == seqid and g["biotype"] == "protein_coding" and g["symbol"]
          and not PLACEHOLDER.match(g["symbol"]) and not LIPO_SYM.match(g["symbol"])]
    left = sorted([g for g in pc if g["end"] < lo], key=lambda g: -g["end"])[:K_FLANK]
    right = sorted([g for g in pc if g["start"] > hi], key=lambda g: g["start"])[:K_FLANK]
    return left, right, (seqid, lo, hi)


def locate(genes: list[dict], left_syms: list[str], right_syms: list[str]) -> list[dict]:
    """Find the anchor region in one species from the flanking-gene symbols. Returns 0-2 region rows."""
    idx = {}
    for g in genes:
        idx.setdefault(g["symbol"].upper(), []).append(g)
    hits = []
    for side, syms in (("L", left_syms), ("R", right_syms)):
        for rank, s in enumerate(syms, 1):
            for g in idx.get(s.upper(), []):
                hits.append(dict(g, side=side, rank=rank, flank=s))
    if not hits:
        return []
    by_seq = {}
    for h in hits:
        by_seq.setdefault(h["seqid"], []).append(h)
    # the seqid holding most flanking genes; a second one if the other side is elsewhere
    order = sorted(by_seq, key=lambda s: -len(by_seq[s]))
    best = by_seq[order[0]]
    sides = {h["side"] for h in best}
    rows = []

    def one_side(hs, side):
        hs = sorted(hs, key=lambda h: h["rank"])
        inner = hs[0]
        far = hs[-1]
        seqlen = inner["seqlen"]
        # the cluster lies on the side of the innermost flank away from the outer ones
        if len(hs) == 1:
            direction = 0
        else:
            direction = 1 if inner["start"] > far["start"] else -1
        if direction >= 0:
            a, b = inner["end"], min(seqlen or inner["end"] + MAX_REGION, inner["end"] + MAX_REGION)
        else:
            a, b = max(1, inner["start"] - MAX_REGION), inner["start"]
        return {"seqid": inner["seqid"], "chrom": inner["chrom"], "seqlen": seqlen, "start": a, "end": b,
                "left_flank": inner["flank"] if side == "L" else "", "right_flank": inner["flank"] if side == "R" else "",
                "n_left": sum(h["side"] == "L" for h in hs), "n_right": sum(h["side"] == "R" for h in hs),
                "open": "toward scaffold end" if direction >= 0 else "toward scaffold start",
                "direction": direction}

    if sides == {"L", "R"}:
        L = sorted([h for h in best if h["side"] == "L"], key=lambda h: h["rank"])
        R = sorted([h for h in best if h["side"] == "R"], key=lambda h: h["rank"])
        li, ri = L[0], R[0]
        a, b = (li["end"], ri["start"]) if li["start"] < ri["start"] else (ri["end"], li["start"])
        rows.append({"seqid": li["seqid"], "chrom": li["chrom"], "seqlen": li["seqlen"], "start": a, "end": b,
                     "left_flank": li["flank"], "right_flank": ri["flank"], "n_left": len(L), "n_right": len(R),
                     "open": "", "direction": 1 if li["start"] < ri["start"] else -1})
    else:
        rows.append(one_side(best, sides.pop()))
        # the other side on another scaffold?
        other_side = "R" if rows[0]["left_flank"] else "L"
        for s in order[1:]:
            hs = [h for h in by_seq[s] if h["side"] == other_side]
            if len(hs) >= 2:
                rows.append(one_side(hs, other_side))
                break
    return rows


def define(cfg, P) -> None:
    build_queries(cfg, P)
    out = P["derived"] / "regions"
    mouse = read_genes(P, "mus_musculus")
    lip = read_tsv(P["derived"] / "annotated" / "lipocalins.tsv")
    lip_keys = {(r["slug"], r["gene_id"]): r for r in lip}
    anchors, flanks, regions, rgenes = [], [], [], []
    for name, a in ANCHORS.items():
        ag = [g for g in mouse if a["match"].match(g["symbol"])]
        seqs = {}
        for g in ag:
            seqs.setdefault(g["seqid"], []).append(g)
        ag = max(seqs.values(), key=len)              # e.g. Mup-ps22 on chr 11 is not part of the cluster
        L, R, (seqid, lo, hi) = flanking(mouse, ag)
        anchors.append({"region": name, "description": a["desc"], "seqid": seqid, "chrom": ag[0]["chrom"],
                        "start": lo, "end": hi, "n_genes": len(ag), "genes": ",".join(g["symbol"] for g in ag)})
        for side, gs in (("L", L), ("R", R)):
            for rank, g in enumerate(gs, 1):
                flanks.append({"region": name, "side": side, "rank": rank, "symbol": g["symbol"],
                               "description": g["description"], "mouse_start": g["start"], "mouse_end": g["end"]})
        ls, rs = [g["symbol"] for g in L], [g["symbol"] for g in R]
        for sp in cfg["refseq"]:
            genes = mouse if sp["slug"] == "mus_musculus" else read_genes(P, sp["slug"])
            if sp["slug"] == "mus_musculus":
                found = [{"seqid": seqid, "chrom": ag[0]["chrom"], "seqlen": ag[0]["seqlen"], "start": L[0]["end"],
                          "end": R[0]["start"], "left_flank": L[0]["symbol"], "right_flank": R[0]["symbol"],
                          "n_left": len(L), "n_right": len(R), "open": "", "direction": 1}]
            else:
                found = locate(genes, ls, rs)
            if not found:
                regions.append({"slug": sp["slug"], "region": name, "part": 0, "found": False})
                continue
            for part, rr in enumerate(found, 1):
                rid = f"{sp['slug']}__{name}__{part}"
                rr = dict(rr, slug=sp["slug"], region=name, part=part, found=True, region_id=rid,
                          length=rr["end"] - rr["start"])
                regions.append(rr)
                for g in genes:
                    if g["seqid"] == rr["seqid"] and g["end"] >= rr["start"] - MARGIN and g["start"] <= rr["end"] + MARGIN:
                        lp = lip_keys.get((sp["slug"], g["gene_id"]))
                        rgenes.append(dict(g, slug=sp["slug"], region=name, region_id=rid,
                                           inside=g["end"] >= rr["start"] and g["start"] <= rr["end"],
                                           lipocalin=bool(lp), flank=g["symbol"] in ls + rs))
        print(f"{name}: mouse {ag[0]['chrom']}:{lo}-{hi} ({len(ag)} genes); flanks L {ls[:3]}... R {rs[:3]}...")
    write_tsv(out / "anchors.tsv", anchors)
    write_tsv(out / "flanks.tsv", flanks)
    write_tsv(out / "regions.tsv", regions, ["region_id", "slug", "region", "part", "found", "seqid", "chrom", "seqlen",
                                             "start", "end", "length", "left_flank", "right_flank", "n_left",
                                             "n_right", "open", "direction"])
    write_tsv(out / "region_genes.tsv", rgenes, ["region_id", "slug", "region", "gene_id", "symbol", "description",
                                                 "biotype", "pseudo", "seqid", "chrom", "start", "end", "strand",
                                                 "inside", "lipocalin", "flank"])


# ----------------------------------------------------------------------------------------------------------------------
# --fetch
# ----------------------------------------------------------------------------------------------------------------------
def fetch(cfg, P) -> None:
    import requests
    s = requests.Session()
    s.headers["User-Agent"] = cfg["user_agent"]
    raw = P["raw"] / "ncbi_regions"
    raw.mkdir(parents=True, exist_ok=True)
    log = []
    for r in read_tsv(P["derived"] / "regions" / "regions.tsv"):
        if r["found"] != "True":
            continue
        a = max(1, int(r["start"]) - MARGIN)
        b = int(r["end"]) + MARGIN
        if int(r["seqlen"] or 0):
            b = min(b, int(r["seqlen"]))
        dest = raw / f"{r['region_id']}.fa"
        if not dest.exists():
            for attempt in range(6):
                time.sleep(0.4)
                try:
                    q = s.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", timeout=(30, 300),
                              params={"db": "nuccore", "id": r["seqid"], "seq_start": a, "seq_stop": b,
                                      "rettype": "fasta", "retmode": "text", "tool": "VoleMUPs"})
                    if q.status_code == 200 and q.text.startswith(">"):
                        part = dest.with_name(dest.name + ".part")
                        part.write_text(q.text, encoding="utf-8")
                        part.replace(dest)
                        break
                except requests.RequestException:
                    pass
                time.sleep(2 ** attempt)
            else:
                raise RuntimeError(f"efetch failed: {r['region_id']}")
        log.append({"region_id": r["region_id"], "seqid": r["seqid"], "fetch_start": a, "fetch_end": b,
                    "file": dest.name})
        print(f"{r['region_id']}: {r['seqid']}:{a}-{b}", flush=True)
    write_tsv(raw / "fetched.tsv", log)
    (raw / "README.txt").write_text("Genome regions around the MUP / OBP1 / LCN2 lipocalin clusters, fetched from NCBI "
                                    "nuccore with E-utilities efetch (rettype=fasta) by VoleMUPs/scripts/"
                                    "step03_regions.py --fetch. fetched.tsv lists accession and coordinates.\n",
                                    encoding="utf-8")


# ----------------------------------------------------------------------------------------------------------------------
# --search
# ----------------------------------------------------------------------------------------------------------------------
def call_status(m: dict) -> str:
    """fragment: under 40% of the query aligned (short matches: old gene debris or chance); disrupted: frameshift or
    premature stop (a pseudogene, an assembly error, or a splice junction miniprot could not model); intact: neither,
    and >= 80% of the query aligned; partial: the rest."""
    if m["qcov"] < 0.4:
        return "fragment"
    if m["frameshifts"] or m["stops"]:
        return "disrupted"
    if m["qcov"] >= 0.8:
        return "intact"
    return "partial"


def search(cfg, P) -> None:
    lower_priority()
    out = P["derived"] / "regions"
    raw = P["raw"] / "ncbi_regions"
    q = P["derived"] / "queries" / "queries.faa"
    fetched = {r["region_id"]: r for r in read_tsv(raw / "fetched.tsv")}
    regions = {r["region_id"]: r for r in read_tsv(out / "regions.tsv") if r["found"] == "True"}
    mp_dir = out / "miniprot"
    mp_dir.mkdir(exist_ok=True)
    gaps, loci = [], []
    for rid, f in fetched.items():
        reg = regions[rid]
        recs = list(fasta_records(raw / f["file"]))
        name, _d, seq = recs[0]
        off = int(f["fetch_start"]) - 1          # region coordinate -> chromosome coordinate
        for a, b in n_runs(seq):
            gaps.append({"region_id": rid, "slug": reg["slug"], "region": reg["region"], "seqid": reg["seqid"],
                         "start": a + off, "end": b + off, "length": b - a + 1,
                         "inside": int(reg["start"]) <= a + off <= int(reg["end"])})
        gff = mp_dir / f"{rid}.gff"
        if not gff.exists():
            with open(str(gff) + ".part", "w") as fh:
                run(["miniprot", "-t", cfg["threads"], "-G", "200k", "-j", "2", "--gff", "--trans", "-N", "60",
                     "--outs=0.3", "-p", "0.3", raw / f["file"], q], stdout=fh)
            Path(str(gff) + ".part").replace(gff)
        models = parse_miniprot_gff(gff)
        for i, locus in enumerate(cluster_loci(models), 1):
            best = locus[0]
            fam = best["query"].split("|")[1]
            best_mup = max((m for m in locus if "|MUP|" in m["query"]), key=lambda m: m["score"], default=None)
            best_obp = max((m for m in locus if "|OBP|" in m["query"]), key=lambda m: m["score"], default=None)
            loci.append({"region_id": rid, "slug": reg["slug"], "region": reg["region"], "locus": i,
                         "seqid": reg["seqid"], "start": best["start"] + off, "end": best["end"] + off,
                         "strand": best["strand"],
                         "inside": int(reg["start"]) <= best["start"] + off <= int(reg["end"]),
                         "family": fam, "best_query": best["query"], "score": best["score"],
                         "identity": round(best["identity"], 3), "qcov": round(best["qcov"], 3),
                         "frameshifts": best["frameshifts"], "stops": best["stops"], "n_exons": len(best["exons"]),
                         "status": call_status(best),
                         "mup_score": best_mup["score"] if best_mup else 0,
                         "mup_identity": round(best_mup["identity"], 3) if best_mup else "",
                         "obp_score": best_obp["score"] if best_obp else 0,
                         "translation": best["translation"]})
        print(f"{rid}: {len(seq):,} bp, {sum(1 for g in gaps if g['region_id'] == rid)} gaps, "
              f"{sum(1 for l in loci if l['region_id'] == rid)} lipocalin loci", flush=True)
    write_tsv(out / "gaps.tsv", gaps)
    write_tsv(out / "region_loci.tsv", loci)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--define", action="store_true")
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--search", action="store_true")
    a = ap.parse_args()
    cfg = load_config()
    P = paths(cfg)
    if a.define:
        define(cfg, P)
    if a.fetch:
        fetch(cfg, P)
    if a.search:
        search(cfg, P)


if __name__ == "__main__":
    main()
