"""Step 04: whole-genome search of the assemblies that have no RefSeq gene annotation (WSL/Mac: miniprot).

Why: the RefSeq annotations (Step 02) cover only some genomes. The UC Davis prairie vole assembly (2026, long reads,
far more complete than MicOch1.0) and the field, common and montane vole assemblies have no gene annotation, so their
lipocalin genes are found here directly from the DNA. Mus spretus (a mouse with a known MUP cluster) is a positive
control: the method should find many MUP genes there.

What it does, per assembly in config 'genomes':
1. miniprot aligns the query proteins (Step 03 queries.faa: MUPs, OBPs, aphrodisin, other lipocalins as decoys) plus
   the proteins of the flanking genes of the three mouse clusters (taken from the closest annotated relative: prairie
   vole for voles, mouse for Mus, desert hamster for the Djungarian hamster) to the whole genome. The genome goes to
   miniprot in pieces of about 400 Mb (whole sequences), which keeps memory near 4 GB.
2. Lipocalin loci: overlapping alignments are grouped into loci; each locus takes the family (MUP / OBP / other) of
   its best-scoring query and is called intact / disrupted / partial as in Step 03.
3. Regions: each flanking gene is placed at its best alignment (rank 1, >= 50% of the protein aligned). The region of
   each cluster is found from the placed flanking genes as in Step 03; assembly gaps (runs of >= 10 N) in it are
   counted from the genome sequence.

Outputs (derived/genomes/<slug>/): miniprot.gff, loci.tsv, flank_hits.tsv, regions.tsv, gaps.tsv
Usage: run_wsl.sh step04_genomes.py [--only slug1,slug2] [--chunk-mb 400]
"""
from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from step03_regions import ANCHORS, call_status
from vm_common import (cluster_loci, fasta_records, genome_fasta, load_config, lower_priority, n_runs,
                       parse_miniprot_gff, paths, read_tsv, refseq_files, run, write_fasta, write_tsv)

FLANK_SOURCE = {"Arvicolinae": "microtus_ochrogaster", "Murinae": "mus_musculus", "Cricetinae": "phodopus_roborovskii"}
MAX_REGION = 4_000_000


def flank_queries(cfg, P, src_slug: str) -> list[tuple[str, str]]:
    import gzip
    flanks = read_tsv(P["derived"] / "regions" / "flanks.tsv")
    want = {f["symbol"].upper(): f for f in flanks}
    lip = {}
    # longest protein per wanted symbol, from the source species' annotation
    genes = {}
    with gzip.open(P["derived"] / "annotated" / f"genes_{src_slug}.tsv.gz", "rt", encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["symbol"].upper() in want:
                genes[r["gene_id"]] = r["symbol"]
    # protein -> gene via the GFF is in step02; re-derive cheaply from the protein FASTA headers is not possible, so
    # use the all_longest_proteins.faa written by step02 (key = slug|gene_id)
    for key, _d, seq in fasta_records(P["derived"] / "annotated" / "all_longest_proteins.faa"):
        slug, gid = key.split("|")
        if slug == src_slug and gid in genes:
            sym = genes[gid]
            if sym not in lip or len(seq) > len(lip[sym]):
                lip[sym] = seq
    out = []
    for sym, seq in lip.items():
        f = want[sym.upper()]
        out.append((f"F|{f['region']}|{f['side']}{f['rank']}|{f['symbol']}", seq))
    return out


def run_miniprot(genome: Path, queries: Path, out_gff: Path, threads, chunk_bp: int) -> None:
    part = Path(str(out_gff) + ".part")
    with tempfile.TemporaryDirectory(dir=str(Path.home())) as tmp, open(part, "w") as out:
        st = {"fh": None, "size": 0, "n": 0}

        def flush():
            if st["fh"] is None:
                return
            st["fh"].close()
            st["n"] += 1
            run(["miniprot", "-t", threads, "-G", "200k", "-j", "2", "--gff", "--trans", "-N", "60", "--outs=0.3",
                 "-p", "0.3", "-P", f"c{st['n']}_", st["fh"].name, queries], stdout=out)
            out.flush()
            Path(st["fh"].name).unlink()
            st["fh"], st["size"] = None, 0
            print(f"  chunk {st['n']} done", flush=True)

        for name, _d, seq in fasta_records(genome):
            if st["fh"] is None:
                st["fh"] = open(Path(tmp) / "chunk.fa", "w")
            st["fh"].write(f">{name}\n")
            for i in range(0, len(seq), 100_000):
                st["fh"].write(seq[i:i + 100_000] + "\n")
            st["size"] += len(seq)
            if st["size"] >= chunk_bp:
                flush()
        flush()
    part.replace(out_gff)


def regions_from_flanks(fh_rows: list[dict], seqlens: dict) -> list[dict]:
    """Same logic as step03_regions.locate(), on miniprot placements of the flanking genes."""
    out = []
    for name in ANCHORS:
        hits = [h for h in fh_rows if h["region"] == name]
        by = {}
        for h in hits:
            by.setdefault(h["contig"], []).append(h)
        if not by:
            out.append({"region": name, "part": 0, "found": False})
            continue
        order = sorted(by, key=lambda c: -len(by[c]))
        best = by[order[0]]
        sides = {h["side"] for h in best}
        parts = []
        if sides == {"L", "R"}:
            L = sorted([h for h in best if h["side"] == "L"], key=lambda h: h["rank"])
            R = sorted([h for h in best if h["side"] == "R"], key=lambda h: h["rank"])
            li, ri = L[0], R[0]
            a, b = (li["end"], ri["start"]) if li["start"] < ri["start"] else (ri["end"], li["start"])
            parts.append({"contig": li["contig"], "start": a, "end": b, "left_flank": li["symbol"],
                          "right_flank": ri["symbol"], "n_left": len(L), "n_right": len(R), "open": ""})
        else:
            for c in order[:2]:
                hs = sorted(by[c], key=lambda h: h["rank"])
                if c != order[0] and len(hs) < 2:
                    continue
                inner, far = hs[0], hs[-1]
                direction = 0 if len(hs) == 1 else (1 if inner["start"] > far["start"] else -1)
                sl = seqlens.get(c, inner["end"] + MAX_REGION)
                if direction >= 0:
                    a, b = inner["end"], min(sl, inner["end"] + MAX_REGION)
                else:
                    a, b = max(1, inner["start"] - MAX_REGION), inner["start"]
                parts.append({"contig": c, "start": a, "end": b,
                              "left_flank": inner["symbol"] if inner["side"] == "L" else "",
                              "right_flank": inner["symbol"] if inner["side"] == "R" else "",
                              "n_left": sum(h["side"] == "L" for h in hs), "n_right": sum(h["side"] == "R" for h in hs),
                              "open": "toward contig end" if direction >= 0 else "toward contig start"})
        for i, p in enumerate(parts, 1):
            out.append(dict(p, region=name, part=i, found=True, length=p["end"] - p["start"]))
    return out


def one_genome(cfg, P, g: dict, chunk_bp: int) -> None:
    out = P["derived"] / "genomes" / g["slug"]
    out.mkdir(parents=True, exist_ok=True)
    genome = genome_fasta(cfg, g)
    src = FLANK_SOURCE.get(g["group"], "microtus_ochrogaster")
    qfile = out / "queries.faa"
    lipq = list(fasta_records(P["derived"] / "queries" / "queries.faa"))
    write_fasta(qfile, [(n, s) for n, _d, s in lipq] + flank_queries(cfg, P, src))
    gff = out / "miniprot.gff"
    if not gff.exists():
        print(f"{g['slug']}: miniprot on {genome.name}", flush=True)
        run_miniprot(genome, qfile, gff, cfg["threads"], chunk_bp)
    models = parse_miniprot_gff(gff)
    lip_models = [m for m in models if m["query"].startswith("Q|")]
    flank_models = [m for m in models if m["query"].startswith("F|")]
    # flanking genes: best placement per query
    fh_rows = []
    best = {}
    for m in flank_models:
        if m["qcov"] < 0.5:
            continue
        if m["query"] not in best or m["score"] > best[m["query"]]["score"]:
            best[m["query"]] = m
    for q, m in best.items():
        _f, region, sr, sym = q.split("|")
        fh_rows.append({"region": region, "side": sr[0], "rank": int(sr[1:]), "symbol": sym, "contig": m["contig"],
                        "start": m["start"], "end": m["end"], "strand": m["strand"], "identity": round(m["identity"], 3),
                        "qcov": round(m["qcov"], 3)})
    write_tsv(out / "flank_hits.tsv", fh_rows)
    # lipocalin loci
    loci = []
    for i, locus in enumerate(cluster_loci(lip_models), 1):
        b = locus[0]
        bm = max((m for m in locus if "|MUP|" in m["query"]), key=lambda m: m["score"], default=None)
        bo = max((m for m in locus if "|OBP|" in m["query"]), key=lambda m: m["score"], default=None)
        loci.append({"slug": g["slug"], "locus": i, "contig": b["contig"], "start": b["start"], "end": b["end"],
                     "strand": b["strand"], "family": b["query"].split("|")[1], "best_query": b["query"],
                     "score": b["score"], "identity": round(b["identity"], 3), "qcov": round(b["qcov"], 3),
                     "frameshifts": b["frameshifts"], "stops": b["stops"], "n_exons": len(b["exons"]),
                     "status": call_status(b), "mup_score": bm["score"] if bm else 0,
                     "mup_identity": round(bm["identity"], 3) if bm else "",
                     "obp_score": bo["score"] if bo else 0, "translation": b["translation"]})
    # sequence lengths + region gaps (second pass over the genome, only the contigs needed)
    contigs = {h["contig"] for h in fh_rows}
    seqlens, seqs = {}, {}
    for name, _d, seq in fasta_records(genome):
        seqlens[name] = len(seq)
        if name in contigs:
            seqs[name] = seq
    regs = regions_from_flanks(fh_rows, seqlens)
    gaps = []
    for r in regs:
        if not r["found"]:
            continue
        r["contig_length"] = seqlens.get(r["contig"], 0)
        seq = seqs.get(r["contig"], "")
        a, b = max(1, r["start"]), min(len(seq), r["end"])
        for x, y in n_runs(seq[a - 1:b]):
            gaps.append({"region": r["region"], "part": r["part"], "contig": r["contig"], "start": x + a - 1,
                         "end": y + a - 1, "length": y - x + 1})
        r["n_gaps"] = sum(1 for z in gaps if z["region"] == r["region"] and z["part"] == r["part"])
        r["gap_bp"] = sum(z["length"] for z in gaps if z["region"] == r["region"] and z["part"] == r["part"])
        r["n_lipocalin_loci"] = sum(1 for l in loci if l["contig"] == r["contig"] and r["start"] <= l["start"] <= r["end"])
    for l in loci:
        l["in_region"] = ",".join(f"{r['region']}" for r in regs if r.get("found") and r["contig"] == l["contig"]
                                  and r["start"] - 50_000 <= l["start"] <= r["end"] + 50_000)
        l["contig_length"] = seqlens.get(l["contig"], 0)
    write_tsv(out / "loci.tsv", loci)
    write_tsv(out / "regions.tsv", regs, ["region", "part", "found", "contig", "contig_length", "start", "end",
                                          "length", "left_flank", "right_flank", "n_left", "n_right", "open",
                                          "n_gaps", "gap_bp", "n_lipocalin_loci"])
    write_tsv(out / "gaps.tsv", gaps, ["region", "part", "contig", "start", "end", "length"])
    tot = sum(seqlens.values())
    write_tsv(out / "assembly.tsv", [{"slug": g["slug"], "accession": g["accession"], "file": genome.name,
                                      "n_sequences": len(seqlens), "total_bp": tot}])
    print(f"{g['slug']}: {len(loci)} lipocalin loci ({sum(l['family'] == 'MUP' for l in loci)} best-matching a MUP), "
          f"{len(fh_rows)} flanking genes placed", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--chunk-mb", type=int, default=400)
    a = ap.parse_args()
    lower_priority()
    cfg = load_config()
    P = paths(cfg)
    only = set(filter(None, a.only.split(",")))
    for g in cfg["genomes"]:
        if only and g["slug"] not in only:
            continue
        one_genome(cfg, P, g, a.chunk_mb * 1_000_000)


if __name__ == "__main__":
    main()
