"""Step 06: are the vole MUP-like and OBP genes switched on in the liver? A quick read-count check on public RNA-seq.

Why: mouse MUPs are made in the liver in enormous amounts (in adult males several per cent of all liver mRNA) and
then filtered into urine. If a vole's MUP-like gene does the same job, its transcripts should be common in vole liver.
NCBI's RefSeq notes (Step 02) say whether a gene model is supported by RNA-seq, but not in which tissue.

Data: male liver RNA-seq from the UC Davis prairie and meadow vole genome project (ENA PRJEB108798: ERR17322872,
ERR17322876), a male C57BL/6 mouse liver as positive control (DRR188168), and bank vole liver (SRR29282985).
Only the first ~120 MB of each read-1 FASTQ file is downloaded (about 1-3 million reads): enough to measure genes
that make up a sizeable share of the liver's mRNA, not enough for rare transcripts.

Counting: every 25-letter word (k-mer) of each target transcript (longest RefSeq transcript of every lipocalin gene of
that species, plus albumin, apolipoprotein A-I and beta-actin as yardsticks) is indexed, on both strands; k-mers shared
by two genes are dropped. A read is counted for a gene when at least 3 of its k-mers (sampled every 4 letters) belong
to that gene and to no other. Results are reads per million (RPM) and per kb of transcript per million (RPKM).

Usage: python scripts/step06_expression.py [--download] [--count]   (--download needs internet; both run anywhere)
Outputs: raw/ena_<run>/<run>_1.head.fastq.gz (+ README.txt); derived/expression/counts.tsv
"""
from __future__ import annotations

import argparse
import gzip
import re
import time
import zlib
from collections import Counter
from pathlib import Path

from vm_common import (fasta_records, gff_attrs, load_config, now, paths, read_tsv, refseq_files, unquote,
                       write_tsv)

# One liver per species (n = 1). Labels and sources are shown on the website; metadata from the ENA/BioSample records.
RUNS = [
    {"run": "ERR17322872", "slug": "microtus_ochrogaster", "label": "prairie vole, adult male liver",
     "source": "UC Davis vole genome project (PRJEB108798; Abuelanin et al. 2026); tissue from Cold Spring Harbor Laboratory, 2017"},
    {"run": "ERR17322876", "slug": "microtus_pennsylvanicus", "label": "meadow vole, adult male liver",
     "source": "UC Davis vole genome project (PRJEB108798; Abuelanin et al. 2026); tissue from Cold Spring Harbor Laboratory, 2017"},
    {"run": "DRR188168", "slug": "mus_musculus", "label": "house mouse C57BL/6, male liver (positive control)",
     "source": "circadian liver study, University of Tokyo (PRJDB7789); wild-type male at ZT0"},
    {"run": "SRR29282985", "slug": "myodes_glareolus", "label": "bank vole, adult female liver",
     "source": "tissue-expression study, Lund University (PRJNA1119683); collected at Revingehed, Sweden, January 2023"},
]
HEAD_BYTES = 120_000_000
K, STRIDE, MIN_HITS = 25, 4, 3
YARDSTICKS = {"Alb": "albumin", "Apoa1": "apolipoprotein A-I", "Actb": "actin, cytoplasmic 1|beta-actin|actin beta"}


def download(cfg, P) -> None:
    import requests
    s = requests.Session()
    s.headers["User-Agent"] = cfg["user_agent"]
    for r in RUNS:
        out = P["raw"] / f"ena_{r['run']}"
        out.mkdir(parents=True, exist_ok=True)
        dest = out / f"{r['run']}_1.head.fastq.gz"
        if dest.exists():
            continue
        rep = s.get("https://www.ebi.ac.uk/ena/portal/api/filereport", timeout=60,
                    params={"accession": r["run"], "result": "read_run",
                            "fields": "run_accession,fastq_ftp,sample_title,library_layout", "format": "tsv"}).text
        row = dict(zip(*[l.split("\t") for l in rep.strip().splitlines()[:2]]))
        url = "https://" + row["fastq_ftp"].split(";")[0]
        part = dest.with_name(dest.name + ".part")
        for attempt in range(8):                 # EBI drops long transfers now and then: resume where it stopped
            have = part.stat().st_size if part.exists() else 0
            if have >= HEAD_BYTES:
                break
            try:
                with s.get(url, stream=True, timeout=(30, 120),
                           headers={"Range": f"bytes={have}-{HEAD_BYTES - 1}"}) as q:
                    q.raise_for_status()
                    if have and q.status_code != 206:
                        have = 0                 # server ignored the range: start again
                    with open(part, "ab" if have else "wb") as fh:
                        for block in q.iter_content(1 << 20):
                            fh.write(block)
            except requests.RequestException as e:
                print(f"  {r['run']}: interrupted at {part.stat().st_size if part.exists() else 0:,} bytes "
                      f"({type(e).__name__}); resuming", flush=True)
                time.sleep(2 ** attempt)
        part.replace(dest)
        (out / "README.txt").write_text(
            f"{r['run']} ({row.get('sample_title', '')}, {row.get('library_layout', '')}): the first {HEAD_BYTES:,} bytes "
            f"of {url}, downloaded {now()} by VoleMUPs/scripts/step06_expression.py with an HTTP Range request. The file "
            "is a truncated gzip stream (the last read is incomplete and is skipped when reading).\n", encoding="utf-8")
        print(f"{r['run']}: {dest.stat().st_size:,} bytes", flush=True)


def reads(path: Path):
    """Sequences from a truncated gzipped FASTQ (stops cleanly at the cut)."""
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    buf = b""
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            try:
                buf += d.decompress(chunk)
            except zlib.error:
                break
            lines = buf.split(b"\n")
            n = (len(lines) - 1) // 4 * 4
            for i in range(0, n, 4):
                yield lines[i + 1].decode()
            buf = b"\n".join(lines[n:])


def transcripts(cfg, P, slug: str) -> dict:
    """gene key -> (label, transcript sequence): lipocalin genes + yardsticks, longest RefSeq transcript each."""
    sp = next(s for s in cfg["refseq"] if s["slug"] == slug)
    f = refseq_files(cfg, sp["accession"])
    lip = [r for r in read_tsv(P["derived"] / "annotated" / "lipocalins.tsv") if r["slug"] == slug]
    lip_ids = {r["gene_id"]: r for r in lip}
    tx_gene, gene_info = {}, {}
    with gzip.open(f["gff"], "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 9:
                continue
            if x[2] in ("gene", "pseudogene"):
                a = gff_attrs(x[8])
                gid = next((d[7:] for d in a.get("Dbxref", "").split(",") if d.startswith("GeneID:")), "")
                gene_info[a["ID"]] = (gid, a.get("Name", ""), unquote(a.get("description", "")))
            elif x[2] in ("mRNA", "transcript", "lnc_RNA") and "transcript_id" in x[8]:
                a = gff_attrs(x[8])
                tx_gene[a["transcript_id"]] = (a.get("Parent", ""), unquote(a.get("product", "")))
    want = {}
    for tid, (par, product) in tx_gene.items():
        gid, sym, desc = gene_info.get(par, ("", "", ""))
        label = None
        if gid in lip_ids:
            label = f"{sym}|{lip_ids[gid]['description']}"
        else:
            for ysym, pat in YARDSTICKS.items():
                if sym.upper() == ysym.upper() or (sym.startswith("LOC") and re.fullmatch(pat, product.split(",")[0], re.I)):
                    label = f"{ysym}|yardstick"
        if label:
            want[tid] = (gid, label)
    best = {}
    for name, _d, seq in fasta_records(f["rna"]):
        if name in want:
            gid, label = want[name]
            if gid not in best or len(seq) > len(best[gid][1]):
                best[gid] = (label, seq.upper())
    return best


def count(cfg, P) -> None:
    from vm_common import revcomp
    rows = []
    for r in RUNS:
        sp = next(s for s in cfg["refseq"] if s["slug"] == r["slug"])
        if refseq_files(cfg, sp["accession"])["rna"] is None:
            print(f"{r['slug']}: no RefSeq transcripts (run step01_download.py --rna {r['slug']})")
            continue
        fq = P["raw"] / f"ena_{r['run']}" / f"{r['run']}_1.head.fastq.gz"
        if not fq.exists():
            continue
        tx = transcripts(cfg, P, r["slug"])
        # k-mers shared by two genes are dropped, except when both are MUP genes: the near-identical mouse MUPs share
        # most of their sequence, so those k-mers count for the MUP family as a whole ('Mup (shared sequence)')
        fam = {gid: ("MUPFAM" if re.search(r"major urinary|alpha-2u", lab, re.I) else gid) for gid, (lab, _s) in tx.items()}
        idx, dup = {}, set()
        for gid, (_lab, seq) in tx.items():
            for s in (seq, revcomp(seq)):
                for i in range(len(s) - K + 1):
                    k = s[i:i + K]
                    old = idx.get(k)
                    if old is not None and old != gid:
                        if fam[gid] == "MUPFAM" and (old == "MUPFAM" or fam.get(old) == "MUPFAM"):
                            idx[k] = "MUPFAM"
                            continue
                        dup.add(k)
                    idx[k] = gid if old is None or old == gid else idx[k]
        for k in dup:
            idx.pop(k, None)
        n, c = 0, Counter()
        for rd in reads(fq):
            n += 1
            hits = Counter()
            for i in range(0, len(rd) - K + 1, STRIDE):
                g = idx.get(rd[i:i + K])
                if g:
                    hits[g] += 1
            if not hits or sum(hits.values()) < MIN_HITS:
                continue
            ids = set(hits)
            if len(ids) == 1:
                c[ids.pop()] += 1
            elif all(fam.get(x, x) == "MUPFAM" for x in ids):
                genes_hit = [x for x in ids if x != "MUPFAM"]
                c[genes_hit[0] if len(genes_hit) == 1 else "MUPFAM"] += 1
        if any(v == "MUPFAM" for v in fam.values()):
            tx = dict(tx, MUPFAM=("Mup (shared sequence)|reads matching sequence shared by several MUP genes", "N" * 1000))
        for gid, (lab, seq) in tx.items():
            sym, desc = lab.split("|", 1)
            rows.append({"run": r["run"], "slug": r["slug"], "sample": r["label"], "gene_id": gid, "symbol": sym,
                         "description": desc, "tx_length": len(seq), "reads": c[gid], "reads_total": n,
                         "rpm": round(c[gid] / n * 1e6, 1) if n else 0,
                         "rpkm": round(c[gid] / n * 1e6 / (len(seq) / 1000), 1) if n else 0})
        top = sorted((x for x in rows if x["run"] == r["run"]), key=lambda x: -x["rpm"])[:5]
        print(f"{r['run']} ({r['slug']}): {n:,} reads; top: " +
              "; ".join(f"{x['symbol']} {x['rpm']:.0f} RPM" for x in top), flush=True)
    write_tsv(P["derived"] / "expression" / "counts.tsv", rows)
    cross_check(cfg, P)


def cross_check(cfg, P) -> None:
    """A check that does not rest on each species' own gene model: every vole liver sample is matched (21-letter
    words, both strands) against the MUP-clade transcripts of all three voles together (prairie, meadow, bank), and
    the reads that match any of them are counted. If the meadow vole model were wrong or incomplete, its reads would
    still match the prairie or bank vole transcript, which are ~90-95% identical."""
    from vm_common import revcomp
    k21 = 21
    voles = ["microtus_ochrogaster", "microtus_pennsylvanicus", "myodes_glareolus"]
    genes = read_tsv(P["derived"] / "summary" / "genes.tsv") if (P["derived"] / "summary" / "genes.tsv").exists() else []
    mup = {(g["slug"], g["id"].split("|")[-1]) for g in genes if g["clade"] == "MUP"}
    idx = {}
    used = []
    for slug in voles:
        tx = transcripts(cfg, P, slug)
        for gid, (lab, seq) in tx.items():
            if (slug, gid) in mup:
                used.append(f"{slug}:{lab.split('|')[0]}")
                for s in (seq, revcomp(seq)):
                    for i in range(len(s) - k21 + 1):
                        idx[s[i:i + k21]] = slug
    out = []
    for r in RUNS:
        if r["slug"] not in voles:
            continue
        fq = P["raw"] / f"ena_{r['run']}" / f"{r['run']}_1.head.fastq.gz"
        if not fq.exists():
            continue
        n = hit = 0
        for rd in reads(fq):
            n += 1
            h = sum(1 for i in range(0, len(rd) - k21 + 1, STRIDE) if rd[i:i + k21] in idx)
            if h >= MIN_HITS:
                hit += 1
        out.append({"run": r["run"], "slug": r["slug"], "sample": r["label"], "reads_total": n, "reads_any_vole_mup": hit,
                    "rpm": round(hit / n * 1e6, 1) if n else 0, "transcripts": ", ".join(used)})
        print(f"cross-check {r['run']} ({r['slug']}): {hit} of {n:,} reads match a vole MUP transcript "
              f"({hit / n * 1e6:.0f} RPM)", flush=True)
    write_tsv(P["derived"] / "expression" / "vole_mup_crosscheck.tsv", out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--count", action="store_true")
    a = ap.parse_args()
    cfg = load_config()
    P = paths(cfg)
    if a.download:
        download(cfg, P)
    if a.count:
        count(cfg, P)


if __name__ == "__main__":
    main()
