"""Shared helpers for the VoleMUPs scripts (not run on its own).

- load_config(): config.yml (+ config.local.yml), with data_root translated to /mnt/z/... when run inside WSL.
- paths(): raw/ and derived/ folders.
- species tables: refseq_species(), genome_species(), genome_fasta().
- small FASTA / TSV / GFF helpers that need only the Python standard library, so the same file works on Windows (the
  download step) and inside WSL (the search steps, micromamba env "genomics").
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJ = HERE.parent


# ----------------------------------------------------------------------------------------------------------------------
# Config and paths
# ----------------------------------------------------------------------------------------------------------------------
def _yaml_load(path: Path) -> dict:
    try:
        import yaml
    except ImportError:                       # tiny fallback is not worth it; both envs have PyYAML
        raise SystemExit("PyYAML is needed (pip install pyyaml)")
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def native_path(p: str) -> str:
    """Z:\\Work\\... -> /mnt/z/Work/... when running under Linux (WSL); unchanged on Windows."""
    if os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", p):
        return "/mnt/" + p[0].lower() + "/" + p[3:].replace("\\", "/")
    return p


def load_config() -> dict:
    cfg = _yaml_load(PROJ / "config.yml")
    local = PROJ / "config.local.yml"
    if local.exists():
        cfg.update(_yaml_load(local))
    cfg["data_root"] = native_path(cfg["data_root"])
    cfg["shared_roots"] = {k: native_path(v) for k, v in cfg.get("shared_roots", {}).items()}
    return cfg


def paths(cfg: dict) -> dict:
    root = Path(cfg["data_root"])
    d = {"root": root, "raw": root / "raw", "derived": root / "derived"}
    for k in ("raw", "derived"):
        d[k].mkdir(parents=True, exist_ok=True)
    return d


def refseq_dir(cfg: dict, acc: str) -> Path:
    return paths(cfg)["raw"] / f"refseq_{acc}"


def genome_fasta(cfg: dict, g: dict) -> Path:
    """Path of the genome FASTA (gz) of an entry of config 'genomes': from its shared_roots folder if that exists on
    this machine, else from <data_root>/raw/genome_<accession>/ (where step01_download.py puts it)."""
    bases = []
    if g.get("shared") and g["shared"] in cfg.get("shared_roots", {}):
        b = Path(cfg["shared_roots"][g["shared"]])
        bases.append(b / g["accession"] if g["shared"] == "monogamy_genomes" else b)
    bases.append(paths(cfg)["raw"] / f"genome_{g['accession']}")
    for base in bases:
        hits = sorted(base.glob(f"{g['accession']}_*_genomic.fna.gz"))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"no genome FASTA for {g['accession']} in {', '.join(map(str, bases))}; "
                            "run scripts/step01_download.py to fetch it")


def refseq_files(cfg: dict, acc: str) -> dict:
    d = refseq_dir(cfg, acc)
    out = {}
    for key, suffix in (("protein", "_protein.faa.gz"), ("gff", "_genomic.gff.gz"), ("report", "_assembly_report.txt"),
                        ("rna", "_rna.fna.gz"), ("genome", "_genomic.fna.gz")):
        hits = sorted(p for p in d.glob(f"{acc}_*{suffix}") if not p.name.endswith(("_rna_from_genomic.fna.gz",
                                                                                     "_cds_from_genomic.fna.gz")))
        if key == "genome":
            hits = [p for p in hits if p.name.endswith("_genomic.fna.gz") and "_rna" not in p.name
                    and "_cds" not in p.name]
        out[key] = hits[0] if hits else None
    return out


# ----------------------------------------------------------------------------------------------------------------------
# Files
# ----------------------------------------------------------------------------------------------------------------------
def opener(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, encoding="utf-8") if str(path).endswith(".gz") else open(path, mode, encoding="utf-8")


def fasta_records(path: Path):
    """Yield (name, description, sequence) from a FASTA file (plain or gz)."""
    name, desc, chunks = None, "", []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(chunks)
                head = line[1:].rstrip("\n")
                name, _, desc = head.partition(" ")
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, desc, "".join(chunks)


def write_fasta(path: Path, records, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
    tmp.replace(path)


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0].keys()) if rows else [])
    tmp = path.with_name(path.name + ".part")
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def read_tsv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1)
    tmp.replace(path)


def read_json(path: Path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


# ----------------------------------------------------------------------------------------------------------------------
# GFF
# ----------------------------------------------------------------------------------------------------------------------
def gff_attrs(s: str) -> dict:
    out = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def unquote(s: str) -> str:
    from urllib.parse import unquote as uq
    return uq(s)


# ----------------------------------------------------------------------------------------------------------------------
# Tools (inside WSL)
# ----------------------------------------------------------------------------------------------------------------------
def tool(name: str) -> str:
    found = shutil.which(name) or shutil.which(name, path=str(Path(sys.executable).parent))
    if not found:
        for env in ("genomics", "vmsearch"):
            cand = Path.home() / "micromamba" / "envs" / env / "bin" / name
            if cand.exists():
                return str(cand)
        raise SystemExit(f"{name} not found (run inside WSL with the micromamba envs genomics / vmsearch)")
    return found


def run(cmd: list, stdout=None, stdin=None, quiet: bool = True) -> str:
    cmd = [tool(cmd[0])] + [str(c) for c in cmd[1:]]
    r = subprocess.run(cmd, stdout=stdout or subprocess.PIPE, stdin=stdin, stderr=subprocess.PIPE)
    if r.returncode != 0:
        raise SystemExit(f"failed ({r.returncode}): {' '.join(cmd)}\n{r.stderr.decode(errors='replace')[-3000:]}")
    return r.stdout.decode(errors="replace") if stdout is None and r.stdout else ""


def lower_priority() -> None:
    try:
        os.nice(10)
    except (AttributeError, OSError):
        pass


COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]


def translate(cds: str) -> str:
    cds = cds.upper()
    return "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M")


# ----------------------------------------------------------------------------------------------------------------------
# miniprot (Li 2023, Bioinformatics 39:btad014) output
# ----------------------------------------------------------------------------------------------------------------------
def parse_miniprot_gff(path: Path) -> list[dict]:
    """One dict per miniprot alignment (adapted from 2027_GenomicCorpusCallosum/scripts/step0_genome_tools.py):
    query, contig, start/end (1-based, inclusive), strand, score (AS), frameshifts (fs), in-frame stops (st),
    identity, query range aligned, exons, and the translation of the genome along the model (--trans)."""
    models, cur = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("##PAF"):
                f = line.rstrip("\n").split("\t")[1:]
                tags = {t[:2]: t[5:] for t in f[12:] if len(t) > 5 and t[2] == ":" and t[4] == ":"}
                cur = {"query": f[0], "qlen": int(f[1]), "qstart": int(f[2]) + 1, "qend": int(f[3]),
                       "strand": f[4], "contig": f[5], "start": int(f[7]) + 1, "end": int(f[8]),
                       "score": int(tags.get("AS", 0)), "frameshifts": int(tags.get("fs", 0)),
                       "stops": int(tags.get("st", 0)), "exons": [], "translation": "", "identity": 0.0,
                       "rank": 0}
                models.append(cur)
            elif line.startswith("##STA") and cur is not None:
                cur["translation"] = line.rstrip("\n").split("\t")[1]
            elif line.startswith("#") or not line.strip():
                continue
            else:
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or cur is None:
                    continue
                a = gff_attrs(f[8])
                if f[2] == "mRNA":
                    cur["rank"] = int(a.get("Rank", 0))
                    cur["identity"] = float(a.get("Identity", 0))
                    cur["positive"] = float(a.get("Positive", 0))
                elif f[2] == "CDS":
                    cur["exons"].append((int(f[3]), int(f[4])))
    for m in models:
        m["exons"].sort()
        m["qcov"] = (m["qend"] - m["qstart"] + 1) / m["qlen"] if m["qlen"] else 0.0
    return models


MAX_LOCUS_SPAN = 30_000   # lipocalin genes span 2-15 kb of DNA from first to last coding exon; longer alignments
                          # join neighbouring genes into one chimera and are dropped before loci are formed


def cluster_loci(models: list[dict], min_overlap: float = 0.3) -> list[list[dict]]:
    """Group alignments of different queries that hit the same gene into loci. Greedy by score: the best alignment
    founds a locus; every other alignment joins the first locus whose founder it overlaps (same contig and strand) by
    at least min_overlap of the shorter of the two, else it founds a new locus. (Founders, not growing spans, so that
    an alignment that runs across two neighbouring genes does not chain them together.) Returns loci, each a list of
    alignments sorted by score (best first), in genome order."""
    loci = []
    by = {}
    models = [m for m in models if m["end"] - m["start"] + 1 <= MAX_LOCUS_SPAN]
    for m in sorted(models, key=lambda m: -m["score"]):
        key = (m["contig"], m["strand"])
        for loc in by.get(key, []):
            f = loc[0]
            ov = min(f["end"], m["end"]) - max(f["start"], m["start"]) + 1
            if ov >= min_overlap * min(f["end"] - f["start"] + 1, m["end"] - m["start"] + 1):
                loc.append(m)
                break
        else:
            loc = [m]
            by.setdefault(key, []).append(loc)
            loci.append(loc)
    loci.sort(key=lambda l: (l[0]["contig"], l[0]["start"]))
    return loci


def n_runs(seq: str, min_len: int = 10) -> list[tuple[int, int]]:
    """Runs of N (assembly gaps) at least min_len long, as 1-based inclusive (start, end) on seq."""
    return [(m.start() + 1, m.end()) for m in re.finditer(r"[Nn]{%d,}" % min_len, seq)]
