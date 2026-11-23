"""Step 01: download what the search needs (run on Windows or the Mac; no Linux tools needed).

1. For every species in config 'refseq': NCBI's RefSeq protein FASTA, gene annotation (GFF3) and assembly report,
   from the NCBI genome FTP site, each checked against NCBI's md5checksums.txt, into raw/refseq_<accession>/ with a
   README.txt. --rna also fetches the RefSeq transcript FASTA (rna.fna.gz), used by the expression check.
2. The reference lipocalins listed in config 'uniprot' (aphrodisin, the bank vole OBPs, arvicolin, the horse, cat and
   dog MUP allergens) from UniProt, into raw/uniprot_reference/.
3. The Pfam lipocalin profile HMM (PF00061) from InterPro, into raw/pfam_PF00061/.

Files that are already complete (md5 matches) are not downloaded again, so the script can simply be rerun.
Usage: python scripts/step01_download.py [--rna] [--only slug1,slug2]
"""
from __future__ import annotations

import argparse
import re
import time
from pathlib import Path

import requests

from vm_common import genome_fasta, load_config, md5, now, paths, refseq_dir

FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"


def session(cfg) -> requests.Session:
    s = requests.Session()
    s.headers["User-Agent"] = cfg["user_agent"]
    return s


def get(s: requests.Session, url: str, tries: int = 6, **kw) -> requests.Response:
    for attempt in range(tries):
        try:
            r = s.get(url, timeout=(30, 180), **kw)
            if r.status_code == 200:
                return r
            if r.status_code not in (429, 500, 502, 503, 504):
                r.raise_for_status()
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    raise RuntimeError(f"failed: {url}")


def ftp_folder(s, acc: str) -> str:
    prefix, digits = acc.split("_")[0], acc.split("_")[1].split(".")[0]
    parent = f"{FTP}/{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/"
    listing = get(s, parent).text
    folders = sorted(set(re.findall(rf'href="({re.escape(acc)}_[^"/]+)/?"', listing)))
    if not folders:
        raise RuntimeError(f"{acc}: no folder in {parent}")
    return parent + folders[0] + "/"


def fetch(s, url: str, dest: Path, want_md5: str | None) -> None:
    if dest.exists() and (want_md5 is None or md5(dest) == want_md5):
        return
    part = dest.with_name(dest.name + ".part")
    for attempt in range(6):
        try:
            with s.get(url, stream=True, timeout=(30, 120)) as r:
                r.raise_for_status()
                with open(part, "wb") as fh:
                    for block in r.iter_content(1 << 22):
                        fh.write(block)
        except requests.RequestException as e:
            print(f"  {dest.name}: interrupted ({type(e).__name__}); retrying", flush=True)
            time.sleep(2 ** attempt)
            continue
        if want_md5 is None or md5(part) == want_md5:
            part.replace(dest)
            return
        time.sleep(2 ** attempt)
    raise RuntimeError(f"{dest.name}: no complete download matching NCBI's md5")


def refseq(cfg, s, sp: dict, with_rna: bool) -> None:
    acc = sp["accession"]
    out = refseq_dir(cfg, acc)
    out.mkdir(parents=True, exist_ok=True)
    base = ftp_folder(s, acc)
    stem = base.rstrip("/").rsplit("/", 1)[1]
    md5s = {}
    for line in get(s, base + "md5checksums.txt").text.splitlines():
        parts = line.split()
        if len(parts) == 2:
            md5s[parts[1].lstrip("./")] = parts[0]
    suffixes = ["_protein.faa.gz", "_genomic.gff.gz", "_assembly_report.txt"] + (["_rna.fna.gz"] if with_rna else [])
    got = []
    for suf in suffixes:
        name = stem + suf
        if name not in md5s:
            print(f"  {acc}: {name} not on the FTP site", flush=True)
            continue
        fetch(s, base + name, out / name, md5s.get(name))
        got.append(name)
    readme = out / "README.txt"
    readme.write_text(
        f"{acc} ({stem.split('_', 2)[-1]}), {sp['name']} ({sp['common']}).\n"
        f"Downloaded {now()} by VoleMUPs/scripts/step01_download.py from {base}\n"
        f"checked against NCBI's md5checksums.txt. Files as received: {', '.join(got)}.\n", encoding="utf-8")
    print(f"{sp['slug']}: {', '.join(got)}", flush=True)


def genome(cfg, s, g: dict) -> None:
    """Whole-genome FASTA of an assembly without a RefSeq annotation (config 'genomes'), unless it is already on this
    machine (shared_roots or raw/genome_<accession>/)."""
    try:
        print(f"{g['slug']}: found {genome_fasta(cfg, g)}", flush=True)
        return
    except FileNotFoundError:
        pass
    acc = g["accession"]
    out = paths(cfg)["raw"] / f"genome_{acc}"
    out.mkdir(parents=True, exist_ok=True)
    base = ftp_folder(s, acc)
    stem = base.rstrip("/").rsplit("/", 1)[-1]
    md5s = {}
    for line in get(s, base + "md5checksums.txt").text.splitlines():
        parts = line.split()
        if len(parts) == 2:
            md5s[parts[1].lstrip("./")] = parts[0]
    name = stem + "_genomic.fna.gz"
    fetch(s, base + name, out / name, md5s.get(name))
    (out / "README.txt").write_text(
        f"{acc}, {g['name']} ({g['common']}).\nDownloaded {now()} by scripts/step01_download.py from {base}{name},\n"
        f"checked against NCBI's md5checksums.txt.\n", encoding="utf-8")
    print(f"{g['slug']}: {name}", flush=True)


def uniprot(cfg, s) -> None:
    out = paths(cfg)["raw"] / "uniprot_reference"
    out.mkdir(parents=True, exist_ok=True)
    recs = []
    for u in cfg["uniprot"]:
        txt = get(s, f"https://rest.uniprot.org/uniprotkb/{u['acc']}.fasta").text
        lines = txt.strip().splitlines()
        recs.append((f"UP|{u['label']}|{u['acc']}", "".join(lines[1:]), lines[0]))
    with open(out / "reference_lipocalins.fasta", "w", encoding="utf-8", newline="\n") as fh:
        for name, seq, _ in recs:
            fh.write(f">{name}\n{seq}\n")
    with open(out / "uniprot_headers.txt", "w", encoding="utf-8", newline="\n") as fh:
        for name, _, head in recs:
            fh.write(f"{name}\t{head}\n")
    (out / "README.txt").write_text(f"UniProt entries {', '.join(u['acc'] for u in cfg['uniprot'])}, downloaded {now()} "
                                    "from https://rest.uniprot.org/uniprotkb/<acc>.fasta by step01_download.py.\n",
                                    encoding="utf-8")
    print(f"uniprot: {len(recs)} proteins", flush=True)


def pfam(cfg, s) -> None:
    acc = cfg["pfam_hmm"]
    out = paths(cfg)["raw"] / f"pfam_{acc}"
    out.mkdir(parents=True, exist_ok=True)
    dest = out / f"{acc}.hmm.gz"
    if not dest.exists():
        r = get(s, f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{acc}?annotation=hmm")
        part = dest.with_name(dest.name + ".part")
        part.write_bytes(r.content)
        part.replace(dest)
    (out / "README.txt").write_text(f"Pfam {acc} profile HMM (gzip), downloaded {now()} from the InterPro API "
                                    f"(https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{acc}?annotation=hmm).\n",
                                    encoding="utf-8")
    print(f"pfam: {dest}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna", default="", help="comma-separated slugs whose RefSeq transcripts to fetch too")
    ap.add_argument("--only", default="")
    ap.add_argument("--skip-refseq", action="store_true")
    ap.add_argument("--skip-genomes", action="store_true",
                    help="do not fetch the whole genomes of the assemblies without an annotation (about 6 GB)")
    a = ap.parse_args()
    cfg = load_config()
    s = session(cfg)
    only = set(filter(None, a.only.split(",")))
    rna = set(filter(None, a.rna.split(",")))
    uniprot(cfg, s)
    pfam(cfg, s)
    if not a.skip_refseq:
        for sp in cfg["refseq"]:
            if only and sp["slug"] not in only:
                continue
            refseq(cfg, s, sp, sp["slug"] in rna)
    if not a.skip_genomes:
        for g in cfg.get("genomes", []):
            if only and g["slug"] not in only:
                continue
            genome(cfg, s, g)


if __name__ == "__main__":
    main()
