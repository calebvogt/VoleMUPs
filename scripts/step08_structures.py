"""Step 08: 3D structures of the lipocalins shown on the project page.

  --alphafold  download AlphaFold database models (config 'structures: alphafold')
  --esmfold    predict structures with ESMFold (config 'structures: esmfold'), for proteins the AlphaFold database
               does not have (e.g. the prairie vole MUP). Runs on the CPU with Hugging Face transformers
               (facebook/esmfold_v1, ~8 GB of weights, cached in <project>/models/huggingface next to data/).

Signal peptides (UniProt 'Signal' feature) are removed. B-factors hold the per-residue confidence (pLDDT, 0-100).
Writes <data_root>/derived/structures/<name>.pdb and structures.tsv.

usage: python scripts/step08_structures.py --alphafold --esmfold [--threads 4]
ESMFold needs torch and transformers (not in environment.yml): pip install torch transformers
"""
from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import requests

from vm_common import load_config, now, paths, write_tsv

AFDB = "https://alphafold.ebi.ac.uk"
UNIPROT = "https://rest.uniprot.org/uniprotkb"


def session(cfg) -> requests.Session:
    s = requests.Session()
    s.headers["User-Agent"] = cfg.get("user_agent", "genome-inventory-script/1.0")
    return s


def uniprot_entry(s, acc: str) -> tuple[str, int]:
    """Sequence and signal-peptide end (0 when none) of a UniProt entry."""
    seq = "".join(s.get(f"{UNIPROT}/{acc}.fasta", timeout=60).text.strip().splitlines()[1:])
    sig = s.get(f"{UNIPROT}/{acc}", params={"fields": "ft_signal", "format": "tsv"}, timeout=60).text
    end = 0
    for tok in sig.replace(";", " ").split():
        if ".." in tok:
            end = int(tok.split("..")[1])
            break
    return seq, end


def trim(pdb_text: str, sig_end: int, scale: float = 1.0) -> tuple[list[str], float, int]:
    """ATOM lines after the signal peptide (B-factor scaled to 0-100), mean pLDDT over CA atoms, residues kept."""
    out, ca = [], []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or int(line[22:26]) <= sig_end:
            continue
        b = float(line[60:66]) * scale
        line = f"{line[:60]}{b:6.2f}{line[66:]}"
        out.append(line)
        if line[12:16].strip() == "CA":
            ca.append(b)
    return out, sum(ca) / len(ca), len(ca)


def write_model(path: Path, remark: str, atoms: list[str]) -> None:
    path.write_text("\n".join([f"REMARK {remark}"] + atoms + ["END", ""]), encoding="utf-8", newline="\n")


def alphafold(cfg, s, out: Path) -> list[dict]:
    rows = []
    for e in cfg["structures"]["alphafold"]:
        seq, sig = uniprot_entry(s, e["acc"])
        meta = s.get(f"{AFDB}/api/prediction/{e['acc']}", timeout=60).json()[0]
        txt = s.get(meta["pdbUrl"], timeout=120).text
        atoms, mean, n = trim(txt, sig)
        ver = meta.get("latestVersion", "")
        write_model(out / f"{e['name']}.pdb", f"AlphaFold DB {meta['entryId']} v{ver} (CC BY 4.0), UniProt {e['acc']}, "
                    f"signal peptide 1-{sig} removed; B-factor = pLDDT; fetched {now()}", atoms)
        rows.append({"name": e["name"], "label": e["label"], "acc": e["acc"], "method": f"AlphaFold DB {meta['entryId']}",
                     "first_res": sig + 1, "last_res": len(seq), "residues": n, "mean_plddt": round(mean, 1)})
        print(f"{e['name']}: {meta['entryId']}, {n} residues, mean pLDDT {mean:.1f}", flush=True)
    return rows


def esmfold(cfg, s, out: Path, threads: int) -> list[dict]:
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding
    torch.set_num_threads(threads)
    cache = Path(cfg["data_root"]).parent / "models" / "huggingface"
    cache.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1", cache_dir=cache)
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", cache_dir=cache, low_cpu_mem_usage=True)
    model.eval()
    model.trunk.set_chunk_size(64)
    print(f"ESMFold loaded in {time.time() - t0:.0f} s", flush=True)
    rows = []
    for e in cfg["structures"]["esmfold"]:
        seq, sig = uniprot_entry(s, e["acc"])
        t1 = time.time()
        with torch.no_grad():
            res = model(**tok([seq], return_tensors="pt", add_special_tokens=False))
        txt = model.output_to_pdb(res)[0]
        bmax = max(float(l[60:66]) for l in txt.splitlines() if l.startswith("ATOM"))
        atoms, mean, n = trim(txt, sig, scale=100.0 if bmax <= 1.0 else 1.0)
        write_model(out / f"{e['name']}.pdb", f"ESMFold (facebook/esmfold_v1) prediction, UniProt {e['acc']}, "
                    f"signal peptide 1-{sig} removed; B-factor = pLDDT; predicted {now()}", atoms)
        rows.append({"name": e["name"], "label": e["label"], "acc": e["acc"], "method": "ESMFold (esmfold_v1)",
                     "first_res": sig + 1, "last_res": len(seq), "residues": n, "mean_plddt": round(mean, 1)})
        print(f"{e['name']}: ESMFold, {len(seq)} aa in {time.time() - t1:.0f} s, {n} residues kept, "
              f"mean pLDDT {mean:.1f}", flush=True)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--alphafold", action="store_true")
    ap.add_argument("--esmfold", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    cfg = load_config()
    s = session(cfg)
    out = paths(cfg)["derived"] / "structures"
    out.mkdir(parents=True, exist_ok=True)
    table = out / "structures.tsv"
    old = {}
    if table.exists():
        lines = table.read_text(encoding="utf-8").splitlines()
        head = lines[0].split("\t")
        old = {r["name"]: r for r in (dict(zip(head, l.split("\t"))) for l in lines[1:])}
    rows = []
    if a.alphafold:
        rows += alphafold(cfg, s, out)
    if a.esmfold:
        rows += esmfold(cfg, s, out, a.threads)
    for r in rows:
        old[r["name"]] = r
    write_tsv(table, list(old.values()))
    print(f"wrote {table}", flush=True)


if __name__ == "__main__":
    os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
    main()
