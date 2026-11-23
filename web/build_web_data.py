"""Export the VoleMUPs results as compact JSON for the website page /projects/vole-mups/.

Reads the derived tables on the Synology (Steps 02-07) and writes, into the site's projects/vole-mups/data/:
- meta.json       species list with counts per clade and region summaries, headline numbers, build date
- tree.json       the focused MUP/OBP gene tree (Newick, rooted) + one record per tip (species, gene, clade, status, where)
- synteny.json    the three lipocalin neighbourhoods (MUP, OBP1, LCN2) per species: region, gaps, genes
- expression.json liver RNA-seq read counts for the lipocalins and three yardstick genes

Files are written to a temporary name and then renamed, so the running Jekyll preview never copies a half-written file.

Usage (any Python 3.10+ with PyYAML; e.g. the fnt env):
  python web/build_web_data.py [--data <data_root>] [--out ../calebvogt.github.io/projects/vole-mups/data]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
from vm_common import load_config  # noqa: E402
from step06_expression import RUNS  # noqa: E402  (labels and sources of the liver samples)

# the website checkout next to this repo (calebvogt.github.io)
SITE_DATA = HERE.parents[1] / "calebvogt.github.io" / "projects" / "vole-mups" / "data"
REGIONS = ["MUP", "OBP1", "LCN2"]
DISPLAY = {"microtus_ochrogaster": "prairie vole (2012)", "microtus_ochrogaster_ucd1": "prairie vole (2026)",
           "microtus_ochrogaster_ucd2": "prairie vole (2026, hap2)", "mus_spretus": "Algerian mouse (control)",
           "castor_canadensis": "beaver"}
NOTES = {"microtus_ochrogaster_ucd2": "The second haplotype of a male holds the autosomes only, so the X-linked OBP "
                                      "cluster is not in it.",
         "microtus_ochrogaster": "MicOch1.0 (2012), a short-read assembly with many gaps.",
         "microtus_ochrogaster_ucd1": "UC Davis 2026 long-read assembly, haplotype 1 (with chromosomes X and Y); "
                                      "no RefSeq annotation, searched with miniprot.",
         "mus_spretus": "Positive control: a mouse species whose genome should hold a large MUP cluster."}


def rd(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_json(path: Path, obj) -> None:
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, separators=(",", ":"), ensure_ascii=False)
    os.replace(tmp, path)


def num(x, default=None):
    try:
        f = float(x)
        return int(f) if f.is_integer() else round(f, 3)
    except (TypeError, ValueError):
        return default


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=None, help="VoleMUPs data folder (default: config.yml data_root)")
    ap.add_argument("--out", default=str(SITE_DATA))
    a = ap.parse_args()
    cfg = load_config()
    D = Path(a.data or cfg["data_root"]) / "derived"
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    species = rd(D / "summary" / "species.tsv")
    genes = rd(D / "summary" / "genes.tsv")
    cls = {r["tip"]: r for r in rd(D / "tree" / "classification.tsv")}
    sp_by = {s["slug"]: s for s in species}
    gene_by_id = {g["id"]: g for g in genes}

    # ---------------------------------------------------------------- tree
    treefile = D / "tree" / "focused" / "mup_obp.rooted.nwk"
    tree_kind = "focused"
    if not treefile.exists():
        treefile, tree_kind = D / "tree" / "full" / "lipocalins.rooted.nwk", "full"
    tips = {}
    newick = treefile.read_text().strip() if treefile.exists() else ""
    if newick:
        import re
        names = re.findall(r"[(,]([^(),:;]+):", newick)
        for nm in names:
            c = cls.get(nm, {})
            g = gene_by_id.get(nm, {})
            slug = c.get("slug", nm.split("|")[0])
            if slug == "uniprot":
                tips[nm] = {"sp": "uniprot", "g": c.get("symbol", "").replace("_", " "), "d": c.get("description", ""),
                            "c": c.get("clade", ""), "st": "reference", "src": "UniProt"}
                continue
            tips[nm] = {"sp": slug, "g": c.get("symbol", ""), "d": g.get("description") or c.get("description", ""),
                        "c": g.get("clade") or c.get("clade", ""), "st": g.get("status") or c.get("status", ""),
                        "loc": f"{c.get('chrom') or c.get('seqid', '')}:{num(c.get('start'), '')}-{num(c.get('end'), '')}",
                        "reg": g.get("region", ""), "src": g.get("source") or c.get("source", ""),
                        "aa": num(c.get("length")), "ref": c.get("best_ref", ""), "pid": num(c.get("best_ref_pid")),
                        "rna": num(g.get("rnaseq_cov")), "rns": num(g.get("rnaseq_intron_samples")),
                        "liv": num(g.get("liver_rpm"))}
    clades = rd(D / "tree" / "clades.tsv")
    write_json(out / "tree.json", {"kind": tree_kind, "newick": newick, "tips": tips,
                                   "clades": clades[0] if clades else {}})

    # ---------------------------------------------------------------- synteny (tables from step07_collect.py)
    flanks = rd(D / "regions" / "flanks.tsv")
    mrows = rd(D / "summary" / "map_rows.tsv")
    mgenes = rd(D / "summary" / "map_genes.tsv")
    mgaps = rd(D / "summary" / "map_gaps.tsv")
    syn = {r: [] for r in REGIONS}
    for r in mrows:
        if r["region"] not in syn:
            continue
        rid = r["region_id"]
        gl = [[int(g["start"]), int(g["end"]), g["strand"], g["symbol"], g["kind"], g["status"],
               g["tip"] if g["tip"] in tips else "", g["desc"]] for g in mgenes if g["region_id"] == rid]
        gp = [[int(x["start"]), int(x["end"])] for x in mgaps if x["region_id"] == rid]
        syn[r["region"]].append({"sp": r["slug"], "part": num(r["part"]), "seq": r["seqid"], "chrom": r["chrom"],
                                 "s": int(r["start"]), "e": int(r["end"]), "lo": int(r["lo"]), "hi": int(r["hi"]),
                                 "flip": r["flip"] == "True", "open": r["open"], "gaps": gp, "genes": gl})
    anchors = {r["region"]: r for r in rd(D / "regions" / "anchors.tsv")}
    write_json(out / "synteny.json", {"regions": syn, "anchors": {k: {"desc": v["description"], "n": num(v["n_genes"]),
                                                                       "chrom": v["chrom"]} for k, v in anchors.items()},
                                      "flanks": [{"r": f["region"], "side": f["side"], "rank": num(f["rank"]),
                                                  "sym": f["symbol"]} for f in flanks]})
    ann_ids = {(g["slug"], g["id"].split("|")[-1]): g for g in genes if g["source"] == "RefSeq annotation"}

    # ---------------------------------------------------------------- expression
    ex = rd(D / "expression" / "counts.tsv")
    samples = {}
    for e in ex:
        meta = next((r for r in RUNS if r["run"] == e["run"]), {})
        sm = samples.setdefault(e["run"], {"run": e["run"], "sp": e["slug"], "label": meta.get("label", e["sample"]),
                                           "source": meta.get("source", ""),
                                           "reads": num(e["reads_total"]), "genes": []})
        info = ann_ids.get((e["slug"], e["gene_id"]), {})
        clade = "yard" if e["description"] == "yardstick" else (e.get("family") or
                                                                 ("MUP" if e["gene_id"] == "MUPFAM" else
                                                                  (info.get("clade") or "lipo")))
        sym = e["symbol"]
        if clade == "MUP" and sym.startswith("LOC"):
            sym = "MUP " + sym                 # e.g. the single vole MUP gene, LOC101994101
        sm["genes"].append([sym, clade if clade in ("MUP", "OBP", "OBP2", "yard") else "lipo",
                            num(e["reads"]), num(e["rpm"]), (info.get("description") or e["description"])[:60]])
    for sm in samples.values():
        sm["genes"].sort(key=lambda x: -(x[3] or 0))
    write_json(out / "expression.json", {"samples": list(samples.values())})

    # ---------------------------------------------------------------- meta
    sp_out = []
    for s in species:
        sp_out.append({"sp": s["slug"], "name": s["name"], "common": DISPLAY.get(s["slug"], s["common"]),
                       "group": s["group"],
                       "src": s["source"], "acc": s["accession"], "note": NOTES.get(s["slug"], ""),
                       "mup": [num(s["MUP_intact"], 0), num(s["MUP_other"], 0)],
                       "obp": [num(s["OBP_intact"], 0), num(s["OBP_other"], 0)],
                       "obp2": [num(s.get("OBP2_intact"), 0), num(s.get("OBP2_other"), 0)],
                       "reg": {r: {"state": s[f"{r}_region"], "kb": num(s[f"{r}_region_kb"]),
                                   "gaps": num(s[f"{r}_gaps"]), "gap_bp": num(s[f"{r}_gap_bp"])} for r in REGIONS}})
    qc = rd(D / "summary" / "qc.tsv")
    write_json(out / "meta.json", {"built": dt.date.today().isoformat(), "species": sp_out,
                                   "qc": {q["check"]: q["value"] for q in qc}})
    print(f"wrote {out}: tree {len(tips)} tips ({tree_kind}), synteny "
          + ", ".join(f"{k} {len(v)} rows" for k, v in syn.items()) + f", expression {len(samples)} samples")

    # ---------------------------------------------------------------- 3D structures (step 08), next to data/
    st = D / "structures"
    if (st / "structures.tsv").exists():
        dest = out.parent / "structures"
        dest.mkdir(parents=True, exist_ok=True)
        rows = rd(st / "structures.tsv")
        for r in rows:
            (dest / f"{r['name']}.pdb").write_bytes((st / f"{r['name']}.pdb").read_bytes())
        write_json(dest / "structures.json", {"models": [{"name": r["name"], "label": r["label"], "acc": r["acc"],
                                                          "method": r["method"], "first": num(r["first_res"]),
                                                          "last": num(r["last_res"]), "plddt": num(r["mean_plddt"])}
                                                         for r in rows]})
        print(f"wrote {dest}: {len(rows)} structures")


if __name__ == "__main__":
    main()
