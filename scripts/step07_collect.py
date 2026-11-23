"""Step 07: one table per question, from Steps 02-06, for the report (reports/01_mup_search.Rmd) and the web page.

Outputs (derived/summary/)
- genes.tsv     every lipocalin of the MUP and OBP clades (and the decoy-family genes inside the three regions): species,
                source (RefSeq annotation / DNA search), clade from the gene tree, status (annotated protein,
                pseudogene, intact / disrupted / partial DNA hit), location, which region it sits in, RNA-seq notes
- species.tsv   per species: genes in the MUP clade and the OBP clade (intact / not), and for each of the three regions
                whether it was found, its length, assembly gaps, and whether it runs off a scaffold end
- regions.tsv   per species x region: the region and every gene in it (for the synteny map)
Run anywhere (pure Python): python scripts/step07_collect.py
"""
from __future__ import annotations

import gzip
from collections import Counter, defaultdict
from pathlib import Path

import re

from vm_common import load_config, paths, read_tsv, write_tsv

re_mup = re.compile(r"^Mup\d+$")

ORDER = ["mus_musculus", "mus_spretus", "rattus_norvegicus", "meriones_unguiculatus", "peromyscus_maniculatus",
         "mesocricetus_auratus", "cricetulus_griseus", "phodopus_roborovskii", "phodopus_sungorus",
         "myodes_glareolus", "arvicola_amphibius", "chionomys_nivalis", "alexandromys_fortis", "microtus_oregoni",
         "microtus_agrestis", "microtus_arvalis", "microtus_pennsylvanicus", "microtus_montanus",
         "microtus_ochrogaster", "microtus_ochrogaster_ucd1", "microtus_ochrogaster_ucd2", "castor_canadensis",
         "equus_caballus", "homo_sapiens"]


def cds_spans(cfg, regions: list[dict]) -> list[dict]:
    """First-to-last coding base of every gene in the regions (RefSeq GFF CDS lines). Gene records can span far more
    DNA than the protein needs (long untranslated or read-through transcripts), so the map draws lipocalins by CDS."""
    import gzip as _gz
    from vm_common import gff_attrs, refseq_files
    want = defaultdict(list)
    for r in regions:
        if r.get("found") and r.get("region_id") and "__" in r["region_id"]:
            want[r["slug"]].append((r["seqid"], r["start"] - 60_000, r["end"] + 60_000))
    out = []
    for sp in cfg["refseq"]:
        if sp["slug"] not in want:
            continue
        seqids = {w[0] for w in want[sp["slug"]]}
        gff = refseq_files(cfg, sp["accession"])["gff"]
        mrna_gene, gene_id, span = {}, {}, {}
        with _gz.open(gff, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split(chr(9), 8)
                if len(f) < 9 or f[0] not in seqids:
                    continue
                a_, b_ = int(f[3]), int(f[4])
                if not any(f[0] == q and a_ <= hi and b_ >= lo for q, lo, hi in want[sp["slug"]]):
                    continue
                at = gff_attrs(f[8].rstrip())
                if f[2] in ("gene", "pseudogene"):
                    gene_id[at["ID"]] = next((x[7:] for x in at.get("Dbxref", "").split(",") if x.startswith("GeneID:")), "")
                elif f[2] in ("mRNA", "transcript"):
                    mrna_gene[at["ID"]] = at.get("Parent", "")
                elif f[2] == "CDS":
                    par = at.get("Parent", "")
                    g = mrna_gene.get(par, par)
                    lo_, hi_ = span.get(g, (a_, b_))
                    span[g] = (min(lo_, a_), max(hi_, b_))
        for g, (a_, b_) in span.items():
            if gene_id.get(g):
                out.append({"slug": sp["slug"], "gene_id": gene_id[g], "cds_start": a_, "cds_end": b_})
    return out


LONG = 3_000_000         # regions longer than this keep only lipocalins, flanks and genes within NEAR of a lipocalin
NEAR = 50_000


def kind_of(clade: str, flank: bool, lipocalin: bool) -> str:
    if clade in ("MUP", "OBP", "OBP2"):
        return clade
    if flank:
        return "flank"
    if lipocalin:
        return "lipo"
    return "gene"


def map_tables(D: Path, regions: list[dict], genes: list[dict], sp_meta: dict, cds: dict, order: list[str]):
    """The rows of the synteny map: for each species x region, the window drawn (the region between the innermost
    flanking genes plus a small pad), whether to draw it reversed (so that the mouse-left flanking genes are on the
    left), the genes in it (lipocalins drawn by their coding span), unannotated DNA hits and assembly gaps."""
    flanks = read_tsv(D / "regions" / "flanks.tsv")
    left_syms = {f["region"]: set() for f in flanks}
    for f in flanks:
        if f["side"] == "L":
            left_syms[f["region"]].add(f["symbol"].upper())
    flank_syms = {f["symbol"].upper() for f in flanks}
    rgenes = read_tsv(D / "regions" / "region_genes.tsv")
    gaps = read_tsv(D / "regions" / "gaps.tsv") if (D / "regions" / "gaps.tsv").exists() else []
    rloci = read_tsv(D / "regions" / "region_loci.tsv") if (D / "regions" / "region_loci.tsv").exists() else []
    by_id = {g["id"]: g for g in genes}
    ann = {(g["slug"], g["id"].split("|")[-1]): g for g in genes if g["source"] == "RefSeq annotation"}
    rows, mg, mgap = [], [], []
    found = [r for r in regions if r["found"]]
    found.sort(key=lambda r: (order.index(r["slug"]) if r["slug"] in order else 99, r["part"], r["region"]))
    for rg in found:
        rid, slug, reg = rg["region_id"], rg["slug"], rg["region"]
        s, e = rg["start"], rg["end"]
        gl, gp = [], []
        if sp_meta.get(slug, {}).get("source") == "refseq":
            for g in rgenes:
                if g["region_id"] != rid:
                    continue
                gs, ge = int(g["start"]), int(g["end"])
                info = ann.get((slug, g["gene_id"]), {})
                clade = info.get("clade", "")
                lipo = g["lipocalin"] == "True"
                if (lipo or clade) and (slug, g["gene_id"]) in cds:
                    gs, ge = cds[(slug, g["gene_id"])]
                is_flank = g["symbol"].upper() in flank_syms
                gl.append({"start": gs, "end": ge, "strand": g["strand"], "symbol": g["symbol"],
                           "kind": kind_of(clade, is_flank, lipo), "status": info.get("status", ""),
                           "tip": f"{slug}|{g['gene_id']}", "desc": (g["description"] or "")[:60]})
            for l in rloci:
                if l["region_id"] != rid or l["status"] == "fragment":
                    continue
                a_, b_ = int(l["start"]), int(l["end"])
                if any(x["start"] <= b_ and a_ <= x["end"] and x["kind"] in ("MUP", "OBP", "OBP2", "lipo") for x in gl):
                    continue
                g2 = by_id.get(f"{rid}|locus{l['locus']}", {})
                q = l["best_query"].split("|")[-1]
                gl.append({"start": a_, "end": b_, "strand": l["strand"], "symbol": "",
                           "kind": kind_of(g2.get("clade", ""), False, True), "status": f"DNA hit, {l['status']}",
                           "tip": "", "desc": f"unannotated; best match {q} ({round(float(l['identity']) * 100)}% identical)"})
            gp = [(int(x["start"]), int(x["end"])) for x in gaps if x["region_id"] == rid]
        else:
            gdir = D / "genomes" / slug
            for f in read_tsv(gdir / "flank_hits.tsv"):
                if f["region"] == reg and f["contig"] == rg["seqid"]:
                    gl.append({"start": int(f["start"]), "end": int(f["end"]), "strand": f["strand"],
                               "symbol": f["symbol"], "kind": "flank", "status": "placed by miniprot", "tip": "",
                               "desc": f"flanking gene ({round(float(f['identity']) * 100)}% identical to the query)"})
            for l in read_tsv(gdir / "loci.tsv"):
                if l["contig"] != rg["seqid"] or l["status"] == "fragment":
                    continue
                tip = f"{slug}|locus{l['locus']}"
                g2 = by_id.get(tip, {})
                q = l["best_query"].split("|")[-1]
                gl.append({"start": int(l["start"]), "end": int(l["end"]), "strand": l["strand"],
                           "symbol": f"locus{l['locus']}", "kind": kind_of(g2.get("clade", ""), False, True),
                           "status": l["status"], "tip": tip,
                           "desc": f"best match {q} ({round(float(l['identity']) * 100)}% identical)"})
            gp = [(int(x["start"]), int(x["end"])) for x in read_tsv(gdir / "gaps.tsv")
                  if x["region"] == reg and x["contig"] == rg["seqid"]]
        fl = [x for x in gl if x["kind"] == "flank"]
        pad = max(8_000, int(0.06 * (e - s)))
        lo, hi = max(1, s - pad), e + pad
        gl = [x for x in gl if x["end"] >= lo and x["start"] <= hi]
        if hi - lo > LONG:
            lip = [x for x in gl if x["kind"] in ("MUP", "OBP", "OBP2", "lipo")]
            gl = [x for x in gl if x["kind"] != "gene" or any(abs(x["start"] - y["start"]) < NEAR for y in lip)]
        L = [x for x in fl if x["symbol"].upper() in left_syms.get(reg, set())]
        R = [x for x in fl if x["symbol"].upper() not in left_syms.get(reg, set())]
        flip = False
        if L and R:
            flip = min(x["start"] for x in L) > min(x["start"] for x in R)
        elif L:
            flip = min(x["start"] for x in L) > (lo + hi) / 2
        elif R:
            flip = min(x["start"] for x in R) < (lo + hi) / 2
        rows.append({"region_id": rid, "region": reg, "slug": slug, "part": rg["part"], "seqid": rg["seqid"],
                     "chrom": rg.get("chrom", ""), "start": s, "end": e, "lo": lo, "hi": hi, "flip": flip,
                     "open": rg.get("open", "")})
        for x in sorted(gl, key=lambda x: x["start"]):
            mg.append(dict(x, region_id=rid))
        for a_, b_ in gp:
            if b_ >= lo and a_ <= hi:
                mgap.append({"region_id": rid, "start": a_, "end": b_})
    return rows, mg, mgap


def main() -> None:
    cfg = load_config()
    P = paths(cfg)
    D = P["derived"]
    out = D / "summary"
    out.mkdir(exist_ok=True)
    sp_meta = {s["slug"]: dict(s, source="refseq") for s in cfg["refseq"]}
    sp_meta.update({g["slug"]: dict(g, source="genome") for g in cfg["genomes"]})
    cls = {r["tip"]: r for r in read_tsv(D / "tree" / "classification.tsv")}
    # clade of each query protein (for DNA hits that are not in a tree: they take their best query's clade)
    q_clade = {r["symbol"]: r["clade"] for r in cls.values()
               if r["slug"] in ("mus_musculus", "rattus_norvegicus", "uniprot") and r["clade"]}
    lip = {r["key"]: r for r in read_tsv(D / "annotated" / "lipocalins.tsv")}

    # --- regions (annotated genomes from step03, others from step04)
    regions = []
    for r in read_tsv(D / "regions" / "regions.tsv"):
        if r["found"] != "True":
            regions.append({"slug": r["slug"], "region": r["region"], "part": 0, "found": False})
            continue
        regions.append({"slug": r["slug"], "region": r["region"], "part": int(r["part"]), "found": True,
                        "seqid": r["seqid"], "chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"]),
                        "length": int(r["length"]), "left_flank": r["left_flank"], "right_flank": r["right_flank"],
                        "open": r["open"], "region_id": r["region_id"]})
    gaps = read_tsv(D / "regions" / "gaps.tsv") if (D / "regions" / "gaps.tsv").exists() else []
    for rg in regions:
        if rg["found"]:
            g = [x for x in gaps if x["region_id"] == rg["region_id"] and x["inside"] == "True"]
            rg["n_gaps"], rg["gap_bp"] = len(g), sum(int(x["length"]) for x in g)
    for g in cfg["genomes"]:
        f = D / "genomes" / g["slug"] / "regions.tsv"
        if not f.exists():
            continue
        for r in read_tsv(f):
            if r["found"] != "True" or (g.get("no_x") and r["region"] == "OBP1"):
                # OBP1 is X-linked; a haplotype assembly of a male without chr X cannot hold it
                regions.append({"slug": g["slug"], "region": r["region"], "part": 0, "found": False,
                                "note": "no X chromosome in this haplotype" if g.get("no_x") else ""})
                continue
            regions.append({"slug": g["slug"], "region": r["region"], "part": int(r["part"]), "found": True,
                            "seqid": r["contig"], "chrom": r["contig"], "start": int(r["start"]), "end": int(r["end"]),
                            "length": int(r["length"]), "left_flank": r["left_flank"], "right_flank": r["right_flank"],
                            "open": r["open"], "n_gaps": int(r["n_gaps"] or 0), "gap_bp": int(r["gap_bp"] or 0),
                            "region_id": f"{g['slug']}__{r['region']}__{r['part']}"})

    def region_of(slug, seqid, start) -> str:
        hits = [rg for rg in regions if rg["slug"] == slug and rg["found"] and rg["seqid"] == seqid
                and rg["start"] - 50_000 <= int(start) <= rg["end"] + 50_000]
        return ",".join(sorted({h["region"] for h in hits}))

    # --- genes: annotated (RefSeq)
    genes = []
    for key, r in lip.items():
        tip = f"{r['slug']}|{r['gene_id']}"
        c = cls.get(tip, {})
        pseudo = r["biotype"] != "protein_coding"
        genes.append({"slug": r["slug"], "id": tip, "source": "RefSeq annotation", "symbol": r["symbol"],
                      "description": r["description"], "seqid": r["seqid"], "chrom": r["chrom"],
                      "start": int(r["start"]), "end": int(r["end"]), "strand": r["strand"],
                      "clade": c.get("clade", "") or ("MUP" if "urinary" in r["description"].lower() and pseudo else
                                                       ("OBP" if pseudo and "odorant" in r["description"].lower() else "")),
                      "clade_from": "tree" if c.get("clade") else ("name (no protein)" if pseudo else ""),
                      "status": "pseudogene (annotated)" if pseudo else "protein-coding (annotated)",
                      "intact": not pseudo, "length_aa": int(r["length"] or 0), "best_ref": c.get("best_ref", ""),
                      "best_ref_pid": c.get("best_ref_pid", ""), "fabp_like": c.get("fabp_like", ""),
                      "rnaseq_cov": r["rnaseq_cov"], "rnaseq_intron_samples": r["rnaseq_intron_samples"],
                      "region": region_of(r["slug"], r["seqid"], r["start"])})
    # --- DNA hits in the regions of annotated genomes that no annotated gene covers (unannotated / pseudogenes)
    ann_by = defaultdict(list)
    for g in genes:
        ann_by[(g["slug"], g["seqid"])].append((g["start"], g["end"]))
    if (D / "regions" / "region_loci.tsv").exists():
        for l in read_tsv(D / "regions" / "region_loci.tsv"):
            s, e = int(l["start"]), int(l["end"])
            if any(a <= e and s <= b for a, b in ann_by[(l["slug"], l["seqid"])]):
                continue
            if l["inside"] != "True":
                continue
            genes.append({"slug": l["slug"], "id": f"{l['region_id']}|locus{l['locus']}", "source": "DNA search (region)",
                          "symbol": "", "description": f"unannotated; best match {l['best_query']}",
                          "seqid": l["seqid"], "chrom": "", "start": s, "end": e, "strand": l["strand"],
                          "clade": q_clade.get(l["best_query"].split("|")[-1], "") or "other",
                          "clade_from": "best-matching query", "status": l["status"],
                          "intact": l["status"] == "intact", "length_aa": len(l["translation"]),
                          "best_ref": l["best_query"].split("|")[-1], "best_ref_pid": round(float(l["identity"]) * 100, 1),
                          "frameshifts": l["frameshifts"], "stops": l["stops"], "qcov": l["qcov"],
                          "region": l["region"]})
    # --- unannotated genomes (step04)
    for g in cfg["genomes"]:
        f = D / "genomes" / g["slug"] / "loci.tsv"
        if not f.exists():
            continue
        for l in read_tsv(f):
            tip = f"{g['slug']}|locus{l['locus']}"
            c = cls.get(tip, {})
            clade = c.get("clade") or q_clade.get(l["best_query"].split("|")[-1], "") or "other"
            genes.append({"slug": g["slug"], "id": tip, "source": "DNA search (whole genome)", "symbol": "",
                          "description": f"best match {l['best_query']}", "seqid": l["contig"], "chrom": l["contig"],
                          "start": int(l["start"]), "end": int(l["end"]), "strand": l["strand"], "clade": clade,
                          "clade_from": "tree" if c.get("clade") else "best-matching query", "status": l["status"],
                          "intact": l["status"] == "intact", "length_aa": len(l["translation"]),
                          "best_ref": c.get("best_ref", "") or l["best_query"].split("|")[-1],
                          "best_ref_pid": c.get("best_ref_pid", "") or round(float(l["identity"]) * 100, 1),
                          "frameshifts": l["frameshifts"], "stops": l["stops"], "qcov": l["qcov"],
                          "region": region_of(g["slug"], l["contig"], l["start"])})
    # expression
    expr = {}
    if (D / "expression" / "counts.tsv").exists():
        for e in read_tsv(D / "expression" / "counts.tsv"):
            expr[(e["slug"], e["gene_id"])] = e
    for g in genes:
        e = expr.get((g["slug"], g["id"].split("|")[-1]))
        g["liver_rpm"] = e["rpm"] if e else ""
    fields = ["slug", "id", "source", "symbol", "description", "seqid", "chrom", "start", "end", "strand", "clade",
              "clade_from", "status", "intact", "length_aa", "best_ref", "best_ref_pid", "frameshifts", "stops", "qcov",
              "fabp_like", "rnaseq_cov", "rnaseq_intron_samples", "liver_rpm", "region"]
    write_tsv(out / "genes.tsv", genes, fields)

    # --- species summary
    rows = []
    for slug in ORDER:
        if slug not in sp_meta:
            continue
        m = sp_meta[slug]
        gs = [g for g in genes if g["slug"] == slug]
        r = {"slug": slug, "name": m["name"], "common": m["common"], "group": m["group"], "source": m["source"],
             "accession": m["accession"]}
        for cl in ("MUP", "OBP", "OBP2"):
            x = [g for g in gs if g["clade"] == cl and g["status"] != "fragment"]
            r[f"{cl}_intact"] = sum(1 for g in x if g["intact"] in (True, "True"))
            r[f"{cl}_other"] = len(x) - r[f"{cl}_intact"]
        for reg in ("MUP", "OBP1", "LCN2"):
            rr = [x for x in regions if x["slug"] == slug and x["region"] == reg]
            f = [x for x in rr if x["found"]]
            r[f"{reg}_region"] = ("open" if any(x.get("open") for x in f) else "closed") if f else "not found"
            r[f"{reg}_region_kb"] = round(sum(x["length"] for x in f) / 1000, 1) if f else ""
            r[f"{reg}_gaps"] = sum(x.get("n_gaps", 0) for x in f) if f else ""
            r[f"{reg}_gap_bp"] = sum(x.get("gap_bp", 0) for x in f) if f else ""
        rows.append(r)
    write_tsv(out / "species.tsv", rows)

    # --- quality checks
    qc = []

    def add(check, value, note=""):
        qc.append({"check": check, "value": value, "note": note})
    mm = [g for g in genes if g["slug"] == "mus_musculus"]
    add("mouse_Mup_genes_annotated", sum(1 for g in mm if re_mup.match(g["symbol"] or "")),
        "protein-coding Mup1-Mup22 in the RefSeq annotation")
    add("mouse_Mup_pseudogenes_annotated", sum(1 for g in mm if (g["symbol"] or "").startswith("Mup-ps")))
    add("mouse_MUP_clade_genes_tree", sum(1 for g in mm if g["clade"] == "MUP" and g["source"] == "RefSeq annotation"
                                          and g["status"].startswith("protein")))
    rl = read_tsv(D / "regions" / "region_loci.tsv") if (D / "regions" / "region_loci.tsv").exists() else []
    mrl = [l for l in rl if l["slug"] == "mus_musculus" and l["region"] == "MUP" and l["family"] == "MUP"
           and l["status"] != "fragment"]
    add("mouse_MUP_region_DNA_loci", len(mrl), "MUP-family miniprot loci in the mouse Mup region (positive control): "
        + ", ".join(f"{k} {v}" for k, v in sorted(Counter(l['status'] for l in mrl).items())))
    for slug in ("mus_spretus", "microtus_ochrogaster_ucd1", "microtus_ochrogaster_ucd2"):
        f = D / "genomes" / slug / "loci.tsv"
        if f.exists():
            ll = [l for l in read_tsv(f) if l["family"] == "MUP" and l["status"] != "fragment"]
            add(f"{slug}_MUP_loci_genome", len(ll), ", ".join(f"{k} {v}" for k, v in sorted(Counter(
                l['status'] for l in ll).items())) + "; whole-genome miniprot, fragments (<40% of query) excluded")
            fr = [l for l in read_tsv(f) if l["family"] == "MUP" and l["status"] == "fragment"]
            add(f"{slug}_MUP_fragments_genome", len(fr), "short MUP-best matches under 40% of the query (not counted)")
    ref = [c for c in cls.values() if c["slug"] == "uniprot"]
    ok = sum(1 for c in ref if (("MUP" in c["symbol"]) and c["clade"] == "MUP") or
             ((c["symbol"].startswith("Obp") or c["symbol"] in ("aphrodisin", "arvicolin")) and c["clade"] == "OBP"))
    add("uniprot_refs_in_expected_clade", f"{ok}/{len(ref)}", "MUP allergens in the MUP clade; aphrodisin, bank vole "
        "OBPs and arvicolin in the OBP clade")
    closed = [r for r in rows if r["MUP_region"] == "closed"]
    add("species_with_closed_MUP_region", f"{len(closed)}/{len(rows)}", "both flanking sides found on one sequence")
    write_tsv(out / "qc.tsv", qc)
    write_tsv(out / "regions.tsv", [x for x in regions], ["slug", "region", "part", "found", "region_id", "seqid",
                                                          "chrom", "start", "end", "length", "left_flank",
                                                          "right_flank", "open", "n_gaps", "gap_bp"])
    cds_rows = cds_spans(cfg, regions)
    write_tsv(out / "region_cds.tsv", cds_rows, ["slug", "gene_id", "cds_start", "cds_end"])
    cds = {(r["slug"], r["gene_id"]): (r["cds_start"], r["cds_end"]) for r in cds_rows}
    mrows, mgenes, mgaps = map_tables(D, regions, genes, sp_meta, cds, ORDER)
    write_tsv(out / "map_rows.tsv", mrows)
    write_tsv(out / "map_genes.tsv", mgenes,
              ["region_id", "start", "end", "strand", "symbol", "kind", "status", "tip", "desc"])
    write_tsv(out / "map_gaps.tsv", mgaps, ["region_id", "start", "end"])
    print(f"{len(genes)} genes, {len(rows)} species, {len(regions)} region rows -> {out}")


if __name__ == "__main__":
    main()
