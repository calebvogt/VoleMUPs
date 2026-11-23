"""Draw the project tile (1200 x 675 PNG, dark, Getz-tile style) from the exported synteny.json.

The tile shows the MUP neighbourhood (the DNA between Slc46a2 and Zfp37) on one shared scale: the 2.2 Mb mouse
cluster of MUP genes and pseudogenes against the single MUP gene of each vole.

Usage: python web/make_tile.py [--site ../calebvogt.github.io/projects/vole-mups]
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

# the website checkout next to this repo (calebvogt.github.io)
SITE = Path(__file__).resolve().parents[2] / "calebvogt.github.io" / "projects" / "vole-mups"
ROWS = [("mus_musculus", "House mouse"), ("rattus_norvegicus", "Norway rat"), ("mesocricetus_auratus", "Golden hamster"),
        ("peromyscus_maniculatus", "Deer mouse"), ("myodes_glareolus", "Bank vole"),
        ("microtus_pennsylvanicus", "Meadow vole"), ("microtus_ochrogaster_ucd1", "Prairie vole")]
BG, INK, MUTED = "#0d0d0d", "#c3c2b7", "#898781"
COL = {"MUP": "#f08a3c", "flank": "#6f6e68", "gene": "#3b3a37", "lipo": "#b9b39a", "OBP": "#3f9bdc", "OBP2": "#9d86e6"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", default=str(SITE))
    a = ap.parse_args()
    site = Path(a.site)
    syn = json.loads((site / "data" / "synteny.json").read_text(encoding="utf-8"))["regions"]["MUP"]
    by = {}
    for r in syn:
        by.setdefault(r["sp"], r)
    fig = plt.figure(figsize=(12, 6.75), dpi=100, facecolor=BG)
    ax = fig.add_axes([0.03, 0.08, 0.95, 0.86])
    ax.set_facecolor(BG)
    ax.axis("off")
    rows = [(s, n) for s, n in ROWS if s in by]
    maxspan = max(by[s]["hi"] - by[s]["lo"] for s, _ in rows)
    x0, x1 = 0.215, 1.0
    for i, (slug, name) in enumerate(rows):
        r = by[slug]
        y = len(rows) - 1 - i
        span = r["hi"] - r["lo"]
        sc = (x1 - x0) / maxspan
        ax.plot([x0, x0 + span * sc], [y, y], color="#4a4945", lw=1.2, solid_capstyle="butt")
        ax.text(0.0, y, name, color=INK if "vole" not in name.lower() else "#ffffff", fontsize=17, va="center",
                ha="left", fontweight="bold" if slug.startswith("microtus_ochrogaster") else "normal")
        for g in sorted(r["genes"], key=lambda g: g[4] != "gene"):
            s, e, strand, kind = g[0], g[1], g[2], g[4]
            a_ = (r["hi"] - e) if r["flip"] else (s - r["lo"])
            xa = x0 + max(0, a_) * sc
            w = max(0.0045 if kind == "MUP" else 0.002, (min(e, r["hi"]) - max(s, r["lo"])) * sc)
            hh = 0.13 if kind == "gene" else 0.30
            fwd = (strand == "+") != r["flip"]
            t = min(0.004, w * 0.45)
            pts = [(xa, y - hh), (xa + w - t, y - hh), (xa + w, y), (xa + w - t, y + hh), (xa, y + hh)] if fwd else \
                [(xa + w, y - hh), (xa + t, y - hh), (xa, y), (xa + t, y + hh), (xa + w, y + hh)]
            pseudo = "pseudo" in (g[5] or "") or "disrupted" in (g[5] or "")
            ax.add_patch(Polygon(pts, closed=True, facecolor=COL.get(kind, COL["gene"]),
                                 alpha=0.45 if pseudo else 1.0, lw=0))
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.8, len(rows) - 0.4)
    # scale bar
    bar = 500_000 * (x1 - x0) / maxspan
    ax.plot([x1 - bar, x1], [-0.55, -0.55], color=MUTED, lw=1.5)
    ax.text(x1 - bar / 2, -0.75, "500 kb", color=MUTED, fontsize=13, ha="center", va="top")
    out = site / "tile.png"
    tmp = out.with_name("tile.tmp.png")
    fig.savefig(tmp, facecolor=BG, dpi=100)
    plt.close(fig)
    os.replace(tmp, out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
