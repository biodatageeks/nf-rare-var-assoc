#!/usr/bin/env python3
"""Grouped bar chart: recall-scaled average precision per dataset, three methods.

One group of bars per benchmark dataset (x-axis), three bars per group:
  - nf-rare-var-assoc      (the reference pipeline; every dataset)
  - RICOPILI + nf-gwas     (every dataset it produced a result for)
  - RICOPILI + STAAR       (Full version; only the datasets it produced a result for)

A method with no result for a dataset simply has no bar there -- STAAR is sparse
by design (autosomes-only, run per dataset separately), and nf-gwas is missing
run_27 (too few cases). The y axis is recall-scaled average precision against the
complete list of causal genes, computed exactly as pairwise_compare.py does, by
reusing its load_arm().

Two output styles (--style):
  standalone (default) -- wide slide/report figure: in-plot title, run_<N> tick labels.
  panel                -- sized to be placed beside other figures on a printed page
                          (full text width of an A4 page with 2cm margins, ~17cm), so it
                          is included at 1:1 and its text stays legible: no in-plot title
                          (the caption carries it), numeric dataset ticks, print-scale fonts.

Usage:
  python three_method_ap_bar.py --runs <runs_dir> --out <out_stem>
  python three_method_ap_bar.py --style panel --out <tex_dir>/aux/three_method_ap_bar
Defaults target this workstation's layout.
"""
from __future__ import annotations

import argparse
import os

import numpy as np

# Reuse the exact recall-scaling + file parsing the pairwise comparison uses.
from pairwise_compare import load_arm, _style_axes, INK, INK_SOFT, GRID, SURFACE

# Categorical slots 1-3 of the dataviz default theme (documented all-pairs CVD-safe).
C_REF   = "#2a78d6"   # blue   -- nf-rare-var-assoc (its identity across the figures)
C_NFGW  = "#eb6834"   # orange -- RICOPILI + nf-gwas
C_STAAR = "#1baf7a"   # aqua   -- RICOPILI + STAAR (Full)

METHODS = [
    ("nf-rare-var-assoc", "nf_rare_var_assoc_eval",   C_REF),
    ("RICOPILI + nf-gwas", "ricopili_nf_gwas_qcmatched_eval",   C_NFGW),
    ("RICOPILI + STAAR",   "ricopili_staar_qcmatched_full_eval", C_STAAR),
]

DEFAULT_RUNS = "/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/tools_comparison/runs"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=DEFAULT_RUNS, help="runs/ directory")
    ap.add_argument("--out", default=None,
                    help="output stem (writes <stem>.png and <stem>.pdf)")
    ap.add_argument("--metric", default="average_precision",
                    choices=["average_precision", "auc_pr", "auc_roc"])
    ap.add_argument("--style", default="standalone", choices=["standalone", "panel"],
                    help="standalone = wide screen figure; panel = print size")
    args = ap.parse_args()
    panel = args.style == "panel"

    out_stem = args.out or os.path.join(args.runs, "pairwise_ricopili_staar",
                                        "three_method_ap_bar")
    os.makedirs(os.path.dirname(out_stem), exist_ok=True)

    # Load each method's recall-scaled measure, indexed by dataset_idx.
    series = {}
    for label, subdir, _ in METHODS:
        df = load_arm(args.runs, subdir)
        series[label] = df[args.metric]
        n_present = int(df[args.metric].notna().sum())
        print(f"{label:22s} {subdir:26s} datasets={n_present}")

    # x-axis = every dataset any method scored, ascending index.
    all_idx = sorted(set().union(*[s.index for s in series.values()]))
    x = np.arange(len(all_idx), dtype=float)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "font.size": 7 if panel else 9,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5 if panel else 0.6,
        "axes.axisbelow": True,
    })

    n_m = len(METHODS)
    bar_w = 0.8 / n_m
    # panel: ~10.9cm x 4.6cm -- 0.64 of an A4 text width, so the PDF is included at
    # ~1:1 and no font shrinks.
    fig, ax = plt.subplots(figsize=(4.35, 1.85) if panel else (12.5, 4.2))

    for j, (label, _subdir, color) in enumerate(METHODS):
        s = series[label]
        vals = np.array([s.get(i, np.nan) for i in all_idx], dtype=float)
        # offset so the group of n_m bars is centred on the tick
        offset = (j - (n_m - 1) / 2) * bar_w
        # Only draw bars that exist (NaN -> no bar, leaving a visible gap).
        mask = ~np.isnan(vals)
        ax.bar(x[mask] + offset, vals[mask], width=bar_w * 0.9,
               color=color, edgecolor=SURFACE, linewidth=0.5,
               label=label, zorder=3)

    ax.set_xticks(x)
    if panel:
        # Bare dataset numbers fit horizontally in the narrower panel; "run_" moves
        # to the axis label, which keeps the tick strip one line high and readable.
        ax.set_xticklabels([str(i) for i in all_idx], fontsize=5)
        ax.set_xlabel("Simulated dataset (run index)", fontsize=6.5)
    else:
        ax.set_xticklabels([f"run_{i}" for i in all_idx], rotation=90, fontsize=7)
    ax.set_xlim(-0.6, len(all_idx) - 0.4)
    ax.set_ylabel("Average precision\n(recall-scaled)" if panel
                  else "Average precision (recall-scaled)",
                  fontsize=6.5 if panel else 9.5)
    if not panel:
        ax.set_title("Per-dataset gene-detection performance across three pipelines "
                     "(recall-scaled average precision)",
                     fontsize=10, color=INK, pad=8)
    _style_axes(ax)
    ax.grid(axis="x", visible=False)
    if panel:
        ax.tick_params(labelsize=5, length=2, pad=1.5)

    handles = [Patch(facecolor=c, label=l) for l, _s, c in METHODS]
    # panel: upper LEFT -- the low-index datasets leave that corner empty, and the
    # right side carries the tallest bars. 5pt matches the paired-difference panel's key.
    ax.legend(handles=handles, frameon=False, fontsize=5 if panel else 8.5,
              loc="upper left" if panel else "upper right", ncol=1, labelcolor=INK,
              handlelength=1.2 if panel else 2.0,
              labelspacing=0.3 if panel else 0.5,
              borderpad=0.2 if panel else 0.4)

    # Footnote: which methods are sparse and why (identity is never colour-alone).
    #fig.text(0.008, -0.02,
    #         "RICOPILI + STAAR shown for its Full (relatedness-modelled) version, "
    #         "produced only for datasets 4, 11, 18 (autosomes-only, run per dataset). "
    #         "RICOPILI + nf-gwas missing run_27 (too few cases).",
    #         fontsize=7, color=INK_SOFT, ha="left")

    written = []
    for ext in ("png", "pdf"):
        p = f"{out_stem}.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        written.append(p)
    plt.close(fig)
    for p in written:
        print("Wrote:", p)


if __name__ == "__main__":
    main()
