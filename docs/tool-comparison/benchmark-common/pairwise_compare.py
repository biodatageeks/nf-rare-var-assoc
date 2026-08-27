#!/usr/bin/env python3
"""
Pairwise comparison of two benchmark arms on the tool-comparison datasets.

TOOL-AGNOSTIC generalisation of nf-gwas-benchmark/pairwise_compare.py: instead of a
hardcoded nf-rare-var-assoc-vs-nf-gwas pair, the two methods are given on the command
line (label + score subdirectory). The defaults reproduce the original nf-gwas
comparison exactly.
Difference is always arm_a - arm_b.

The two arms were run on the SAME simulated phenotype datasets but the datasets
differ widely in difficulty (number and effect size of the causal genes), so pooling
scores is dominated by between-dataset variance. The correct analysis is PAIRWISE:
for each dataset, take metric(A) - metric(B) and test whether the mean paired
difference is non-zero. This removes the dataset-difficulty variance.

Primary metric: average precision (AP). We also report auc_pr / auc_roc.

RECALL SCALING: each ranking metric is measured only over the causal genes that
actually appear in the results. To fold gene-level recall in, every metric is
multiplied by  n_causal_genes_in_results / n_causal_genes  BEFORE the comparison.

MISSING DATASETS (--missing):
  One method can fail to produce a result for a dataset the other handled -- for example
  a combined method whose quality control leaves too few cases, on a dataset the
  reference pipeline's own filtering kept analysable. Two ways to treat that dataset:
    drop  (default) -- pair only on datasets BOTH arms produced (intersection). Lenient:
            it silently excludes the dataset, which rewards a method for failing to run.
    zero  -- a method that produced NOTHING scores 0 (zero causal-gene recall, which is
            exactly what the recall-scaling formula yields for an empty result set), and
            the pair is kept. This penalises the failure instead of hiding it, and is the
            honest choice when one method genuinely could not analyse a valid dataset.
  Either way the COVERAGE (how many datasets each method produced a result for) is
  reported, because "ran 30/30 vs 29/30" is itself a result.

Inputs (per method, one file per dataset_idx):
  <runs>/<eval_subdir>/results/compute_score/<prefix>_dataset_idx_<N>_compute_score_auc_summary.csv

The old single-dataset run (filenames containing 'run18') is ignored.

Usage:
  python pairwise_compare.py                       # the default pair
  python pairwise_compare.py \
      --arm-a nf_rare_var_assoc nf_rare_var_assoc_eval \
      --arm-b ricopili_staar   chain_ricopili_staar_eval \
      --missing zero --runs <dir> --out <dir>
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

DEFAULT_RUNS = (
    "/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/"
    "tools_comparison/runs"
)

# Metrics to report; AP is primary.
METRICS = ["average_precision", "auc_pr", "auc_roc"]
METRIC_LABELS = {
    "average_precision": "Average precision",
    "auc_pr": "AUC-PR",
    "auc_roc": "AUC-ROC",
}

IDX_RE = re.compile(r"_dataset_idx_(\d+)_compute_score_auc_summary\.csv$")

# --- Figure palette (validated for colour-vision deficiency) ---
# Diverging by sign of the paired difference: first method better vs second.
C_A      = "#2a78d6"   # first method better  (diff > 0)
C_B      = "#e34948"   # second method better (diff < 0)
C_POINT  = "#2a78d6"   # single-series scatter marks
INK      = "#0b0b0b"   # text-primary
INK_SOFT = "#52514e"   # text-secondary
GRID     = "#e6e5e1"   # recessive grid / zero reference
SURFACE  = "#ffffff"   # figures sit on white


def load_arm(runs_dir: str, eval_subdir: str) -> pd.DataFrame:
    """Return a DataFrame indexed by dataset_idx with one recall-scaled column per metric."""
    cs_dir = os.path.join(runs_dir, eval_subdir, "results", "compute_score")
    pattern = os.path.join(cs_dir, "*_dataset_idx_*_compute_score_auc_summary.csv")
    rows = []
    for path in sorted(glob.glob(pattern)):
        base = os.path.basename(path)
        # Skip the stale single-dataset run (e.g. tools_comparison_run18_...).
        if "run18" in base:
            continue
        m = IDX_RE.search(base)
        if not m:
            continue
        idx = int(m.group(1))
        df = pd.read_csv(path)
        # proxy_scoring_mode='none' => a single 'exact' row; be defensive anyway.
        if "r2_threshold" in df.columns and (df["r2_threshold"] == "exact").any():
            row = df.loc[df["r2_threshold"] == "exact"].iloc[0]
        else:
            row = df.iloc[0]
        # Recall scaling factor: fraction of true causal genes present in results.
        n_in = float(row["n_causal_genes_in_results"])
        n_tot = float(row["n_causal_genes"])
        recall = (n_in / n_tot) if n_tot > 0 else np.nan
        rec = {"dataset_idx": idx, "recall_scale": recall,
               "n_causal_genes_in_results": n_in, "n_causal_genes": n_tot}
        for met in METRICS:
            rec[met] = (float(row[met]) * recall) if met in df.columns else np.nan
        rows.append(rec)
    if not rows:
        sys.exit(f"ERROR: no auc_summary files found under {cs_dir}")
    out = pd.DataFrame(rows).set_index("dataset_idx").sort_index()
    out = out[~out.index.duplicated(keep="first")]
    return out


def paired_stats(diff: np.ndarray) -> dict:
    """Mean/SD of paired differences + tests that the mean difference is non-zero."""
    diff = np.asarray(diff, dtype=float)
    diff = diff[~np.isnan(diff)]
    n = diff.size
    mean = float(np.mean(diff))
    sd = float(np.std(diff, ddof=1)) if n > 1 else float("nan")
    se = sd / np.sqrt(n) if n > 1 else float("nan")
    t_res = stats.ttest_1samp(diff, 0.0)
    try:
        w_res = stats.wilcoxon(diff)
        w_p = float(w_res.pvalue)
    except ValueError:
        w_p = float("nan")
    ci_low = ci_high = float("nan")
    if n > 1 and se == se:  # se not NaN
        tcrit = stats.t.ppf(0.975, df=n - 1)
        ci_low, ci_high = mean - tcrit * se, mean + tcrit * se
    return {
        "n": n,
        "mean_diff": mean,
        "sd_diff": sd,
        "se_diff": se,
        "ci95_low": ci_low,
        "ci95_high": ci_high,
        "t_stat": float(t_res.statistic),
        "t_pvalue": float(t_res.pvalue),
        "wilcoxon_pvalue": w_p,
        "n_a_wins": int(np.sum(diff > 0)),
        "n_b_wins": int(np.sum(diff < 0)),
        "n_ties": int(np.sum(diff == 0)),
    }


def _style_axes(ax) -> None:
    """Thin marks, recessive axes, no top/right spines."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_SOFT)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_SOFT, labelcolor=INK, length=3, width=0.8)
    ax.title.set_color(INK)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def _sig_star(p: float) -> str:
    return ("***" if p < 0.001 else "**" if p < 0.01 else
            "*" if p < 0.05 else "n.s.")


def _bootstrap_pvalue_fan(diffs: np.ndarray, n_boot: int, seed: int):
    """P-value trajectory over many RANDOM accumulation orders.

    For n_boot random permutations of the paired differences, recompute the paired
    t-test p-value on every growing prefix (n = 2..N). Returns, per n, the median and
    the 5-95 / 25-75 percentile envelope across orders. Every order ends on the SAME
    full set, so the fan pinches to the single full-sample p at n = N. Fully
    vectorised (no per-permutation scipy calls)."""
    d = diffs[~np.isnan(diffs)]
    N = d.size
    if N < 3:
        return None
    rng = np.random.default_rng(seed)
    perms = np.stack([rng.permutation(d) for _ in range(n_boot)])   # (B, N)
    ns = np.arange(1, N + 1, dtype=float)
    csum = np.cumsum(perms, axis=1)
    csum2 = np.cumsum(perms ** 2, axis=1)
    mean = csum / ns
    with np.errstate(invalid="ignore", divide="ignore"):
        var = (csum2 - csum ** 2 / ns) / (ns - 1)          # NaN/inf at n=1 (df 0)
        tstat = mean / np.sqrt(var / ns)
        p = 2.0 * stats.t.sf(np.abs(tstat), df=(ns - 1)[None, :])   # (B, N)
    xs = np.arange(2, N + 1)
    cols = p[:, 1:N]                                        # n = 2..N
    med = np.nanmedian(cols, axis=0)
    lo5, hi95 = np.nanpercentile(cols, [5, 95], axis=0)
    lo25, hi75 = np.nanpercentile(cols, [25, 75], axis=0)
    frac_sig = np.nanmean(cols < 0.05, axis=0)             # share of orders with p<0.05
    return xs, med, lo5, hi95, lo25, hi75, frac_sig


def make_plots(per_ds: pd.DataFrame, summ: pd.DataFrame,
               a_label: str, b_label: str, out_dir: str,
               n_boot: int = 2000, boot_seed: int = 0,
               style: str = "standalone",
               display_a: str | None = None, display_b: str | None = None) -> list[str]:
    """Three figures per comparison. Returns the paths written.

    1. Paired per-dataset difference in average precision (the headline), one bar
       per dataset, coloured by which method wins, with the mean +/- 95% CI overlaid.
    2. Paired scatter of one method against the other, with the y=x identity line.
    3. Forest plot: mean paired difference +/- 95% CI for all three metrics.
    Each is written as both .png (300 dpi, for previewing) and .pdf (vector, for print).
    Difference is always A - B, matching the tables.

    style='panel' re-sizes figure 1 to sit beside other figures on a printed page:
    ~6.6cm wide, print-scale fonts, and
    the in-plot title / per-dataset y labels dropped -- at that size they would set
    at 3-5pt, so the caption carries the win counts, the test and the sorting
    instead. Figures 2-3 keep their standalone size.
    """
    panel = style == "panel"
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
    except Exception as exc:  # pragma: no cover - matplotlib optional
        print(f"[plots] skipped -- matplotlib unavailable ({exc})")
        return []

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "font.size": 9,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.axisbelow": True,
    })
    # Figure text may spell the arms properly (display_*); filenames always use the
    # machine labels, so the output paths do not depend on the display names.
    a_name = display_a or a_label.replace("_", " ")
    b_name = display_b or b_label.replace("_", " ")
    stem = os.path.join(out_dir, f"{a_label}__vs__{b_label}")
    written: list[str] = []

    def _save(fig, suffix: str) -> None:
        for ext in ("png", "pdf"):
            p = f"{stem}__{suffix}.{ext}"
            fig.savefig(p, dpi=300, bbox_inches="tight")
            written.append(p)
        plt.close(fig)

    # ---- Figure 1: headline paired difference in average precision --------------
    d = per_ds["average_precision__diff"].dropna().sort_values()
    if len(d):
        s = summ.loc["average_precision"]
        colors = [C_A if v > 0 else C_B for v in d.to_numpy()]
        y = np.arange(len(d))
        h = max(2.4, 0.26 * len(d) + 1.0)
        # panel: ~5.4cm x 5.2cm -- 0.35 of an A4 text width, included at ~1:1.
        fig, ax = plt.subplots(figsize=(2.32, 2.05) if panel else (6.6, h))
        ax.barh(y, d.to_numpy(), color=colors, height=0.72, zorder=3)
        ax.axvline(0, color=INK_SOFT, lw=1.0, zorder=2)
        # mean +/- 95% CI band
        ax.axvspan(s["ci95_low"], s["ci95_high"], color=INK_SOFT, alpha=0.10, zorder=1)
        ax.axvline(s["mean_diff"], color=INK, lw=1.0 if panel else 1.4, ls="--", zorder=4,
                   label=f"mean {s['mean_diff']:+.3f}")
        star = _sig_star(float(s["t_pvalue"]))
        if panel:
            # Per-dataset ticks would set at ~4.5pt here; the panel's message is the
            # sorted distribution and where the mean/CI sits, so the ticks go and the
            # caption states the win counts, the paired t test and the sorting.
            ax.set_yticks([])
            ax.set_ylabel(f"{int(s['n'])} datasets, sorted by Δ", fontsize=6)
            ax.set_xlabel("Δ average precision (recall-scaled)", fontsize=6)
            # 4 ticks max: the default locator's labels collide at this width.
            ax.xaxis.set_major_locator(plt.MaxNLocator(4))
            ax.tick_params(labelsize=5)
            ax.xaxis.get_offset_text().set_fontsize(6)
        else:
            ax.set_yticks(y)
            ax.set_yticklabels([str(int(i)) for i in d.index])
            ax.set_ylabel("dataset")
            ax.set_xlabel(f"Δ average precision  ({a_name} − {b_name}, recall-scaled)")
            ax.set_title(
                f"Per-dataset paired difference\n"
                f"{a_name} wins {int(s['n_a_wins'])} / {b_name} wins {int(s['n_b_wins'])}"
                f"   ·   paired t p={s['t_pvalue']:.3g} ({star}), n={int(s['n'])}",
                fontsize=9.5, loc="left")
        ax.grid(axis="y", visible=False)
        ax.margins(y=0.01)
        legend = [
            Patch(facecolor=C_A, label=f"{a_name} better"),
            Patch(facecolor=C_B, label=f"{b_name} better"),
        ]
        if panel:
            # The mean/CI marks carry the statistical claim, so name them in the key.
            legend.append(Line2D([0], [0], color=INK, lw=1.0, ls="--",
                                 label=f"mean {s['mean_diff']:+.3f} (95% CI)"))
        if panel:
            # Mid-height rows are the near-zero differences, so a vertically centred key
            # pushed to the right sits over empty plot area whatever the data does.
            ax.legend(handles=legend, frameon=False, fontsize=5,
                      loc="center right", bbox_to_anchor=(1.04, 0.5),
                      handlelength=1.2, labelspacing=0.3, borderpad=0.2)
        else:
            ax.legend(handles=legend, frameon=False, fontsize=8, loc="lower right")
        _style_axes(ax)
        _save(fig, "paired_diff_ap")

    # ---- Figure 2: paired scatter A vs B ---------------------------------------
    sc = per_ds[["average_precision__a", "average_precision__b"]].dropna()
    if len(sc):
        xa = sc["average_precision__a"].to_numpy()
        xb = sc["average_precision__b"].to_numpy()
        lo = float(min(xa.min(), xb.min()))
        hi = float(max(xa.max(), xb.max()))
        pad = 0.04 * (hi - lo if hi > lo else 1.0)
        lo, hi = lo - pad, hi + pad
        n_a_better = int(np.sum(xa > xb))
        n_b_better = int(np.sum(xb > xa))
        fig, ax = plt.subplots(figsize=(4.8, 4.8))
        ax.plot([lo, hi], [lo, hi], color=INK_SOFT, lw=1.0, ls="--", zorder=2)
        ax.scatter(xa, xb, s=46, color=C_POINT, edgecolor=SURFACE, linewidth=1.2,
                   alpha=0.9, zorder=3)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"average precision  —  {a_name}")
        ax.set_ylabel(f"average precision  —  {b_name}")
        ax.set_title("Paired per-dataset scores (recall-scaled)", fontsize=9.5,
                     loc="left")
        # which side of the identity line each region favours
        ax.text(0.04, 0.96, f"above line: {b_name} better ({n_b_better})",
                transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
                color=INK_SOFT)
        ax.text(0.96, 0.04, f"below line: {a_name} better ({n_a_better})",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
                color=INK_SOFT)
        _style_axes(ax)
        _save(fig, "scatter_ap")

    # ---- Figure 3: forest plot of mean differences across metrics --------------
    mets = [m for m in METRICS if m in summ.index]
    if mets:
        y = np.arange(len(mets))[::-1]
        means = summ.loc[mets, "mean_diff"].to_numpy(dtype=float)
        los = summ.loc[mets, "ci95_low"].to_numpy(dtype=float)
        his = summ.loc[mets, "ci95_high"].to_numpy(dtype=float)
        colors = [C_A if m > 0 else C_B for m in means]
        fig, ax = plt.subplots(figsize=(6.2, 2.4))
        ax.axvline(0, color=INK_SOFT, lw=1.0, zorder=2)
        for yi, m, l, h, c in zip(y, means, los, his, colors):
            ax.plot([l, h], [yi, yi], color=c, lw=2.0, zorder=3,
                    solid_capstyle="round")
            ax.scatter([m], [yi], s=42, color=c, edgecolor=SURFACE, linewidth=1.2,
                       zorder=4)
        for yi, met, m in zip(y, mets, means):
            p = float(summ.loc[met, "t_pvalue"])
            ax.text(m, yi + 0.16, f"{m:+.3f} ({_sig_star(p)})", ha="center",
                    va="bottom", fontsize=7.5, color=INK)
        ax.set_yticks(y)
        ax.set_yticklabels([METRIC_LABELS.get(m, m) for m in mets])
        ax.set_ylim(-0.6, len(mets) - 0.4)
        ax.set_xlabel(f"mean paired difference  ({a_name} − {b_name})   with 95% CI")
        ax.set_title("Effect size across metrics", fontsize=9.5, loc="left")
        ax.grid(axis="y", visible=False)
        _style_axes(ax)
        _save(fig, "forest_metrics")

    # ---- Figure 4: p-value vs number of datasets (sequential, index order) -----
    # Datasets are already randomly indexed, so adding them in ascending index order
    # is a random accretion. For each prefix of n datasets we recompute the paired
    # test on average precision and plot how the p-value settles as n grows -- the
    # evidence that the result is not an artefact of stopping at a particular n.
    d_all = per_ds["average_precision__diff"].to_numpy(dtype=float)
    idx_order = list(per_ds.index)
    if len(d_all) >= 2:
        ns, t_p, w_p = [], [], []
        for n in range(2, len(d_all) + 1):
            d = d_all[:n]
            d = d[~np.isnan(d)]
            if d.size < 2:
                continue
            ns.append(n)
            t_p.append(float(stats.ttest_1samp(d, 0.0).pvalue))
            try:
                w_p.append(float(stats.wilcoxon(d).pvalue))
            except ValueError:
                w_p.append(np.nan)
        fig, ax = plt.subplots(figsize=(6.6, 3.6))
        ax.axhline(0.05, color=INK_SOFT, lw=1.0, ls="--", zorder=2)
        ax.text(ns[0], 0.05, " p = 0.05", color=INK_SOFT, fontsize=7.5,
                va="bottom", ha="left")
        ax.plot(ns, t_p, color=C_A, lw=2.0, marker="o", ms=4,
                solid_capstyle="round", label="paired t-test", zorder=4)
        ax.plot(ns, w_p, color="#eb6834", lw=2.0, ls="--", marker="s", ms=3.5,
                solid_capstyle="round", label="Wilcoxon", zorder=3)
        # direct-label the final values
        ax.annotate(f"{t_p[-1]:.2f}", (ns[-1], t_p[-1]), color=C_A, fontsize=8,
                    xytext=(4, 0), textcoords="offset points", va="center")
        ax.set_ylim(0, 1.0)
        ax.set_xlim(ns[0] - 0.5, ns[-1] + 1.5)
        ax.set_xlabel("number of datasets (added in ascending index order)")
        ax.set_ylabel("p-value  —  average precision")
        ax.set_title("P-value vs sample size (sequential)", fontsize=9.5, loc="left")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
        _style_axes(ax)
        _save(fig, "pvalue_vs_n")

    # ---- Figure 5: bootstrapped p-value envelope over random orders ------------
    # Fig 4 is one accumulation order (index order); this is the distribution over
    # MANY random orders, which separates "the conclusion is stable" from "we happened
    # to add datasets in a lucky order". The fan pinches to the full-sample p at n=N.
    if n_boot and len(d_all) >= 3:
        fan = _bootstrap_pvalue_fan(d_all, n_boot, boot_seed)
        if fan is not None:
            xs, med, lo5, hi95, lo25, hi75, frac_sig = fan
            from matplotlib.ticker import MaxNLocator
            fig, ax = plt.subplots(figsize=(6.6, 3.6))
            ax.axhline(0.05, color=INK_SOFT, lw=1.0, ls="--", zorder=2)
            ax.text(xs[0], 0.05, " p = 0.05", color=INK_SOFT, fontsize=7.5,
                    va="bottom", ha="left")
            ax.fill_between(xs, lo5, hi95, color=C_A, alpha=0.14, lw=0,
                            label="5-95th percentile", zorder=2)
            ax.fill_between(xs, lo25, hi75, color=C_A, alpha=0.26, lw=0,
                            label="25-75th percentile", zorder=3)
            ax.plot(xs, med, color=C_A, lw=2.0, zorder=4, label="median", solid_capstyle="round")
            ax.scatter([xs[-1]], [med[-1]], s=34, color=C_A, edgecolor=SURFACE,
                       linewidth=1.2, zorder=5)
            ax.annotate(f"{med[-1]:.2f}", (xs[-1], med[-1]), color=C_A, fontsize=8,
                        xytext=(4, 0), textcoords="offset points", va="center")
            ax.set_ylim(0, 1.0)
            ax.set_xlim(xs[0] - 0.5, xs[-1] + 1.5)
            ax.set_xlabel("number of datasets (random draw without replacement)")
            ax.set_ylabel("p-value  —  average precision")
            ax.set_title(f"P-value envelope over {n_boot} random accumulation orders",
                         fontsize=9.5, loc="left")
            ax.legend(frameon=False, fontsize=8, loc="upper right")
            ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
            _style_axes(ax)
            _save(fig, "pvalue_vs_n_bootstrap")
            n_peak = int(xs[int(np.argmax(frac_sig))])
            print(f"[bootstrap] over {n_boot} random orders, the largest share reaching "
                  f"p<0.05 at any n is {frac_sig.max():.1%} (at n={n_peak}); "
                  f"at the full n={int(xs[-1])} it is {frac_sig[-1]:.1%}.")

    for p in written:
        print(f"Wrote: {p}")
    return written


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=DEFAULT_RUNS, help="runs/ directory")
    ap.add_argument("--arm-a", nargs=2, metavar=("LABEL", "EVAL_SUBDIR"),
                    default=["nf_rare_var_assoc", "nf_rare_var_assoc_eval"],
                    help="first method: label and its score subdirectory under <runs>/")
    ap.add_argument("--arm-b", nargs=2, metavar=("LABEL", "EVAL_SUBDIR"),
                    default=["nf_gwas", "nf_gwas_eval"],
                    help="arm B label and its eval subdir under <runs>/")
    ap.add_argument("--out", default=None,
                    help="output directory for CSVs + figures (default: <runs>/pairwise_comparison)")
    ap.add_argument("--missing", choices=["drop", "zero"], default="drop",
                    help="datasets one arm did not produce: 'drop' pairs only on the "
                         "intersection (lenient); 'zero' keeps them and scores the "
                         "missing arm 0 (penalises the failure). Default: drop.")
    ap.add_argument("--no-plots", action="store_true",
                    help="skip the figures (write only the CSVs)")
    ap.add_argument("--display-a", default=None, metavar="NAME",
                    help="how arm A is spelled in figure text (default: its label with "
                         "underscores as spaces); output filenames are unaffected")
    ap.add_argument("--display-b", default=None, metavar="NAME",
                    help="how arm B is spelled in figure text")
    ap.add_argument("--style", default="standalone", choices=["standalone", "panel"],
                    help="'panel' re-sizes the paired-difference figure to sit "
                         "beside others on a printed page; see make_plots()")
    ap.add_argument("--bootstrap", type=int, default=2000, metavar="N",
                    help="random accumulation orders for the p-value-envelope figure "
                         "(0 disables it). Default: 2000.")
    ap.add_argument("--bootstrap-seed", type=int, default=0,
                    help="RNG seed for the bootstrap figure (reproducibility). Default: 0.")
    args = ap.parse_args()

    a_label, a_subdir = args.arm_a
    b_label, b_subdir = args.arm_b
    out_dir = args.out or os.path.join(args.runs, "pairwise_comparison")
    os.makedirs(out_dir, exist_ok=True)

    a = load_arm(args.runs, a_subdir)
    b = load_arm(args.runs, b_subdir)

    union = a.index.union(b.index)
    only_a = a.index.difference(b.index)   # datasets only A produced -> B failed them
    only_b = b.index.difference(a.index)   # datasets only B produced -> A failed them

    print("=" * 70)
    print(f"Pairwise comparison: {a_label} vs {b_label}   [missing={args.missing}]")
    print(f"  (paired by dataset_idx; difference = {a_label} - {b_label})")
    print("=" * 70)
    print(f"coverage: {a_label} produced {len(a)}/{len(union)} datasets, "
          f"{b_label} produced {len(b)}/{len(union)}")
    if len(only_a):
        print(f"  {b_label} FAILED to produce (only {a_label} has): {sorted(only_a)}")
    if len(only_b):
        print(f"  {a_label} FAILED to produce (only {b_label} has): {sorted(only_b)}")

    if args.missing == "zero":
        # An arm that produced no result for a dataset scores 0 on every recall-scaled
        # metric (empty result set => zero causal-gene recall). Keep the dataset paired.
        idx = union.sort_values()
        a_use = a.reindex(idx)
        b_use = b.reindex(idx)
        for df in (a_use, b_use):
            for met in METRICS:
                df[met] = df[met].fillna(0.0)
        n_paired = len(idx)
        if len(only_a) or len(only_b):
            print(f"  -> penalising: {len(only_b)} dataset(s) scored 0 for {a_label}, "
                  f"{len(only_a)} scored 0 for {b_label} (kept in the pairing).")
    else:
        idx = a.index.intersection(b.index).sort_values()
        a_use, b_use = a, b
        n_paired = len(idx)
        if len(only_a) or len(only_b):
            print(f"  -> dropping {len(only_a.union(only_b))} unpaired dataset(s) "
                  f"(lenient; not counted against either arm).")
    print(f"paired (scored) datasets: {n_paired}")
    print()

    per_ds = pd.DataFrame(index=idx)
    for met in METRICS:
        per_ds[f"{met}__a"] = a_use.loc[idx, met]
        per_ds[f"{met}__b"] = b_use.loc[idx, met]
        per_ds[f"{met}__diff"] = per_ds[f"{met}__a"] - per_ds[f"{met}__b"]
    per_ds.index.name = "dataset_idx"
    per_ds_path = os.path.join(out_dir, "pairwise_per_dataset.csv")
    per_ds.to_csv(per_ds_path)

    summ_rows = []
    for met in METRICS:
        s = paired_stats(per_ds[f"{met}__diff"].to_numpy())
        s["metric"] = met
        s["mean_a"] = float(per_ds[f"{met}__a"].mean())
        s["mean_b"] = float(per_ds[f"{met}__b"].mean())
        summ_rows.append(s)
    summ = pd.DataFrame(summ_rows).set_index("metric")
    cols = ["n", "mean_a", "mean_b", "mean_diff", "sd_diff", "se_diff",
            "ci95_low", "ci95_high", "t_stat", "t_pvalue", "wilcoxon_pvalue",
            "n_a_wins", "n_b_wins", "n_ties"]
    summ = summ[cols]
    summ_path = os.path.join(out_dir, "pairwise_summary.csv")
    summ.to_csv(summ_path)

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 50)
    print("Metrics are RECALL-SCALED: value * n_causal_genes_in_results / n_causal_genes.")
    print()
    print("Per-dataset average precision (primary metric, recall-scaled):")
    show = per_ds[["average_precision__a", "average_precision__b",
                   "average_precision__diff"]].rename(columns={
        "average_precision__a": f"AP_{a_label}",
        "average_precision__b": f"AP_{b_label}",
        "average_precision__diff": "diff",
    })
    print(show.round(4).to_string())
    print()
    print(f"Paired difference summary ({a_label} - {b_label}):")
    print(summ.round(5).to_string())
    print()

    ap_row = summ.loc["average_precision"]
    star = ("***" if ap_row["t_pvalue"] < 0.001 else
            "**" if ap_row["t_pvalue"] < 0.01 else
            "*" if ap_row["t_pvalue"] < 0.05 else "n.s.")
    print("-" * 70)
    print("HEADLINE (recall-scaled average precision):")
    print(f"  mean difference = {ap_row['mean_diff']:+.4f}  "
          f"(SD {ap_row['sd_diff']:.4f}, n={int(ap_row['n'])})")
    print(f"  95% CI          = [{ap_row['ci95_low']:+.4f}, {ap_row['ci95_high']:+.4f}]")
    print(f"  paired t-test   = t({int(ap_row['n'])-1})={ap_row['t_stat']:.3f}, "
          f"p={ap_row['t_pvalue']:.3g} ({star})")
    print(f"  Wilcoxon        = p={ap_row['wilcoxon_pvalue']:.3g}")
    print(f"  wins            = {a_label} {int(ap_row['n_a_wins'])} / "
          f"{b_label} {int(ap_row['n_b_wins'])} / ties {int(ap_row['n_ties'])}")
    print("-" * 70)
    print()
    print(f"Wrote: {per_ds_path}")
    print(f"Wrote: {summ_path}")
    if not args.no_plots:
        make_plots(per_ds, summ, a_label, b_label, out_dir,
                   n_boot=args.bootstrap, boot_seed=args.bootstrap_seed,
                   style=args.style,
                   display_a=args.display_a, display_b=args.display_b)


if __name__ == "__main__":
    main()
