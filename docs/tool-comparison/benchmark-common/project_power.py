#!/usr/bin/env python3
"""
Projected statistical power vs number of datasets, for a paired benchmark comparison.

Answers a planning question the descriptive comparison cannot: "if we GENERATED more
datasets from the same simulator, how likely is the paired test to reach p < 0.05?"
This is an EXTRAPOLATION beyond the data we have, so it rests on an assumption that
must be stated plainly:

  It assumes the TRUE per-dataset effect equals the CURRENT point estimate
  (Cohen's d_z = mean_diff / sd_diff, measured on the datasets we already ran).

That is the optimistic-but-standard basis for a power projection. The current result's
95% CI includes 0, so the true effect may well be null -- in which case no amount of
extra data crosses 0.05 (the rejection rate stays at the 5% false-positive level). Read
the curve as "power IF the effect is as large as we currently measure", i.e. a best case.

Two independent estimates are drawn, and they should agree:
  * analytic  -- noncentral-t power of the one-sample (paired) t-test at effect d_z.
  * bootstrap -- resample the observed paired differences WITH replacement to size N
                 and count how often p < 0.05 (captures the skew the analytic curve
                 assumes away).

Input is the `pairwise_per_dataset.csv` that pairwise_compare.py writes (column
`<metric>__diff`, default average_precision). Output: a figure (.png + .pdf) and a
printed table.

Usage:
  python project_power.py --per-dataset-csv <.../pairwise_per_dataset.csv> \
      --out <dir> [--proposed 50 60] [--target-power 0.8] [--max-n 160]
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

# Palette (matches pairwise_compare.py).
C_ANALYTIC = "#2a78d6"   # blue
C_BOOT     = "#eb6834"   # orange
INK        = "#0b0b0b"
INK_SOFT   = "#52514e"
GRID       = "#e6e5e1"
SURFACE    = "#ffffff"


def analytic_power(n: np.ndarray, dz: float, alpha: float = 0.05) -> np.ndarray:
    """Two-sided one-sample t-test power at standardized effect dz for sample size n."""
    n = np.asarray(n, dtype=float)
    df = n - 1
    tcrit = stats.t.ppf(1 - alpha / 2, df)
    ncp = dz * np.sqrt(n)
    return stats.nct.sf(tcrit, df, ncp) + stats.nct.cdf(-tcrit, df, ncp)


def bootstrap_power(diffs: np.ndarray, sizes, n_boot: int, seed: int, alpha: float = 0.05):
    """P(p<0.05) and median p when N differences are resampled with replacement."""
    rng = np.random.default_rng(seed)
    powers, med_p = [], []
    for N in sizes:
        samp = rng.choice(diffs, size=(n_boot, N), replace=True)
        mean = samp.mean(axis=1)
        sd = samp.std(axis=1, ddof=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            t = mean / (sd / np.sqrt(N))
            p = 2.0 * stats.t.sf(np.abs(t), df=N - 1)
        powers.append(float(np.mean(p < alpha)))
        med_p.append(float(np.median(p)))
    return np.array(powers), np.array(med_p)


def n_for_power(dz: float, target: float, alpha: float = 0.05, nmax: int = 5000) -> int | None:
    for n in range(3, nmax + 1):
        if analytic_power(n, dz, alpha) >= target:
            return n
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-dataset-csv", required=True,
                    help="pairwise_per_dataset.csv from pairwise_compare.py")
    ap.add_argument("--metric", default="average_precision",
                    help="metric prefix; uses column '<metric>__diff'. Default: average_precision")
    ap.add_argument("--out", default=None, help="output dir (default: alongside the CSV)")
    ap.add_argument("--proposed", type=int, nargs="+", default=[50, 60],
                    help="proposed total dataset counts to highlight. Default: 50 60")
    ap.add_argument("--target-power", type=float, default=0.8,
                    help="power threshold to mark and solve for. Default: 0.8")
    ap.add_argument("--max-n", type=int, default=160, help="right edge of the curve. Default: 160")
    ap.add_argument("--bootstrap", type=int, default=3000, help="bootstrap replicates. Default: 3000")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()

    col = f"{args.metric}__diff"
    df = pd.read_csv(args.per_dataset_csv)
    if col not in df.columns:
        raise SystemExit(f"ERROR: column {col!r} not in {args.per_dataset_csv}")
    d = df[col].dropna().to_numpy(dtype=float)
    n0 = d.size
    mean = float(d.mean())
    sd = float(d.std(ddof=1))
    dz = mean / sd
    p_now = float(stats.ttest_1samp(d, 0.0).pvalue)
    out_dir = args.out or os.path.dirname(os.path.abspath(args.per_dataset_csv))
    os.makedirs(out_dir, exist_ok=True)

    n_target = n_for_power(dz, args.target_power, args.alpha)

    print("=" * 70)
    print(f"Projected power for {args.metric} (paired t-test, alpha={args.alpha})")
    print("=" * 70)
    print(f"observed: n={n0}, mean diff={mean:+.4f}, sd={sd:.4f}, Cohen d_z={dz:.3f}, "
          f"current p={p_now:.3f}")
    print("ASSUMES the true effect equals this estimate; the current 95% CI includes 0,")
    print("so if the effect is really null, power stays at the 5% false-positive rate.")
    print()

    report_ns = sorted(set([n0] + list(args.proposed) + [80, 100] +
                           ([n_target] if n_target else [])))
    boot_pow, boot_medp = bootstrap_power(d, report_ns, args.bootstrap, args.seed, args.alpha)
    print(f"{'N':>5}  {'analytic power':>14}  {'bootstrap power':>15}  {'median p (boot)':>15}")
    for N, bp, mp in zip(report_ns, boot_pow, boot_medp):
        tag = ""
        if N == n0:
            tag = "  <- current"
        elif N in args.proposed:
            tag = "  <- proposed"
        elif N == n_target:
            tag = f"  <- {int(args.target_power*100)}% power"
        print(f"{N:>5}  {analytic_power(N, dz):>14.2f}  {bp:>15.2f}  {mp:>15.3f}{tag}")
    if n_target:
        print(f"\n=> {int(args.target_power*100)}% power needs about n={n_target} datasets "
              f"(analytic). At the proposed {args.proposed}, power is only "
              f"{', '.join(f'{analytic_power(N, dz):.0%}' for N in args.proposed)}.")

    # ---- Figure -------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
    except Exception as exc:  # pragma: no cover
        print(f"[plot] skipped -- matplotlib unavailable ({exc})")
        return

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "font.size": 9,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.axisbelow": True,
    })
    ns = np.arange(n0, args.max_n + 1)
    curve = analytic_power(ns, dz, args.alpha)
    dense_sizes = list(range(n0, args.max_n + 1, 10))
    dense_pow, _ = bootstrap_power(d, dense_sizes, args.bootstrap, args.seed, args.alpha)

    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    # proposed range band
    if len(args.proposed) >= 2:
        ax.axvspan(min(args.proposed), max(args.proposed), color=INK_SOFT, alpha=0.08,
                   zorder=1, label=f"proposed {min(args.proposed)}-{max(args.proposed)}")
    ax.axhline(args.target_power, color=INK_SOFT, lw=1.0, ls="--", zorder=2)
    ax.text(ns[-1], args.target_power, f" {int(args.target_power*100)}% power",
            color=INK_SOFT, fontsize=7.5, va="bottom", ha="right")
    ax.plot(ns, curve, color=C_ANALYTIC, lw=2.2, zorder=4,
            label="analytic (noncentral-t)", solid_capstyle="round")
    ax.scatter(dense_sizes, dense_pow, s=26, color=C_BOOT, edgecolor=SURFACE,
               linewidth=1.0, zorder=5, label="bootstrap resample")
    # current n marker
    ax.scatter([n0], [analytic_power(n0, dz)], s=48, color=INK, zorder=6)
    ax.annotate(f"now: n={n0}, {analytic_power(n0, dz):.0%}",
                (n0, analytic_power(n0, dz)), color=INK, fontsize=8,
                xytext=(6, -2), textcoords="offset points", va="top")
    if n_target and n_target <= args.max_n:
        ax.axvline(n_target, color=INK_SOFT, lw=0.8, ls=":", zorder=2)
        ax.annotate(f"n={n_target}", (n_target, 0.02), color=INK_SOFT, fontsize=7.5,
                    xytext=(3, 0), textcoords="offset points", va="bottom")
    ax.set_ylim(0, 1.0)
    ax.set_xlim(n0, args.max_n)
    ax.set_xlabel("total number of datasets")
    ax.set_ylabel("probability of reaching p < 0.05  (power)")
    ax.set_title(
        "Projected power vs number of datasets\n"
        f"if the true effect equals the current estimate "
        f"(mean Δ={mean:+.3f}, Cohen d_z={dz:.2f})",
        fontsize=9.5, loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_SOFT); ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_SOFT, labelcolor=INK, length=3, width=0.8)

    stem = os.path.join(out_dir, f"power_projection_{args.metric}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches="tight")
        print(f"Wrote: {stem}.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    main()
