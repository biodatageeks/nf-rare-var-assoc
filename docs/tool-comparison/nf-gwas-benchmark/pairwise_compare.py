#!/usr/bin/env python3
"""
Pairwise comparison of nf-rare-var-assoc vs nf-gwas on the tool-comparison benchmark.

The two arms were run on the SAME set of simulated phenotype datasets but the
datasets differ wildly in difficulty (number/effect of causal genes), so pooling
the per-arm scores and comparing means is dominated by between-dataset variance
(huge SD) -- see runs/nf_*_eval/results/aggregate_scores. The correct analysis is
PAIRWISE: for each dataset, take metric(nf-rare-var-assoc) - metric(nf-gwas) and
test whether the mean of those paired differences is non-zero. This removes the
dataset-difficulty variance.

Primary metric: average precision (AP), the sklearn-equivalent column in each
per-dataset `*_compute_score_auc_summary.csv`. We also report auc_pr / auc_roc.

RECALL SCALING: each ranking metric (auc_roc, auc_pr, average_precision) is a
ranking quality measured only over the causal genes that actually appear in the
results. A pipeline can score well on the few causal genes it surfaced while
missing most of them. To fold gene-level recall into the metric, every metric is
multiplied by  n_causal_genes_in_results / n_causal_genes  (the fraction of true
causal genes present in the result set) BEFORE the pairwise comparison.

Inputs (per arm, one file per dataset_idx):
  runs/<arm>_eval/results/compute_score/<prefix>_dataset_idx_<N>_compute_score_auc_summary.csv

The old single-dataset run (filenames containing 'run18') is ignored: only files
matching '<prefix>_eval_dataset_idx_<N>_...' (no 'run18') are used.

Usage:
  python pairwise_compare.py                 # default paths under $RUNS
  python pairwise_compare.py --runs <dir> --out <dir>
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

# arm label -> (eval subdir, filename glob). The glob deliberately excludes the
# stale 'run18' single-dataset files (they do not match '*_eval_dataset_idx_*'
# because their token is 'run18_eval' / 'run18_nfgwas_eval' -- we filter below).
ARMS = {
    "nf_rare_var_assoc": "nf_rare_var_assoc_eval",
    "nf_gwas": "nf_gwas_eval",
}

# Metrics to report; AP is primary.
METRICS = ["average_precision", "auc_pr", "auc_roc"]

IDX_RE = re.compile(r"_dataset_idx_(\d+)_compute_score_auc_summary\.csv$")


def load_arm(runs_dir: str, eval_subdir: str) -> pd.DataFrame:
    """Return a DataFrame indexed by dataset_idx with one column per metric."""
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
    # If an idx somehow appears twice (it should not after run18 filtering), keep first.
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
    # Paired (one-sample on the differences) t-test, two-sided.
    t_res = stats.ttest_1samp(diff, 0.0)
    # Wilcoxon signed-rank as a distribution-free robustness check.
    try:
        w_res = stats.wilcoxon(diff)
        w_p = float(w_res.pvalue)
    except ValueError:
        # wilcoxon errors if all diffs are zero.
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
        "n_rva_wins": int(np.sum(diff > 0)),
        "n_gwas_wins": int(np.sum(diff < 0)),
        "n_ties": int(np.sum(diff == 0)),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=DEFAULT_RUNS, help="runs/ directory")
    ap.add_argument("--out", default=None,
                    help="output directory for CSVs (default: <runs>/pairwise_comparison)")
    args = ap.parse_args()

    out_dir = args.out or os.path.join(args.runs, "pairwise_comparison")
    os.makedirs(out_dir, exist_ok=True)

    rva = load_arm(args.runs, ARMS["nf_rare_var_assoc"])
    gwas = load_arm(args.runs, ARMS["nf_gwas"])

    common = rva.index.intersection(gwas.index)
    only_rva = rva.index.difference(gwas.index)
    only_gwas = gwas.index.difference(rva.index)

    print("=" * 70)
    print("Pairwise comparison: nf-rare-var-assoc vs nf-gwas")
    print("  (paired by dataset_idx; difference = nf-rare-var-assoc - nf-gwas)")
    print("=" * 70)
    print(f"datasets: nf-rare-var-assoc={len(rva)}, nf-gwas={len(gwas)}, "
          f"paired (common)={len(common)}")
    if len(only_rva):
        print(f"  only in nf-rare-var-assoc: {sorted(only_rva)}")
    if len(only_gwas):
        print(f"  only in nf-gwas: {sorted(only_gwas)}")
    print()

    # Per-dataset paired table for the primary metric (+ all metrics' diffs saved).
    per_ds = pd.DataFrame(index=common.sort_values())
    for met in METRICS:
        per_ds[f"{met}__rva"] = rva.loc[common, met]
        per_ds[f"{met}__gwas"] = gwas.loc[common, met]
        per_ds[f"{met}__diff"] = per_ds[f"{met}__rva"] - per_ds[f"{met}__gwas"]
    per_ds.index.name = "dataset_idx"
    per_ds_path = os.path.join(out_dir, "pairwise_per_dataset.csv")
    per_ds.to_csv(per_ds_path)

    # Summary stats per metric.
    summ_rows = []
    for met in METRICS:
        s = paired_stats(per_ds[f"{met}__diff"].to_numpy())
        s["metric"] = met
        s["mean_rva"] = float(per_ds[f"{met}__rva"].mean())
        s["mean_gwas"] = float(per_ds[f"{met}__gwas"].mean())
        summ_rows.append(s)
    summ = pd.DataFrame(summ_rows).set_index("metric")
    cols = ["n", "mean_rva", "mean_gwas", "mean_diff", "sd_diff", "se_diff",
            "ci95_low", "ci95_high", "t_stat", "t_pvalue", "wilcoxon_pvalue",
            "n_rva_wins", "n_gwas_wins", "n_ties"]
    summ = summ[cols]
    summ_path = os.path.join(out_dir, "pairwise_summary.csv")
    summ.to_csv(summ_path)

    # Pretty print.
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 50)
    print("Metrics are RECALL-SCALED: value * n_causal_genes_in_results / n_causal_genes.")
    print()
    print("Per-dataset average precision (primary metric, recall-scaled):")
    show = per_ds[["average_precision__rva", "average_precision__gwas",
                   "average_precision__diff"]].rename(columns={
        "average_precision__rva": "AP_rva",
        "average_precision__gwas": "AP_gwas",
        "average_precision__diff": "diff",
    })
    print(show.round(4).to_string())
    print()
    print("Paired difference summary (nf-rare-var-assoc - nf-gwas):")
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
    print(f"  wins            = nf-rare-var-assoc {int(ap_row['n_rva_wins'])} / "
          f"nf-gwas {int(ap_row['n_gwas_wins'])} / ties {int(ap_row['n_ties'])}")
    print("-" * 70)
    print()
    print(f"Wrote: {per_ds_path}")
    print(f"Wrote: {summ_path}")


if __name__ == "__main__":
    main()
