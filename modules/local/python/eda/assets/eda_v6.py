#!/usr/bin/env python3
"""Exploratory data analysis for a multisample VCF (v6).

Single-pass cyvcf2 scan (replaces v5's two-pass pysam scan):
- Per-variant stats + per-sample aggregates.
- Heatmaps built in the same pass.  v5 needed a second pass only to learn the
  per-group maxima that size the heatmap bins; here those maxima are supplied
  up front (currently hardcoded, see _TEMP_MAX_BY_GROUP) so one pass suffices.

FORMAT fields (DP/GQ/DS) are read as vectorized numpy arrays via
``v.format(field)`` instead of v5's per-sample Python loop -- this is the
main speedup. Genotype-derived stats (het rate, missingness, the GT
encoding for the GT-vs-DS heatmap) come from ``v.gt_types`` /
``v.genotype.array()`` and are equivalent to v5's per-sample encoding.
"""
import os
from pathlib import Path
import sys
import argparse
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
from cyvcf2 import VCF
import time
import csv as _csv


tmp_matplotlib_dir = Path("matplotlib_tmp")
tmp_matplotlib_dir.mkdir(exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(tmp_matplotlib_dir.resolve())

plt.rcParams.update({'font.size': 20})

stats_dir = "eda_stats"
Path(stats_dir).mkdir(exist_ok=True)

output_dir = "plots"
Path(output_dir).mkdir(exist_ok=True)


HEATMAP_PAIRS = [('GQ', 'DS'), ('DP', 'DS'), ('GQ', 'DP'), ('GT', 'DS')]

# Fixed heatmap axis maxima, identical for every group (overall and per
# phenotype).  v5 derived these from per-group data maxima in a separate pass;
# fixed values let us size the bins up front and build the heatmaps in a single
# pass.  Values beyond a max fall outside the 2D histogram and are dropped --
# acceptable for these coarse diagnostic plots.
HEATMAP_STAT_MAX = {'DP': 300, 'GQ': 100, 'GT': 2, 'DS': 2}


def emit_stat(name, data):
    arr = np.asarray(data, dtype=float)
    with open(str(Path(stats_dir) / f"{name}.csv"), 'w', newline='') as fh:
        w = _csv.writer(fh)
        if arr.ndim == 1:
            for v in arr:
                w.writerow(['' if np.isnan(v) else v])
        else:
            for row in arr:
                w.writerow(['' if np.isnan(v) else v for v in row])


def current_milli_time():
    return round(time.time() * 1000)


def load_phenotype(phenotype_file):
    print(f"load_phenotype({phenotype_file})  ts = {current_milli_time()}")
    return pl.read_csv(phenotype_file, separator='\t')


def build_groups(samples, pheno_df):
    phenotypes = pheno_df['Y1'].unique().to_list()
    pheno_map = {row[0]: row[1] for row in pheno_df.select(['IID', 'Y1']).iter_rows()}
    group_indices = {p: [] for p in phenotypes}
    for idx, sample in enumerate(samples):
        pheno = pheno_map.get(sample)
        if pheno in group_indices:
            group_indices[pheno].append(idx)
    for p in group_indices:
        group_indices[p] = np.asarray(group_indices[p], dtype=int)
    return phenotypes, group_indices


def _init_vstats(stats, group_labels, percentiles):
    vstats = {}
    for stat in stats:
        vstats[stat] = {}
        for group in group_labels:
            vstats[stat][group] = {'mean': []}
            for p in percentiles:
                vstats[stat][group][p] = []
    return vstats


def _append_variant_stats(vstats, stat, values, group_labels, group_indices, percentiles):
    for group in group_labels:
        if group == 'all':
            group_vals = values
        else:
            idxs = group_indices[group]
            group_vals = values[idxs]
        if np.isnan(group_vals).all():
            mean_val = np.nan
            perc_vals = [np.nan for _ in percentiles]
        else:
            mean_val = np.nanmean(group_vals)
            perc_vals = [np.quantile(group_vals, p / 100.0) for p in percentiles]
        vstats[stat][group]['mean'].append(mean_val)
        for p, v in zip(percentiles, perc_vals):
            vstats[stat][group][p].append(v)


def _append_missingness(missingness_variants, gt_vals, group_labels, group_indices):
    is_missing = np.isnan(gt_vals)
    for group in group_labels:
        if group == 'all':
            miss_rate = float(np.mean(is_missing))
        else:
            idxs = group_indices[group]
            if idxs.size == 0:
                miss_rate = np.nan
            else:
                miss_rate = float(np.mean(is_missing[idxs]))
        missingness_variants[group].append(miss_rate)


def _stat_bins_arange(max_val):
    if max_val is None or not np.isfinite(max_val) or max_val < 0:
        return None
    return np.arange(0, max_val + 1, 1)


def _stat_bins_linspace(max_val):
    if max_val is None or not np.isfinite(max_val) or max_val < 0:
        return None
    return np.linspace(0, max_val, 201)


def _hist2d(values1, values2, xbins, ybins):
    mask = np.isfinite(values1) & np.isfinite(values2)
    if not np.any(mask):
        return np.zeros((len(xbins) - 1, len(ybins) - 1), dtype=int)
    hist, _, _ = np.histogram2d(values1[mask], values2[mask], bins=[xbins, ybins])
    return hist.astype(int)


def _format_field(v, field, n_samples):
    """Return a length-n_samples float array for a per-sample FORMAT field.

    Missing values become NaN. Int fields (DP/GQ) use cyvcf2's negative
    sentinels for missing; floats (DS) come back as NaN already. All three
    fields are non-negative, so masking ``< 0`` covers the int sentinels
    without touching valid data.
    """
    arr = v.format(field)
    if arr is None:
        return np.full(n_samples, np.nan)
    arr = np.asarray(arr, dtype=float).reshape(n_samples, -1)[:, 0]
    arr[arr < 0] = np.nan
    return arr


def _encode_gt(v, gt_types, n_samples):
    """Vectorized equivalent of v5's encode_genotype over all samples.

    f(j,k) = (k*(k+1)/2) + j with j=min allele, k=max allele, computed over
    valid (non-negative) alleles only so haploid calls (e.g. male chrX, stored
    as ``[a, -1]``) encode like v5 rather than reading as missing. Samples whose
    genotype has any missing allele are flagged via gt_types == 3 (UNKNOWN,
    gts012 mode) -- equivalent to v5's ``None in gt`` -- and set to NaN.
    """
    g = v.genotype.array()
    alleles = g[:, :-1]  # drop trailing phase column
    valid = alleles >= 0
    big = np.iinfo(alleles.dtype).max
    amin = np.where(valid, alleles, big).min(axis=1).astype(np.int64)
    amax = np.where(valid, alleles, -1).max(axis=1).astype(np.int64)
    enc = (amax * (amax + 1) // 2 + amin).astype(float)
    enc[gt_types == 3] = np.nan
    return enc


def scan(vcf_file, phenotypes, group_indices, use_dosage, percentiles, chunk_size):
    print(f"scan(vcf_file={vcf_file}, phenotypes={phenotypes}, use_dosage={use_dosage}, percentiles={percentiles}, chunk_size={chunk_size})  ts = {current_milli_time()}")
    vcf = VCF(vcf_file, gts012=True)
    samples = list(vcf.samples)
    n_samples = len(samples)

    # Per-phenotype heatmaps need the real per-group indices even in the >5
    # branch (where group_indices below is collapsed to 'all'), mirroring v5
    # where main() passed the original indices to the pass-2 heatmap builder.
    orig_group_indices = group_indices

    if len(phenotypes) > 5:
        group_labels = ['all']
        group_indices = {'all': np.arange(n_samples, dtype=int)}
    else:
        group_labels = list(phenotypes)

    stats = ['DP', 'GQ'] + (['DS'] if use_dosage else [])
    vstats = _init_vstats(stats, group_labels, percentiles)

    missingness_variants = {g: [] for g in group_labels}
    # Overall per-variant mean across all samples -- always computed (feeds the
    # "Mean {stat} Across Variants" plots 1_/1b_, which v1 emitted unconditionally).
    overall_means = {stat: [] for stat in stats}
    dp_diff = []

    allele_freq = []
    snp_count = 0
    indel_count = 0
    chrom_counts = {}

    sum_dp = np.zeros(n_samples, dtype=float)
    sum_gq = np.zeros(n_samples, dtype=float)
    sum_ds = np.zeros(n_samples, dtype=float) if use_dosage else None
    cnt_dp = np.zeros(n_samples, dtype=int)
    cnt_gq = np.zeros(n_samples, dtype=int)
    cnt_ds = np.zeros(n_samples, dtype=int) if use_dosage else None
    missing_gt_count = np.zeros(n_samples, dtype=int)
    het_count = np.zeros(n_samples, dtype=int)

    total_variants = 0

    case_indices = None
    control_indices = None
    if len(phenotypes) == 2:
        case_indices = group_indices[phenotypes[0]]
        control_indices = group_indices[phenotypes[1]]

    # Heatmaps (only when dosage is available, matching v5's pass-2 condition).
    # Bins are sized from fixed maxima (HEATMAP_STAT_MAX), identical for the
    # overall and per-phenotype heatmaps, and accumulated over chunks of variants
    # within this single pass.
    heatmaps = {}
    buf = {'DP': [], 'GQ': [], 'DS': [], 'GT': []}

    def _bins(s1, s2):
        return _stat_bins_arange(HEATMAP_STAT_MAX.get(s1)), _stat_bins_linspace(HEATMAP_STAT_MAX.get(s2))

    if use_dosage:
        for s1, s2 in HEATMAP_PAIRS:
            xbins, ybins = _bins(s1, s2)
            if xbins is None or ybins is None:
                continue
            heatmaps[(s1, s2, None)] = np.zeros((len(xbins) - 1, len(ybins) - 1), dtype=int)
            for pheno in phenotypes:
                heatmaps[(s1, s2, pheno)] = np.zeros((len(xbins) - 1, len(ybins) - 1), dtype=int)

    def flush_chunk():
        if not buf['DP']:
            return
        mats = {k: np.vstack(buf[k]) for k in ('DP', 'GQ', 'DS', 'GT')}
        for s1, s2 in HEATMAP_PAIRS:
            key = (s1, s2, None)
            if key not in heatmaps:
                continue
            xbins, ybins = _bins(s1, s2)
            heatmaps[key] += _hist2d(mats[s1].ravel(), mats[s2].ravel(), xbins, ybins)
            for pheno in phenotypes:
                idxs = orig_group_indices[pheno]
                if idxs.size == 0:
                    continue
                v1 = mats[s1][:, idxs].ravel()
                v2 = mats[s2][:, idxs].ravel()
                heatmaps[(s1, s2, pheno)] += _hist2d(v1, v2, xbins, ybins)
        for k in buf:
            buf[k].clear()

    for v in vcf:
        total_variants += 1
        if total_variants % 10000 == 0:
            print('.', end='', flush=True)

        alts = v.ALT or ['.']
        is_snp = len(v.REF) == 1 and all(len(alt) == 1 for alt in alts)
        if is_snp:
            snp_count += 1
        else:
            indel_count += 1

        chrom_counts[v.CHROM] = chrom_counts.get(v.CHROM, 0) + 1
        af = v.INFO.get('AF')
        if af is None:
            allele_freq.append(np.nan)
        elif isinstance(af, (tuple, list)):
            allele_freq.append(float(af[0]))
        else:
            allele_freq.append(float(af))

        dp_vals = _format_field(v, 'DP', n_samples)
        gq_vals = _format_field(v, 'GQ', n_samples)
        ds_vals = _format_field(v, 'DS', n_samples) if use_dosage else None

        gt_types = v.gt_types
        gt_vals = _encode_gt(v, gt_types, n_samples)

        # Per-sample aggregates (vectorized).
        valid_dp = ~np.isnan(dp_vals)
        sum_dp += np.where(valid_dp, dp_vals, 0.0)
        cnt_dp += valid_dp
        valid_gq = ~np.isnan(gq_vals)
        sum_gq += np.where(valid_gq, gq_vals, 0.0)
        cnt_gq += valid_gq
        if use_dosage:
            valid_ds = ~np.isnan(ds_vals)
            sum_ds += np.where(valid_ds, ds_vals, 0.0)
            cnt_ds += valid_ds

        missing_gt_count += (gt_types == 3)
        het_count += (gt_types == 1)

        _append_variant_stats(vstats, 'DP', dp_vals, group_labels, group_indices, percentiles)
        _append_variant_stats(vstats, 'GQ', gq_vals, group_labels, group_indices, percentiles)
        if use_dosage:
            _append_variant_stats(vstats, 'DS', ds_vals, group_labels, group_indices, percentiles)

        overall_means['DP'].append(np.nan if np.isnan(dp_vals).all() else np.nanmean(dp_vals))
        overall_means['GQ'].append(np.nan if np.isnan(gq_vals).all() else np.nanmean(gq_vals))
        if use_dosage:
            overall_means['DS'].append(np.nan if np.isnan(ds_vals).all() else np.nanmean(ds_vals))

        _append_missingness(missingness_variants, gt_vals, group_labels, group_indices)

        if case_indices is not None and control_indices is not None:
            if case_indices.size > 0 and control_indices.size > 0:
                case_mean = np.nanmean(dp_vals[case_indices])
                control_mean = np.nanmean(dp_vals[control_indices])
                if np.isnan(case_mean) or np.isnan(control_mean):
                    dp_diff.append(np.nan)
                else:
                    dp_diff.append(abs(case_mean - control_mean))

        if use_dosage:
            buf['DP'].append(dp_vals)
            buf['GQ'].append(gq_vals)
            buf['DS'].append(ds_vals)
            buf['GT'].append(gt_vals)
            if len(buf['DP']) >= chunk_size:
                flush_chunk()

    flush_chunk()
    vcf.close()
    print("")

    mean_dp = np.where(cnt_dp > 0, sum_dp / cnt_dp, np.nan)
    mean_gq = np.where(cnt_gq > 0, sum_gq / cnt_gq, np.nan)
    if use_dosage:
        mean_ds = np.where(cnt_ds > 0, sum_ds / cnt_ds, np.nan)
    else:
        mean_ds = None

    if total_variants > 0:
        missingness_samples = missing_gt_count / float(total_variants)
        het_rates = het_count / float(total_variants)
    else:
        missingness_samples = np.full(n_samples, np.nan, dtype=float)
        het_rates = np.full(n_samples, np.nan, dtype=float)

    return {
        'samples': samples,
        'group_labels': group_labels,
        'group_indices': group_indices,
        'vstats': vstats,
        'overall_means': overall_means,
        'missingness_variants': missingness_variants,
        'dp_diff': np.asarray(dp_diff, dtype=float),
        'allele_freq': np.asarray(allele_freq, dtype=float),
        'variant_types': (snp_count, indel_count),
        'chrom_counts': chrom_counts,
        'mean_dp': mean_dp,
        'mean_gq': mean_gq,
        'mean_ds': mean_ds,
        'missingness_samples': missingness_samples,
        'het_rates': het_rates,
        'heatmaps': heatmaps,
        'phenotypes': phenotypes,
        'use_dosage': use_dosage,
        'percentiles': percentiles,
    }


# Plotting functions

def plot_variant_stats(vstats, overall_mean, phenotypes, stat, percentiles, stat_label='Statistic', include_log_scale=True):
    print(f"plot_variant_stats(stat={stat})  ts = {current_milli_time()}")

    # Overall mean across variants (all samples) -- always produced, regardless of the
    # phenotype-count branch, matching v1's unconditional plots 1_/1b_.
    mean_values = np.asarray(overall_mean, dtype=float)
    mean_values = mean_values[~np.isnan(mean_values)]
    plt.figure(figsize=(9, 6))
    sns.histplot(mean_values, bins=50)
    plt.title(f'Mean {stat} Across Variants')
    plt.xlabel(stat_label)
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/1_{stat}_mean_variants.png')
    plt.savefig(f'{output_dir}/1_{stat}_mean_variants.svg', format="svg")
    plt.close()

    if include_log_scale:
        plt.figure(figsize=(9, 6))
        sns.histplot(mean_values, bins=50, log_scale=True)
        plt.title(f'Mean {stat} Across Variants - log scale')
        plt.xlabel(f'{stat_label} (log)')
        plt.ylabel('Variant Count')
        plt.savefig(f'{output_dir}/1b_{stat}_mean_variants.png')
        plt.savefig(f'{output_dir}/1b_{stat}_mean_variants.svg', format="svg")
        plt.close()

    if len(phenotypes) > 5:
        for perc in percentiles:
            perc_values = np.asarray(vstats[stat]['all'][perc], dtype=float)
            plt.figure(figsize=(9, 6))
            sns.histplot(perc_values[~np.isnan(perc_values)], bins=50)
            plt.title(f'{stat} ({perc}th Percentile) Across Variants')
            plt.xlabel(stat_label)
            plt.ylabel('Variant Count')
            plt.savefig(f'{output_dir}/2_{stat}_percentile_{perc}_variants.png')
            plt.savefig(f'{output_dir}/2_{stat}_percentile_{perc}_variants.svg', format="svg")
            plt.close()

            if include_log_scale:
                plt.figure(figsize=(9, 6))
                sns.histplot(perc_values[~np.isnan(perc_values)], bins=50, log_scale=True)
                plt.title(f'{stat} ({perc}th Percentile) Across Variants - log scale')
                plt.xlabel(f'{stat_label} (log)')
                plt.ylabel('Variant Count')
                plt.savefig(f'{output_dir}/2b_{stat}_percentile_{perc}_variants.png')
                plt.savefig(f'{output_dir}/2b_{stat}_percentile_{perc}_variants.svg', format="svg")
                plt.close()
    else:
        for perc in percentiles + ['mean']:
            key = 'mean' if perc == 'mean' else perc
            plt.figure(figsize=(10, 6))
            for pheno in phenotypes:
                pheno_values = np.asarray(vstats[stat][pheno][key], dtype=float)
                pheno_values = pheno_values[~np.isnan(pheno_values)]
                emit_stat(f'variant_stats_{stat}_p{perc}_pheno{pheno}', pheno_values)
                sns.histplot(pheno_values, bins=50, label=f'{pheno} ({perc})', alpha=0.5)
            plt.title(f'Mean {stat} by Phenotype' if perc == 'mean' else f'{stat} ({perc}th Percentile) by Phenotype')
            plt.xlabel(stat_label)
            plt.ylabel('Variant Count')
            plt.legend()
            plt.savefig(f'{output_dir}/3_{stat}_percentile_{perc}_by_phenotype_variants.png')
            plt.savefig(f'{output_dir}/3_{stat}_percentile_{perc}_by_phenotype_variants.svg', format="svg")
            plt.close()

            if include_log_scale:
                plt.figure(figsize=(10, 6))
                for pheno in phenotypes:
                    pheno_values = np.asarray(vstats[stat][pheno][key], dtype=float)
                    pheno_values = pheno_values[~np.isnan(pheno_values)]
                    sns.histplot(pheno_values, bins=50, label=f'{pheno} ({perc})', alpha=0.5, log_scale=True)
                plt.title(f'Mean {stat} by Phenotype - log scale' if perc == 'mean' else f'{stat} ({perc}th Percentile) by Phenotype - log scale')
                plt.xlabel(f'{stat_label} (log)')
                plt.ylabel('Variant Count')
                plt.legend()
                plt.savefig(f'{output_dir}/3b_{stat}_percentile_{perc}_by_phenotype_variants.png')
                plt.savefig(f'{output_dir}/3b_{stat}_percentile_{perc}_by_phenotype_variants.svg', format="svg")
                plt.close()


def plot_sample_stats(sample_means, samples, pheno_df, phenotypes, stat, stat_label='Statistic'):
    print(f"plot_sample_stats(stat={stat})  ts = {current_milli_time()}")
    sample_stats = pl.DataFrame({'IID': samples, 'value': sample_means}).join(
        pheno_df.select(['IID', 'Y1']), on='IID', how='left')

    if len(phenotypes) > 5:
        plt.figure(figsize=(9, 6))
        sns.histplot(sample_stats['value'].drop_nulls().to_numpy(), bins=50)
        plt.title(f'Mean {stat} Across Samples')
        plt.xlabel(stat_label)
        plt.ylabel('Sample Count')
        plt.savefig(f'{output_dir}/4_{stat}_samples.png')
        plt.savefig(f'{output_dir}/4_{stat}_samples.svg', format="svg")
        plt.close()
    else:
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_data = sample_stats.filter(pl.col('Y1') == pheno)['value'].drop_nulls().to_numpy()
            if len(pheno_data) > 0:
                emit_stat(f'sample_stats_{stat}_pheno{pheno}', pheno_data)
                sns.histplot(pheno_data, bins=50, label=pheno, alpha=0.5)
        plt.title(f'Mean {stat} by Phenotype (Samples)')
        plt.xlabel(stat_label)
        plt.ylabel('Sample Count')
        plt.legend()
        plt.savefig(f'{output_dir}/5_{stat}_by_phenotype_samples.png')
        plt.savefig(f'{output_dir}/5_{stat}_by_phenotype_samples.svg', format="svg")
        plt.close()


def plot_missingness(missingness_variants, missingness_samples, samples, pheno_df, phenotypes):
    print(f"plot_missingness()  ts = {current_milli_time()}")

    if len(phenotypes) > 5:
        mv = np.asarray(missingness_variants['all'], dtype=float)
        plt.figure(figsize=(9, 6))
        sns.histplot(mv, bins=50)
        plt.title('Missingness Rate Across Variants')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Variant Count')
        plt.savefig(f'{output_dir}/6_missingness_variants.png')
        plt.savefig(f'{output_dir}/6_missingness_variants.svg', format="svg")
        plt.close()
    else:
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_missing = np.asarray(missingness_variants[pheno], dtype=float)
            emit_stat(f'missingness_variants_pheno{pheno}', pheno_missing)
            sns.histplot(pheno_missing, bins=50, label=pheno, alpha=0.5)
        plt.title('Missingness Rate by Phenotype (Variants)')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Variant Count')
        plt.legend()
        plt.savefig(f'{output_dir}/7_missingness_by_phenotype_variants.png')
        plt.savefig(f'{output_dir}/7_missingness_by_phenotype_variants.svg', format="svg")
        plt.close()

    missingness_samples_df = pl.DataFrame({'IID': samples, 'missing_rate': missingness_samples}).join(
        pheno_df.select(['IID', 'Y1']), on='IID', how='left')

    if len(phenotypes) > 5:
        plt.figure(figsize=(9, 6))
        sns.histplot(missingness_samples_df['missing_rate'].drop_nulls().to_numpy(), bins=50)
        plt.title('Missingness Rate Across Samples')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Sample Count')
        plt.savefig(f'{output_dir}/8_missingness_samples.png')
        plt.savefig(f'{output_dir}/8_missingness_samples.svg', format="svg")
        plt.close()
    else:
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_data = missingness_samples_df.filter(pl.col('Y1') == pheno)['missing_rate'].drop_nulls().to_numpy()
            if len(pheno_data) > 0:
                emit_stat(f'missingness_samples_pheno{pheno}', pheno_data)
                sns.histplot(pheno_data, bins=50, label=pheno, alpha=0.5)
        plt.title('Missingness Rate by Phenotype (Samples)')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Sample Count')
        plt.legend()
        plt.savefig(f'{output_dir}/9_missingness_by_phenotype_samples.png')
        plt.savefig(f'{output_dir}/9_missingness_by_phenotype_samples.svg', format="svg")
        plt.close()


def plot_dp_differences(dp_diff):
    print(f"plot_dp_differences()  ts = {current_milli_time()}")
    if dp_diff.size == 0:
        return
    dp_diff = dp_diff[~np.isnan(dp_diff)]
    if dp_diff.size == 0:
        return
    emit_stat('dp_diff', dp_diff)
    plt.figure(figsize=(9, 6))
    sns.histplot(dp_diff, bins=50)
    plt.title('Absolute DP Differences (Cases vs Controls)')
    plt.xlabel('Absolute DP Difference')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/10_dp_abs_diff_cases_controls.png')
    plt.savefig(f'{output_dir}/10_dp_abs_diff_cases_controls.svg', format="svg")
    plt.close()

    plt.figure(figsize=(9, 6))
    sns.histplot(dp_diff, bins=50, log_scale=True)
    plt.title('Absolute DP Differences (Cases vs Controls) - log scale')
    plt.xlabel('Absolute DP Difference (log)')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/10b_dp_abs_diff_cases_controls.png')
    plt.savefig(f'{output_dir}/10b_dp_abs_diff_cases_controls.svg', format="svg")
    plt.close()


def plot_allele_frequency(allele_freq):
    print(f"plot_allele_frequency()  ts = {current_milli_time()}")
    af = allele_freq[~np.isnan(allele_freq)]
    emit_stat('allele_freq', af)
    plt.figure(figsize=(9, 6))
    sns.histplot(af, bins=100, log_scale=True)
    plt.title('Alternate Allele Frequency Distribution - log scale')
    plt.xlabel('Allele Frequency')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/11_allele_frequency.png')
    plt.savefig(f'{output_dir}/11_allele_frequency.svg', format="svg")
    plt.close()


def plot_variant_types(snp_count, indel_count):
    print(f"plot_variant_types()  ts = {current_milli_time()}")
    emit_stat('variant_types', np.array([snp_count, indel_count]))
    df = pl.DataFrame({'Variant_Type': [0, 1], 'count': [snp_count, indel_count]}).to_pandas()
    plt.figure(figsize=(9, 6))
    sns.barplot(data=df, x='Variant_Type', y='count')
    plt.title('Variant Type Distribution')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/12_variant_types.png')
    plt.savefig(f'{output_dir}/12_variant_types.svg', format="svg")
    plt.close()


def plot_chrom_density(chrom_counts):
    print(f"plot_chrom_density()  ts = {current_milli_time()}")
    chroms = sorted(chrom_counts.keys())
    counts = [chrom_counts[c] for c in chroms]
    chrom_df = pl.DataFrame({'CHROM': chroms, 'len': counts}).to_pandas()
    chrom_df[['CHROM', 'len']].to_csv(str(Path(stats_dir) / 'chrom_density.csv'), index=False)
    plt.figure(figsize=(12, 6))
    sns.barplot(x='CHROM', y='len', data=chrom_df)
    plt.title('Variant Count by Chromosome')
    plt.xlabel('Chromosome')
    plt.ylabel('Variant Count')
    plt.xticks(rotation=45)
    plt.savefig(f'{output_dir}/13_chrom_variant_density.png')
    plt.savefig(f'{output_dir}/13_chrom_variant_density.svg', format="svg")
    plt.close()


def plot_heterozygosity(het_rates, samples, pheno_df, phenotypes):
    print(f"plot_heterozygosity()  ts = {current_milli_time()}")
    het_df = pl.DataFrame({'IID': samples, 'Heterozygosity': het_rates}).join(
        pheno_df.select(['IID', 'Y1']), on='IID', how='left')

    if len(phenotypes) > 5:
        plt.figure(figsize=(9, 6))
        sns.histplot(het_df['Heterozygosity'].to_numpy(), bins=50)
        plt.title('Heterozygosity Rate Across Samples')
        plt.xlabel('Heterozygosity Rate')
        plt.ylabel('Sample Count')
        plt.savefig(f'{output_dir}/14_heterozygosity_samples.png')
        plt.savefig(f'{output_dir}/14_heterozygosity_samples.svg', format="svg")
        plt.close()
    else:
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_data = het_df.filter(pl.col('Y1') == pheno)['Heterozygosity'].to_numpy()
            if len(pheno_data) > 0:
                emit_stat(f'het_pheno{pheno}', pheno_data)
                sns.histplot(pheno_data, bins=50, label=pheno, alpha=0.5)
        plt.title('Heterozygosity Rate by Phenotype')
        plt.xlabel('Heterozygosity Rate')
        plt.ylabel('Sample Count')
        plt.legend()
        plt.savefig(f'{output_dir}/15_heterozygosity_by_phenotype_samples.png')
        plt.savefig(f'{output_dir}/15_heterozygosity_by_phenotype_samples.svg', format="svg")
        plt.close()


def plot_boxplots(sample_means, samples, pheno_df, stat, stat_label):
    print(f"plot_boxplots(stat={stat})  ts = {current_milli_time()}")
    sample_stats = pl.DataFrame({'IID': samples, 'value': sample_means}).join(
        pheno_df.select(['IID', 'Y1']), on='IID', how='left').to_pandas()

    for _pheno in sorted(sample_stats['Y1'].dropna().unique()):
        emit_stat(f'boxplot_{stat}_pheno{_pheno}', sample_stats[sample_stats['Y1'] == _pheno]['value'].to_numpy())
    plt.figure(figsize=(10, 6))
    try:
        sns.boxplot(data=sample_stats, x='Y1', y='value')
    except ValueError as e:
        print(f"ValueError {e}   sample_stats:\n{sample_stats}")
    plt.title(f'{stat_label} by Phenotype (Samples)')
    plt.xlabel('Phenotype')
    plt.ylabel(stat_label)
    plt.savefig(f'{output_dir}/16_{stat}_boxplot_by_phenotype.png')
    plt.savefig(f'{output_dir}/16_{stat}_boxplot_by_phenotype.svg', format="svg")
    plt.close()


def plot_stat_vs_stat(heatmaps, phenotypes, stat1, stat2):
    print(f"plot_stat_vs_stat()  ts = {current_milli_time()}")
    key = (stat1, stat2, None)
    if key in heatmaps:
        heatmap = heatmaps[key]
        emit_stat(f'heatmap_{stat2}_vs_{stat1}', heatmap)
        plt.figure(figsize=(12, 10))
        sns.heatmap(heatmap.T, cmap='viridis', norm='log', cbar_kws={'label': 'Count (log scale)'})
        plt.xlabel(f'{stat1}')
        plt.ylabel(f'{stat2}')
        plt.title(f'Heatmap of {stat2} distribution per {stat1}')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/17_{stat2}_vs_{stat1}.png')
        plt.savefig(f'{output_dir}/17_{stat2}_vs_{stat1}.svg', format="svg")
        plt.close()

    if len(phenotypes) <= 5:
        for pheno in phenotypes:
            key_p = (stat1, stat2, pheno)
            if key_p not in heatmaps:
                continue
            heatmap = heatmaps[key_p]
            emit_stat(f'heatmap_{stat2}_vs_{stat1}_pheno{pheno}', heatmap)
            plt.figure(figsize=(12, 10))
            sns.heatmap(heatmap.T, cmap='viridis', norm='log', cbar_kws={'label': 'Count (log scale)'})
            plt.xlabel(f'{stat1}')
            plt.ylabel(f'{stat2}')
            plt.title(f'Heatmap of {stat2} distribution per {stat1} by Phenotype {pheno}')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/17_{stat2}_vs_{stat1}_by_phenotype_{pheno}.png')
            plt.savefig(f'{output_dir}/17_{stat2}_vs_{stat1}_by_phenotype_{pheno}.svg', format="svg")
            plt.close()


# Main execution

def main(vcf_file, phenotype_file, percentiles, use_dosage, chunk_size):
    print(f"main()  percentiles = {percentiles}  use_dosage = {use_dosage}  ts = {current_milli_time()}", flush=True)

    pheno_df = load_phenotype(phenotype_file)
    vcf = VCF(vcf_file)
    samples = list(vcf.samples)
    vcf.close()
    phenotypes, group_indices = build_groups(samples, pheno_df)

    res = scan(vcf_file, phenotypes, group_indices, use_dosage, percentiles, chunk_size)

    plot_variant_stats(res['vstats'], res['overall_means']['DP'], phenotypes, 'DP', percentiles, stat_label='Depth of Coverage')
    plot_sample_stats(res['mean_dp'], samples, pheno_df, phenotypes, 'DP', stat_label='Depth of Coverage')

    plot_variant_stats(res['vstats'], res['overall_means']['GQ'], phenotypes, 'GQ', percentiles, stat_label='Genotype Quality')
    plot_sample_stats(res['mean_gq'], samples, pheno_df, phenotypes, 'GQ', stat_label='Genotype Quality')

    if use_dosage:
        plot_variant_stats(res['vstats'], res['overall_means']['DS'], phenotypes, 'DS', percentiles, stat_label='Genotype Dosage', include_log_scale=False)
        plot_sample_stats(res['mean_ds'], samples, pheno_df, phenotypes, 'DS', stat_label='Genotype Dosage')

    plot_missingness(res['missingness_variants'], res['missingness_samples'], samples, pheno_df, phenotypes)
    plot_dp_differences(res['dp_diff'])
    plot_allele_frequency(res['allele_freq'])
    plot_variant_types(res['variant_types'][0], res['variant_types'][1])
    plot_chrom_density(res['chrom_counts'])
    plot_heterozygosity(res['het_rates'], samples, pheno_df, phenotypes)

    plot_boxplots(res['mean_dp'], samples, pheno_df, 'DP', 'Depth of Coverage')
    plot_boxplots(res['mean_gq'], samples, pheno_df, 'GQ', 'Genotype Quality')
    if use_dosage:
        plot_boxplots(res['mean_ds'], samples, pheno_df, 'DS', 'Genotype Dosage')

        heatmaps = res['heatmaps']
        plot_stat_vs_stat(heatmaps, phenotypes, stat1='GQ', stat2='DS')
        plot_stat_vs_stat(heatmaps, phenotypes, stat1='DP', stat2='DS')
        plot_stat_vs_stat(heatmaps, phenotypes, stat1='GQ', stat2='DP')
        plot_stat_vs_stat(heatmaps, phenotypes, stat1='GT', stat2='DS')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--use-dosage", required=True)
    parser.add_argument("--process-name", required=True)
    parser.add_argument("--chunk-size", type=int, default=500)
    args = parser.parse_args()

    percentiles = [1, 50]
    use_dosage = str(args.use_dosage).lower() == "true"
    main(args.vcf, args.phenotype, percentiles, use_dosage, args.chunk_size)

    # Write versions.yml
    with open('versions.yml', 'w') as f:
        f.write(f"{args.process_name}:\n")
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write("    d3js: v7\n")
