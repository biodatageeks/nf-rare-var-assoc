#!/usr/bin/env python3
"""Exploratory data analysis for a multisample VCF (v5).

Two-pass, chunked pysam scan:
- Pass 1: per-variant stats + per-sample aggregates, track heatmap maxima.
- Pass 2 (use_dosage=true): build exact heatmaps with np.histogram2d.
"""
import os
from pathlib import Path
import sys
import argparse
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
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


HOMOZYGOUS_VALUES = {0, 2, 5, 9}


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


# Function to encode genotype as integer using f(j,k) = (k*(k+1)/2) + j
def encode_genotype(gt):
    if gt is None or None in gt:
        return None  # Missing genotype
    alleles = sorted(gt)
    j = alleles[0]
    k = alleles[-1]  # For haploid (e.g. male chrX), k == j; for diploid, k is second allele
    return (k * (k + 1) // 2) + j  # Integer encoding


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


def _update_max(max_by_group, values, group_labels, group_indices, stat):
    for group in group_labels:
        if group == 'all':
            vals = values
        else:
            vals = values[group_indices[group]]
        if vals.size == 0:
            continue
        if np.isnan(vals).all():
            continue
        m = float(np.nanmax(vals))
        prev = max_by_group[group].get(stat)
        if prev is None or m > prev:
            max_by_group[group][stat] = m


def scan_pass1(vcf_file, phenotypes, group_indices, use_dosage, percentiles, chunk_size):
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    n_samples = len(samples)

    if len(phenotypes) > 5:
        group_labels = ['all']
        group_indices = {'all': np.arange(n_samples, dtype=int)}
    else:
        group_labels = list(phenotypes)

    max_group_labels = ['all'] + (list(phenotypes) if len(phenotypes) <= 5 else [])
    max_group_indices = {'all': np.arange(n_samples, dtype=int)}
    for pheno in phenotypes:
        if pheno in group_indices:
            max_group_indices[pheno] = group_indices[pheno]

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

    max_by_group = {g: {} for g in max_group_labels}

    case_indices = None
    control_indices = None
    if len(phenotypes) == 2:
        case_indices = group_indices[phenotypes[0]]
        control_indices = group_indices[phenotypes[1]]

    for record in vcf.fetch():
        total_variants += 1
        if total_variants % 100000 == 0:
            print('.', end='', flush=True)

        is_snp = all(len(record.ref) == 1 and len(alt) == 1 for alt in (record.alts or ['.']))
        if is_snp:
            snp_count += 1
        else:
            indel_count += 1

        chrom_counts[record.chrom] = chrom_counts.get(record.chrom, 0) + 1
        if 'AF' in record.info:
            allele_freq.append(float(record.info['AF'][0]))
        else:
            allele_freq.append(np.nan)

        dp_vals = np.empty(n_samples, dtype=float)
        gq_vals = np.empty(n_samples, dtype=float)
        ds_vals = np.empty(n_samples, dtype=float) if use_dosage else None
        gt_vals = np.empty(n_samples, dtype=float)

        for i, sample in enumerate(samples):
            sample_data = record.samples[sample]

            dp = sample_data.get('DP', np.nan)
            if dp is None:
                dp = np.nan
            dp_vals[i] = dp
            if not np.isnan(dp):
                sum_dp[i] += dp
                cnt_dp[i] += 1

            gq = sample_data.get('GQ', np.nan)
            if gq is None:
                gq = np.nan
            gq_vals[i] = gq
            if not np.isnan(gq):
                sum_gq[i] += gq
                cnt_gq[i] += 1

            if use_dosage:
                ds = sample_data.get('DS', np.nan)
                if ds is None:
                    ds = np.nan
                ds_vals[i] = ds
                if not np.isnan(ds):
                    sum_ds[i] += ds
                    cnt_ds[i] += 1

            gt = sample_data.get('GT')
            gt_enc = encode_genotype(gt)
            if gt_enc is None:
                gt_vals[i] = np.nan
                missing_gt_count[i] += 1
            else:
                gt_vals[i] = float(gt_enc)
                if gt_enc not in HOMOZYGOUS_VALUES:
                    het_count[i] += 1

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
            _update_max(max_by_group, dp_vals, max_group_labels, max_group_indices, 'DP')
            _update_max(max_by_group, gq_vals, max_group_labels, max_group_indices, 'GQ')
            _update_max(max_by_group, gt_vals, max_group_labels, max_group_indices, 'GT')
            _update_max(max_by_group, ds_vals, max_group_labels, max_group_indices, 'DS')

    vcf.close()

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
        'max_by_group': max_by_group,
        'phenotypes': phenotypes,
        'use_dosage': use_dosage,
        'percentiles': percentiles,
    }


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


def scan_pass2_heatmaps(vcf_file, samples, phenotypes, group_indices, max_by_group, chunk_size):
    n_samples = len(samples)

    heatmaps = {}
    overall_pairs = [('GQ', 'DS'), ('DP', 'DS'), ('GQ', 'DP'), ('GT', 'DS')]

    overall_max = max_by_group.get('all', {})
    for s1, s2 in overall_pairs:
        xbins = _stat_bins_arange(overall_max.get(s1))
        ybins = np.linspace(0, 2, 201)
        if xbins is not None:
            heatmaps[(s1, s2, None)] = np.zeros((len(xbins) - 1, len(ybins) - 1), dtype=int)

    for pheno in phenotypes:
        maxes = max_by_group.get(pheno, {})
        for s1, s2 in overall_pairs:
            xbins = _stat_bins_arange(maxes.get(s1))
            ybins = _stat_bins_linspace(maxes.get(s2))
            if xbins is not None and ybins is not None:
                heatmaps[(s1, s2, pheno)] = np.zeros((len(xbins) - 1, len(ybins) - 1), dtype=int)

    vcf = pysam.VariantFile(vcf_file)
    chunk_dp = []
    chunk_gq = []
    chunk_ds = []
    chunk_gt = []

    def flush_chunk():
        if not chunk_dp:
            return
        dp_mat = np.vstack(chunk_dp)
        gq_mat = np.vstack(chunk_gq)
        ds_mat = np.vstack(chunk_ds)
        gt_mat = np.vstack(chunk_gt)

        stat_map = {
            'DP': dp_mat,
            'GQ': gq_mat,
            'DS': ds_mat,
            'GT': gt_mat,
        }

        for s1, s2 in overall_pairs:
            key = (s1, s2, None)
            if key in heatmaps:
                xbins = _stat_bins_arange(overall_max.get(s1))
                ybins = np.linspace(0, 2, 201)
                v1 = stat_map[s1].ravel()
                v2 = stat_map[s2].ravel()
                heatmaps[key] += _hist2d(v1, v2, xbins, ybins)

            for pheno in phenotypes:
                key_p = (s1, s2, pheno)
                if key_p in heatmaps:
                    idxs = group_indices[pheno]
                    if idxs.size == 0:
                        continue
                    maxes = max_by_group.get(pheno, {})
                    xbins = _stat_bins_arange(maxes.get(s1))
                    ybins = _stat_bins_linspace(maxes.get(s2))
                    v1 = stat_map[s1][:, idxs].ravel()
                    v2 = stat_map[s2][:, idxs].ravel()
                    heatmaps[key_p] += _hist2d(v1, v2, xbins, ybins)

        chunk_dp.clear()
        chunk_gq.clear()
        chunk_ds.clear()
        chunk_gt.clear()

    for record in vcf.fetch():
        dp_vals = np.empty(n_samples, dtype=float)
        gq_vals = np.empty(n_samples, dtype=float)
        ds_vals = np.empty(n_samples, dtype=float)
        gt_vals = np.empty(n_samples, dtype=float)

        for i, sample in enumerate(samples):
            sample_data = record.samples[sample]

            dp = sample_data.get('DP', np.nan)
            if dp is None:
                dp = np.nan
            dp_vals[i] = dp

            gq = sample_data.get('GQ', np.nan)
            if gq is None:
                gq = np.nan
            gq_vals[i] = gq

            ds = sample_data.get('DS', np.nan)
            if ds is None:
                ds = np.nan
            ds_vals[i] = ds

            gt = sample_data.get('GT')
            gt_enc = encode_genotype(gt)
            if gt_enc is None:
                gt_vals[i] = np.nan
            else:
                gt_vals[i] = float(gt_enc)

        chunk_dp.append(dp_vals)
        chunk_gq.append(gq_vals)
        chunk_ds.append(ds_vals)
        chunk_gt.append(gt_vals)

        if len(chunk_dp) >= chunk_size:
            flush_chunk()

    flush_chunk()
    vcf.close()

    return heatmaps


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
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    vcf.close()
    phenotypes, group_indices = build_groups(samples, pheno_df)

    pass1 = scan_pass1(vcf_file, phenotypes, group_indices, use_dosage, percentiles, chunk_size)

    plot_variant_stats(pass1['vstats'], pass1['overall_means']['DP'], phenotypes, 'DP', percentiles, stat_label='Depth of Coverage')
    plot_sample_stats(pass1['mean_dp'], samples, pheno_df, phenotypes, 'DP', stat_label='Depth of Coverage')

    plot_variant_stats(pass1['vstats'], pass1['overall_means']['GQ'], phenotypes, 'GQ', percentiles, stat_label='Genotype Quality')
    plot_sample_stats(pass1['mean_gq'], samples, pheno_df, phenotypes, 'GQ', stat_label='Genotype Quality')

    if use_dosage:
        plot_variant_stats(pass1['vstats'], pass1['overall_means']['DS'], phenotypes, 'DS', percentiles, stat_label='Genotype Dosage', include_log_scale=False)
        plot_sample_stats(pass1['mean_ds'], samples, pheno_df, phenotypes, 'DS', stat_label='Genotype Dosage')

    plot_missingness(pass1['missingness_variants'], pass1['missingness_samples'], samples, pheno_df, phenotypes)
    plot_dp_differences(pass1['dp_diff'])
    plot_allele_frequency(pass1['allele_freq'])
    plot_variant_types(pass1['variant_types'][0], pass1['variant_types'][1])
    plot_chrom_density(pass1['chrom_counts'])
    plot_heterozygosity(pass1['het_rates'], samples, pheno_df, phenotypes)

    plot_boxplots(pass1['mean_dp'], samples, pheno_df, 'DP', 'Depth of Coverage')
    plot_boxplots(pass1['mean_gq'], samples, pheno_df, 'GQ', 'Genotype Quality')
    if use_dosage:
        plot_boxplots(pass1['mean_ds'], samples, pheno_df, 'DS', 'Genotype Dosage')

        heatmaps = scan_pass2_heatmaps(vcf_file, samples, phenotypes, group_indices, pass1['max_by_group'], chunk_size)
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
