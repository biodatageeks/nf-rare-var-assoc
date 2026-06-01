#!/usr/bin/env python3
"""Exploratory data analysis for a multisample VCF.

Memory-bounded rewrite (T18b): instead of materializing a variants x samples
matrix in Python, every per-variant / per-sample / per-bin aggregate is computed
as a streaming DataFusion GROUP BY via polars-bio.  The VCF genotypes struct is
unnested into a long (variant_idx, sample_idx, GT, DP, GQ, DS) stream with a
CROSS JOIN against generate_series; aggregates flow through DataFusion's default
streaming hash-aggregate (state bounded by group count), so the full matrix is
never buffered.  Plotting and the emitted-stats schema are unchanged from the
previous implementation.
"""
import argparse
import os
import sys
import time
import csv as _csv
from pathlib import Path

tmp_matplotlib_dir = Path("matplotlib_tmp")
tmp_matplotlib_dir.mkdir(exist_ok=True)
os.environ['MPLCONFIGDIR'] = str(tmp_matplotlib_dir.resolve())

import numpy as np
import polars as pl
import polars_bio as pb
import matplotlib.pyplot as plt
import seaborn as sns
import pysam

plt.rcParams.update({'font.size': 20})

output_dir = "plots"
Path(output_dir).mkdir(exist_ok=True)
stats_dir = "eda_stats"
Path(stats_dir).mkdir(exist_ok=True)


def current_milli_time():
    return round(time.time() * 1000)


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


# ---------------------------------------------------------------------------
# SQL building blocks
# ---------------------------------------------------------------------------

# GT arrives as strings: '0/0', '0|1', '0', '1', './.', '.'.  Encode with the
# diploid index f(j,k) = k*(k+1)/2 + j (j<=k); haploid 'a' encodes as f(a,a).
# Missing (any '.' allele) -> NULL.  Het = the two alleles differ (haploid -> 0).
def _gt_split_sql(col, sep):
    a = f"CAST(split_part({col}, '{sep}', 1) AS BIGINT)"
    b = f"CAST(split_part({col}, '{sep}', 2) AS BIGINT)"
    hi = f"GREATEST({a}, {b})"
    lo = f"LEAST({a}, {b})"
    return f"({hi} * ({hi} + 1) / 2 + {lo})"


def _gt_encode_sql(col):
    hap = f"(CAST({col} AS BIGINT) * (CAST({col} AS BIGINT) + 3) / 2)"
    return f"""CASE
        WHEN {col} IS NULL OR {col} LIKE '.%' OR {col} LIKE '%/.' OR {col} LIKE '%|.' THEN NULL
        WHEN {col} LIKE '%/%' THEN {_gt_split_sql(col, '/')}
        WHEN {col} LIKE '%|%' THEN {_gt_split_sql(col, '|')}
        ELSE {hap}
    END"""


def _gt_is_het_sql(col):
    het_slash = f"CASE WHEN split_part({col}, '/', 1) != split_part({col}, '/', 2) THEN 1.0 ELSE 0.0 END"
    het_pipe = f"CASE WHEN split_part({col}, '|', 1) != split_part({col}, '|', 2) THEN 1.0 ELSE 0.0 END"
    return f"""CASE
        WHEN {col} IS NULL OR {col} LIKE '.%' OR {col} LIKE '%/.' OR {col} LIKE '%|.' THEN 0.0
        WHEN {col} LIKE '%/%' THEN {het_slash}
        WHEN {col} LIKE '%|%' THEN {het_pipe}
        ELSE 0.0
    END"""


def _idx_list_sql(indices):
    return ','.join(str(int(i)) for i in indices)


# UNNEST the genotypes struct into a long per-(variant, sample) stream.  The
# CROSS JOIN against generate_series gives a deterministic sample_idx aligned
# with the VCF header sample order.  This is the single source the GROUP BY
# queries reduce over; it is never collected to Python.
def _unnested_cte(n_samples):
    return f"""
indexed AS (
    SELECT
        ROW_NUMBER() OVER () AS variant_idx,
        chrom,
        start + 1 AS pos,
        ref,
        alt,
        "AF"[1] AS af,
        genotypes."GT" AS gt_arr,
        genotypes."DP" AS dp_arr,
        genotypes."GQ" AS gq_arr,
        genotypes."DS" AS ds_arr
    FROM vcf_table
),
unnested AS (
    SELECT
        i.variant_idx,
        i.chrom,
        CAST(s.value AS BIGINT) AS sample_idx,
        i.gt_arr[CAST(s.value AS BIGINT) + 1] AS gt,
        i.dp_arr[CAST(s.value AS BIGINT) + 1] AS dp,
        i.gq_arr[CAST(s.value AS BIGINT) + 1] AS gq,
        i.ds_arr[CAST(s.value AS BIGINT) + 1] AS ds
    FROM indexed i
    CROSS JOIN generate_series(0, {n_samples - 1}) AS s
)
"""


def _filter_clause(indices):
    """FILTER predicate restricting an aggregate to a phenotype group, or '' for all samples."""
    if indices is None:
        return ""
    return f"FILTER (WHERE sample_idx IN ({_idx_list_sql(indices)}))"


def _null_guard_sql(col, indices, perc):
    """Reproduce np.quantile NaN semantics: percentile is NULL if any group value is NULL."""
    if indices is None:
        any_null = f"SUM(CASE WHEN {col} IS NULL THEN 1 ELSE 0 END)"
        filt = ""
    else:
        idx_sql = _idx_list_sql(indices)
        any_null = f"SUM(CASE WHEN sample_idx IN ({idx_sql}) AND {col} IS NULL THEN 1 ELSE 0 END)"
        filt = f"FILTER (WHERE sample_idx IN ({idx_sql}))"
    return (f"CASE WHEN {any_null} > 0 THEN NULL "
            f"ELSE approx_percentile_cont({col}, {perc / 100.0}) {filt} END")


# ---------------------------------------------------------------------------
# Data layer: streaming GROUP BY reductions (no matrix materialized)
# ---------------------------------------------------------------------------

def load_phenotype(phenotype_file):
    print(f"load_phenotype({phenotype_file})  ts = {current_milli_time()}", flush=True)
    return pl.read_csv(phenotype_file, separator='\t')


def read_samples(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    vcf.close()
    return samples


def build_groups(samples, pheno_df):
    """Map each phenotype value to the list of sample_idx (VCF header order) that carry it."""
    idx_df = pl.DataFrame({'IID': samples, 'sample_idx': list(range(len(samples)))})
    joined = idx_df.join(pheno_df.select(['IID', 'Y1']), on='IID', how='left')
    phenotypes = pheno_df['Y1'].unique().to_list()
    groups = {}
    for p in phenotypes:
        groups[p] = joined.filter(pl.col('Y1') == p)['sample_idx'].to_list()
    return groups, phenotypes


def compute_variant_stats(stat, group_indices, percentiles):
    """Per-variant AVG and percentile arrays, per group label.

    Returns {label: {'mean': np.ndarray, <perc>: np.ndarray}}, each array aligned
    to variant order (1..N).  NaN marks NULL (AVG over no values / guarded percentile).
    """
    print(f"compute_variant_stats(stat={stat})  ts = {current_milli_time()}", flush=True)
    col = stat.lower()
    selects = []
    for label, idxs in group_indices.items():
        filt = _filter_clause(idxs)
        selects.append(f"AVG({col}) {filt} AS mean_{label}")
        for p in percentiles:
            selects.append(f"{_null_guard_sql(col, idxs, p)} AS p{p}_{label}")
    sql = (f"WITH {_unnested_cte(N_SAMPLES)} "
           f"SELECT variant_idx, {', '.join(selects)} "
           f"FROM unnested GROUP BY variant_idx ORDER BY variant_idx")
    df = pb.sql(sql).collect()
    result = {}
    for label in group_indices:
        result[label] = {'mean': df[f'mean_{label}'].to_numpy().astype(float)}
        for p in percentiles:
            result[label][p] = df[f'p{p}_{label}'].to_numpy().astype(float)
    return result


def compute_sample_means(stat):
    """Per-sample mean of a stat across all variants, aligned to sample_idx order."""
    print(f"compute_sample_means(stat={stat})  ts = {current_milli_time()}", flush=True)
    col = stat.lower()
    sql = (f"WITH {_unnested_cte(N_SAMPLES)} "
           f"SELECT sample_idx, AVG({col}) AS m "
           f"FROM unnested GROUP BY sample_idx ORDER BY sample_idx")
    df = pb.sql(sql).collect()
    return df['m'].to_numpy().astype(float)


def compute_missingness_variants(group_indices):
    """Per-variant missingness rate per group label, aligned to variant order."""
    print(f"compute_missingness_variants()  ts = {current_milli_time()}", flush=True)
    selects = []
    for label, idxs in group_indices.items():
        filt = _filter_clause(idxs)
        selects.append(f"AVG(CASE WHEN gt IS NULL THEN 1.0 ELSE 0.0 END) {filt} AS miss_{label}")
    sql = (f"WITH {_unnested_cte(N_SAMPLES)} "
           f"SELECT variant_idx, {', '.join(selects)} "
           f"FROM unnested GROUP BY variant_idx ORDER BY variant_idx")
    df = pb.sql(sql).collect()
    return {label: df[f'miss_{label}'].to_numpy().astype(float) for label in group_indices}


def compute_missingness_samples():
    """Per-sample GT missingness rate, aligned to sample_idx order."""
    sql = (f"WITH {_unnested_cte(N_SAMPLES)} "
           f"SELECT sample_idx, AVG(CASE WHEN gt IS NULL THEN 1.0 ELSE 0.0 END) AS miss "
           f"FROM unnested GROUP BY sample_idx ORDER BY sample_idx")
    df = pb.sql(sql).collect()
    return df['miss'].to_numpy().astype(float)


def compute_heterozygosity():
    """Per-sample heterozygosity rate (fraction of variants where the sample is het)."""
    print(f"compute_heterozygosity()  ts = {current_milli_time()}", flush=True)
    sql = (f"WITH {_unnested_cte(N_SAMPLES)} "
           f"SELECT sample_idx, AVG({_gt_is_het_sql('gt')}) AS het "
           f"FROM unnested GROUP BY sample_idx ORDER BY sample_idx")
    df = pb.sql(sql).collect()
    return df['het'].to_numpy().astype(float)


def compute_dp_diff(case_indices, control_indices):
    """Per-variant |mean(DP) over cases - mean(DP) over controls|, aligned to variant order."""
    print(f"compute_dp_diff()  ts = {current_milli_time()}", flush=True)
    sql = (f"WITH {_unnested_cte(N_SAMPLES)} "
           f"SELECT variant_idx, "
           f"AVG(dp) {_filter_clause(case_indices)} AS case_dp, "
           f"AVG(dp) {_filter_clause(control_indices)} AS control_dp "
           f"FROM unnested GROUP BY variant_idx ORDER BY variant_idx")
    df = pb.sql(sql).collect()
    case_dp = df['case_dp'].to_numpy().astype(float)
    control_dp = df['control_dp'].to_numpy().astype(float)
    return np.abs(case_dp - control_dp)


def compute_variant_level():
    """Per-variant AF, variant type (0=SNP,1=indel), chrom -- one pass, no unnest."""
    print(f"compute_variant_level()  ts = {current_milli_time()}", flush=True)
    sql = ("SELECT chrom, "
           "\"AF\"[1] AS af, "
           "CASE WHEN length(ref) = 1 AND length(alt) = 1 THEN 0 ELSE 1 END AS variant_type "
           "FROM vcf_table")
    return pb.sql(sql).collect()


def compute_cell_pair(stat1, stat2, indices):
    """All (stat1, stat2) cell values where both are non-null, for one group (or all).

    Returns two aligned numpy arrays.  This is the only query that collects a per-cell
    (long) result; two columns of ~variants*samples floats is O(100 MB), far below the
    process budget, and feeds the identical np.histogram2d the previous code used --
    keeping the exact-tier heatmap match free of binning round-off risk.
    """
    print(f"compute_cell_pair({stat1} vs {stat2})  ts = {current_milli_time()}", flush=True)
    expr1 = _gt_encode_sql('gt') if stat1 == 'GT' else stat1.lower()
    expr2 = _gt_encode_sql('gt') if stat2 == 'GT' else stat2.lower()
    where = f"WHERE ({expr1}) IS NOT NULL AND ({expr2}) IS NOT NULL"
    if indices is not None:
        where += f" AND sample_idx IN ({_idx_list_sql(indices)})"
    sql = (f"WITH {_unnested_cte(N_SAMPLES)} "
           f"SELECT CAST(({expr1}) AS DOUBLE) AS s1, CAST(({expr2}) AS DOUBLE) AS s2 "
           f"FROM unnested {where}")
    df = pb.sql(sql).collect()
    return df['s1'].to_numpy().astype(float), df['s2'].to_numpy().astype(float)


# ---------------------------------------------------------------------------
# Plotting (plt/sns calls and emitted-stats keys unchanged from prior version)
# ---------------------------------------------------------------------------

def plot_variant_stats(vstats, phenotypes, stat, percentiles=[1, 10, 50],
                       stat_label='Statistic', include_log_scale=True):
    print(f"plot_variant_stats(stat={stat})  ts = {current_milli_time()}", flush=True)
    mean_values = vstats['all']['mean']

    plt.figure(figsize=(9, 6))
    sns.histplot(mean_values[~np.isnan(mean_values)], bins=50)
    plt.title(f'Mean {stat} Across Variants')
    plt.xlabel(stat_label)
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/1_{stat}_mean_variants.png')
    plt.savefig(f'{output_dir}/1_{stat}_mean_variants.svg', format="svg")
    plt.close()

    if include_log_scale:
        plt.figure(figsize=(9, 6))
        sns.histplot(mean_values[~np.isnan(mean_values)], bins=50, log_scale=True)
        plt.title(f'Mean {stat} Across Variants - log scale')
        plt.xlabel(f'{stat_label} (log)')
        plt.ylabel('Variant Count')
        plt.savefig(f'{output_dir}/1b_{stat}_mean_variants.png')
        plt.savefig(f'{output_dir}/1b_{stat}_mean_variants.svg', format="svg")
        plt.close()

    if len(phenotypes) > 5:
        for perc in percentiles:
            perc_values = vstats['all'][perc]
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
                pheno_values = vstats[pheno][key]
                pheno_values = pheno_values[~np.isnan(pheno_values)]
                emit_stat(f'variant_stats_{stat}_p{perc}_pheno{pheno}', pheno_values)
                sns.histplot(pheno_values, bins=50, label=f'{pheno} ({perc})', alpha=0.5)
            plt.title(f'Mean {stat} by Phenotype' if perc == 'mean'
                      else f'{stat} ({perc}th Percentile) by Phenotype')
            plt.xlabel(stat_label)
            plt.ylabel('Variant Count')
            plt.legend()
            plt.savefig(f'{output_dir}/3_{stat}_percentile_{perc}_by_phenotype_variants.png')
            plt.savefig(f'{output_dir}/3_{stat}_percentile_{perc}_by_phenotype_variants.svg', format="svg")
            plt.close()

            if include_log_scale:
                plt.figure(figsize=(10, 6))
                for pheno in phenotypes:
                    pheno_values = vstats[pheno][key]
                    pheno_values = pheno_values[~np.isnan(pheno_values)]
                    sns.histplot(pheno_values, bins=50, label=f'{pheno} ({perc})', alpha=0.5, log_scale=True)
                plt.title(f'Mean {stat} by Phenotype - log scale' if perc == 'mean'
                          else f'{stat} ({perc}th Percentile) by Phenotype - log scale')
                plt.xlabel(f'{stat_label} (log)')
                plt.ylabel('Variant Count')
                plt.legend()
                plt.savefig(f'{output_dir}/3b_{stat}_percentile_{perc}_by_phenotype_variants.png')
                plt.savefig(f'{output_dir}/3b_{stat}_percentile_{perc}_by_phenotype_variants.svg', format="svg")
                plt.close()


def plot_sample_stats(sample_means, samples, pheno_df, phenotypes, stat, stat_label='Statistic'):
    print(f"plot_sample_stats(stat={stat})  ts = {current_milli_time()}", flush=True)
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


def plot_missingness(miss_variants, miss_samples, samples, pheno_df, phenotypes):
    print(f"plot_missingness()  ts = {current_milli_time()}", flush=True)
    if len(phenotypes) > 5:
        plt.figure(figsize=(9, 6))
        sns.histplot(miss_variants['all'], bins=50)
        plt.title('Missingness Rate Across Variants')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Variant Count')
        plt.savefig(f'{output_dir}/6_missingness_variants.png')
        plt.savefig(f'{output_dir}/6_missingness_variants.svg', format="svg")
        plt.close()
    else:
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_missing = miss_variants[pheno]
            emit_stat(f'missingness_variants_pheno{pheno}', pheno_missing)
            sns.histplot(pheno_missing, bins=50, label=pheno, alpha=0.5)
        plt.title('Missingness Rate by Phenotype (Variants)')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Variant Count')
        plt.legend()
        plt.savefig(f'{output_dir}/7_missingness_by_phenotype_variants.png')
        plt.savefig(f'{output_dir}/7_missingness_by_phenotype_variants.svg', format="svg")
        plt.close()

    sample_miss = pl.DataFrame({'IID': samples, 'missing_rate': miss_samples}).join(
        pheno_df.select(['IID', 'Y1']), on='IID', how='left')
    if len(phenotypes) > 5:
        plt.figure(figsize=(9, 6))
        sns.histplot(sample_miss['missing_rate'].drop_nulls().to_numpy(), bins=50)
        plt.title('Missingness Rate Across Samples')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Sample Count')
        plt.savefig(f'{output_dir}/8_missingness_samples.png')
        plt.savefig(f'{output_dir}/8_missingness_samples.svg', format="svg")
        plt.close()
    else:
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_data = sample_miss.filter(pl.col('Y1') == pheno)['missing_rate'].drop_nulls().to_numpy()
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


def plot_dp_differences(abs_diff, phenotypes):
    print(f"plot_dp_differences()  ts = {current_milli_time()}", flush=True)
    if len(phenotypes) != 2:
        return
    abs_diff = abs_diff[~np.isnan(abs_diff)]
    emit_stat('dp_diff', abs_diff)
    plt.figure(figsize=(9, 6))
    sns.histplot(abs_diff, bins=50)
    plt.title('Absolute DP Differences (Cases vs Controls)')
    plt.xlabel('Absolute DP Difference')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/10_dp_abs_diff_cases_controls.png')
    plt.savefig(f'{output_dir}/10_dp_abs_diff_cases_controls.svg', format="svg")
    plt.close()

    plt.figure(figsize=(9, 6))
    sns.histplot(abs_diff, bins=50, log_scale=True)
    plt.title('Absolute DP Differences (Cases vs Controls) - log scale')
    plt.xlabel('Absolute DP Difference (log)')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/10b_dp_abs_diff_cases_controls.png')
    plt.savefig(f'{output_dir}/10b_dp_abs_diff_cases_controls.svg', format="svg")
    plt.close()


def plot_allele_frequency(af):
    print(f"plot_allele_frequency()  ts = {current_milli_time()}", flush=True)
    af = af[~np.isnan(af)]
    emit_stat('allele_freq', af)
    plt.figure(figsize=(9, 6))
    sns.histplot(af, bins=100, log_scale=True)
    plt.title('Alternate Allele Frequency Distribution - log scale')
    plt.xlabel('Allele Frequency')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/11_allele_frequency.png')
    plt.savefig(f'{output_dir}/11_allele_frequency.svg', format="svg")
    plt.close()


def plot_variant_types(variant_df):
    print(f"plot_variant_types()  ts = {current_milli_time()}", flush=True)
    snp_count = int((variant_df['variant_type'] == 0).sum())
    indel_count = int((variant_df['variant_type'] == 1).sum())
    emit_stat('variant_types', np.array([snp_count, indel_count]))
    plt.figure(figsize=(9, 6))
    sns.countplot(data=variant_df.rename({'variant_type': 'Variant_Type'}).to_pandas(), x='Variant_Type')
    plt.title('Variant Type Distribution')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/12_variant_types.png')
    plt.savefig(f'{output_dir}/12_variant_types.svg', format="svg")
    plt.close()


def plot_chrom_density(variant_df):
    print(f"plot_chrom_density()  ts = {current_milli_time()}", flush=True)
    chrom_counts = variant_df.group_by('chrom').len().sort('chrom').to_pandas()
    chrom_counts = chrom_counts.rename(columns={'chrom': 'CHROM', 'len': 'len'})
    chrom_counts[['CHROM', 'len']].to_csv(str(Path(stats_dir) / 'chrom_density.csv'), index=False)
    plt.figure(figsize=(12, 6))
    sns.barplot(x='CHROM', y='len', data=chrom_counts)
    plt.title('Variant Count by Chromosome')
    plt.xlabel('Chromosome')
    plt.ylabel('Variant Count')
    plt.xticks(rotation=45)
    plt.savefig(f'{output_dir}/13_chrom_variant_density.png')
    plt.savefig(f'{output_dir}/13_chrom_variant_density.svg', format="svg")
    plt.close()


def plot_heterozygosity(het_rates, samples, pheno_df, phenotypes):
    print(f"plot_heterozygosity()  ts = {current_milli_time()}", flush=True)
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
    print(f"plot_boxplots(stat={stat})  ts = {current_milli_time()}", flush=True)
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


def plot_stat_vs_stat(group_indices, phenotypes, stat1, stat2):
    print(f"plot_stat_vs_stat({stat1} vs {stat2})  ts = {current_milli_time()}", flush=True)
    stat1_values, stat2_values = compute_cell_pair(stat1, stat2, None)
    if len(stat1_values) == 0 or max(stat1_values) <= 0:
        return

    stat1_bins = np.arange(0, max(stat1_values) + 1, 1)
    stat2_bins = np.linspace(0, 2, 201)
    heatmap, xedges, yedges = np.histogram2d(stat1_values, stat2_values, bins=[stat1_bins, stat2_bins])
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
            s1, s2 = compute_cell_pair(stat1, stat2, group_indices[pheno])
            if len(s1) == 0:
                continue
            stat1_bins = np.arange(0, np.nanmax(s1) + 1, 1)
            stat2_bins = np.linspace(0, np.nanmax(s2), 201)
            heatmap, xedges, yedges = np.histogram2d(s1, s2, bins=[stat1_bins, stat2_bins])
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(vcf_file, phenotype_file, percentiles, use_dosage):
    global N_SAMPLES
    print(f"main()  percentiles = {percentiles}  use_dosage = {use_dosage}  "
          f"ts = {current_milli_time()}", flush=True)

    samples = read_samples(vcf_file)
    N_SAMPLES = len(samples)
    pheno_df = load_phenotype(phenotype_file)
    pb.register_vcf(vcf_file, name='vcf_table',
                    format_fields=['GT', 'DP', 'GQ', 'DS'], info_fields=['AF'])

    groups, phenotypes = build_groups(samples, pheno_df)
    # group_indices used for per-variant aggregates: each phenotype plus 'all'
    variant_groups = dict(groups)
    variant_groups['all'] = None  # None -> unfiltered aggregate over all samples

    for stat, label in [('DP', 'Depth of Coverage'), ('GQ', 'Genotype Quality')]:
        vstats = compute_variant_stats(stat, variant_groups, percentiles)
        plot_variant_stats(vstats, phenotypes, stat, percentiles=percentiles, stat_label=label)
        plot_sample_stats(compute_sample_means(stat), samples, pheno_df, phenotypes, stat, stat_label=label)

    if use_dosage:
        vstats = compute_variant_stats('DS', variant_groups, percentiles)
        plot_variant_stats(vstats, phenotypes, 'DS', percentiles=percentiles,
                           stat_label='Genotype Dosage', include_log_scale=False)
        plot_sample_stats(compute_sample_means('DS'), samples, pheno_df, phenotypes, 'DS',
                          stat_label='Genotype Dosage')

    plot_missingness(compute_missingness_variants(variant_groups), compute_missingness_samples(),
                     samples, pheno_df, phenotypes)

    if len(phenotypes) == 2:
        plot_dp_differences(compute_dp_diff(groups[phenotypes[0]], groups[phenotypes[1]]), phenotypes)

    variant_df = compute_variant_level()
    plot_allele_frequency(variant_df['af'].to_numpy().astype(float))
    plot_variant_types(variant_df)
    plot_chrom_density(variant_df)
    plot_heterozygosity(compute_heterozygosity(), samples, pheno_df, phenotypes)

    plot_boxplots(compute_sample_means('DP'), samples, pheno_df, 'DP', 'Depth of Coverage')
    plot_boxplots(compute_sample_means('GQ'), samples, pheno_df, 'GQ', 'Genotype Quality')
    if use_dosage:
        plot_boxplots(compute_sample_means('DS'), samples, pheno_df, 'DS', 'Genotype Dosage')
        plot_stat_vs_stat(groups, phenotypes, stat1='GQ', stat2='DS')
        plot_stat_vs_stat(groups, phenotypes, stat1='DP', stat2='DS')
        plot_stat_vs_stat(groups, phenotypes, stat1='GQ', stat2='DP')
        plot_stat_vs_stat(groups, phenotypes, stat1='GT', stat2='DS')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--phenotype', required=True)
    parser.add_argument('--use-dosage', default='false')
    parser.add_argument('--process-name', default='EXPLORATORY_DATA_ANALYSIS')
    args = parser.parse_args()

    main(args.vcf, args.phenotype, [1, 50], args.use_dosage == 'true')

    with open('versions.yml', 'w') as f:
        f.write(f'{args.process_name}:\n')
        f.write(f'    python: {sys.version.split()[0]}\n')
        f.write(f'    polars-bio: {pb.__version__}\n')
