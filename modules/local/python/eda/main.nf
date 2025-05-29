process EXPLORATORY_DATA_ANALYSIS {

    tag "$meta.id"
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/python_tools:1.0.1'

    input:
    tuple val(meta), path(vcf), path(tbi), path(phenotype_file)
    
    output:
    tuple val(meta), path("plots/*.png"), emit: plots
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
from pathlib import Path

# Create output directory for plots
output_dir = "plots"
Path(output_dir).mkdir(exist_ok=True)

# Function to load VCF.gz and extract DP, GQ, and missingness
def load_vcf(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    data = {'CHROM': [], 'POS': [], 'REF': [], 'ALT': [], 'AF': []}
    for sample in samples:
        data[f'DP_{sample}'] = []
        data[f'GQ_{sample}'] = []
        data[f'GT_{sample}'] = []  # To track missingness
    
    for record in vcf.fetch():
        data['CHROM'].append(record.chrom)
        data['POS'].append(record.pos)
        data['REF'].append(record.ref)
        data['ALT'].append(','.join(record.alts or ['.']))
        if 'AF' in record.info:
            data['AF'].append(float(record.info['AF'][0]))
        for sample in samples:
            sample_data = record.samples[sample]

            # Depth of coverage (DP)
            dp = sample_data.get('DP', np.nan)
            data[f'DP_{sample}'].append(dp)
            
            # Genotype quality (GQ)
            gq = sample_data.get('GQ', np.nan)
            data[f'GQ_{sample}'].append(gq)
            
            # Genotype (to check for missingness)
            gt = sample_data['GT']
            gt_value = np.nan if gt is None or None in gt else '/'.join(map(str, gt))
            data[f'GT_{sample}'].append(gt_value)
    
    vcf_df = pd.DataFrame(data)
    return vcf_df, samples

# Load phenotype file
def load_phenotype(phenotype_file):
    pheno_df = pd.read_csv(phenotype_file, sep='\t')
    return pheno_df

# Plotting function for variant-level statistics
def plot_variant_stats(df, pheno_df, samples, stat, percentiles=[1, 10, 50], stat_label='Statistic'):
    # Compute statistics per variant
    stat_cols = [f'{stat}_{sample}' for sample in samples]
    stats = df[stat_cols]
    phenotypes = pheno_df['Y1'].unique()
    
    # Mean
    mean_values = stats.mean(axis=1)
    plt.figure(figsize=(8, 6))
    sns.histplot(mean_values.dropna(), bins=50)
    plt.title(f'Mean {stat} Across Variants')
    plt.xlabel(stat_label)
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/1_{stat}_mean_variants.png')
    plt.close()

    plt.figure(figsize=(8, 6))
    sns.histplot(mean_values.dropna(), bins=50, log_scale=True)
    plt.title(f'Mean {stat} Across Variants (log scale)')
    plt.xlabel(stat_label + " (log)")
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/1b_{stat}_mean_variants.png')
    plt.close()

    if len(phenotypes) > 5:
        # Overall histograms
        for perc in percentiles:
            perc_values = stats.quantile(perc / 100, axis=1)
            plt.figure(figsize=(8, 6))
            sns.histplot(perc_values.dropna(), bins=50)
            plt.title(f'{stat} ({perc}th Percentile) Across Variants')
            plt.xlabel(stat_label)
            plt.ylabel('Variant Count')
            plt.savefig(f'{output_dir}/2_{stat}_percentile_{perc}_variants.png')
            plt.close()
    else:
        # Per-phenotype histograms
        for perc in percentiles + ['mean']:
            plt.figure(figsize=(10, 6))
            for pheno in phenotypes:
                pheno_samples = pheno_df[pheno_df['Y1'] == pheno]['IID'].values
                pheno_cols = [f'{stat}_{sample}' for sample in pheno_samples if f'{stat}_{sample}' in df.columns]
                if pheno_cols:
                    if perc == 'mean':
                        pheno_values = df[pheno_cols].mean(axis=1)
                    else:
                        pheno_values = df[pheno_cols].quantile(perc / 100, axis=1)
                    sns.histplot(pheno_values.dropna(), bins=50, label=f'{pheno} ({perc})', alpha=0.5)
            plt.title(f'{stat} ({perc}th Percentile) by Phenotype')
            plt.xlabel(stat_label)
            plt.ylabel('Variant Count')
            plt.legend()
            plt.savefig(f'{output_dir}/3_{stat}_percentile_{perc}_by_phenotype_variants.png')
            plt.close()

# Plotting function for sample-level statistics
def plot_sample_stats(df, pheno_df, samples, stat, stat_label='Statistic'):
    stat_cols = [f'{stat}_{sample}' for sample in samples]
    sample_stats = df[stat_cols].mean().reset_index()
    sample_stats['IID'] = sample_stats['index'].str.replace(f'{stat}_', '')
    sample_stats['value'] = sample_stats[0]
    
    # Merge with phenotype
    sample_stats = sample_stats.merge(pheno_df[['IID', 'Y1']], on='IID', how='left')
    phenotypes = pheno_df['Y1'].unique()
    
    if len(phenotypes) > 5:
        # Overall histogram
        plt.figure(figsize=(8, 6))
        sns.histplot(sample_stats['value'].dropna(), bins=50)
        plt.title(f'Mean {stat} Across Samples')
        plt.xlabel(stat_label)
        plt.ylabel('Sample Count')
        plt.savefig(f'{output_dir}/4_{stat}_samples.png')
        plt.close()
    else:
        # Per-phenotype histogram
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_data = sample_stats[sample_stats['Y1'] == pheno]['value']
            if not pheno_data.empty:
                sns.histplot(pheno_data.dropna(), bins=50, label=pheno, alpha=0.5)
        plt.title(f'Mean {stat} by Phenotype (Samples)')
        plt.xlabel(stat_label)
        plt.ylabel('Sample Count')
        plt.legend()
        plt.savefig(f'{output_dir}/5_{stat}_by_phenotype_samples.png')
        plt.close()

# Plot missingness per variant and sample
def plot_missingness(df, pheno_df, samples):
    gt_cols = [f'GT_{sample}' for sample in samples]
    missingness = df[gt_cols].isna().mean(axis=0)  # Per sample
    missingness_variants = df[gt_cols].isna().mean(axis=1)  # Per variant
    phenotypes = pheno_df['Y1'].unique()
    
    if len(phenotypes) > 5:
        # Variant-level missingness
        plt.figure(figsize=(8, 6))
        sns.histplot(missingness_variants, bins=50)
        plt.title('Missingness Rate Across Variants')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Variant Count')
        plt.savefig(f'{output_dir}/6_missingness_variants.png')
        plt.close()
    else:
        # Per-phenotype variant-level missingness
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_samples = pheno_df[pheno_df['Y1'] == pheno]['IID'].values
            pheno_cols = [f'GT_{sample}' for sample in pheno_samples if f'GT_{sample}' in df.columns]
            if pheno_cols:
                pheno_missing = df[pheno_cols].isna().mean(axis=1)
                sns.histplot(pheno_missing, bins=50, label=pheno, alpha=0.5)
        plt.title('Missingness Rate by Phenotype (Variants)')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Variant Count')
        plt.legend()
        plt.savefig(f'{output_dir}/7_missingness_by_phenotype_variants.png')
        plt.close()
    
    # Sample-level missingness
    missingness_samples = missingness.reset_index()
    missingness_samples['IID'] = missingness_samples['index'].str.replace('GT_', '')
    missingness_samples['missing_rate'] = missingness_samples[0]
    missingness_samples = missingness_samples.merge(pheno_df[['IID', 'Y1']], on='IID', how='left')
    
    if len(phenotypes) > 5:
        plt.figure(figsize=(8, 6))
        sns.histplot(missingness_samples['missing_rate'], bins=50)
        plt.title('Missingness Rate Across Samples')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Sample Count')
        plt.savefig(f'{output_dir}/8_missingness_samples.png')
        plt.close()
    else:
        # Per-phenotype sample-level missingness
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_data = missingness_samples[missingness_samples['Y1'] == pheno]['missing_rate']
            if not pheno_data.empty:
                sns.histplot(pheno_data, bins=50, label=pheno, alpha=0.5)
        plt.title('Missingness Rate by Phenotype (Samples)')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Sample Count')
        plt.legend()
        plt.savefig(f'{output_dir}/9_missingness_by_phenotype_samples.png')
        plt.close()

# Plot absolute DP differences for cases vs controls
def plot_dp_differences(df, pheno_df, samples):
    phenotypes = pheno_df['Y1'].unique()
    if len(phenotypes) == 2:
        case_samples = pheno_df[pheno_df['Y1'] == phenotypes[0]]['IID'].values
        control_samples = pheno_df[pheno_df['Y1'] == phenotypes[1]]['IID'].values
        case_cols = [f'DP_{sample}' for sample in case_samples if f'DP_{sample}' in df.columns]
        control_cols = [f'DP_{sample}' for sample in control_samples if f'DP_{sample}' in df.columns]
        
        if case_cols and control_cols:
            case_dp = df[case_cols].mean(axis=1)
            control_dp = df[control_cols].mean(axis=1)
            abs_diff = np.abs(case_dp - control_dp)
            plt.figure(figsize=(8, 6))
            sns.histplot(abs_diff.dropna(), bins=50)
            plt.title('Absolute DP Differences (Cases vs Controls)')
            plt.xlabel('Absolute DP Difference')
            plt.ylabel('Variant Count')
            plt.savefig(f'{output_dir}/10_dp_abs_diff_cases_controls.png')
            plt.close()

            plt.figure(figsize=(8, 6))
            sns.histplot(abs_diff.dropna(), bins=50, log_scale=True)
            plt.title('Absolute DP Differences (Cases vs Controls) - log scale')
            plt.xlabel('Absolute DP Difference (log)')
            plt.ylabel('Variant Count')
            plt.savefig(f'{output_dir}/10b_dp_abs_diff_cases_controls.png')
            plt.close()

def plot_allele_frequency(vcf_df, pheno_df):
    phenotypes = pheno_df['Y1'].unique()
    # TODO: add plots per phenotype

    af = vcf_df['AF']
    plt.figure(figsize=(8, 6))
    sns.histplot(af, bins=100)
    plt.title('Alternate Allele Frequency Distribution')
    plt.xlabel('Allele Frequency')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/11_allele_frequency.png')
    plt.close()

def plot_variant_types(vcf_df, pheno_df, samples):
    vcf_df['Variant_Type'] = vcf_df.apply(
        lambda x: 'SNP' if len(x['REF']) == 1 and len(x['ALT'].split(',')[0]) == 1 else 'Indel', axis=1)
    plt.figure(figsize=(8, 6))
    sns.countplot(data=vcf_df, x='Variant_Type')
    plt.title('Variant Type Distribution')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/12_variant_types.png')
    plt.close()

def plot_chrom_density(vcf_df):
    chrom_counts = vcf_df['CHROM'].value_counts()
    plt.figure(figsize=(12, 6))
    sns.barplot(x=chrom_counts.index, y=chrom_counts.values)
    plt.title('Variant Count by Chromosome')
    plt.xlabel('Chromosome')
    plt.ylabel('Variant Count')
    plt.xticks(rotation=45)
    plt.savefig(f'{output_dir}/13_chrom_variant_density.png')
    plt.close()

def plot_heterozygosity(vcf_df, pheno_df, samples):
    het_rates = []
    for sample in samples:
        gt_col = f'GT_{sample}'
        het_rate = (vcf_df[gt_col].str.contains('0/1|1/0', na=False)).mean()
        het_rates.append({'IID': sample, 'Heterozygosity': het_rate})
    het_df = pd.DataFrame(het_rates).merge(pheno_df[['IID', 'Y1']], on='IID', how='left')
    phenotypes = pheno_df['Y1'].unique()

    if len(phenotypes) > 5:
        plt.figure(figsize=(8, 6))
        sns.histplot(het_df['Heterozygosity'], bins=50)
        plt.title('Heterozygosity Rate Across Samples')
        plt.xlabel('Heterozygosity Rate')
        plt.ylabel('Sample Count')
        plt.savefig(f'{output_dir}/14_heterozygosity_samples.png')
        plt.close()
    else:
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_data = het_df[het_df['Y1'] == pheno]['Heterozygosity']
            if not pheno_data.empty:
                sns.histplot(pheno_data, bins=50, label=pheno, alpha=0.5)
        plt.title('Heterozygosity Rate by Phenotype')
        plt.xlabel('Heterozygosity Rate')
        plt.ylabel('Sample Count')
        plt.legend()
        plt.savefig(f'{output_dir}/15_heterozygosity_by_phenotype_samples.png')
        plt.close()

def plot_boxplots(vcf_df, pheno_df, samples, stat, stat_label):
    stat_cols = [f'{stat}_{sample}' for sample in samples]
    sample_stats = vcf_df[stat_cols].mean().reset_index()
    sample_stats['IID'] = sample_stats['index'].str.replace(f'{stat}_', '')
    sample_stats = sample_stats.merge(pheno_df[['IID', 'Y1']], on='IID', how='left')
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=sample_stats, x='Y1', y=0)
    plt.title(f'{stat_label} by Phenotype (Samples)')
    plt.xlabel('Phenotype')
    plt.ylabel(stat_label)
    plt.savefig(f'{output_dir}/16_{stat}_boxplot_by_phenotype.png')
    plt.close()

# Main execution
def main(vcf_file, phenotype_file, percentiles):
    # Load data
    vcf_df, samples = load_vcf(vcf_file)
    pheno_df = load_phenotype(phenotype_file)
    
    # Plot DP (variant and sample level)
    plot_variant_stats(vcf_df, pheno_df, samples, 'DP', percentiles=percentiles, stat_label='Depth of Coverage')
    plot_sample_stats(vcf_df, pheno_df, samples, 'DP', stat_label='Depth of Coverage')
    
    # Plot GQ (variant and sample level)
    plot_variant_stats(vcf_df, pheno_df, samples, 'GQ', percentiles=percentiles, stat_label='Genotype Quality')
    plot_sample_stats(vcf_df, pheno_df, samples, 'GQ', stat_label='Genotype Quality')
    
    plot_missingness(vcf_df, pheno_df, samples)
    
    # Plot DP differences if exactly two phenotypes
    plot_dp_differences(vcf_df, pheno_df, samples)

    plot_allele_frequency(vcf_df, pheno_df)
    plot_variant_types(vcf_df, pheno_df, samples)
    plot_chrom_density(vcf_df)
    plot_heterozygosity(vcf_df, pheno_df, samples)
    
    plot_boxplots(vcf_df, pheno_df, samples, 'DP', 'Depth of Coverage')
    plot_boxplots(vcf_df, pheno_df, samples, 'GQ', 'Genotype Quality')

# Run the analysis
vcf_file = "${vcf}"
phenotype_file = "${phenotype_file}"
percentiles = [1, 10, 50]
main(vcf_file, phenotype_file, percentiles)


# Write versions.yml
with open('versions.yml', 'w') as f:
    f.write('${task.process}:\\n')
    f.write(f'    python: {sys.version.split()[0]}\\n')
    f.write(f'    d3js: v7\\n')
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir plots
    touch plots/heterozygosity_by_phenotype_samples.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}