process EXPLORATORY_DATA_ANALYSIS {

    tag "$meta.id"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/python_tools:1.0.4'

    input:
    tuple val(meta), path(vcf), path(tbi), path(phenotype_file)
    
    output:
    tuple val(meta), path("plots/*.png"), emit: plots
    tuple val(meta), path("plots/*.svg"), emit: plots_svg
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
import os
from pathlib import Path
tmp_matplotlib_dir = Path("matplotlib_tmp")
tmp_matplotlib_dir.mkdir(exist_ok=True)
os.environ['MPLCONFIGDIR'] = str(tmp_matplotlib_dir.resolve())
import sys
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import time


plt.rcParams.update({'font.size': 20})


def current_milli_time():
    return round(time.time() * 1000)

# Create output directory for plots
output_dir = "plots"
Path(output_dir).mkdir(exist_ok=True)

# Function to encode genotype as integer using f(j,k) = (k*(k+1)/2) + j
def encode_genotype(gt):
    if gt is None or None in gt:
        return None  # Missing genotype
    j, k = sorted(gt)  # Sort to treat j/k and k/j the same (unphased or phased)
    return (k * (k + 1) // 2) + j  # Integer encoding

# Function to load VCF.gz and extract DP, GQ, and missingness
def load_vcf(vcf_file):
    print(f"load_vcf(\${vcf_file})  ts = \${current_milli_time()}", flush=True)
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    data = {'CHROM': [], 'POS': [], 'Variant_Type': [], 'AF': []}
    for sample in samples:
        data[f'DP_{sample}'] = []
        data[f'GQ_{sample}'] = []
        data[f'GT_{sample}'] = []

    i = 0
    for record in vcf.fetch():
        i += 1
        if i % 100000 == 0:
            print('.', end='', flush=True)

        data['CHROM'].append(record.chrom)
        data['POS'].append(record.pos)
        # Determine variant type (0 for SNP, 1 for Indel)
        is_snp = all(len(record.ref) == 1 and len(alt) == 1 for alt in (record.alts or ['.']))
        data['Variant_Type'].append(0 if is_snp else 1)
        # Allele frequency
        data['AF'].append(float(record.info['AF'][0]) if 'AF' in record.info else np.nan)
        
        for sample in samples:
            sample_data = record.samples[sample]
            # Depth of coverage (DP)
            dp = sample_data.get('DP', np.nan)
            data[f'DP_{sample}'].append(dp)
            # Genotype quality (GQ)
            gq = sample_data.get('GQ', np.nan)
            data[f'GQ_{sample}'].append(gq)
            # Genotype (encoded as integer)
            gt = sample_data['GT']
            data[f'GT_{sample}'].append(encode_genotype(gt))
 
    print(f"Loaded \${i} records from the vcf file  ts = \${current_milli_time()}", flush=True)
    
    # Convert to Polars DataFrame with appropriate dtypes
    dtypes = {'CHROM': pl.Utf8, 'POS': pl.Int32, 'Variant_Type': pl.Int8, 'AF': pl.Float32}
    for sample in samples:
        dtypes[f'DP_{sample}'] = pl.Float32
        dtypes[f'GQ_{sample}'] = pl.Float32
        dtypes[f'GT_{sample}'] = pl.Int16  # Int16 for multiallelic support
    vcf_df = pl.DataFrame(data, schema=dtypes)
    return vcf_df, samples

# Load phenotype file
def load_phenotype(phenotype_file):
    print(f"load_phenotype(\${phenotype_file})  ts = \${current_milli_time()}")
    pheno_df = pl.read_csv(phenotype_file, separator='\t')
    return pheno_df

# Plotting function for variant-level statistics
def plot_variant_stats(df, pheno_df, samples, stat, percentiles=[1, 10, 50], stat_label='Statistic'):
    print(f"plot_variant_stats(stat=\${stat})  ts = \${current_milli_time()}")
    stat_cols = [f'{stat}_{sample}' for sample in samples]
    phenotypes = pheno_df['Y1'].unique().to_list()
    
    # Mean
    mean_values = df.select(stat_cols).mean_horizontal().to_numpy()
    plt.figure(figsize=(9, 6))
    sns.histplot(mean_values[~np.isnan(mean_values)], bins=50)
    plt.title(f'Mean {stat} Across Variants')
    plt.xlabel(stat_label)
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/1_{stat}_mean_variants.png')
    plt.savefig(f'{output_dir}/1_{stat}_mean_variants.svg', format="svg")
    plt.close()

    plt.figure(figsize=(9, 6))
    sns.histplot(mean_values[~np.isnan(mean_values)], bins=50, log_scale=True)
    plt.title(f'Mean {stat} Across Variants - log scale')
    plt.xlabel(f'{stat_label} (log)')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/1b_{stat}_mean_variants.png')
    plt.savefig(f'{output_dir}/1b_{stat}_mean_variants.svg', format="svg")
    plt.close()

    if len(phenotypes) > 5:
        # Overall histograms
        for perc in percentiles:
            perc_values = np.quantile(df.select(stat_cols).to_numpy(), perc / 100, axis=1)
            
            plt.figure(figsize=(9, 6))
            sns.histplot(perc_values[~np.isnan(perc_values)], bins=50)
            plt.title(f'{stat} ({perc}th Percentile) Across Variants')
            plt.xlabel(stat_label)
            plt.ylabel('Variant Count')
            plt.savefig(f'{output_dir}/2_{stat}_percentile_{perc}_variants.png')
            plt.savefig(f'{output_dir}/2_{stat}_percentile_{perc}_variants.svg', format="svg")
            plt.close()

            plt.figure(figsize=(9, 6))
            sns.histplot(perc_values[~np.isnan(perc_values)], bins=50, log_scale=True)
            plt.title(f'{stat} ({perc}th Percentile) Across Variants - log scale')
            plt.xlabel(f'{stat_label} (log)')
            plt.ylabel('Variant Count')
            plt.savefig(f'{output_dir}/2b_{stat}_percentile_{perc}_variants.png')
            plt.savefig(f'{output_dir}/2b_{stat}_percentile_{perc}_variants.svg', format="svg")
            plt.close()
    else:
        # Per-phenotype histograms
        for perc in percentiles + ['mean']:
            plt.figure(figsize=(10, 6))
            for pheno in phenotypes:
                pheno_samples = pheno_df.filter(pl.col('Y1') == pheno)['IID'].to_list()
                pheno_cols = [f'{stat}_{sample}' for sample in pheno_samples if f'{stat}_{sample}' in df]
                if pheno_cols:
                    if perc == 'mean':
                        pheno_values = df.select(pheno_cols).mean_horizontal().to_numpy()
                    else:
                        pheno_values = np.quantile(df.select(pheno_cols).to_numpy(), perc / 100, axis=1)
                    sns.histplot(pheno_values[~np.isnan(pheno_values)], bins=50, label=f'{pheno} ({perc})', alpha=0.5)
            plt.title(f'Mean {stat} by Phenotype' if perc == 'mean' else f'{stat} ({perc}th Percentile) by Phenotype')
            plt.xlabel(stat_label)
            plt.ylabel('Variant Count')
            plt.legend()
            plt.savefig(f'{output_dir}/3_{stat}_percentile_{perc}_by_phenotype_variants.png')
            plt.savefig(f'{output_dir}/3_{stat}_percentile_{perc}_by_phenotype_variants.svg', format="svg")
            plt.close()

            plt.figure(figsize=(10, 6))
            for pheno in phenotypes:
                pheno_samples = pheno_df.filter(pl.col('Y1') == pheno)['IID'].to_list()
                pheno_cols = [f'{stat}_{sample}' for sample in pheno_samples if f'{stat}_{sample}' in df]
                if pheno_cols:
                    if perc == 'mean':
                        pheno_values = df.select(pheno_cols).mean_horizontal().to_numpy()
                    else:
                        pheno_values = df.select(pheno_cols).quantile(perc / 100).to_numpy()
                    sns.histplot(pheno_values[~np.isnan(pheno_values)], bins=50, label=f'{pheno} ({perc})', alpha=0.5, log_scale=True)
            plt.title(f'Mean {stat} by Phenotype - log scale' if perc == 'mean' else f'{stat} ({perc}th Percentile) by Phenotype - log scale')
            plt.xlabel(f'{stat_label} (log)')
            plt.ylabel('Variant Count')
            plt.legend()
            plt.savefig(f'{output_dir}/3b_{stat}_percentile_{perc}_by_phenotype_variants.png')
            plt.savefig(f'{output_dir}/3b_{stat}_percentile_{perc}_by_phenotype_variants.svg', format="svg")
            plt.close()

# Plotting function for sample-level statistics
def plot_sample_stats(df, pheno_df, samples, stat, stat_label='Statistic'):
    print(f"plot_sample_stats(stat=\${stat})  ts = \${current_milli_time()}")
    stat_cols = [f'{stat}_{sample}' for sample in samples]
    sample_stats = df.select(stat_cols).to_pandas().mean().reset_index()
    sample_stats['IID'] = sample_stats['index'].str.replace(f'{stat}_', '')
    sample_stats['value'] = sample_stats[0]
    sample_stats = pl.from_pandas(sample_stats[['IID', 'value']]).join(pheno_df.select(['IID', 'Y1']), on='IID', how='left')
    phenotypes = pheno_df['Y1'].unique().to_list()
    
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
                sns.histplot(pheno_data, bins=50, label=pheno, alpha=0.5)
        plt.title(f'Mean {stat} by Phenotype (Samples)')
        plt.xlabel(stat_label)
        plt.ylabel('Sample Count')
        plt.legend()
        plt.savefig(f'{output_dir}/5_{stat}_by_phenotype_samples.png')
        plt.savefig(f'{output_dir}/5_{stat}_by_phenotype_samples.svg', format="svg")
        plt.close()

# Plot missingness per variant and sample
def plot_missingness(df, pheno_df, samples):
    print(f"plot_missingness()  ts = \${current_milli_time()}")
    gt_cols = [f'GT_{sample}' for sample in samples]
    missingness_samples = (
        df.select(pl.col(gt_cols).is_null().mean())  # Calculate mean of nulls per column
        .unpivot()  # transpose
        .select(pl.col("variable").alias("index"), pl.col("value").alias("missing_rate"))
        .to_pandas()
    )
    missingness_samples['IID'] = missingness_samples['index'].str.replace('GT_', '')
    missingness_samples = missingness_samples.merge(pheno_df.to_pandas()[['IID', 'Y1']], on='IID', how='left')
    missingness_variants = df.select(pl.mean_horizontal(pl.col(gt_cols).is_null())).to_series().to_pandas()  # Per variant
    
    phenotypes = pheno_df['Y1'].unique().to_list()
    
    if len(phenotypes) > 5:
        plt.figure(figsize=(9, 6))
        sns.histplot(missingness_variants, bins=50)
        plt.title('Missingness Rate Across Variants')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Variant Count')
        plt.savefig(f'{output_dir}/6_missingness_variants.png')
        plt.savefig(f'{output_dir}/6_missingness_variants.svg', format="svg")
        plt.close()
    else:
        plt.figure(figsize=(10, 6))
        for pheno in phenotypes:
            pheno_samples = pheno_df.filter(pl.col('Y1') == pheno)['IID'].to_list()
            pheno_cols = [f'GT_{sample}' for sample in pheno_samples if f'GT_{sample}' in df]
            if pheno_cols:
                pheno_missing = df.select(pl.mean_horizontal(pl.col(pheno_cols).is_null())).to_series().to_numpy()
                sns.histplot(pheno_missing, bins=50, label=pheno, alpha=0.5)
        plt.title('Missingness Rate by Phenotype (Variants)')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Variant Count')
        plt.legend()
        plt.savefig(f'{output_dir}/7_missingness_by_phenotype_variants.png')
        plt.savefig(f'{output_dir}/7_missingness_by_phenotype_variants.svg', format="svg")
        plt.close()
    
    if len(phenotypes) > 5:
        plt.figure(figsize=(9, 6))
        sns.histplot(missingness_samples['missing_rate'].drop_nulls().to_numpy(), bins=50)
        plt.title('Missingness Rate Across Samples')
        plt.xlabel('Missingness Rate')
        plt.ylabel('Sample Count')
        plt.savefig(f'{output_dir}/8_missingness_samples.png')
        plt.savefig(f'{output_dir}/8_missingness_samples.svg', format="svg")
        plt.close()
    else:
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
        plt.savefig(f'{output_dir}/9_missingness_by_phenotype_samples.svg', format="svg")
        plt.close()

# Plot absolute DP differences for cases vs controls
def plot_dp_differences(df, pheno_df, samples):
    print(f"plot_dp_differences()  ts = \${current_milli_time()}")
    phenotypes = pheno_df['Y1'].unique().to_list()
    if len(phenotypes) == 2:
        case_samples = pheno_df.filter(pl.col('Y1') == phenotypes[0])['IID'].to_list()
        control_samples = pheno_df.filter(pl.col('Y1') == phenotypes[1])['IID'].to_list()
        case_cols = [f'DP_{sample}' for sample in case_samples if f'DP_{sample}' in df]
        control_cols = [f'DP_{sample}' for sample in control_samples if f'DP_{sample}' in df]
        
        if case_cols and control_cols:
            case_dp = df.select(case_cols).mean_horizontal().to_numpy()
            control_dp = df.select(control_cols).mean_horizontal().to_numpy()
            abs_diff = np.abs(case_dp - control_dp)
            plt.figure(figsize=(9, 6))
            sns.histplot(abs_diff[~np.isnan(abs_diff)], bins=50)
            plt.title('Absolute DP Differences (Cases vs Controls)')
            plt.xlabel('Absolute DP Difference')
            plt.ylabel('Variant Count')
            plt.savefig(f'{output_dir}/10_dp_abs_diff_cases_controls.png')
            plt.savefig(f'{output_dir}/10_dp_abs_diff_cases_controls.svg', format="svg")
            plt.close()

            plt.figure(figsize=(9, 6))
            sns.histplot(abs_diff[~np.isnan(abs_diff)], bins=50, log_scale=True)
            plt.title('Absolute DP Differences (Cases vs Controls) - log scale')
            plt.xlabel('Absolute DP Difference (log)')
            plt.ylabel('Variant Count')
            plt.savefig(f'{output_dir}/10b_dp_abs_diff_cases_controls.png')
            plt.savefig(f'{output_dir}/10b_dp_abs_diff_cases_controls.svg', format="svg")
            plt.close()

def plot_allele_frequency(vcf_df):
    print(f"plot_allele_frequency()  ts = \${current_milli_time()}")
    af = vcf_df['AF'].to_numpy()
    plt.figure(figsize=(9, 6))
    sns.histplot(af[~np.isnan(af)], bins=100, log_scale=True)
    plt.title('Alternate Allele Frequency Distribution - log scale')
    plt.xlabel('Allele Frequency')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/11_allele_frequency.png')
    plt.savefig(f'{output_dir}/11_allele_frequency.svg', format="svg")
    plt.close()

def plot_variant_types(vcf_df, samples):
    print(f"plot_variant_types()  ts = \${current_milli_time()}")
    plt.figure(figsize=(9, 6))
    sns.countplot(data=vcf_df.to_pandas(), x='Variant_Type')
    plt.title('Variant Type Distribution')
    plt.ylabel('Variant Count')
    plt.savefig(f'{output_dir}/12_variant_types.png')
    plt.savefig(f'{output_dir}/12_variant_types.svg', format="svg")
    plt.close()

def plot_chrom_density(vcf_df):
    print(f"plot_chrom_density()  ts = \${current_milli_time()}")
    chrom_counts = vcf_df.group_by('CHROM').len().sort('CHROM').to_pandas()
    plt.figure(figsize=(12, 6))
    sns.barplot(x='CHROM', y='len', data=chrom_counts)
    plt.title('Variant Count by Chromosome')
    plt.xlabel('Chromosome')
    plt.ylabel('Variant Count')
    plt.xticks(rotation=45)
    plt.savefig(f'{output_dir}/13_chrom_variant_density.png')
    plt.savefig(f'{output_dir}/13_chrom_variant_density.svg', format="svg")
    plt.close()

def plot_heterozygosity(vcf_df, pheno_df, samples):
    print(f"plot_heterozygosity()  ts = \${current_milli_time()}")
    het_rates = []
    for sample in samples:
        gt_col = f'GT_{sample}'
        # Homozygous genotypes: 0/0 (0), 1/1 (2), 2/2 (5), 3/3 (9) - ignore 4/4 and higher as they are extremely rare
        homozygous_values = [0, 2, 5, 9]
        # Heterozygous if GT is not null and not in homozygous_values
        het_rate = vcf_df.select(
            pl.col(gt_col).is_not_null() & (~pl.col(gt_col).is_in(homozygous_values))
        ).mean()[gt_col][0]
        het_rates.append({'IID': sample, 'Heterozygosity': het_rate})
    het_df = pl.from_dicts(het_rates).join(pheno_df.select(['IID', 'Y1']), on='IID', how='left')
    phenotypes = pheno_df['Y1'].unique().to_list()

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
                sns.histplot(pheno_data, bins=50, label=pheno, alpha=0.5)
        plt.title('Heterozygosity Rate by Phenotype')
        plt.xlabel('Heterozygosity Rate')
        plt.ylabel('Sample Count')
        plt.legend()
        plt.savefig(f'{output_dir}/15_heterozygosity_by_phenotype_samples.png')
        plt.savefig(f'{output_dir}/15_heterozygosity_by_phenotype_samples.svg', format="svg")
        plt.close()

def plot_boxplots(vcf_df, pheno_df, samples, stat, stat_label):
    print(f"plot_boxplots(stat=\${stat})  ts = \${current_milli_time()}")
    stat_cols = [f'{stat}_{sample}' for sample in samples]
    sample_stats = vcf_df.select(stat_cols).to_pandas().mean().reset_index()
    sample_stats['IID'] = sample_stats['index'].str.replace(f'{stat}_', '')
    sample_stats['value'] = sample_stats[0]
    sample_stats = pl.from_pandas(sample_stats[['IID', 'value']]).join(pheno_df.select(['IID', 'Y1']), on='IID', how='left')
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=sample_stats.to_pandas(), x='Y1', y='value')
    plt.title(f'{stat_label} by Phenotype (Samples)')
    plt.xlabel('Phenotype')
    plt.ylabel(stat_label)
    plt.savefig(f'{output_dir}/16_{stat}_boxplot_by_phenotype.png')
    plt.savefig(f'{output_dir}/16_{stat}_boxplot_by_phenotype.svg', format="svg")
    plt.close()

# Main execution
def main(vcf_file, phenotype_file, percentiles):
    print(f"main()  ts = \${current_milli_time()}", flush=True)

    vcf_df, samples = load_vcf(vcf_file)
    pheno_df = load_phenotype(phenotype_file)
    
    plot_variant_stats(vcf_df, pheno_df, samples, 'DP', percentiles=percentiles, stat_label='Depth of Coverage')
    plot_sample_stats(vcf_df, pheno_df, samples, 'DP', stat_label='Depth of Coverage')
    
    plot_variant_stats(vcf_df, pheno_df, samples, 'GQ', percentiles=percentiles, stat_label='Genotype Quality')
    plot_sample_stats(vcf_df, pheno_df, samples, 'GQ', stat_label='Genotype Quality')
    
    plot_missingness(vcf_df, pheno_df, samples)
    
    plot_dp_differences(vcf_df, pheno_df, samples)
    plot_allele_frequency(vcf_df)
    plot_variant_types(vcf_df, samples)
    plot_chrom_density(vcf_df)
    plot_heterozygosity(vcf_df, pheno_df, samples)
    
    plot_boxplots(vcf_df, pheno_df, samples, 'DP', 'Depth of Coverage')
    plot_boxplots(vcf_df, pheno_df, samples, 'GQ', 'Genotype Quality')

# Run the analysis
vcf_file = "${vcf}"
phenotype_file = "${phenotype_file}"
percentiles = [1, 50]
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
