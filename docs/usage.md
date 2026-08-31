# Usage

How to prepare the input, run the pipeline, and set every parameter.

- [Requirements](#requirements)
- [Input files](#input-files)
- [Running the pipeline](#running-the-pipeline)
- [Profiles](#profiles)
- [Running several datasets in one go](#running-several-datasets-in-one-go)
- [Running on a cluster](#running-on-a-cluster)
- [Parameter reference](#parameter-reference)
- [Troubleshooting](#troubleshooting)

## Requirements

- [Nextflow](https://www.nextflow.io/) 25.10.2 or newer, and Java 17 or newer.
- A container engine: Docker, Singularity, Apptainer, Podman, Shifter or Charliecloud.
  Conda and Mamba also work, but are slower and give weaker reproducibility.

Nothing has to be downloaded in advance. The VEP cache and the reference genome are
fetched automatically the first time the preparation step runs.

## Input files

### The VCF

`--input_vcf` is a single multi-sample VCF, optionally compressed. A multi-sample VCF may
be produced from gvcf files for example by [GLNexus](https://github.com/dnanexus-rnd/GLnexus). The VCF file should also contain an allele frequency field in the INFO column.
By default the pipeline "prepares" the VCF: normalises, deduplicates and annotates the
VCF, so this doesn't have to be done beforehand by the user.

For genotype dosages to be computed, the FORMAT field must contain genotype likelihoods
(`PL`) alongside `GT`, `GQ` and `DP`. Without them the pipeline still runs, but on hard
genotype calls only (in such case `GT`, `GQ` and `DP` are mandatory).

Underscores in sample names are replaced with hyphens, because some downstream tools
treat the underscore as a field separator. This applies to the phenotype files too, so
the two stay consistent. Change the character with
`--bcftools_replace_sample_names_sed_arg` and `--rscript_build_phenotypes_options`.

### Case and control lists, or a phenotype file

You must supply one of these two, and they are mutually exclusive.

**Option 1: case and control lists.** Two plain text files, one sample identifier per
line:

```text title="cases.txt"
HG03925
HG03926
HG03927
```

Pass them with `--input_cases` and `--input_controls`. The pipeline builds the phenotype
file for you: every sample listed in the case file gets 1, every sample listed in the
control file gets 0.

A case file may optionally have a second, tab-separated column of numeric values, in
which case the phenotype is treated as quantitative rather than case/control.

**Option 2: a phenotype file.** A tab-separated file with a header, in the format REGENIE
expects:

```text title="phenotype.tsv"
FID	IID	Y1
HG00096	HG00096	0
HG00097	HG00097	1
```

Pass it with `--input_phenotype`. `FID` and `IID` are the family and individual
identifiers, and may be the same value. `Y1` is the phenotype: 1 for a case, 0 for a
control.

### The variant grouping file

`--input_masks` defines which VEP consequences belong to which impact group. The default,
[`assets/default.masks`](../assets/default.masks), defines three groups and is used unless
you supply your own. It is tab-separated, with a group name and then a comma-separated
list of VEP consequence terms:

```text
Mask_High	stop_gained,stop_lost,start_lost,frameshift_variant,splice_donor_variant,...
Mask_Mod	missense_variant,inframe_insertion,inframe_deletion,...
Mask_HighMod	stop_gained,stop_lost,...,missense_variant,inframe_insertion,...
```

`Mask_High` covers variants predicted to disrupt the protein, `Mask_Mod` covers those
predicted to alter it moderately, and `Mask_HighMod` is the union of the two. Association
tests are run for each group separately.

## Running the pipeline

```bash
nextflow run main.nf -profile docker \
    --input_vcf     /path/to/input.vcf.gz \
    --input_cases   /path/to/cases.txt \
    --input_controls /path/to/controls.txt \
    --project_name  myproject \
    --outdir        results
```

`--project_name` is a short identifier that appears in output filenames. `--outdir` is
where results are written.

### Two ways of running

**The full pipeline** (`--skip_preparation false --skip_reporting false`, the default)
runs everything, including a nested call to
[nf-prepare-vcf](https://github.com/biodatageeks/nf-prepare-vcf) that normalises, annotates
and computes dosages for the input VCF. This is the mode to use for a real analysis.

**Skipping preparation and reporting** (`--skip_preparation true --skip_reporting true`)
expects an already-prepared VCF and stops after the association tests. This is what
`nextflow-gene-assoc-tuner` uses when tuning the pipeline's parameters on simulated data:
`nf-prepare-vcf` is run once, and this pipeline is then run many times with different
settings. The two flags are independent - you can skip only one of them - but this
combination was the one that mattered in practice, because preparation and reporting
could be avoided for parameter optimization runs.

If you pass `--skip_preparation true`, the VCF must already be normalised, have
multi-allelic sites split, carry `CHROM_POS_REF_ALT` variant identifiers, and contain VEP
consequence annotations in a `CSQ` field. Otherwise the association step will match no
variants.

## Profiles

Combine profiles with commas, for example `-profile docker,low_resources`.

| Profile | What it does |
|---|---|
| `docker`, `podman`, `singularity`, `apptainer` | container engine to use |
| `conda`, `mamba` | run from conda environments instead of containers |
| `low_resources` | reduce the memory and CPU each process requests, so the pipeline fits on a workstation. Required for the test suite. |
| `medium_resources` | a middle setting between `low_resources` and the defaults |
| `slurm`, `hq` | submit work to a Slurm cluster, or to HyperQueue |
| `nocache` | disable Nextflow's task cache |
| `debug` | keep working directories and print more detail |
| `test`, `test_full`, `test_sim_chr22`, `test_skip_preparation_and_reporting` | small built-in datasets for checking that an installation works |

## Running several datasets in one go

`--input_phenotype` accepts a comma-separated list of phenotype files. Each is combined
with the same VCF and analysed independently within a single run. In such case it is best to
run VCF preparation once and then take advantage of `--skip_preparation true`, which avoids
repeating the expensive preparation step.

For results to be kept apart, each filename must follow the pattern

```
<anything>_dataset_idx_<N>_<anything>.phenotype.txt
```

where `<N>` is a number. The pipeline reads `<N>` out of the filename and labels that
dataset's outputs `<project_name>_dataset_idx_<N>`. A file that does not match the
pattern is still processed, but is labelled with `--project_name` alone -- so if you pass
several such files, their results will collide.

## Running on a cluster

Use `-profile slurm` or `-profile hq` together with a container profile, for example
`-profile apptainer,slurm`.

Two parameters are often needed on shared systems:

- `--tmpdir` sets the temporary directory for every process. Useful when the node's
  default temporary space is too small for the intermediate files.
- `--errorStrategy` overrides how failures are handled. In addition to Nextflow's own
  strategies (`retry`, `ignore`, `terminate`, `finish`), this pipeline accepts
  `retryThenIgnore`, which retries a failing task and then continues the run without it
  if it still fails. This is necessary is some of the datasets might fail but we
  would still be interested in results for other datasets.

`--cpu_support_avx2 false` switches to a PLINK2 build that does not require AVX2
instructions. Set it if the pipeline fails with an illegal-instruction error on older
processors.

## Parameter reference

### Required

| Parameter | Description |
|---|---|
| `--input_vcf` | Path to the input multi-sample VCF |
| `--outdir` | Directory to write results to |

Plus **either** `--input_phenotype`, **or** both `--input_cases` and `--input_controls`.

### Input

| Parameter | Default | Description |
|---|---|---|
| `--input_cases` | -- | File listing case sample identifiers, one per line |
| `--input_controls` | -- | File listing control sample identifiers, one per line |
| `--input_phenotype` | -- | Tab-separated phenotype file, or a comma-separated list of them |
| `--input_masks` | `assets/default.masks` | Definition of the VEP consequence impact groups |
| `--project_name` | -- | Short identifier used in output filenames |
| `--hild_path` | `assets/hg38_hild.txt` | Regions of high or unusual linkage disequilibrium, excluded from pruning |

### Behaviour

| Parameter | Default | Description |
|---|---|---|
| `--skip_preparation` | `false` | Skip VCF preparation and expect an already-prepared VCF |
| `--skip_reporting` | `false` | Skip all reporting steps |
| `--use_dosage` | `false` | Use the DS dosage field rather than hard genotype calls in the association tests |
| `--publish_intermediate` | `false` | Copy intermediate files into `--outdir` as well as the final results |
| `--regenie_step1_kinship_filtering` | `false` | Apply relatedness-based sample filtering to the REGENIE step 1 input |
| `--cpu_support_avx2` | `true` | Use the AVX2-optimised PLINK2 build; set to `false` on older processors |
| `--tmpdir` | Nextflow default | Temporary directory for all processes |
| `--errorStrategy` | Nextflow default | Process error strategy; also accepts `retryThenIgnore` |

### Variant and genotype quality filters

Applied to the VCF before any PLINK2 step. The first four are per-variant, computed
across all samples; the last three are per-genotype, and a genotype failing them is set
to missing rather than the variant being removed.

| Parameter | Default | Description |
|---|---|---|
| `--filter_vcf_qual_min` | `25` | Minimum variant QUAL score |
| `--filter_vcf_avg_gq_min` | `25` | Minimum average genotype quality across samples |
| `--filter_vcf_avg_dp_min` | `25` | Minimum average read depth across samples |
| `--filter_vcf_avg_dp_max` | `200` | Maximum average read depth across samples |
| `--filter_vcf_sample_gq_min` | `20` | Genotype quality below which a genotype is set to missing |
| `--filter_vcf_sample_dp_min` | `20` | Minimum read depth for a genotype |
| `--filter_vcf_sample_dp_max` | `250` | Maximum read depth for a genotype |

### PLINK2 quality control

After the shared initial steps the data follows four paths, each with its own thresholds:
inbreeding-coefficient filtering (in the table below identified by "ICF"), 
principal-component analysis ("PCA"), REGENIE step 1 ("R1"), and REGENIE step 2 ("R2").
The reasoning behind the separate settings is in [pipeline.md](pipeline.md).

| Parameter | Default | Description | Data flow paths |
|---|---|---|---|
| `--plink2_makepgen_1_options` | `--double-id --vcf-half-call missing --split-par b38 --1` | Initial VCF-to-pgen import | ICF, PCA, R1, R2 |
| `--plink2_makepgen_2_options` | `--impute-sex max-female-xf=0.2 min-male-xf=0.8` | Sex imputation thresholds | ICF, PCA, R1, R2 |
| `--plink2_missing_per_pheno_options` | `--geno 0.2` | Per-phenotype variant missingness filter | ICF, PCA, R1, R2 |
| `--plink2_makepgen_3_options` | `--geno 0.1 --hwe 1e-13 0.001 --mac 70 --maf 0.01` | Common-variant filtering | ICF, PCA, R1 |
| `--inbreeding_outliers_range_stds` | `3` | Standard deviations from the mean inbreeding coefficient beyond which a sample is an outlier | ICF |
| `--plink2_indep_pairwise_options` | `--mind 0.1` | Samples missingness filtering | ICF, PCA |
| `--plink2_indep_pairwise_window` | `50 5 0.2` | Linkage-disequilibrium pruning window | ICF |
| `--plink2_indep_pairwise_window_pca` | `500 50 0.2` | Linkage-disequilibrium pruning window | PCA |
| `--plink2_write_snplist_qc_options` | `--mind 0.1` | Samples missingness filtering | R1 |
| `--plink2_king_cutoff_threshold_pca` | `0.0884` | Relatedness cutoff used to choose the samples principal components are computed from | PCA |
| `--plink2_pca_settings` | `allele-wts 10` | Principal-component settings, including how many components to compute | PCA |
| `--plink2_makepgen_4_options` | `--geno 0.2` | Variant missingness filtering | R2 |
| `--plink2_makepgen_5_options` | `--mind 0.2` | Samples missingness filtering | R2 |
| `--plink2_write_snplist_step2_options` | `--mind 0.2` | Samples missingness filtering (this filtering is actually not needed, after one of the code refactorings it now duplicates the work done in plink2_makepgen_5_options) | R2 |
| `--plink2_export_other_options` | `dosage=DS --double-id --vcf-half-call missing --split-par b38 --1 --export Av` | Dosage export settings | R2 |
| `--plink2_import_dosage_options` | `skip0=1 skip1=2 id-delim=_ chr-col-num=1 pos-col-num=4 ref-first --make-pgen` | Dosage re-import settings | R2 |

The relatedness cutoff only decides which samples the principal components are computed
from. **Related samples are not removed from the association test** -- see
[pipeline.md](pipeline.md).

### Association testing

| Parameter | Default | Description |
|---|---|---|
| `--regenie_step1_options` | `--bt --bsize 100 --lowmem --covarColList PC1_AVG,PC2_AVG` | REGENIE step 1 (whole-genome model) |
| `--regenie_step2_options` | `--bt --minMAC 1 --ref-first --firth --approx --bsize 200 --lowmem --aaf-bins 0.01,0.05,0.1,1 --write-mask --write-mask-snplist --vc-tests skato --covarColList PC1_AVG,PC2_AVG` | REGENIE step 2 (association tests) |
| `--rscript_annotate_options` | `--min_top_annotations 30 --max_annotations 62 --quantile_threshold 0.25 --include-intergenic FALSE` | Variant annotation and grouping |
| `--rscript_vcf2aaf_options` | `AF_nfe_gnomad AF` | Which INFO fields to take alternative allele frequencies from, in order of preference |

`--bt` selects a case/control analysis; remove it for a quantitative phenotype.
`--covarColList` must name principal components that actually exist - with the default
`allele-wts 10`, that is `PC1_AVG` through `PC10_AVG`, and also `SEX` can be used, which is
the imputed sex.

### Reporting

| Parameter | Default | Description |
|---|---|---|
| `--manhattan_annotation_enabled` | `true` | Label significant genes on the Manhattan plot |
| `--annotation_min_log10p` | `3` | Least -log10(p) a gene needs before it is labelled |
| `--plot_ylimit` | `0` | Upper limit of the Manhattan plot's y axis; 0 means fit to the data |
| `--phenotypes_apply_rint` | `false` | Apply a rank-based inverse normal transformation to the phenotype before reporting |
| `--multiqc_config` | -- | Custom MultiQC configuration file |
| `--multiqc_title` | -- | Title for the MultiQC report |
| `--multiqc_logo` | -- | Logo for the MultiQC report |

### VEP annotation

Used only when `--skip_preparation` is `false`.

| Parameter | Default | Description |
|---|---|---|
| `--vep_cache_url` | Ensembl release 113, GRCh38 | Where to download the VEP cache from |
| `--ref_fasta_url` | Ensembl release 113, GRCh38 | Where to download the reference genome from |
| `--vep_annotate_options` | see `nextflow.config` | VEP command-line options |
| `--vep_updatecache_options` | `--AUTO acf --ASSEMBLY GRCh38` | VEP cache download options |

To use a different genome build, change `--vep_cache_url`, `--ref_fasta_url` and the
`--split-par` value inside `--plink2_makepgen_1_options` together. The default
`assets/hg38_hild.txt` is also build-specific.

## Troubleshooting

**The association step reports no variants, or no genes.** Might be a mismatch in
variant identifiers or chromosome names between the VCF and the grouping files. With
`--skip_preparation true`, check that the VCF has `CHROM_POS_REF_ALT` identifiers and a
`CSQ` field.

**An illegal-instruction error in a PLINK2 step.** The processor does not support AVX2.
Set `--cpu_support_avx2 false`.

**A process runs out of memory.** Try `-profile medium_resources`, or
raise the limits in a custom configuration file passed with `-c`.

**A process runs out of temporary space.** Set `--tmpdir` to a directory on a larger
filesystem.

**REGENIE stops with "all phenotypes have less than 10 cases".** There are genuinely too
few cases left after quality control. Loosening the missingness thresholds may help, but
below roughly a dozen cases no burden test can be run.

**Sex imputation gives implausible results, or many samples come out ambiguous.** The
default thresholds (`max-female-xf=0.2 min-male-xf=0.8`) assume reasonable coverage of
chromosome X. Exome data covers chromosome X sparsely, so the inbreeding estimate is
noisier and more samples land between the two thresholds. Widen them in
`--plink2_makepgen_2_options` if needed.

**Checking what was filtered out and where.** Most steps records how many samples and
variants it received and returned. Open
`<outdir>/generate_tracking_report/*_sankey_report.html` or the txt reports - they show 
the data paths through the pipeline, and may help in finding the step responsible 
for an unexpected loss of data.
