# nf-rare-var-assoc

[![GitHub Actions CI Status](https://github.com/psuszyns/rare-var-assoc-nf/actions/workflows/ci.yml/badge.svg)](https://github.com/psuszyns/rare-var-assoc-nf/actions/workflows/ci.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.10.2-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**nf-rare-var-assoc** is a Nextflow pipeline for rare-variant association studies.
It ingests a multi-sample VCF and case/control sample lists, delegates VCF preparation
(normalisation, VEP annotation, variant filtering) to
[nf-prepare-vcf](https://github.com/psuszyns/nf-prepare-vcf), then runs PCA,
per-phenotype missingness and F-coefficient filtering, Regenie step1/step2 GWAS, and
produces association results alongside an HTML tracking and reporting artifact.

## Usage

### Required parameters

| Parameter | Description |
|---|---|
| `--input_vcf` | Path to the input multi-sample VCF (`.vcf.gz` with `.tbi` index) |
| `--outdir` | Output directory (absolute path) |

### Optional input parameters

| Parameter | Description |
|---|---|
| `--input_cases` | Directory containing per-phenotype case sample-ID files |
| `--input_controls` | Directory containing per-phenotype control sample-ID files |
| `--input_phenotype` | Path to a tab-separated phenotype file (alternative to cases/controls directories) |
| `--input_masks` | Path to a masks file defining variant sets for burden/SKAT tests |
| `--project_name` | Short project identifier used in output filenames |

### Run command

```bash
nextflow run psuszyns/nf-rare-var-assoc \
    --input_vcf /path/to/input.vcf.gz \
    --input_cases /path/to/cases/ \
    --input_controls /path/to/controls/ \
    --outdir results \
    --project_name myproject
```

### Pipeline behaviour

| Parameter | Default | Description |
|---|---|---|
| `--skip_preparation` | `false` | Skip VCF preparation (see production modes below) |
| `--skip_reporting` | `false` | Skip the reporting subworkflow |
| `--use_dosage` | `false` | Use DS (dosage) field instead of hard genotype calls in Regenie |
| `--publish_intermediate` | `false` | Copy intermediate files to `--outdir` |
| `--regenie_step1_kinship_filtering` | `false` | Enable kinship-based sample filtering in Regenie step 1 |
| `--cpu_support_avx2` | `true` | Use AVX2-optimised PLINK2 binary; set to `false` on older CPUs |
| `--tmpdir` | _(Nextflow default)_ | Override temporary directory for all processes |
| `--errorStrategy` | _(Nextflow default)_ | Override process error strategy (e.g. `retryThenIgnore`, `terminate`) |

### Production modes

**Main path (`--skip_preparation false`, default)** — runs the full pipeline including a
nested call to `nf-prepare-vcf` for VCF preparation (bcftools normalisation, VEP
annotation, variant filtering). Use this for clinical-centre runs where raw VCF input
is provided.

**HPC param-tuning path (`--skip_preparation true`)** — skips preparation and expects
a pre-prepared VCF (normalised, VEP-annotated, biallelic, unique IDs). Use this on
PLGrid HPC when tuning Regenie parameters with data already prepared out-of-band.

### VCF quality filters

These parameters apply only when `--skip_preparation false`. They control the
per-variant and per-genotype quality thresholds applied during VCF preparation.

| Parameter | Default | Description |
|---|---|---|
| `--filter_and_enhance_vcf_qual_min` | `25` | Minimum variant QUAL score |
| `--filter_and_enhance_vcf_avg_gq_min` | `25` | Minimum average GQ across samples |
| `--filter_and_enhance_vcf_avg_dp_min` | `25` | Minimum average DP across samples |
| `--filter_and_enhance_vcf_avg_dp_max` | `200` | Maximum average DP across samples |
| `--filter_and_enhance_vcf_sample_gq_min` | `20` | Per-sample GQ below which GT is set to missing |
| `--filter_and_enhance_vcf_sample_dp_min` | `20` | Per-sample minimum DP |
| `--filter_and_enhance_vcf_sample_dp_max` | `250` | Per-sample maximum DP |

### PLINK2 QC parameters

| Parameter | Default | Description |
|---|---|---|
| `--plink2_makepgen_1_options` | `--double-id --vcf-half-call missing --split-par b38 --1` | Options for initial VCF-to-pgen import (preparation path) |
| `--plink2_makepgen_2_options` | `--impute-sex max-female-xf=0.2 min-male-xf=0.8` | Sex imputation thresholds |
| `--plink2_makepgen_3_options` | `--geno 0.1 --hwe 1e-13 0.001 --mac 70 --maf 0.01` | Common-variant QC filtering |
| `--plink2_makepgen_4_options` | `--geno 0.2` | Post-QC genotype missingness filter |
| `--plink2_makepgen_5_options` | `--mind 0.2` | Per-sample missingness filter |
| `--plink2_write_snplist_qc_options` | `--maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 --write-samples` | SNP list QC for Regenie step 1 |
| `--plink2_write_snplist_step2_options` | `--mind 0.2` | Additional filters for Regenie step 2 variant list |
| `--plink2_indep_pairwise_window_pca` | `500 50 0.2` | LD pruning window settings for PCA |
| `--plink2_king_cutoff_threshold_pca` | `0.0884` | KING relatedness cutoff for PCA sample selection |
| `--plink2_pca_settings` | `allele-wts 10` | PCA computation settings (number of PCs) |
| `--inbreeding_outliers_range_stds` | `6` | Standard deviations from mean F-coefficient for outlier detection |

### Regenie parameters

| Parameter | Default | Description |
|---|---|---|
| `--regenie_step1_options` | `--bt --bsize 100 --lowmem` | Regenie step 1 options (ridge regression) |
| `--regenie_step2_options` | `--bt --ref-first --firth --approx --bsize 200 --lowmem --aaf-bins 0.01,0.05,0.1,1 --write-mask --write-mask-snplist --vc-tests skato` | Regenie step 2 options (association testing) |
| `--phenotypes_apply_rint` | `false` | Apply rank-inverse normal transformation to phenotypes |

### Reporting and annotation parameters

| Parameter | Default | Description |
|---|---|---|
| `--rscript_annotate_options` | `--min_top_annotations 50 --quantile_threshold 0.25 --include-intergenic FALSE` | Options for the variant annotation R script |
| `--manhattan_annotation_enabled` | `true` | Annotate top hits on Manhattan plots |
| `--annotation_min_log10p` | `1` | Minimum -log10(p) threshold for Manhattan plot annotation |
| `--plot_ylimit` | `0` | Fixed y-axis limit for Manhattan plots (0 = auto) |

## Testing

Requires [nf-test](https://www.nf-test.com/) and Nextflow >=25.10.2.

```bash
# Fast CI suite (tagged "ci"; runs in minutes on prepared fixtures)
nf-test test --profile podman,low_resources --tag ci

# Full suite (includes slow tests: VEP annotation, full workflow with preparation)
nf-test test --profile podman,low_resources
```

The `low_resources` profile caps memory at 7 GB; required for running on development
machines and CI runners.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the
[`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the
[nf-core](https://nf-co.re) community, reused here under the
[MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## Copyright

Copyright (c) 2026 Piotr Suszyński, Tomasz Gambin. All rights reserved.

This software is proprietary. No part of it may be reproduced, distributed, or used
in any form without explicit written permission from the copyright holders.

A licence to use this software free of charge for scientific and non-profit purposes
may be obtained from the authors upon request.
