# Comparison against nf-gwas

Scripts that run this pipeline and [genepi/nf-gwas](https://github.com/genepi/nf-gwas)
on the same simulated datasets and compare the results.

In this comparison nf-gwas is given the gene groupings and principal components computed
by this pipeline, because nf-gwas cannot produce them itself. The measured difference
therefore reflects only the earlier stages -- quality control and the use of dosages
derived from genotype likelihoods -- and not the two tools as a whole. The comparisons in
[../chained-benchmark/](../chained-benchmark/) and
[../chained-benchmark-nf-gwas/](../chained-benchmark-nf-gwas/) exist because of this
limitation: there, each tool is given only the raw data and does everything else itself.

These scripts are kept unchanged so their result stays reproducible. The
generalised versions are in [../benchmark-common/](../benchmark-common/).

## Input data

- **Raw VCF**: `…/tools_comparison/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz`
  -- exome data for chromosomes 12, 22 and X, genome build GRCh38, 3202 samples. The
  FORMAT field contains GT, GQ, DP, PL and AD, but no DS. This is what nf-gwas reads.
- **Prepared VCF**: `…/tools_comparison/prepared.vcf.gz` -- the output of nf-prepare-vcf
  (multi-allelic sites split, VEP consequences added, dosages computed). This is what
  this pipeline reads, run with `--skip_preparation true` so that the slow preparation
  step is not repeated for every dataset.
- **Simulated datasets**: `…/tools_comparison/datasets/run_<N>/`, where `N` is the
  dataset number. Each contains `gcta_simu/*_gcta_simu.phenotype.txt` (columns
  `FID IID Y1`, case/control) and the known-answer files in `select_genes/`.

To run only some datasets, set `DATASET_IDXS="4 18 27"`. The default is every `run_<N>`
directory found.

## Why nf-gwas reads the raw VCF

The point of the comparison is that nf-gwas does not normalise the data and does not
compute dosages from genotype likelihoods. It receives only the multi-allelic split and
variant-ID assignment that it cannot perform itself and that the shared gene groupings
depend on. Because the raw VCF has no DS field, plink2's `dosage=DS` option falls back to
hard genotype calls without warning -- so nf-gwas works from hard calls, while this
pipeline uses dosages. That is the main real difference between the two, not an artefact
of these scripts. Giving nf-gwas the prepared VCF would hand it this pipeline's
normalisation and dosages for free.

## How multiple datasets are handled

The two pipelines accept multiple datasets differently, so the scripts treat them
differently:

- **This pipeline** accepts a comma-separated list of phenotype files and combines each
  with the VCF, so all datasets run in **one** `nextflow run`, each identified as
  `tools_comparison_dataset_idx_<N>`.
- **nf-gwas** has no equivalent option -- `phenotypes_columns` refers to columns of a
  single file -- and produces a single `Y1.regenie`. It is therefore run **once per
  dataset** in a loop, and each result file is renamed to include `_dataset_idx_<N>_`.
- **nf-eval-gene-assoc**, which scores the results, reads the dataset number from
  filenames and matches the known answers to the results by that number. It is run once
  per pipeline, using file patterns that match all datasets.

## Scripts

| Script | What it does |
|---|---|
| `run_nf_rare_var_assoc.sh` | Runs this pipeline over all datasets in one `nextflow run`, using the prepared VCF and `--skip_preparation`. `--publish_intermediate true` makes each dataset's gene groupings and principal-component covariates available to the nf-gwas run. Deletes large intermediate directories afterwards, then scores the results. **Run this first.** |
| `run_nf_gwas.sh` | For each dataset in turn: prepares input with bcftools and plink2 in a container (sample filtering, multi-allelic split, `CHROM_POS_REF_ALT` variant IDs, bed/bim/fam files, covariates), runs nf-gwas in gene-based mode reusing that dataset's gene groupings and covariates, and renames the REGENIE result. Cleans up each dataset's working directory, then scores all results. **Run this second.** |
| `pairwise_compare.py` | Compares the two sets of scores. For each dataset it reads `*_auc_summary.csv`, scales each measure by the fraction of causal genes that appear in the results at all, pairs the two pipelines by dataset number, and reports the mean difference in average precision (and in AUC-PR and AUC-ROC) with standard deviation, 95% confidence interval, and paired t-test and Wilcoxon p-values. Writes `runs/pairwise_comparison/pairwise_per_dataset.csv` and `pairwise_summary.csv`. Requires pandas and scipy. **Run this third.** |

```bash
./run_nf_rare_var_assoc.sh    # all datasets in one run, plus reusable groupings and covariates
./run_nf_gwas.sh              # per-dataset preparation and nf-gwas, reusing those files
python pairwise_compare.py    # compare the two sets of scores

# only some datasets:      DATASET_IDXS="4 18 27" ./run_nf_rare_var_assoc.sh
# keep the intermediates:  CLEANUP=false ./run_nf_rare_var_assoc.sh
```

Results are written under `…/tools_comparison/runs/{nf_rare_var_assoc,nf_gwas}/`
(`results/`, `work/` and Nextflow trace files), and the scores under
`…/runs/{…}_eval/results/compute_score/`, one set of files per dataset.

## Settings chosen to keep the comparison fair

- **Statistical test**: both pipelines run REGENIE's gene-based SKAT-O test on the *same*
  gene sets. nf-gwas cannot build gene groupings from VEP consequences, so it reuses this
  pipeline's `.annotations` and `.setlist` files; the group definitions come from the
  shared `assets/default.masks`.
- **Covariates**: nf-gwas reuses this pipeline's `PC1_AVG..PC4_AVG`, since it has no
  principal-component analysis of its own. Sex is deliberately excluded, so that both
  pipelines use exactly the same covariates.
- **Quality-control settings**: the mapping is documented in the comments in
  `run_nf_gwas.sh`. nf-gwas has a single quality-control step, configured from this
  pipeline's `plink2_makepgen_3_options` plus `--mind` from
  `plink2_write_snplist_qc_options`. Steps that nf-gwas has no equivalent for --
  inbreeding-coefficient filtering, relatedness estimation, sex imputation,
  per-phenotype missingness filtering, separate settings for the two REGENIE steps, and a
  custom allele-frequency file -- are left at nf-gwas defaults.
- **Allele-frequency groups**: both use a single `--aaf-bins 0.2` (in nf-gwas,
  `regenie_gene_aaf: 0.2`), so the gene-based test is set up identically in both.

## Practical notes

- `plink2 --make-pgen` keeps the ID column from the VCF. Setting `CHROM_POS_REF_ALT` IDs
  during preparation is therefore what makes the shared gene groupings usable. The
  multi-allelic split and ID assignment are **required**: without them REGENIE matches no
  variants at all. Removing duplicate records remains optional (`REMOVE_DUPS`).
- With `--skip_preparation` this pipeline does not run VEP, because the prepared VCF
  already contains the consequence annotations. The scoring step also runs no VEP, since
  it works only from the known-answer files and the REGENIE results. No VEP cache is
  needed for this comparison.
- `prepared.vcf.gz` must be the nf-prepare-vcf output **for the same raw VCF**, so that
  the variant IDs match the ones the gene groupings refer to. Otherwise REGENIE matches
  no variants.
- The nf-gwas source at
  `/data/git/doktorat_pw/nf_rare_var_assoc__tools_comparison/genepi/nf-gwas` is v1.0.11,
  but only the v1.0.9 container is available locally, which is what the script uses by
  default. To use v1.0.11, pull that container and set `NFGWAS_CONTAINER`.
- nf-gwas has no podman profile, so the script enables podman through an inline `-c`
  configuration file.
- Rootless podman on this machine needs `--userns=keep-id` and `:z` on mounts before
  containers can write into host-owned directories. Both the manual preparation commands
  and the Nextflow configuration files set this.
