# Output

Everything the pipeline writes, and how to read it. All paths are relative to
`--outdir`.

Most directories appear only when `--publish_intermediate true` is set. The tables below
say which is which. When several phenotype files are analysed in one run, each dataset's
files are prefixed with `<project_name>_dataset_idx_<N>_`; otherwise the prefix is
`<project_name>_`.

## Written by default

### `rscript_manhattan_qq_plots/` -- the main report

A single self-contained HTML file per phenotype. This is the place to start. It contains:

- a Manhattan plot of the gene-level results, with significant genes labelled (controlled
  by `--manhattan_annotation_enabled` and `--annotation_min_log10p`),
- a QQ plot with the genomic inflation factor, which shows whether the test is well
  calibrated,
- a table of the strongest associations,
- summary statistics for the phenotype,
- the principal-component plot,
- the exploratory data-quality figures,
- the parameters the run used.

### `regenie_step2/` -- the association results

| File | Contents |
|---|---|
| `*_step2_<phenotype>.regenie` | the association results |
| `*_step2_masks.snplist` | which variants ended up in each gene and impact group |
| `*_step2.log` | REGENIE's own log |

The results file is space-separated, with a header preceded by comment lines recording
the impact-group definitions. One row per gene, impact group, allele-frequency threshold
and test:

```
CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
12 125654 IQSEC3.Mask_High.0.2 ref Mask_High.0.2 0.0148267 1699 ADD 0.00478484 0.372464 0.000165031 0.00447436 NA
```

| Column | Meaning |
|---|---|
| `ID` | `<gene>.<impact group>.<allele-frequency threshold>` |
| `A1FREQ` | combined frequency of the variants in that group |
| `N` | samples analysed |
| `TEST` | `ADD` for the burden test, `ADD-SKATO` for SKAT-O, and so on -- one row each |
| `BETA`, `SE` | effect size and its standard error (burden tests only; SKAT-O has no single effect size) |
| `LOG10P` | **-log10 of the p-value**. This is the column to sort on -- larger means more significant. |

A gene appears many times, once per combination of impact group, frequency threshold and
test, so filter to the combination you care about before interpreting anything.

### `rscript_buildreports/` -- per-variant carrier tables

For every variant that entered a gene group, the allele counts and frequencies in cases
and in controls:

| File | Contents |
|---|---|
| `*_annotated_snps.csv` | per variant: allele count, allele frequency and number of missing genotypes, separately for cases and controls |
| `*_annotated_snps_with_sample_ids.csv` | the same, plus the identifiers of the samples carrying each variant, as semicolon-separated lists of heterozygous and homozygous carriers |
| `*_res_log10p_1_annotated.csv` | for every gene above -log10(p) = 1, its association result joined to the per-variant table and to each variant's VEP consequence |

The last file is the one to open after finding an interesting gene in the results: it
shows which variants drove the signal, how often each occurs in cases and in controls,
and what VEP predicted them to do.

### `generate_tracking_report/` -- how the data flowed through the pipeline

| File | Contents |
|---|---|
| `*_sankey_report.html` | a diagram of samples and variants entering and leaving each step |
| `*_pipeline_report.txt` | the same as plain text, one block per step |

Every step records how many samples and variants it received and returned, together with
the parameters it used. This is the fastest way to find where data was lost:

```
Process: BCFTOOLS_VIEW_AND_FILTER2_..._viewfilter2
Input Variants: 141706
Output Variants: 110086
Input Samples: 3202
Output Samples: 1368
Parameters:   --samples-file samples_transformed.txt  --output-type z --write-index=tbi
```

A count of `-1` means the step does not report that number, not that everything was
removed.

### `multiqc/`

`multiqc_report.html` combines the quality output of the individual tools into one page,
with the parsed numbers in `multiqc_data/`.

### `pipeline_info/`

Nextflow's own records: `execution_report_*.html` (per-task runtime and memory),
`execution_timeline_*.html`, `execution_trace_*.txt` and `pipeline_dag_*.html`. Use the
trace file to find which steps are slow or run out of memory.

## Written only with `--publish_intermediate true`

### `exploratory_data_analysis/plots/`

Roughly fifty figures, as both PNG and SVG, describing the data before association
testing. They are also embedded in the main HTML report, so you only need this directory
if you want the image files themselves. They cover:

- **genotype quality, read depth and dosage** -- distributions per variant and per sample,
  plotted separately for cases and controls (files beginning `1_`, `3_`, `5_`, `16_`),
- **missingness** -- per variant and per sample, by phenotype (`7_`, `9_`),
- **case/control differences** in read depth (`10_`), which is a common source of false
  positives,
- **allele frequency, variant types and variant density per chromosome** (`11_`, `12_`,
  `13_`),
- **heterozygosity per sample** (`15_`), where outliers indicate contamination or poor
  quality,
- **relationships between measures** -- depth against quality, dosage against depth,
  dosage against genotype (`17_`), which show whether the derived dosages behave sensibly.

Files with a `b` suffix (`1b_`, `3b_`, `10b_`) are the same figure on a different scale.

### Association inputs and population structure

| Directory | Contents |
|---|---|
| `regenie_step1/` | the whole-genome model REGENIE fits before the association tests |
| `plink2_pca/` | principal components |
| `plink2_projection_score/` | every sample's score on those components, used as covariates |
| `draw_pc_plot/` | the principal-component plot |
| `plink2_king_cutoff/` | which samples were used to compute the components |
| `plink2_indep_pairwise/` | the linkage-disequilibrium-pruned variant sets |
| `calculate_f_outliers/` | inbreeding coefficients and the samples flagged as outliers |
| `merge_results/` | association results merged across chunks, one file per phenotype |

### Gene groupings

`bcftools/` holds the per-dataset `.annotations`, `.setlist` and `.aaf` files: which
variants belong to which gene, what VEP predicted for each, and their allele frequencies.
These define what the association test actually tested.

### Preparation and quality control

| Directory | Contents |
|---|---|
| `prepare/`, `bcftools_norm/`, `bcftools_annotate/`, `vep_annotate/`, `vep_updatecache/`, `fix_zero_pl/` | VCF preparation: normalisation, variant identifiers, VEP annotation, dosage computation |
| `bcftools_view_and_filter2/` | the quality-filtered VCF |
| `bcftools_replace_sample_names/`, `bcftools_index/`, `bcftools_vcf2frq/` | sample renaming, indexing, allele-frequency extraction |
| `plink2_makepgen/`, `plink2_makebed/`, `plink2_write_snplist/`, `plink2_export_other/`, `plink2_import_dosage/` | genotype format conversions and the quality-controlled variant sets |
| `plink19_makeset/`, `plink19_makebed/` | PLINK 1.9 conversions used by the principal-component step |
| `rscript_build_phenotypes/`, `phenotype/` | the phenotype file built from the case and control lists |
| `rscript_assign_annotations/`, `rscript_vcftoaaf/` | intermediate annotation and allele-frequency files |
| `check_x_chrom_present/`, `extract/`, `rename/`, `download_file/` | small bookkeeping steps |

Every directory also contains a `*_tracking.json` file with the sample and variant counts
that the tracking report is built from.

## A note on interpreting results

The Manhattan and QQ plots in the main report should be read together. A QQ plot whose
points lift away from the diagonal along their whole length, and a genomic inflation
factor well above 1, mean the test statistics are inflated -- usually because of
uncorrected population structure -- and the p-values cannot be taken at face value. Check
the principal-component plot for unexpected structure, and consider including more
components through `--covarColList` in `--regenie_step2_options`.

Very small case counts have the opposite effect: burden statistics become discrete and
conservative, and the inflation factor drops below 1. That is expected rather than a
problem, but it does mean the analysis has little power.
