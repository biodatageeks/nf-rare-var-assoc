# Test Config Profiles + nextflow_schema.json Updates (T15)

See [../test-quality-and-cleanup-plan.md](../test-quality-and-cleanup-plan.md) for goals,
conventions, and the full task list.

---

## §2 — Test configuration profiles

| Profile | Purpose |
|---|---|
| `test.config` | Medium-data VCF, cases/controls input, `use_dosage=true`, full pipeline smoke test. |
| `test_full.config` | Outdated; never used in practice. Delete or document a concrete difference vs `test.config`. |
| `test_sim_chr22.config` | Simulated chr22, phenotype-file input. Superseded by `three_chr_unprepared`; delete once IT tests are in place. |
| `test_sim_chr22_2.config` | Simulated chr22 v2, phenotype wildcard. Same disposition as above. |
| `test_skip_preparation_and_reporting.config` | Pre-prepared VCF (`skip_preparation=true`, `skip_reporting=true`, `use_dosage=true`). The main parameter-optimization scenario: `nf-prepare-vcf` runs first, `nf-rare-var-assoc` is then invoked by `nextflow-gene-assoc-tuner`. |
| `benchmark.config` | Only used with `workflows/benchmark_implementations.nf`. Remove together with that workflow (T14). |

The canonical fixture for integration tests is defined in §0.4 of the main plan
(`assets/three_chr_unprepared/`). Briefly: 3 chromosomes including X, 3202 samples, and
500/1k/2k/5k/10k variant variants. The X chromosome and >1-chromosome shape are required
because Regenie uses LOCO and the pipeline has custom X handling (sex imputation, sex as
covariate). Start with `unprepared_rand_500.vcf.gz`; escalate to larger sizes only if a
tool complains about low variance.

For subworkflow tests that need pgen/pvar/psam as input, T8 records the prepare-subworkflow
output into `assets/three_chr_unprepared/prepared_500/`.

---

## §8 — `nextflow_schema.json` updates

When the workflow runs, Nextflow emits
`WARN: The following invalid input values have been detected` for every param defined in
`nextflow.config` but missing from the schema. The full list of missing params is below.

**Pipeline behaviour flags** (add to a "Pipeline options" group):

| Param | Default | Type | Description |
|---|---|---|---|
| `skip_preparation` | `false` | boolean | Skip VCF preparation subworkflow (assumes input is already prepared) |
| `skip_reporting` | `false` | boolean | Skip HTML report generation |
| `use_dosage` | `false` | boolean | Use dosage (DS field) instead of hard genotype calls in Regenie |
| `publish_intermediate` | `false` | boolean | Publish intermediate files to outdir |
| `regenie_step1_kinship_filtering` | `false` | boolean | Enable kinship-based sample filtering in Regenie step 1 |
| `tmpdir` | `null` | string | Temporary directory override |
| `errorStrategy` | *(none)* | string | Override process error strategy (e.g. `retryThenIgnore`, `terminate`) |
| `help_full` | `false` | boolean | Show full help including hidden params |
| `show_hidden` | `false` | boolean | Show hidden params in help output |

**VCF filtering parameters** (add to a "VCF preparation" group):

| Param | Default | Type | Description |
|---|---|---|---|
| `bcftools_view_1_options` | `"--output-type z --write-index=tbi"` | string | bcftools view options for initial sample-subsetting step |
| `bcftools_replace_sample_names_sed_arg` | `"s/_/-/g"` | string | sed expression for sample-name normalisation |
| `filter_and_enhance_vcf_qual_min` | `"25"` | string | Minimum QUAL score |
| `filter_and_enhance_vcf_avg_gq_min` | `"25"` | string | Minimum average GQ across samples |
| `filter_and_enhance_vcf_avg_dp_min` | `"25"` | string | Minimum average DP across samples |
| `filter_and_enhance_vcf_avg_dp_max` | `"200"` | string | Maximum average DP across samples |
| `filter_and_enhance_vcf_sample_gq_min` | `"20"` | string | Minimum per-sample GQ (set GT to missing below threshold) |
| `filter_and_enhance_vcf_sample_dp_min` | `"20"` | string | Minimum per-sample DP |
| `filter_and_enhance_vcf_sample_dp_max` | `"250"` | string | Maximum per-sample DP |
| `filter_and_enhance_vcf_calc_ds_min_gq` | `"1"` | string | Minimum GQ for dosage calculation |

**Tool script paths — phase-2 deletion candidates** (mark `hidden: true`; remove when modules deleted):

| Param | Default | Type | Description |
|---|---|---|---|
| `view_and_filter2_polarsbio_script` | `"${projectDir}/modules/local/python/view_and_filter2_polarsbio/assets/filter.py"` | string | Path to polars-based view+filter script |
| `filter_and_enhance_vcf_polarsbio_script` | `"${projectDir}/modules/local/python/filter_and_enhance_vcf_polarsbio/assets/filter.py"` | string | Path to polars-based filter+enhance script |

**PLINK2 makepgen options** (add under "PLINK2 options"):

| Param | Default | Type | Description |
|---|---|---|---|
| `plink2_makepgen_1_options` | `"--double-id --vcf-half-call missing --split-par b38 --1"` | string | Options for initial VCF→pgen import |
| `plink2_makepgen_1_vcf_input_options` | `""` | string | Additional VCF input options for pgen import step 1 |
| `plink2_makepgen_2_options` | `"--impute-sex max-female-xf=0.2 min-male-xf=0.8"` | string | Options for sex imputation step |
| `plink2_makepgen_3_options` | `"--geno 0.1 --hwe 1e-13 0.001 --mac 70 --maf 0.01"` | string | Options for common-variant QC filtering |
| `plink2_makepgen_4_options` | `"--geno 0.2"` | string | Options for step-4 pgen filter |
| `plink2_makepgen_5_options` | `"--mind 0.2"` | string | Options for per-sample missingness filter |
| `plink19_makeset_options` | `"--allow-extra-chr"` | string | Options for PLINK1.9 makeset step in PCA subworkflow |

---

## T15 — Update `nextflow_schema.json`

- Add every param above to `nextflow_schema.json` under the appropriate group
  (`pipeline_options`, `vcf_preparation`, `plink2_options`).
- For each param, set: `type`, `default`, `description`, and `hidden: true` for the
  deletion-candidate script paths.
- **Verify**:
  `nextflow run main.nf -profile test_skip_preparation_and_reporting,podman --outdir /tmp/_schema_check -preview`
  (or any quick invocation) and confirm there is no
  `WARN: The following invalid input values have been detected:` line in `.nextflow.log`.
- **Done-when**: zero schema warnings on a standard run.
