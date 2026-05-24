# Test Quality and Cleanup Plan (Design)

Date: 2026-05-21

## Goals

1. Stabilize and fix all failing nf-test unit tests.
2. Add missing unit tests for modules with non-trivial logic (Python/R); skip tests for dead code
   and modules that merely wrap an external tool with no custom logic.
3. Design and implement integration tests for subworkflows and the main workflow — these are the
   most valuable tests in this codebase.
4. Remove dead code (unused modules, benchmark workflow, orphan test files).
5. Update `nextflow_schema.json` for all params defined in `nextflow.config` but absent from schema.
6. Remove nf-core template boilerplate comments throughout the codebase.

**Order of work**: Fix failing tests → remove dead code imports + rename dead dirs → write missing
unit tests for complex logic → write integration tests → update schema → (after all integration
tests pass) delete the renamed `__to_delete` directories → remove comments.
Do not write tests for modules confirmed for deletion.

**Testing philosophy**: Module-level tests that simply invoke an external tool (plink, bcftools)
have low value — we are not testing whether those tools work. The high-value tests are:
(a) modules containing complex Python or R logic that we wrote, and
(b) subworkflow and workflow integration tests that verify channel wiring and business logic.
Tests for wrapper-only modules should be marked low-priority and run only in the full suite.
All tests must run with `nf-test test --profile podman`. Stub tests are not used — they verify
nothing meaningful.

---

## 1. Dead Code — Phase 1 (Confirm Now, Remove Before Writing New Tests)

### 1a. Modules with zero references outside the benchmark workflow

Verified by grep — none appear in any workflow or subworkflow other than
`workflows/benchmark_implementations.nf`:

| Module path | Reason |
|---|---|
| `modules/local/vcftools/filter_qual_dp` | No references anywhere |
| `modules/local/bcftools/filter_qual_dp` | No references anywhere |
| `modules/local/combo/filter2` | No references anywhere |
| `modules/local/combo/prepare` | No references anywhere |
| `modules/local/plink19/makebed` | No references anywhere |
| `modules/local/rscript/vcf2aaf` | Benchmark-only |
| `modules/local/bcftools/vcf2aaf` | Benchmark-only |
| `modules/local/rscript/assign_annotations` | Benchmark-only |
| `modules/local/bcftools/assign_annotations/main_shell.nf` | Benchmark-only |

### 1b. Modules imported somewhere but never actually called anywhere

Verified by code inspection — `include` statements exist but no call site is present in any
workflow or subworkflow:

| Module | Imported in | What to do |
|---|---|---|
| `modules/local/plink2/export_bgen` | main workflow | Remove import + rename dir to `__to_delete` suffix |
| `modules/local/bgenix` | main workflow | Remove import + rename dir |
| `modules/local/qctool` | main workflow | Remove import + rename dir |
| `modules/local/cmds/download_file` | main workflow | Remove import + rename dir |
| `modules/local/bcftools/tag2tag` | main workflow + prepare subworkflow | Remove both imports + rename dir (module is never called anywhere) |

### 1c. Stale or dead imports in `subworkflows/local/prepare/main.nf`

Additional imports in the prepare subworkflow that are never called there:

| Import | Situation | What to do |
|---|---|---|
| `BCFTOOLS_VCF2PSAM` | Never called in prepare; module IS called in the main workflow | Remove the import line from prepare/main.nf only — keep the module |
| `BCFTOOLS_VCF2FRQ` | Same | Remove the import line from prepare/main.nf only |
| `VIEW_AND_FILTER2_POLARSBIO` | Never called anywhere; replaced by `BCFTOOLS_VIEW_AND_FILTER2` | Remove import + rename dir |
| `FILTER_AND_ENHANCE_VCF_POLARSBIO` | Never called in prepare; phase-2 deletion candidate | Remove import line now; rename dir in phase 2 |

### 1d. Benchmark workflow

`workflows/benchmark_implementations.nf` is not referenced from `main.nf`. Remove it together
with `conf/benchmark.config` once modules in 1a are deleted.

### 1e. Two-phase removal strategy

Because we might have incorrectly identified some modules as dead, we will not delete them
immediately. The safe approach:

1. **Phase 1 (now)**: Remove all `include` statements for dead modules. Rename the module
   directories by appending `__to_delete` (e.g., `modules/local/bgenix__to_delete`). Delete
   associated test files for those directories (they are definitively broken and reference
   dead code).
2. **Phase 2 (after all integration tests pass)**: Once integration tests confirm nothing is
   missing, permanently delete the `__to_delete` directories and commit.

**Test files to delete in phase 1** (along with removing imports):

| Test file | Reason |
|---|---|
| `modules/local/plink19/makebed/tests/main.nf.test` | Module is dead code (§1a) |
| `modules/local/rscript/vcf2aaf/tests/main.nf.test` | Module is dead code (§1a) |
| `modules/local/rscript/assign_annotations/tests/main.nf.test` | Module is dead code (§1a) |
| `modules/local/plink2/export_bgen/tests/main.nf.test` | Module is dead code (§1b) |
| `modules/local/bgenix/tests/main.nf.test` | Module is dead code (§1b) |
| `modules/local/qctool/tests/main.nf.test` | Module is dead code (§1b) |

### 1f. Dead code — Phase 2 (VCF preparation refactoring, plan but do not remove yet)

Some preparation logic currently in `nf-rare-var-assoc` is being migrated to the upstream
`nf-prepare-vcf` pipeline to avoid duplication. The following modules are candidates for removal
in phase 2, once that refactoring is defined and agreed upon:

- `modules/local/python/filter_and_enhance_vcf_polarsbio`
- `modules/local/python/fix_zero_PL`
- `modules/local/combo/filter_and_enhance_vcf`

Do **not** write tests for these modules — they are likely to be removed.

### 1g. nf-core template comments cleanup (last step)

Many modules still have `// TODO nf-core:` blocks copied from the nf-core template. Examples:
`modules/local/plink2/export_bgen`, `modules/local/qctool`, `modules/local/plink2/makebed`.
Remove all of these in a final cleanup pass after tests are stable.

---

## 2. Test Configuration Profiles

### Purpose of each `conf/*.config` profile

| Profile | Purpose |
|---|---|
| `test.config` | Medium-data VCF, cases/controls input, `use_dosage=true`, full pipeline smoke test. |
| `test_full.config` | Outdated; never used in practice. Candidate for deletion or repurposing. |
| `test_sim_chr22.config` | Simulated chr22 dataset, phenotype file input instead of cases/controls. |
| `test_sim_chr22_2.config` | Simulated chr22 v2, phenotype wildcard pattern, fewer overrides. |
| `test_skip_preparation_and_reporting.config` | Pre-prepared VCF (`skip_preparation=true`, `skip_reporting=true`, `use_dosage=true`). Mirrors the main optimization use-case: `nf-prepare-vcf` runs first, `nf-rare-var-assoc` is invoked by `nextflow-gene-assoc-tuner` for parameter optimization. |
| `benchmark.config` | Only used with `workflows/benchmark_implementations.nf`. Remove together with that workflow. |

### Test data fixture

The canonical fixture for integration tests is in `assets/three_chr_unprepared/`. These files
contain variants from 3 chromosomes (including the X chromosome) with 3202 samples and varying
variant counts:

| File | Variants |
|---|---|
| `unprepared_rand_500.vcf.gz` | ~500 |
| `unprepared_rand_1k.vcf.gz` | ~1 000 |
| `unprepared_rand_2k.vcf.gz` | ~2 000 |
| `unprepared_rand_5k.vcf.gz` | ~5 000 |
| `unprepared_rand_10k.vcf.gz` | ~10 000 |

These are superior to the old chr22-only files for two reasons:
1. Regenie uses LOCO (leave-one-chromosome-out) and requires more than one chromosome.
2. The pipeline has custom X chromosome handling (sex imputation, sex as covariate) that cannot
   be tested with chr22-only data.

Start with `unprepared_rand_500.vcf.gz` for integration tests. If a test tool (particularly
Regenie) complains about low-variance variants or insufficient data, escalate to
`unprepared_rand_1k.vcf.gz` or larger. 3202 samples is fixed; do not create sample-filtered
subsets unless there is a compelling runtime reason.

For subworkflow tests that require pgen/pvar/psam as input (e.g., PCA), record the output of
the prepare subworkflow on the 500-variant file and commit the result to `assets/`.

### Deduplication

- `test_sim_chr22` vs `test_sim_chr22_2`: superseded by `assets/three_chr_unprepared`. Both chr22
  profiles are candidates for deletion once integration tests are in place.
- `test` vs `test_full`: Both use medium_data with cases/controls. Keep `test`; delete or document
  `test_full` with a concrete difference.

---

## 3. Failing nf-test Unit Tests

**Second run** (with `--profile podman`): `nf_test_2026_05_21_14_33_console_out.txt`
93 tests, **67 failed**, 26 passed.

Compared to first run (74 failed): running with podman fixed the 7 "missing tool" failures
(`bcftools/index` 4 tests + `multiqc` 3 tests now pass).

### Category A — Module signature mismatch (test calls wrong number of arguments)

Module code was updated to add or consolidate inputs; tests were not updated.

| Test file | Symptom | Module's actual signature | Fix |
|---|---|---|---|
| `bcftools/annotate/tests` | "declares 2 inputs, called with 3" | `tuple(meta, input, index, annotations, annotations_index, header_lines, rename_chrs)` + `val(out_name_part)` | Update test: merge 3rd arg into tuple |
| `bcftools/norm/tests` | "declares 2 inputs, called with 3" | `tuple(meta, vcf, tbi, fasta, tracking_in)` + `val(out_name_part)` | Update test: `tracking_in` moved into tuple |
| `bcftools/view/tests` | "declares 3 inputs, called with 6" | `tuple(meta, vcf, index, regions, targets, samples, snplist, tracking_in)` + `val(input_args)` + `val(out_name_part)` | Update test: many paths now in the tuple |
| `cmds/merge_results/tests` | "declares 1 input, called with 2" | `tuple(meta, phenotype, regenie_chromosomes, csv_concat_py_script)` | Update test: call with 1 tuple |
| `plink2/makebed/tests` | "declares 5 inputs, called with 3" | `tuple(meta, pgen, pvar, psam, vcf, frq, samples_filtering_file, variants_filtering_file, tracking_in)` + 4× `val(...)` | Update test: call with 5 args |
| `rscript/manhattan_qq_plots/tests` | "declares 1 input, called with 6" | `tuple(meta, phenotype_file, pc_plot_file, eda_plots, setlist_file, phenotype, regenie_merged, mask_file, gwas_report_template, r_functions_file, rmd_pheno_stats_file)` | Update test: call with 1 large tuple |
| `vep/annotate/tests` | "declares 4 inputs, called with 5" | `tuple(meta, input_vcf, input_vcf_tbi, vep_cache)` + `val(species)` + `val(fasta_path)` + `val(input_args)` | Update test: call with 4 args |
| `qctool/tests` | "declares 3 inputs, called with 4" | `tuple(meta, bgen_in, sample_in)` + `val(out_name_part)` + `val(input_args)` | **Delete this test** — module is dead code (never called, see §1b) |

### Category B — Wrong module path or wrong process name in test

| Test file | Symptom | Root cause | Fix |
|---|---|---|---|
| `bcftools/filter/tests` | "Cannot find BCFTOOLS_VIEW in filter/main.nf" | Test imports `BCFTOOLS_VIEW` from the filter module, but that process is named `BCFTOOLS_FILTER` | Change test to import `BCFTOOLS_FILTER` and test the filter module; or delete and consolidate into `bcftools/view/tests` |
| `rscript/assign_annotations/tests` | "Cannot find RSCRIPT_ANNOTATE in assign_annotations/main.nf" | Process is `RSCRIPT_ASSIGN_ANNOTATIONS`; also module is dead code | **Delete this test** — module is dead code (§1a) |
| `rscript/buildreports/tests` | "Cannot find RSCRIPT_ANNOTATE in assign_annotations/main.nf" | Test includes `rscript/assign_annotations` looking for `RSCRIPT_ANNOTATE` — a broken transitive dependency on dead code | Remove the broken include; update test to exercise only `RSCRIPT_BUILDREPORTS` |
| `regenie/step1/tests` | "Can't find rscript/annotate/main.nf" | Test includes `rscript/annotate/main.nf` (old name); module was renamed to `rscript/assign_annotations` but then replaced by `bcftools/assign_annotations` in production | Update test include to use `bcftools/assign_annotations` — after benchmarks, `BCFTOOLS_ASSIGN_ANNOTATIONS` replaced `RSCRIPT_ASSIGN_ANNOTATIONS` in all production code |
| `regenie/step2/tests` | Same as step1 | Same root cause | Same fix |

### Category C — Missing executables (fixed by running with `--profile podman`)

Running `nf-test test --profile podman` is **required** for all tests. These failures were
resolved in the second run by switching to podman:

| Test file | Was failing because | Status |
|---|---|---|
| `bcftools/index/tests` | `bcftools` not in PATH | Now passes (4/4) |
| `bcftools/vcf2frq/tests` | `bcftools` not in PATH | Still failing — snapshot mismatch (see Cat E) |
| `nf-core/multiqc/tests` | `multiqc` not in PATH | Now passes (3/3) |

### Category D — Null path / dead code (delete tests, not fix)

| Test file | Symptom | Decision |
|---|---|---|
| `plink19/makebed/tests` | "Path value cannot be null" | **Delete** — module is dead code (§1a) |
| `plink2/export_bgen/tests` | "Path value cannot be null" | **Delete** — module is dead code (§1b) |
| `bgenix/tests` | Transitive failure from PLINK2_EXPORT_BGEN | **Delete** — module is dead code (§1b) |
| `qctool/tests` | "declares 3 inputs, called with 4" | **Delete** — module is dead code (§1b) |
| `plink2/write_snplist/tests` | Transitive failure from PLINK2_MAKEBED null path | Fix PLINK2_MAKEBED test first; once that passes, this should pass too |

### Category E — Snapshot md5 mismatch

| Test file | What changed | Fix |
|---|---|---|
| `bcftools/vcf2frq/tests` | `versions.yml` md5 and output file content differ from snapshot | Re-record: `nf-test test --profile podman --update-snapshot modules/local/bcftools/vcf2frq/tests/main.nf.test` |

---

## 4. Passing Tests After Second Run (26/93)

- `bcftools/index` — 4 tests (2 real, 2 stub): now pass with podman
- `python/merge_sex_covar` — 1 test: passes
- `nf-core/multiqc` — 3 tests (2 real, 1 stub): now pass with podman
- All 18 nf-core utility subworkflow tests: pass

---

## 5. Missing Module Tests

### 5a. High-priority: modules with complex Python/R logic (write these)

| Module | What to test |
|---|---|
| `python/eda` | EDA plots: assert HTML/PNG files emitted; edge case: empty phenotype input |
| `python/draw_pc_plot` | PCA scatter plot: assert output plot file exists and is non-empty |
| `python/generate_tracking_report` | Tracking JSON → HTML: assert output file exists, verify a known key in the report |
| `python/calc_f_outliers` | F-coefficient outliers: verify output has correct columns; edge: no outliers present |
| `python/vcf2aaf` | VCF to AAF: verify output TSV columns and row counts against known input |
| `bcftools/assign_annotations` | Annotation from a small INFO TSV: assert annotated VCF has new INFO fields |

### 5b. Low-priority: wrapper-only modules (only run in full suite; write stubs if at all)

These modules contain little or no custom logic — they invoke an external tool with configured
arguments. Testing them has low value; if tests exist they should be marked with `tag "full"`.

- `plink2/makepgen`, `plink2/indep_pairwise`, `plink2/king_cutoff`, `plink2/pca`,
  `plink2/projection_score`, `plink2/het`, `plink2/import_dosage`, `plink2/export_other`,
  `plink19/makeset`, `bcftools/replace_sample_names`, `bcftools/vcf2psam`,
  `bcftools/view_and_filter2`, `cmds/check_x_chrom_present`, `cmds/extract_phenotypes_and_samples`,
  `cmds/rename`, `rscript/build_phenotypes`

### 5c. Skip writing tests

- `python/filter_and_enhance_vcf_polarsbio` — likely removed in phase 2 refactoring
- `python/fix_zero_PL` — likely removed in phase 2 refactoring
- `combo/filter_and_enhance_vcf` — likely removed in phase 2 refactoring
- `python/view_and_filter2_polarsbio` — dead code (replaced by bcftools version, §1c)
- `vep/updatecache` — downloads external data only, no custom logic; tested manually

---

## 6. Integration Test Design

**This is the most important section.** Subworkflow and workflow tests exercise channel wiring,
business logic, and data transformations that are impossible to cover with module-level tests.

All integration tests must:
- Run with `--profile podman`
- Use real (non-stub) execution
- Use `assets/three_chr_unprepared/unprepared_rand_500.vcf.gz` as the primary fixture
  (3202 samples, ~500 variants across 3 chromosomes including X)
- Have assertions focused on business outcomes, not just `process.success`

### IT-1: `subworkflows/local/prepare`

Input: `unprepared_rand_500.vcf.gz` + optional sample rename file

Assertions:
- Output VCF has no multiallelic sites (no comma in ALT field)
- Variant IDs in the ID column are unique
- VEP annotation fields (CSQ) are present in INFO
- Output VCF is bgzipped + tabix-indexed
- `tracking_out` JSON contains preparation step keys with non-zero variant counts

### IT-2: `subworkflows/local/pca`

Input: pgen/pvar/psam recorded from IT-1 output (commit to `assets/`) or generated from the
500-variant VCF.

Assertions:
- `.eigenvec` has expected column count (sample ID + N PC columns per `params.plink2_pca_settings`)
- `.eigenval` values are in descending order
- PCA plot HTML file is non-empty
- `tracking_out` captures sample counts before/after kinship filtering (when enabled)

### IT-3: `subworkflows/local/filter_missing_per_pheno`

Input: pgen/pvar/psam (from IT-2 assets) + phenotype file with two phenotypes, some samples
present in only one.

Assertions:
- Each phenotype produces its own filtered pgen/pvar/psam
- Sample counts differ between phenotypes (phenotype-specific exclusion works)
- `tracking_out` captures per-phenotype missingness

### IT-4: `subworkflows/local/f_coefficient_filtering`

Input: pgen/pvar/psam with at least 2 artificial outlier samples inserted.

Assertions:
- Outlier samples appear in the exclusion output file
- Output pgen/pvar/psam has fewer samples than input

### IT-5: `subworkflows/local/reporting`

Input: pre-computed Regenie output files + small phenotype/annotation files.

Assertions:
- HTML report file emitted and non-empty
- Manhattan/QQ plot HTML contains expected section headers

### IT-6: `workflows/rare-var-assoc` — fast path

Config: `skip_preparation=true`, `skip_reporting=true`, `use_dosage=true`
Input: pre-prepared VCF + psam + phenotype file + masks

Assertions:
- Regenie step 1 produces ridge model files
- Regenie step 2 produces `.regenie.gz` per phenotype
- `MERGE_RESULTS` output has expected columns (CHROM, GENPOS, ID, …)
- MultiQC HTML report emitted
- `tracking_out` JSON files present in `outdir/pipeline_info`

Note: this is the main parameter-optimization scenario (mirrors
`test_skip_preparation_and_reporting.config`).

### IT-7: `workflows/rare-var-assoc` — full reporting path

Config: `skip_preparation=false`, `skip_reporting=false`, `use_dosage=false`
Input: `unprepared_rand_500.vcf.gz` + phenotype file

Assertions:
- All IT-6 assertions apply
- VEP annotation applied (CSQ fields in VCF before Regenie)
- EDA plots present in output
- GWAS HTML reports present

Note: this is the "production" use case and will be the slowest test. Run in full suite only.

---

## 7. Runtime and Test Organisation

- All tests run with `nf-test test --profile podman`. This is mandatory.
- No stub tests. Stubs verify nothing meaningful.
- No `--shard` parallelism (caused problems in prior attempts).
- Tag tests: `tag "ci"` for fast integration tests (IT-1 through IT-6) and high-value module
  tests; `tag "full"` for slow tests (IT-7) and low-priority wrapper module tests. CI runs only
  `tag "ci"`.
- Shared small fixture: use `assets/three_chr_unprepared/unprepared_rand_500.vcf.gz` as the
  single canonical test input across all tests that need a VCF. Avoid test-specific fixture files.

---

## 8. nextflow_schema.json Updates

When the workflow runs, Nextflow emits `WARN: The following invalid input values have been
detected` for all params not registered in the schema. The full list of missing params
(confirmed by cross-referencing `nextflow.config` against the schema) is below.

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
| `bcftools_view_1_options` | `"--output-type z --write-index=tbi"` | string | bcftools view options for initial sample subsetting step |
| `bcftools_replace_sample_names_sed_arg` | `"s/_/-/g"` | string | sed expression for sample name normalisation |
| `filter_and_enhance_vcf_qual_min` | `"25"` | string | Minimum QUAL score |
| `filter_and_enhance_vcf_avg_gq_min` | `"25"` | string | Minimum average GQ across samples |
| `filter_and_enhance_vcf_avg_dp_min` | `"25"` | string | Minimum average DP across samples |
| `filter_and_enhance_vcf_avg_dp_max` | `"200"` | string | Maximum average DP across samples |
| `filter_and_enhance_vcf_sample_gq_min` | `"20"` | string | Minimum per-sample GQ (set GT to missing below threshold) |
| `filter_and_enhance_vcf_sample_dp_min` | `"20"` | string | Minimum per-sample DP |
| `filter_and_enhance_vcf_sample_dp_max` | `"250"` | string | Maximum per-sample DP |
| `filter_and_enhance_vcf_calc_ds_min_gq` | `"1"` | string | Minimum GQ for dosage calculation |

**Tool script paths — phase-2 deletion candidates** (add temporarily; remove when modules deleted):

| Param | Default | Type | Description |
|---|---|---|---|
| `view_and_filter2_polarsbio_script` | `"${projectDir}/modules/local/python/view_and_filter2_polarsbio/assets/filter.py"` | string | Path to polars-based view+filter script |
| `filter_and_enhance_vcf_polarsbio_script` | `"${projectDir}/modules/local/python/filter_and_enhance_vcf_polarsbio/assets/filter.py"` | string | Path to polars-based filter+enhance script |

**PLINK2 makepgen options** (add to an existing "PLINK2 options" group or create one):

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

## 9. Execution Plan

| # | Task | Deliverable / done-when |
|---|---|---|
| T1 | Confirm dead code list (§1a–§1c); get sign-off | This doc approved |
| T2 | Phase 1 removal: remove all dead `include` lines; rename dead module dirs to `__to_delete`; delete associated test files (§1e) | No more dead includes; `nf-test` count drops; pipeline still runs |
| T3 | Fix Cat-B: wrong process name / path in tests | Cat-B tests pass |
| T4 | Fix Cat-A: update test call signatures to match current module signatures | Cat-A tests pass |
| T5 | Fix Cat-D (plink2/write_snplist transitive): provide required inputs | Passes |
| T6 | Fix Cat-E (bcftools/vcf2frq snapshot): re-record snapshot with podman | Passes |
| T7 | Write missing unit tests for high-priority python/R modules (§5a) | Each new test passes on real run with podman |
| T8 | Implement IT-1 and IT-2 (prepare + PCA subworkflows); commit fixture pgen assets | Both pass |
| T9 | Implement IT-3 and IT-4 (filter_missing_per_pheno + f_coefficient_filtering) | Both pass |
| T10 | Implement IT-5 (reporting subworkflow) | Passes |
| T11 | Implement IT-6 (full workflow, fast path) | Passes |
| T12 | Implement IT-7 (full workflow, reporting path) | Passes |
| T13 | Phase 2 removal: permanently delete `__to_delete` directories | No `__to_delete` dirs remain; all tests still pass |
| T14 | Update `nextflow_schema.json` (§8) | No `WARN: invalid input values` from Nextflow on a standard run |
| T15 | Remove nf-core `// TODO` template comments from `modules/local/` | No TODO blocks remain |
