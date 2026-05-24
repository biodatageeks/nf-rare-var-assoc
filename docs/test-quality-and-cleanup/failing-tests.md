# Failing nf-test Unit Tests (T3, T4, T5, T6)

See [../test-quality-and-cleanup-plan.md](../test-quality-and-cleanup-plan.md) for goals,
conventions (including the **quality bar for assertions**), and the full task list.

**Second run** (with `--profile podman`): `nf_test_2026_05_21_14_33_console_out.txt` —
93 tests, **67 failed**, 26 passed. Switching from local to podman fixed the 7 "missing tool"
failures (`bcftools/index` 4 + `multiqc` 3).

**Important when fixing**: signature/path fixes that get a test to compile are not the whole
job. Once the test runs, audit its existing assertions against the **quality bar for
assertions** in the main plan. If the only assertion is `process.success` or "output file
exists", add at least one business-meaningful assertion before declaring the task done.

---

## §3 Category A — Module signature mismatch

Module code was updated to add or consolidate inputs; tests were not updated.

| Test file | Symptom | Module's actual signature | Fix |
|---|---|---|---|
| `modules/local/bcftools/annotate/tests/main.nf.test` | "declares 2 inputs, called with 3" | `tuple(meta, input, index, annotations, annotations_index, header_lines, rename_chrs)` + `val(out_name_part)` | Merge 3rd arg into tuple |
| `modules/local/bcftools/norm/tests/main.nf.test` | "declares 2 inputs, called with 3" | `tuple(meta, vcf, tbi, fasta, tracking_in)` + `val(out_name_part)` | `tracking_in` moved into tuple |
| `modules/local/bcftools/view/tests/main.nf.test` | "declares 3 inputs, called with 6" | `tuple(meta, vcf, index, regions, targets, samples, snplist, tracking_in)` + `val(input_args)` + `val(out_name_part)` | Many paths now in the tuple |
| `modules/local/cmds/merge_results/tests/main.nf.test` | "declares 1 input, called with 2" | `tuple(meta, phenotype, regenie_chromosomes, csv_concat_py_script)` | Call with 1 tuple |
| `modules/local/plink2/makebed/tests/main.nf.test` | "declares 5 inputs, called with 3" | `tuple(meta, pgen, pvar, psam, vcf, frq, samples_filtering_file, variants_filtering_file, tracking_in)` + 4× `val(...)` | Call with 5 args |
| `modules/local/rscript/manhattan_qq_plots/tests/main.nf.test` | "declares 1 input, called with 6" | `tuple(meta, phenotype_file, pc_plot_file, eda_plots, setlist_file, phenotype, regenie_merged, mask_file, gwas_report_template, r_functions_file, rmd_pheno_stats_file)` | Call with 1 large tuple |
| `modules/local/vep/annotate/tests/main.nf.test` | "declares 4 inputs, called with 5" | `tuple(meta, input_vcf, input_vcf_tbi, vep_cache)` + `val(species)` + `val(fasta_path)` + `val(input_args)` | Call with 4 args |
| `modules/local/qctool/tests/main.nf.test` | "declares 3 inputs, called with 4" | `tuple(meta, bgen_in, sample_in)` + `val(out_name_part)` + `val(input_args)` | **Delete this test file** — module is dead code (§1b) |

Note: the source of truth for each signature is the `input:` block in
`modules/local/<group>/<name>/main.nf`, not this table.

## §3 Category B — Wrong module path / process name

| Test file | Symptom | Root cause | Fix |
|---|---|---|---|
| `modules/local/bcftools/filter/tests/main.nf.test` | "Cannot find BCFTOOLS_VIEW in filter/main.nf" | Test imports `BCFTOOLS_VIEW` but process is named `BCFTOOLS_FILTER` | Change test to import `BCFTOOLS_FILTER` |
| `modules/local/rscript/assign_annotations/tests/main.nf.test` | "Cannot find RSCRIPT_ANNOTATE in assign_annotations/main.nf" | Process is `RSCRIPT_ASSIGN_ANNOTATIONS`; module is dead | **Delete this test file** |
| `modules/local/rscript/buildreports/tests/main.nf.test` | "Cannot find RSCRIPT_ANNOTATE in assign_annotations/main.nf" | Broken transitive `include` of dead code | Remove broken include; exercise only `RSCRIPT_BUILDREPORTS` |
| `modules/local/regenie/step1/tests/main.nf.test` | "Can't find rscript/annotate/main.nf" | After benchmarks, `BCFTOOLS_ASSIGN_ANNOTATIONS` replaced `RSCRIPT_ASSIGN_ANNOTATIONS` in production | Update include to `bcftools/assign_annotations/main.nf`; update process name |
| `modules/local/regenie/step2/tests/main.nf.test` | Same as step1 | Same | Same fix |

## §3 Category C — Missing executables (fixed by `--profile podman`)

Resolved in the second run; no further action required beyond mandating podman everywhere.
The only remaining Cat-C test is `bcftools/vcf2frq`, which is now a snapshot mismatch (Cat E).

## §3 Category D — Null path / dead code (delete, do not fix)

| Test file | Decision |
|---|---|
| `modules/local/plink19/makebed/tests/main.nf.test` | **Delete** — module is dead (§1a) |
| `modules/local/plink2/export_bgen/tests/main.nf.test` | **Delete** — module is dead (§1b) |
| `modules/local/bgenix/tests/main.nf.test` | **Delete** — module is dead (§1b) |
| `modules/local/qctool/tests/main.nf.test` | **Delete** — module is dead (§1b) |
| `modules/local/plink2/write_snplist/tests/main.nf.test` | Transitive failure from PLINK2_MAKEBED null path — fix MAKEBED first (T4); this should then pass (T5) |

## §3 Category E — Snapshot md5 mismatch

| Test file | Fix |
|---|---|
| `modules/local/bcftools/vcf2frq/tests/main.nf.test` | Re-record snapshot (T6) and audit the diff before committing |

## §4 Passing tests (after second run) — for the record

26/93 pass: `bcftools/index` (4), `python/merge_sex_covar` (1), `nf-core/multiqc` (3),
all 18 nf-core utility subworkflow tests.

---

## Tasks

### T3 — Fix Cat-B failures (wrong process name / path)

- **Files to edit**:
  - `modules/local/bcftools/filter/tests/main.nf.test` — change `BCFTOOLS_VIEW` to
    `BCFTOOLS_FILTER`; update include path to filter's `main.nf`
  - `modules/local/rscript/buildreports/tests/main.nf.test` — remove the include pointing at
    `rscript/assign_annotations`; exercise only `RSCRIPT_BUILDREPORTS`
  - `modules/local/regenie/step1/tests/main.nf.test` — change `rscript/annotate/main.nf` to
    `bcftools/assign_annotations/main.nf`; update process name to `BCFTOOLS_ASSIGN_ANNOTATIONS`
  - `modules/local/regenie/step2/tests/main.nf.test` — same as step1
- **Already removed by T2**: `modules/local/rscript/assign_annotations__to_delete/tests/main.nf.test`.
- **After fix compiles**, audit existing assertions against the quality bar in the main plan.
- **Verify**: `nf-test test --profile podman <each file>`.
- **Done-when**: all four named tests pass with at least one meaningful assertion each.

### T4 — Fix Cat-A failures (signature mismatches)

For each Cat-A test (excluding the dead `qctool` row), rewrite the workflow call to match the
current module signature. **The signature in `modules/local/<group>/<name>/main.nf` is the
source of truth, not the table above.**

- **Files to edit**:
  - `modules/local/bcftools/annotate/tests/main.nf.test`
  - `modules/local/bcftools/norm/tests/main.nf.test`
  - `modules/local/bcftools/view/tests/main.nf.test`
  - `modules/local/cmds/merge_results/tests/main.nf.test`
  - `modules/local/plink2/makebed/tests/main.nf.test`
  - `modules/local/rscript/manhattan_qq_plots/tests/main.nf.test`
  - `modules/local/vep/annotate/tests/main.nf.test`
- **Heads-up**: several of these (e.g. `manhattan_qq_plots`, `plink2/makebed`) take many
  input files. Just changing the call shape without providing realistic inputs will not work
  — you will need to author or borrow fixtures. For `manhattan_qq_plots`, lean on the IT-5
  reporting fixture (see [integration-tests.md](integration-tests.md)).
- **After each fix runs**, audit existing assertions against the quality bar. For modules
  with real custom logic (e.g. `cmds/merge_results` concatenates per-chromosome regenie
  outputs and re-sorts), add a row-count or column-presence assertion; for pure wrappers
  (e.g. `bcftools/annotate` with rename_chrs), one assertion checking that the rename
  actually took effect is enough.
- **Verify** (per file): `nf-test test --profile podman <file>`.
- **Done-when**: all seven tests pass individually and
  `nf-test test --profile podman --tag ci` shows zero Cat-A failures.

### T5 — Fix Cat-D transitive failure (`plink2/write_snplist`)

- **Prerequisite**: T4 done (PLINK2_MAKEBED test passes).
- **File to edit**: `modules/local/plink2/write_snplist/tests/main.nf.test` — provide
  required input files now that `PLINK2_MAKEBED` is invokable.
- **Verify**: `nf-test test --profile podman modules/local/plink2/write_snplist/tests/main.nf.test`.
- **Done-when**: passes with at least one meaningful assertion.

### T6 — Re-record Cat-E snapshot (`bcftools/vcf2frq`)

- **Command**:
  `nf-test test --profile podman --update-snapshot modules/local/bcftools/vcf2frq/tests/main.nf.test`.
- **Then**: `git diff modules/local/bcftools/vcf2frq/tests/main.nf.test.snap` — eyeball the
  diff. Only commit if changes are explained by (a) a `versions.yml` bcftools version bump,
  or (b) a deliberate change in module behaviour you can name.
- **Done-when**: snapshot committed; the test passes without `--update-snapshot`.
