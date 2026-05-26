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

Originally a single task; on 2026-05-25 it was split into four sub-tasks because the
"rename the include" fix described in the original write-up was insufficient. In every case
the test's `setup` chain or input shapes also drift from what the current production
workflow uses (signatures were consolidated; several upstream modules in the chain are now
`__to_delete`). Each sub-task is sized to a single test file and is independent of the
others, so they can be parallelised.

The current production callers — `workflows/rare-var-assoc.nf` and the prepare/PCA/etc.
subworkflows — are the source of truth for how each module is invoked now (input tuple
shape, what feeds it). Match the test setup to those callers rather than recreating the old
chain.

#### T3a — Fix `bcftools/filter` test (process rename + signature reshape)

- **File**: `modules/local/bcftools/filter/tests/main.nf.test`.
- **Why it fails today**:
  1. References `BCFTOOLS_VIEW` but the module process is `BCFTOOLS_FILTER`.
  2. Uses the old 6-input shape (`tuple(meta, vcf, tbi)` + 5 separate vals for
     regions/targets/samples/snplist/out_name_part). The current signature is one
     consolidated tuple `(meta, vcf, index, regions, targets, samples, snplist, tracking_in)`
     + `val(input_args)` + `val(out_name_part)` — i.e. only 3 inputs.
- **Fix**:
  - Rename process and `name` string to `BCFTOOLS_FILTER`.
  - Reshape every test's `when { process { ... } }` block to the new 3-input form. Use the
    production call in `workflows/rare-var-assoc.nf` (lines 140–160) as the template — empty
    lists `[]` for unused path slots, `[]` for `tracking_in` when not chaining.
  - Delete the four `- stub` test variants entirely. Stub tests violate §0.2 of the main
    plan ("never write stub tests — they verify nothing meaningful").
  - The four non-stub tests that exercise different `--write-index` configs (default, csi,
    tbi) and the regions/targets variant are worth keeping; consolidate to ≤2 tests if their
    only difference is the resulting index extension.
- **Assertions** (every test must meet the quality bar):
  - The `*_tracking.json` file's `outputs.variants` value matches the number of records in
    the output VCF (parse with `bcftools view -H | wc -l` in a `path.text` check or
    `path(...).readLines().size()`), OR assert that `variants_out ≤ variants_in` when a
    filter argument is exercised. `versions.yml` snapshot is fine for the YAML metadata.
- **Verify**: `nf-test test --profile podman modules/local/bcftools/filter/tests/main.nf.test`.
- **Done-when**: test passes; each remaining test case carries at least one assertion that
  names a concrete invariant (per §"Quality bar for assertions" in the main plan).

#### T3b — Fix `rscript/buildreports` test (strip dead-code setup chain)

- **File**: `modules/local/rscript/buildreports/tests/main.nf.test`.
- **Why it fails today**: the `setup` block runs PLINK2_EXPORT_BGEN, BGENIX,
  RSCRIPT_ANNOTATE (`rscript/assign_annotations` — now dead), RSCRIPT_VCFTOAAF
  (`rscript/vcf2aaf` — now dead), and the old-signature REGENIE_STEP1/STEP2. The first
  compile error is "Cannot find RSCRIPT_ANNOTATE in assign_annotations/main.nf"; even after
  fixing that, the chain is unrecoverable because three upstream modules in it are
  `__to_delete`.
- **Fix**:
  - Delete the entire `setup` block.
  - Feed `RSCRIPT_BUILDREPORTS` directly from canonical fixtures. The module's six inputs
    (`regenie_step2_masks_snplist`, `regenie_step2_Y1_regenie`, `vcf`, `phenotype`,
    `annotations`, `r_script_ch`) are all small text/CSV files; commit a minimal hand-rolled
    fixture set under `modules/local/rscript/buildreports/tests/fixtures/`. Use values
    consistent with the medium_data 1kGP chr22 dataset (a regenie output stub with
    ~10 rows, an annotations file with the same gene IDs, a phenotype file matching the
    sample IDs in `1kGP_cases.txt` + `1kGP_controls.txt`).
  - The R script itself lives at
    `modules/local/rscript/buildreports/assets/build_reports.R` — keep that reference.
- **Assertions** (reporting carve-out applies — smoke level OK per main plan):
  - Assert all three output CSVs exist and are non-empty (≥1 data row past header).
  - Assert the column header of `*_res_log10p_1_annotated.csv` includes the columns the
    downstream reporting templates rely on (look at `assets/build_reports.R` or grep the
    HTML/Rmd templates for which columns are read). One named column is enough.
  - Snapshot `versions.yml`.
- **Verify**: `nf-test test --profile podman modules/local/rscript/buildreports/tests/main.nf.test`.
- **Done-when**: test passes; smoke-level assertions in place.

#### T3c — Rewrite `regenie/step1` test against the current production shape

- **File**: `modules/local/regenie/step1/tests/main.nf.test`.
- **Why it fails today**:
  1. `setup` chain includes PLINK2_EXPORT_BGEN, BGENIX, RSCRIPT_VCFTOAAF — all dead code —
     and `RSCRIPT_ANNOTATE` (imported from `rscript/annotate/main.nf`, a path that does not
     exist). `RSCRIPT_ANNOTATE`'s job was to build the regenie mask files and filter
     variant annotations — that responsibility now belongs to `BCFTOOLS_ASSIGN_ANNOTATIONS`
     in production. (Phenotype generation is a separate concern handled by
     `RSCRIPT_BUILD_PHENOTYPES` in `subworkflows/local/utils_nfcore_rare-var-assoc_pipeline`,
     not by `RSCRIPT_ANNOTATE`.)
  2. The setup also uses `PLINK2_MAKEBED`, which is no longer the production code path —
     genotype-data conversion has moved to `PLINK2_MAKEPGEN` everywhere in
     `workflows/rare-var-assoc.nf`.
  3. Current `REGENIE_STEP1` signature is a single consolidated 9-element tuple
     `(meta, pgen, pvar, psam, qc_pass_id, qc_pass_snplist, phenotype, covar_file, tracking_in)`
     + `val(input_args)`. The test passes 5 separate inputs.
- **Fix** (use `workflows/rare-var-assoc.nf` lines 323–409 as the template):
  - Rebuild `setup` using only live modules:
    - `PLINK2_MAKEPGEN` on the medium_data 1kGP chr22 VCF. Signature is in
      `modules/local/plink2/makepgen/main.nf` — input[0] is a 10-element tuple, then 5
      `val()` args (samples-filter-type, variants-filter-type, vcf-input-options,
      out-name-part, input-args). Mirror the production call at lines 334–343 (use empty
      `[]` for filter files when not exercised).
    - `PLINK2_WRITE_SNPLIST` on PLINK2_MAKEPGEN's `out_pgen_pvar_psam`. Mirror the
      production call at lines 323–327. Use `Channel.value('writesnp_pass')` for
      `out_name_part` and a minimal QC option string.
  - Phenotype: commit a minimal phenotype TSV under
    `modules/local/regenie/step1/tests/fixtures/phenotype.tsv` matching the 1kGP chr22
    sample IDs and one binary trait column (header `FID IID Y1`). Do not regenerate
    phenotypes via any module in the test — keep this as a static fixture.
  - Build the consolidated tuple input as the production workflow does (lines 390–403):
    `tuple(meta, pgen, pvar, psam, qc_pass_id, qc_pass_snplist, phenotype, [], [])` —
    empty `covar_file` and `tracking_in` when not exercising those paths.
- **Assertions**:
  - `*_pred.list` is non-empty and references the `.loco` file by name.
  - `*.loco` is non-empty.
  - Parse the tracking JSON: `inputs.samples` ≥ 1, `inputs.variants` ≥ 1.
  - Snapshot `versions.yml`.
- **Verify**: `nf-test test --profile podman modules/local/regenie/step1/tests/main.nf.test`.
- **Done-when**: test passes; assertions in place.
- **No test-level dependency**: `PLINK2_MAKEPGEN` has no test directory and the
  `PLINK2_WRITE_SNPLIST` module itself works (its own test failure in §3 Cat-D is a
  pre-existing harness issue, not a code defect). T3c can be implemented independently of
  T4 and T5.

#### T3d — Rewrite `regenie/step2` test against the current production shape

- **File**: `modules/local/regenie/step2/tests/main.nf.test`.
- **Why it fails today**: same root causes as T3c plus dependencies on `qctool` and
  `bgenix` (both dead). Current `REGENIE_STEP2` signature is a single 12-element tuple
  + `val(input_args)` + `val(return_snplist)`.
- **Fix** (use `workflows/rare-var-assoc.nf` lines 405–488 as the template):
  - Reuse T3c's setup up to and including REGENIE_STEP1.
  - Add `BCFTOOLS_ASSIGN_ANNOTATIONS` and `PYTHON_VCFTOAAF` to the chain. These are the
    current production replacements for the dead `RSCRIPT_ANNOTATE` /  `RSCRIPT_VCFTOAAF`
    R-based modules: `BCFTOOLS_ASSIGN_ANNOTATIONS` builds the mask annotations + setlist;
    `PYTHON_VCFTOAAF` produces the alt-allele-frequency file. Both take small tuple
    inputs; signatures in their respective `main.nf`. For `BCFTOOLS_ASSIGN_ANNOTATIONS`,
    the python script asset lives at
    `modules/local/bcftools/assign_annotations/assets/assign_annotations.py`. For
    `PYTHON_VCFTOAAF`, the python script asset lives at
    `modules/local/python/vcf2aaf/assets/vcf2aaf.py`.
  - Build the consolidated 12-element tuple as in lines 466–479:
    `tuple(meta, pgen, pvar, psam, phenotype, annotations, setlist, aaf, step1_pred_list,
    [], masks, [])` — empty `covar_file` and `tracking_in`, masks from
    `assets/default.masks`.
  - Pass `Channel.value('--bt --ref-first --firth --approx --bsize 200 --lowmem --aaf-bins 0.01,0.05,0.1,1 --write-mask --vc-tests skato')`
    and `Channel.value(true)` for `return_snplist` so `masks_snplist` is emitted.
- **Assertions**:
  - `*_step2.regenie` exists and contains a header row with the expected regenie columns
    (`CHROM POS ID ALLELE0 ALLELE1 ... TEST BETA SE CHISQ LOG10P` or whatever current
    regenie emits — read the file head with `path(...).readLines()[0]`).
  - `*_step2.regenie` has at least one data row.
  - `*_masks.snplist` is non-empty.
  - Snapshot `versions.yml`.
- **Verify**: `nf-test test --profile podman modules/local/regenie/step2/tests/main.nf.test`.
- **Done-when**: test passes; assertions in place.
- **Heads-up**: depends on T3c (reuse its setup). Plan for ~10–15 min wall clock per test
  run.

#### T3 sub-task status

| Sub-task | Status | Depends on |
|---|---|---|
| T3a — `bcftools/filter` | ✅ Done 2026-05-25 — 3 tests pass | none |
| T3b — `rscript/buildreports` | ✅ Done 2026-05-25 — 1 test passes; fixtures under `tests/fixtures/` | none |
| T3c — `regenie/step1` | ✅ Done 2026-05-25 — 1 test passes; uses `prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz` for genuine LOCO; Welch t-test (p<0.01) confirms model captures YRI vs CHB population structure | none |
| T3d — `regenie/step2` | ✅ Done 2026-05-25 — 1 test passes; ##MASKS comment-stripping needed; WWC3 highest LOG10P assertion (p<0.01); ~30 s wall clock; tagged `"full"` | T3c (reuses its setup chain) |

- **Combined done-when**: all four sub-tasks done; `nf-test test --profile podman --tag ci`
  shows zero Cat-B failures.

### T4 — Fix Cat-A failures (signature mismatches)

Originally a single task; on 2026-05-25 it was split into seven sub-tasks. The seven test
files vary widely in size and required setup: some just need the call tuple reshaped
(`cmds/merge_results`), others need fresh fixtures (`plink2/makebed`, `vep/annotate`), and
one (`rscript/manhattan_qq_plots`) is blocked on the IT-5 reporting fixture from T11 and is
deferred. Sub-tasks T4a..T4f are independent and can be parallelised.

**For every sub-task**: the signature in `modules/local/<group>/<name>/main.nf` is the
source of truth (not the table in §3 Cat-A). The production caller — the subworkflow that
already uses the module — is the template for how to shape inputs. After the test compiles
and runs, audit existing assertions against the **quality bar for assertions** in the main
plan before declaring done.

| Module | Production caller (template) |
|---|---|
| `bcftools/annotate` | [subworkflows/local/prepare/main.nf:93](../../subworkflows/local/prepare/main.nf#L93) |
| `bcftools/norm` | [subworkflows/local/prepare/main.nf:79](../../subworkflows/local/prepare/main.nf#L79) |
| `bcftools/view` | [subworkflows/local/prepare/main.nf:66](../../subworkflows/local/prepare/main.nf#L66), [workflows/rare-var-assoc.nf:375](../../workflows/rare-var-assoc.nf#L375) |
| `cmds/merge_results` | [workflows/rare-var-assoc.nf:515](../../workflows/rare-var-assoc.nf#L515) |
| `plink2/makebed` | [subworkflows/local/pca/main.nf:25](../../subworkflows/local/pca/main.nf#L25) |
| `rscript/manhattan_qq_plots` | [subworkflows/local/reporting/main.nf:23](../../subworkflows/local/reporting/main.nf#L23) |
| `vep/annotate` | [workflows/rare-var-assoc.nf:171](../../workflows/rare-var-assoc.nf#L171) |

#### T4a — Fix `cmds/merge_results` test (smallest)

- **File**: `modules/local/cmds/merge_results/tests/main.nf.test`.
- **Why it fails today**: current signature is one tuple
  `(meta, phenotype, regenie_chromosomes, csv_concat_py_script)` — the test passes the
  python script as a separate `input[0]` and the regenie channel as `input[1]`.
- **Fix**: collapse both arguments into a single tuple per the production call. Reuse the
  two existing test cases ("merging 4 regenie files", "merging fileA and fileB") — assertions
  there already meet the quality bar (column header + row count). Drop the `// TODO nf-core:`
  comment block at the top.
- **Verify**: `nf-test test --profile podman modules/local/cmds/merge_results/tests/main.nf.test`.
- **Done-when**: both tests pass; existing assertions retained.

#### T4b — Fix `bcftools/annotate` test

- **File**: `modules/local/bcftools/annotate/tests/main.nf.test` (373 lines, many variants).
- **Why it fails today**: signature is now one consolidated 7-element tuple
  `(meta, input, index, annotations, annotations_index, header_lines, rename_chrs)` +
  `val(out_name_part)`. Tests pass 3 separate inputs in the old shape.
- **Fix**:
  - Reshape every test's input block to the consolidated tuple form. Use empty `[]` for the
    path slots a given test does not exercise (e.g. `header_lines`, `rename_chrs`).
  - Audit the per-test `.config` files (`vcf.config`, `vcf_gz_index.config`, etc.) — they
    should still drive `task.ext.args` and `task.ext.prefix`. No changes there expected.
  - Consolidate near-duplicate test variants if the only difference is the index extension.
- **Assertions** (every test must meet the quality bar):
  - For tests exercising `--rename-chrs`: assert that a known chromosome name in the input
    VCF appears renamed in the output (`bcftools view -H | head -1 | cut -f1`).
  - For tests exercising `--annotations`: assert that an INFO field added by the annotation
    file appears in `bcftools view -h` output.
  - `versions.yml` snapshot is fine for the YAML metadata.
- **Verify**: `nf-test test --profile podman modules/local/bcftools/annotate/tests/main.nf.test`.
- **Done-when**: every retained test passes with at least one business-meaningful assertion.

#### T4c — Fix `bcftools/norm` test

- **File**: `modules/local/bcftools/norm/tests/main.nf.test` (580 lines, largest in T4).
- **Why it fails today**: signature is one 5-element tuple
  `(meta, vcf, tbi, fasta, tracking_in)` + `val(out_name_part)`. Tests pass 3 separate
  inputs.
- **Fix**:
  - Reshape to the consolidated tuple form per the production call. Pass `[]` for
    `tracking_in` when not chaining.
  - Consolidate near-duplicate test variants — at most 3 retained (split-multiallelic, join,
    and a no-op/`--check-ref` variant) if the existing 580-line file repeats the same shape
    many times.
- **Assertions**:
  - Splitting tests: variant count in output > variant count in input when a multiallelic
    test input is used; tracking JSON `outputs.variants > inputs.variants`.
  - Joining tests: opposite inequality.
- **Verify**: `nf-test test --profile podman modules/local/bcftools/norm/tests/main.nf.test`.
- **Done-when**: retained tests pass; assertions in place.

#### T4d — Fix `bcftools/view` test

- **File**: `modules/local/bcftools/view/tests/main.nf.test` (315 lines).
- **Why it fails today**: signature is one 8-element tuple
  `(meta, vcf, index, regions, targets, samples, snplist, tracking_in)` + `val(input_args)`
  + `val(out_name_part)`. Tests pass 6 separate inputs.
- **Fix**:
  - Reshape to the consolidated tuple form. The production call at
    `workflows/rare-var-assoc.nf:375` (BCFTOOLS_VIEW_2 — sample subsetting) and
    `subworkflows/local/prepare/main.nf:66` (BCFTOOLS_VIEW_1 — region/target filtering)
    cover the two main shapes.
  - Keep two regions/targets/samples variants — one that exercises `regions`, one that
    exercises `samples`. Drop stub variants if present.
- **Assertions**:
  - Region test: assert output VCF contains only records on the requested chromosome
    (`bcftools view -H | awk '{print $1}' | sort -u`).
  - Sample test: assert the output VCF header lists exactly the requested sample(s)
    (`bcftools view -h | tail -1`).
  - Tracking JSON consistency: `outputs.variants ≤ inputs.variants` when a filter is applied.
- **Verify**: `nf-test test --profile podman modules/local/bcftools/view/tests/main.nf.test`.
- **Done-when**: retained tests pass; assertions in place.

#### T4e — Fix `plink2/makebed` test (**unblocks T5**)

- **File**: `modules/local/plink2/makebed/tests/main.nf.test`.
- **Why it fails today**: signature is one 9-element tuple
  `(meta, pgen, pvar, psam, vcf, frq, samples_filtering_file, variants_filtering_file, tracking_in)`
  + 4× `val(samples_filtering_type, variants_filtering_type, out_name_part, input_args)`.
  The current test passes 3 inputs (old shape: 5-path tuple + 2 vals).
- **Fix**:
  - Reshape the single active test ("small vcf - bed,bim,fam") to the new 5-input form.
    Use the production call at [subworkflows/local/pca/main.nf:25](../../subworkflows/local/pca/main.nf#L25) as the template.
  - For the VCF input path, use the canonical fixture
    `assets/three_chr_unprepared/unprepared_rand_500.vcf.gz` (§0.4) — the existing
    `genomics/homo_sapiens/illumina/plink/test.rnaseq.{bed,bim,fam}` data path can be kept
    if the test exercises the bed/bim/fam input branch; otherwise switch to VCF input via
    `vcf` slot with empty `pgen/pvar/psam`.
  - Delete the two commented-out block tests at the bottom (`split-par`, `impute-sex`) —
    they reference the dead PLINK19_MAKEBED and a since-removed signature.
- **Assertions**:
  - The emitted `.bed`/`.bim`/`.fam` trio is non-empty.
  - The `.fam` row count equals the sample count expected for the input VCF (3202 for the
    full fixture; document the expected count if a smaller subset is used).
  - The tracking JSON's `outputs.variants` is ≥ 1 and `outputs.samples` matches the `.fam`
    row count.
  - Snapshot `versions.yml`.
- **Verify**: `nf-test test --profile podman modules/local/plink2/makebed/tests/main.nf.test`.
- **Done-when**: test passes; assertions in place. **T5 unblocked.**

#### T4f — Fix `vep/annotate` test

- **File**: `modules/local/vep/annotate/tests/main.nf.test`.
- **Why it fails today**: signature is one 4-element tuple
  `(meta, input_vcf, input_vcf_tbi, vep_cache)` + 3× `val(species, fasta_path, input_args)`.
  The current test passes the VCF as a 3-element tuple and `vep_cache` as a separate
  `input[1]` (old 5-input shape).
- **Fix**:
  - Move `vep_cache` into the tuple, leaving 4 inputs.
  - The existing test references `${projectDir}/../vep_cachedir/homo_sapiens_refseq` — that
    path exists on this machine but is not portable. Document this in a short comment and
    tag the test `tag "full"` since the cache is large and machine-specific. CI will skip it.
  - Drop the `// TODO nf-core:` comment block at the top.
- **Assertions** (smoke level — VEP output is annotation-heavy, hard to assert business
  invariants on):
  - The output VCF header (`bcftools view -h`) contains a `##INFO=<ID=CSQ,...>` line.
  - At least one record in the output carries a non-empty `CSQ=` field.
  - Snapshot `versions.yml`.
- **Verify**: `nf-test test --profile podman --tag full modules/local/vep/annotate/tests/main.nf.test`.
- **Done-when**: test passes locally; tagged `"full"` because the VEP cache is
  machine-specific.

#### T4g — Fix `rscript/manhattan_qq_plots` test (**deferred**)

- **File**: `modules/local/rscript/manhattan_qq_plots/tests/main.nf.test`.
- **Deferred reason**: the module takes an 11-element tuple including `phenotype_file`,
  `pc_plot_file`, `eda_plots`, `setlist_file`, `regenie_merged`, `mask_file`,
  `gwas_report_template`, `r_functions_file`, and `rmd_pheno_stats_file`. Hand-crafting all
  of these for an isolated module test is high effort and low value when the IT-5 reporting
  subworkflow integration test (T11) will exercise the same path with realistic upstream
  fixtures.
- **Plan**: implement this as part of T11 (IT-5), or immediately after, by reusing the IT-5
  fixture set — keeping the per-module test as a thin re-entry point.
- **Assertions** (reporting carve-out — smoke level OK):
  - The output `*.html` exists and is non-empty.
  - `versions.yml` snapshot.
- **Verify**: `nf-test test --profile podman modules/local/rscript/manhattan_qq_plots/tests/main.nf.test`.
- **Done-when**: T11 IT-5 fixtures exist and this test reuses them with smoke-level
  assertions.

#### T4 sub-task status

| Sub-task | Status | Depends on |
|---|---|---|
| T4a — `cmds/merge_results` | ✅ Done 2026-05-26 — 2 tests pass; meta added to tuple; `.get(2)` for file path | none |
| T4b — `bcftools/annotate` | ✅ Done 2026-05-26 — 5 tests pass (11→5, 4 stubs + 2 near-dups dropped); bcf.config→-Ov; fixture files for header_lines/rename_chrs | none |
| T4c — `bcftools/norm` | ✅ Done 2026-05-26 — 3 tests pass (16→3, all stubs dropped); split test uses `tests/fixtures/multiallelic_sarscov2.vcf.gz` (8 records, 3 multiallelic → 12 after split), join uses `tests/fixtures/split_sarscov2.vcf.gz` (12 → 8), basic norm uses nf-core sarscov2 test.vcf.gz (9 = 9); tracking JSON asserts strict `>`, `<`, `==` respectively | none |
| T4d — `bcftools/view` | ✅ Done 2026-05-26 — 2 tests pass (9 old→2 new; 4 stubs + 3 write-index-only variants dropped); regions BED fixture (2/9 variants pass); samples fixture (header asserts 'test'); tracking JSON strict inequalities | none |
| T4e — `plink2/makebed` (unblocks T5) | ✅ Done 2026-05-26 — 2 tests pass; (1) VCF branch + psam sex info + keep_3_samples.txt fixture; outputs.samples (3) + outputs.variants >= 1 asserted; (2) pgen branch (setup: PLINK2_MAKEPGEN) mirroring production pca/main.nf:25 exactly; fam count (3202) + inputs.samples (3202) asserted; --split-par b38 + --max-alleles 2 needed; plink2 omits "remaining" lines without filters → outputs.{samples,variants} = -1 in that case | none |
| T4f — `vep/annotate` | ✅ Done 2026-05-26 — 1 test passes; vep_cache moved into 4-element tuple; linesGzip asserts CSQ header + at least one CSQ= data record; tagged "full" (machine-local vep_cachedir) | none |
| T4g — `rscript/manhattan_qq_plots` | deferred | T11 (IT-5 fixtures) |

- **Combined done-when**: T4a..T4f all done; `nf-test test --profile podman --tag ci` shows
  zero Cat-A failures from these six modules. T4g closes once T11 lands.

### T5 — Fix Cat-D transitive failure (`plink2/write_snplist`)

✅ Done 2026-05-26 — 1 test passes.

Setup uses PLINK2_MAKEPGEN on `tests/fixtures/prepared_100.vcf.gz` (first 100 chr12
variants extracted from the prepared medium_data fixture; unique IDs, no multiallelic sites
— the unprepared_rand_500 fixture was unusable here because `--write-snplist` rejects
duplicate variant IDs). Assertions: snplist == 100 lines; .id == 3203 lines (header +
3202 samples); tracking.inputs.{samples,variants} == 3202/100;
tracking.outputs.variants == snplist row count.

### T6 — Re-record Cat-E snapshot (`bcftools/vcf2frq`)

- **Command**:
  `nf-test test --profile podman --update-snapshot modules/local/bcftools/vcf2frq/tests/main.nf.test`.
- **Then**: `git diff modules/local/bcftools/vcf2frq/tests/main.nf.test.snap` — eyeball the
  diff. Only commit if changes are explained by (a) a `versions.yml` bcftools version bump,
  or (b) a deliberate change in module behaviour you can name.
- **Done-when**: snapshot committed; the test passes without `--update-snapshot`.
