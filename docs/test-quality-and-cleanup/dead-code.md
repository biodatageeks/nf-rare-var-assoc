# Dead Code Cleanup (T2, T2.5, T14, T16)

See [../test-quality-and-cleanup-plan.md](../test-quality-and-cleanup-plan.md) for goals,
conventions, and the full task list.

Two-phase strategy: phase 1 removes `include` statements and renames module directories to
`__to_delete` (recoverable if we discover the module is actually used); phase 2 deletes the
renamed directories permanently after integration tests confirm nothing is missing.

---

## §1a — Modules with zero references outside the benchmark workflow

None appear in any workflow or subworkflow other than `workflows/benchmark_implementations.nf`:

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

## §1b — Modules imported but never called

`include` statements exist but no call site is present in any workflow or subworkflow:

| Module | Imported in | What to do |
|---|---|---|
| `modules/local/plink2/export_bgen` | main workflow | Remove import + rename dir to `__to_delete` |
| `modules/local/bgenix` | main workflow | Remove import + rename dir |
| `modules/local/qctool` | main workflow | Remove import + rename dir |
| `modules/local/cmds/download_file` | main workflow | Remove import + rename dir |
| `modules/local/bcftools/tag2tag` | main workflow + prepare subworkflow | Remove both imports + rename dir |

## §1c — Stale imports in `subworkflows/local/prepare/main.nf`

| Import | Situation | What to do |
|---|---|---|
| `BCFTOOLS_VCF2PSAM` | Never called in prepare; module IS called in the main workflow | Remove import line from prepare only — keep the module |
| `BCFTOOLS_VCF2FRQ` | Same | Remove import line from prepare only |
| `VIEW_AND_FILTER2_POLARSBIO` | Never called anywhere; replaced by `BCFTOOLS_VIEW_AND_FILTER2` | Remove import + rename dir |
| `FILTER_AND_ENHANCE_VCF_POLARSBIO` | Never called in prepare; phase-2 deletion candidate | Remove import now; rename dir in phase 2 |

## §1d — Benchmark workflow

`workflows/benchmark_implementations.nf` is not referenced from `main.nf`. Remove it together
with `conf/benchmark.config` once §1a modules are deleted.

## §1e — Test files deleted in phase 1 (along with imports)

| Test file | Reason |
|---|---|
| `modules/local/plink19/makebed/tests/main.nf.test` | Module is dead code (§1a) |
| `modules/local/rscript/vcf2aaf/tests/main.nf.test` | Module is dead code (§1a) |
| `modules/local/rscript/assign_annotations/tests/main.nf.test` | Module is dead code (§1a) |
| `modules/local/plink2/export_bgen/tests/main.nf.test` | Module is dead code (§1b) |
| `modules/local/bgenix/tests/main.nf.test` | Module is dead code (§1b) |
| `modules/local/qctool/tests/main.nf.test` | Module is dead code (§1b) |

## §1f — Phase 2 candidates from the VCF-preparation refactor

Preparation logic was migrated to the upstream `nf-prepare-vcf` pipeline in Phase 3 (PB2,
2026-05-29). As of that change **none of these modules are referenced by any prepare-path code**
— they are confirmed-orphaned and slated for removal in PB4:

- `modules/local/python/filter_and_enhance_vcf_polarsbio`
- `modules/local/python/fix_zero_PL`
- `modules/local/combo/filter_and_enhance_vcf`

Do **not** write tests for these modules — they are being removed. The old IT-1 `use_dosage=true`
path that used to exercise `combo/filter_and_enhance_vcf` is gone (IT-1/IT-1b were retired in PB3;
see [integration-tests.md](integration-tests.md)). PB4 should `rg`-check and delete these three
plus this repo's now-unused prep copies (`bcftools/norm`, `bcftools/annotate`, `vep/annotate`,
`vep/updatecache`, `bcftools/filter`) iff unreferenced.

## §1g — Abandoned `tests/` scaffold (delete in T2.5)

The whole top-level `tests/` directory (`tests/workflows/*.nf.test.template`,
`tests/workflows/run_initial_tests.sh`, `tests/workflows/README.md`,
`tests/workflows/expected/`, `tests/fixtures/*`, `tests/results/`, `tests/nextflow.config`,
`tests/README.md`) is a leftover from an earlier aborted test design
(see `docs/skip_preparation_skip_reporting_test_design__old.md`). Pipeline code has been
substantially updated since; none of the `.template` files were ever converted to runnable
tests; the `tests/fixtures/` subdirs contain only `README.md` + empty `.gitkeep` placeholders.

Decision: **delete the entire `tests/` directory** in T2.5. The old design doc itself is
kept on disk (the `__old` suffix already marks it as superseded).

## §1h — nf-core template comments cleanup (T16)

Many modules still have `// TODO nf-core:` blocks copied from the nf-core template (examples:
`modules/local/plink2/export_bgen`, `modules/local/qctool`, `modules/local/plink2/makebed`).
Remove them in a final cleanup pass after tests are stable.

---

## Tasks

### T1 — Confirm dead-code list (✅ Done 2026-05-24)

Confirmed with the developer that §1a–§1c modules are truly unused.

### T2 — Phase 1 dead-code removal (✅ Done 2026-05-24)

Removed 5 includes from `workflows/rare-var-assoc.nf` and 5 from
`subworkflows/local/prepare/main.nf`; renamed 15 module dirs/files to `__to_delete`; deleted
6 obsolete test files; 3 471 lines removed.

### T2.5 — Delete abandoned `tests/` scaffold (✅ Done 2026-05-25, partial revert)

`git rm` of 24 files + `rm -r` of untracked `tests/results/`. `tests/nextflow.config`
restored: it is nf-test's default `configFile` (set in repo-root `nf-test.config`) and
removing it aborts the harness. Contents are stale — see T2.6.

### T2.6 — Prune `tests/nextflow.config` contents (✅ Done 2026-05-25)

Reduced to the `process { resourceLimits = [...] }` block; all the workflow-level params
it carried (`filter_and_enhance_vcf_*`, `use_dosage / skip_preparation / skip_reporting`,
`medium_data` `input_vcf / input_controls / input_cases`, `project_name`, `outdir`,
`rscript_vcf2aaf_options`) were only consumed by the deleted `tests/workflows/*.template`
scaffolds. Workflow ITs (T8+) will set their own params via a sibling `nextflow.config`
next to each `.nf.test`, pointing at the canonical `assets/three_chr_unprepared/` fixtures
(§0.4). Verified: `nf-test test --profile podman modules/local/bcftools/filter/tests/main.nf.test`
— 3/3 pass.

### T14 — Phase 2 dead-code removal

- **Prerequisite**: all tests in T3..T13 pass.
- **Action**: `git rm -r` every `*__to_delete*` path:
  - `modules/local/bcftools/{assign_annotations/main_shell__to_delete.nf, filter_qual_dp__to_delete, tag2tag__to_delete, vcf2aaf__to_delete}`
  - `modules/local/bgenix__to_delete`
  - `modules/local/cmds/download_file__to_delete`
  - `modules/local/combo/{filter2__to_delete, prepare__to_delete}`
  - `modules/local/plink19/makebed__to_delete`
  - `modules/local/plink2/export_bgen__to_delete`
  - `modules/local/python/view_and_filter2_polarsbio__to_delete`
  - `modules/local/qctool__to_delete`
  - `modules/local/rscript/{assign_annotations__to_delete, vcf2aaf__to_delete}`
  - `modules/local/vcftools/filter_qual_dp__to_delete`
  - `workflows/benchmark_implementations.nf`, `conf/benchmark.config`
- **Verify**: `nf-test test --profile podman --tag ci` still passes;
  `find . -name '*__to_delete*'` returns empty.
- **Done-when**: no `__to_delete` paths remain; all CI tests still pass.

### T16 — Remove nf-core template comments

- **Targets**: every `// TODO nf-core:` block in `modules/local/`. Search:
  `rg -n '// TODO nf-core:' modules/local/`.
- **Action**: delete the comment line(s). Do not touch any non-template TODO.
- **Verify**: `rg -n '// TODO nf-core:' modules/local/` returns no results.
- **Done-when**: no nf-core TODO blocks remain in `modules/local/`.
