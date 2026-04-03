# Workflow-level tests (Phase 1 scaffold)

This directory contains workflow-level semantic test templates for `RARE_VAR_ASSOC`.

## Current state

- `*.nf.test.template` files are prepared for scenarios `S01`, `S03`, `S04`, `S05`, `S06`, `S09`, and smoke `S11`.
- Each template includes scenario mapping and TODO points for fixture-specific semantic assertions.

## Runnable now

- `run_initial_tests.sh` executes an initial smoke+differential chunk:
   - run 1: `use_dosage=true`
   - run 2: `use_dosage=false`
   - checks pipeline info artifacts for both runs
   - verifies required processes completed in Nextflow trace and that `MERGE_RESULTS` is absent when `skip_reporting=true`
   - validates every `regenie_step2/*.regenie` file has numeric `P` values within `[0,1]`
   - performs semantic comparison of REGENIE rows between dosage modes and fails if no row-level statistic changes

Run command:

```bash
tests/workflows/run_initial_tests.sh
```

## Why templates (not runnable tests yet)

Semantic tests require engineered fixture inputs and finalized `expected/` files.
Without those artifacts, runnable tests would provide only shallow success checks and miss the business-logic regressions this suite is designed to catch.

## To finalize a scenario

1. Materialize fixture input files under `tests/fixtures/<fixture_name>/`.
2. Materialize expected truth files under `tests/fixtures/<fixture_name>/expected/`.
3. Rename a template from `*.nf.test.template` to `*.nf.test`.
4. Replace placeholder assertions with semantic checks against the fixture contract.
5. Run:
   - `nf-test test tests/workflows/<file>.nf.test`

## Recommended implementation order

1. `rare_var_assoc_semantic_dosage.nf.test.template` (`S01`, `S02`)
2. `rare_var_assoc_semantic_pca.nf.test.template` (`S03`, `S04`)
3. `rare_var_assoc_semantic_filters.nf.test.template` (`S05`, `S06`)
4. `rare_var_assoc_semantic_masks.nf.test.template` (`S09`)
5. `skip_preparation_skip_reporting.nf.test.template` (`S11`)
