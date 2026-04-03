# Workflow-level test assets

This directory hosts workflow-level semantic regression assets for `RARE_VAR_ASSOC`.

## Current contents

- `workflows/`: test specifications and nf-test templates.
- `fixtures/`: engineered truth fixtures for semantic assertions.

## Phase 1 scaffold status

- Fixture directories for `S01`, `S03`, `S04`, `S05`, `S06`, `S09`, and smoke `S11` have been created.
- Each fixture has an `expected/` directory for deterministic expected outputs.
- nf-test files are added as templates (`*.nf.test.template`) and should be converted to runnable `*.nf.test` once synthetic fixture inputs are finalized.

## Planned runnable workflow test files

- `workflows/rare_var_assoc_semantic_dosage.nf.test`
- `workflows/rare_var_assoc_semantic_pca.nf.test`
- `workflows/rare_var_assoc_semantic_filters.nf.test`
- `workflows/rare_var_assoc_semantic_masks.nf.test`
- `workflows/rare_var_assoc_semantic_signal.nf.test`
- `workflows/skip_preparation_skip_reporting.nf.test`
