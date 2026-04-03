# Fixture: pca_truth_influence

Scenario:
- `S04` covariate influence from PCs

## Truth design

- Input should induce at least one association result shift when PC covariates are included.
- Two deterministic runs are required: with and without `PC*_AVG` covariates.

## Required expected artifacts

- `expected/pc_covariate_delta.tsv`: expected rows with non-zero deltas.
- `expected/influence_thresholds.yml`: epsilon values and rationale.

## Assertion contract

- At least one effect estimate or p-value changes beyond epsilon.
- "No change" is a failure.
