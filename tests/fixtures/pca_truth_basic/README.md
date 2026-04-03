# Fixture: pca_truth_basic

Scenario:
- `S03` PCA numerical correctness

## Truth design

- Include a deterministic subset with known PCA projection results.
- Keep filtering and sample ordering identical to independent recomputation pipeline.

## Required expected artifacts

- `expected/pca_reference.tsv`: independent PCA reference values (`IID`, `PC1`, `PC2`, ...).
- `expected/pca_tolerance.yml`: numeric tolerances used in assertions.

## Assertion contract

- Pipeline PCs match independent recomputation within configured tolerance.
