# Fixture: missingness_truth_edges

Scenario:
- `S06` missingness boundary filtering per phenotype

## Truth design

- Include samples just below, at, and above missingness thresholds.
- Ensure deterministic per-phenotype expected sample retention.

## Required expected artifacts

- `expected/missingness_expected_removed.txt`
- `expected/missingness_expected_retained.txt`
- `expected/missingness_expected_counts.tsv`

## Assertion contract

- IDs and counts before/after filtering match exactly.
