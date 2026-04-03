# Fixture: dosage_truth_basic

Scenarios:
- `S01` dosage influence (`use_dosage=true` vs `false`)
- `S02` hard genotype behavior when dosage is disabled

## Truth design

- Construct at least one variant where hard `GT` is non-informative but `DS` differs between cases and controls.
- Keep deterministic sample order and fixed phenotype assignment.

## Required expected artifacts

- `expected/dosage_true_vs_false_delta.tsv`: expected changed association rows and delta thresholds.
- `expected/causal_variant_expectation.tsv`: expected causal variant behavior for both runs.

## Assertion contract

- At least one p-value or effect size differs beyond epsilon.
- Expected causal row is stronger/present with `use_dosage=true` and absent/weaker with `use_dosage=false`.
