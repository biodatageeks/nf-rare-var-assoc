# Fixture: fcoef_truth_edges

Scenario:
- `S05` F-coefficient outlier filtering

## Truth design

- Include known high and low F outliers plus non-outlier controls.
- Engineer values close to threshold edges.

## Required expected artifacts

- `expected/fcoef_outliers_expected_removed.txt`
- `expected/fcoef_controls_expected_retained.txt`

## Assertion contract

- Known outliers are absent after `F_COEFFICIENT_FILTERING`.
- Known inliers are retained.
