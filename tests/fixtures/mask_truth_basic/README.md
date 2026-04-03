# Fixture: mask_truth_basic

Scenario:
- `S09` annotation mask enforcement

## Truth design

- Include variants with mixed allowed/disallowed annotations under selected masks.
- Place allowed and disallowed rows in the same loci set where possible.

## Required expected artifacts

- `expected/allowed_annotations_expected.txt`
- `expected/disallowed_annotations_expected_absent.txt`

## Assertion contract

- Disallowed annotations are absent from setlist/analysis outputs.
- Allowed annotations remain present.
