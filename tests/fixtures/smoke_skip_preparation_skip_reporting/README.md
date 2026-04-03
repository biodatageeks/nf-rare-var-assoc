# Fixture: smoke_skip_preparation_skip_reporting

Scenario:
- `S11` structural smoke checks for `skip_preparation=true` and `skip_reporting=true`

## Input basis

- Existing medium data assets can be reused:
  - `assets/medium_data/prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz`
  - `assets/medium_data/1kGP_cases_200.txt`
  - `assets/medium_data/1kGP_controls_400.txt`

## Required expected artifacts

- `expected/key_outputs_present.txt`: key files expected after successful run.

## Assertion contract

- Workflow completes successfully.
- Key outputs and pipeline info artifacts exist.
