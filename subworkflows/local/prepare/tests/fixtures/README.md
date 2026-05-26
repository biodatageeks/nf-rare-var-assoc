# PREPARE subworkflow test fixtures

## `all_samples.txt`

Concatenation of `assets/three_chr_unprepared/cases.txt` (1202 samples) and
`assets/three_chr_unprepared/controls.txt` (2000 samples). One sample IID per line, no header.

Mirrors the `all.samples` file produced at runtime by `JOIN_CASES_AND_CONTROLS` in
`subworkflows/local/utils_nfcore_rare-var-assoc_pipeline/main.nf` (cases first, then controls).
PREPARE expects this single "all samples" file rather than the cases/controls pair, because
PIPELINE_INITIALISATION joins them before invoking PREPARE.

Regenerate:

```bash
cat assets/three_chr_unprepared/cases.txt assets/three_chr_unprepared/controls.txt \
    > subworkflows/local/prepare/tests/fixtures/all_samples.txt
```
