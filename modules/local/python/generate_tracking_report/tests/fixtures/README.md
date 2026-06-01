# generate_tracking_report test fixtures

Two chained tracking JSON files that represent a two-step pipeline run.

## step1.tracking.json

First pipeline step (`WORKFLOW:PREPARE_STEP`), no predecessor.
- Input: 500 variants, 100 samples
- Output: 450 variants, 98 samples

## step2.tracking.json

Second pipeline step (`WORKFLOW:FILTER_STEP`), predecessor is `WORKFLOW:PREPARE_STEP`
(the `process_name` from `step1.tracking.json`).
- Input: 450 variants, 98 samples
- Output: 400 variants, 95 samples

The chaining via `predecessor` exercises the DAG construction path in the report script.
The text report (`*_pipeline_report.txt`) emits `Process: PREPARE_STEP` and
`Process: FILTER_STEP` (short names from `full_name.split(":")[-1]`), which the test
asserts are both present.
