# Design Doc: Business-Logic Testing for `RARE_VAR_ASSOC`

## 1. Context and Priority

This design is focused on validating business logic correctness, not only whether a process ran.

Top priority is to catch silent scientific/analytic regressions, such as:

- dosage (`DS`) computed but not actually used,
- principal components computed but not influencing regression,
- filtering steps executing but not actually filtering the intended records/samples,
- variant annotation/mask constraints not enforced.

Here we focus on workflow testing, not individual process testing. We have some process tests but they are outdated and need to be corrected, but this will be done in a separate step.

## 2. Scope

In scope:

- Workflow-level correctness for `workflows/rare-var-assoc.nf`.
- Test scenarios that validate numerical and logical behavior of key analysis steps.
- Focus on branches reached when `skip_preparation=true` and `skip_reporting=true`.

Out of scope:

- Full module unit-test refresh.
- Performance benchmarking.

## 3. Core Business Questions to Verify

Questions:

- Is the workflow correctly "wired up", are the outputs of processes used for input to downstream processes (if not final outputs)?
- Are the outputs of filtering steps (including vcf.gz and plink formats filtering, variants and samples filtering, heterozygosity filtering, kinship filtering, missingness filtering) actually filtered?
- Are high/low F-coefficient outlier samples removed according to configured thresholds?
- Are variants with annotations not permitted by masks excluded from setlist/analysis?
- Are the parameters expected to be passed to certain processes actually passed to them and used?
- Is downstream PLINK and then REGENIE input using dosage values (soft genotypes) when `use_dosage=true`, not implicitly falling back to hard `GT`?
- Can the whole pipeline find signal injected to the input vcf?
- Were principal components computed correctly?
- Are PCs actually used by Regenie, and are influencing results, when `covarColList` includes them?

Required checks:

- Differential check: run the same mini dataset with `use_dosage=true` vs `false`; assert at least one association statistic differs (same seed and deterministic settings) and/or construct a fixture containing variants with in-between dosage (so not strictly 0.0, 1.0 or 2.0) the should have caused the variant to be found as causal if dosage is used and not if it isn't. For example: all the samples have GT 0/0 (so with hard genotypes no signal will be found) but in the cases group we have a lot of samples with DS around 0.4 while in the controls all the samples have DS close to 0.0
- Cross-check a small deterministic subset against independent PCA recomputation (same filtered genotype matrix).
- Influence test: run with and without PC covariates; assert at least one phenotype/variant p-value or effect estimate changes beyond epsilon.
- Construct fixture containing known F outliers. Assert those sample IDs are present before filtering and absent after `F_COEFFICIENT_FILTERING`. Assert non-outlier controls remain.
- Construct fixture with known missingness edge cases. Assert counts and specific IDs before/after `FILTER_MISSING_PER_PHENO`.
- Construct fixture with mixed allowed/disallowed annotations and assert allowed/disallowed are present or removed from the output.

## 4. Test Strategy

## 4.1 Two-layer strategy

Layer A: deterministic semantic regression tests (primary)

- Small synthetic or tightly controlled fixtures with known expected outcomes.
- Assertions on exact IDs, values, and deltas.

Layer B: structural workflow smoke checks (secondary)

- Basic success/output presence.
- Process/branch sanity checks for debugging.

## 4.2 Why synthetic fixtures are necessary

For business logic checks, random real-data subsets are insufficient because expected truth is unknown. We need fixtures where expected outcomes are precomputed and stable.

Recommended fixture sets:

- `dosage_truth_*` for DS validation.
- `pca_truth_*` for PC validation and influence.
- `fcoef_truth_*` for F outlier removal.
- `kinship_truth_*` for relatedness removal.
- `missingness_truth_*` for missingness boundaries.
- `depth_and_quality_truth_*` for filtering based on DP and GQ.
- `mask_truth_*` for annotation inclusion/exclusion.

## 5. Scenario Matrix

Legend:

- Priority: P0 (must-have regression guard), P1 (high value), P2 (nice-to-have in this phase).
- Type: `semantic` means value/ID assertions, `smoke` means structural sanity only.

| ID | Scenario | Fixture | Run profile / key params | Assertions (pass criteria) | Type | Priority |
|---|---|---|---|---|---|---|
| S01 | DS actually influences association when enabled | `dosage_truth_basic` | Two runs on same input: `use_dosage=true` vs `false`; same seed/resources | At least one association statistic differs beyond epsilon; expected causal signal appears only (or substantially stronger) with `use_dosage=true` | semantic | P0 |
| S02 | GT fallback behavior is explicit when dosage disabled | `dosage_truth_basic` | `use_dosage=false` | No dosage-only expected signal detected; output consistent with hard genotype expectation | semantic | P0 |
| S03 | PCA numerical correctness on deterministic subset | `pca_truth_basic` | `skip_preparation=true`, deterministic subset | Pipeline PC scores match independent recomputation within tolerance | semantic | P0 |
| S04 | PC covariates influence Regenie results | `pca_truth_influence` | Two runs: with and without PCs in `covarColList` | At least one p-value or effect estimate changes beyond epsilon | semantic | P0 |
| S05 | F-coefficient outlier filtering removes expected IDs only | `fcoef_truth_edges` | thresholds set to include known outliers | Known outlier IDs absent after `F_COEFFICIENT_FILTERING`; known inliers retained | semantic | P0 |
| S06 | Missingness per phenotype filtering on boundary values | `missingness_truth_edges` | boundary threshold values (just below/at/above) | Exact expected sample IDs and counts before/after filtering | semantic | P0 |
| S07 | Kinship filtering removes related samples as configured | `kinship_truth_basic` | deterministic related pairs and cutoff | Related samples removed according to policy; unrelated controls retained | semantic | P1 |
| S08 | Depth/quality filters remove intended records | `depth_and_quality_truth_edges` | edge-case DP/GQ records | Variants/samples with failing DP/GQ removed; passing edge records retained | semantic | P1 |
| S09 | Annotation mask constraints enforced | `mask_truth_basic` | mixed allowed/disallowed annotations in same input | Disallowed annotations absent from setlist/analysis; allowed annotations present | semantic | P0 |
| S10 | End-to-end injected-signal detectability | `signal_injection_truth` | full mini-pipeline run, deterministic | Known injected signal detected (expected top hit or rank window) | semantic | P1 |
| S11 | Wiring/regression smoke for skip-prep/skip-report branch | existing small real subset | `skip_preparation=true`, `skip_reporting=true` | Run succeeds and produces expected key artifacts for debug triage | smoke | P1 |

Suggested initial subset to implement first (minimum viable semantic suite): `S01`, `S03`, `S04`, `S05`, `S06`, `S09` plus smoke `S11`.

## 6. Implementation Plan and Current Status

### 6.1 Test file layout

- Template files created:
	- `tests/workflows/rare_var_assoc_semantic_dosage.nf.test.template` -> `S01`, `S02`
	- `tests/workflows/rare_var_assoc_semantic_pca.nf.test.template` -> `S03`, `S04`
	- `tests/workflows/rare_var_assoc_semantic_filters.nf.test.template` -> `S05`, `S06` (Phase 1)
	- `tests/workflows/rare_var_assoc_semantic_masks.nf.test.template` -> `S09`
	- `tests/workflows/rare_var_assoc_semantic_signal.nf.test.template` -> `S10` (Phase 2)
	- `tests/workflows/skip_preparation_skip_reporting.nf.test.template` -> `S11` (smoke)
- Planned runnable filenames remain `*.nf.test` after fixture finalization.

### 6.2 Fixture organization

- `tests/fixtures/dosage_truth_basic/`
- `tests/fixtures/pca_truth_basic/`
- `tests/fixtures/pca_truth_influence/`
- `tests/fixtures/fcoef_truth_edges/`
- `tests/fixtures/missingness_truth_edges/`
- `tests/fixtures/mask_truth_basic/`
- `tests/fixtures/smoke_skip_preparation_skip_reporting/`
- Phase 2 planned additions:
	- `tests/fixtures/kinship_truth_basic/`
	- `tests/fixtures/depth_and_quality_truth_edges/`
	- `tests/fixtures/signal_injection_truth/`

Each fixture directory should include:

- input files (vcf/plink, phenotype/covariate, masks as needed),
- `expected/` with normalized expected outputs,
- `README.md` describing the engineered truth and exact expected deltas.

### 6.3 Assertion contracts (to keep tests stable)

- Compare normalized tabular outputs (sorted by key columns, stable numeric formatting).
- Use explicit tolerances for float comparisons (e.g. PC values, p-values, betas).
- Assert specific IDs/counts for sample/variant filters.
- For influence tests, fail on "no change" when change is expected.
- Keep one assertion file per scenario to simplify failure diagnosis.

### 6.4 Determinism requirements

- Fix random seeds where supported by tools.
- Pin tool/container versions in test config.
- Normalize locale-dependent formatting in post-processing assertions.
- Keep fixture size tiny to reduce runtime and CI variability.

### 6.5 Delivery phases

Phase 1 (P0 core, immediate):

- Implement `S01`, `S03`, `S04`, `S05`, `S06`, `S09`, `S11`.
- Goal: protect highest-risk scientific regressions.

Phase 2 (P1 expansion):

- Implement `S07`, `S08`, `S10`.

Phase 3 (maintenance hardening):

- Add stricter diagnostics in assertion failures.
- Review flaky behavior and tighten determinism further.

### 6.6 Current status

- Design intent and business questions: defined.
- Scenario matrix and priorities: defined in this document.
- Fixture scaffolding (Phase 1): initialized.
	- Created fixture directories and `expected/` placeholders for `dosage_truth_basic`, `pca_truth_basic`, `pca_truth_influence`, `fcoef_truth_edges`, `missingness_truth_edges`, `mask_truth_basic`, and smoke `smoke_skip_preparation_skip_reporting`.
	- Added fixture-specific `README.md` files with assertion contracts and required expected artifacts.
- Semantic test implementation: started as templates.
	- Added workflow test templates under `tests/workflows/*.nf.test.template` for Phase 1 and Phase 2 scenarios.
	- Templates include scenario mapping, parameter stubs, and TODO assertion points.
- Test harness configuration: refreshed.
	- Updated `tests/nextflow.config` from outdated defaults to a usable baseline aligned with `skip_preparation=true` and `skip_reporting=true` branch testing.
	- Added `tests/workflows/README.md` with finalization steps for converting templates into runnable semantic tests.
- Initial runnable chunk: implemented.
	- Added executable harness `tests/workflows/run_initial_tests.sh`.
	- Added `podman` test profile usage in harness execution.
	- Harness runs two end-to-end workflow executions (`use_dosage=true` and `use_dosage=false`) with `skip_preparation=true` and `skip_reporting=true`.
	- Harness asserts both runs produce pipeline info artifacts, required process completions in trace (header-based parsing of trace columns), and valid REGENIE table semantics (`LOG10P` present, numeric, and non-negative).
	- Harness performs row-level semantic comparison between dosage modes and fails if no REGENIE statistics change.
	- Current observed run result: semantic comparison succeeded with many changed rows between dosage modes, confirming sensitivity to `use_dosage` in this branch.
- Runnable nf-test specs (`*.nf.test`) with final semantic assertions: pending fixture input finalization and expected-output materialization.
- Smoke test refresh for skip-prep/skip-report: template prepared; runnable assertions pending.
- Module-level process test refresh: intentionally deferred (out of scope here).

## 7. Important Implementation Notes

- Prefer assertions on data semantics (IDs, numeric values, deltas) over process-presence assertions.
- Keep fixtures tiny but engineered to trigger each behavior unambiguously.
- Avoid non-determinism: fixed seeds, sorted outputs, normalized formatting before comparisons.
- For value-based checks, include tolerances and explicit rationale.
- Treat "no change when change is expected" as failure for influence tests (dosage and PCs).
- `assets/medium_data/1kGP_cases_200.txt`, `assets/medium_data/1kGP_controls_400.txt` and `assets/medium_data/prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz` file can be used as input for some tests directly or after selecting a subset. `prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz` contains 5000 variants for 3200 samples (small subset from 1000 genomes dataset)

## 8. Suggested Execution Approach

Run a focused semantic test suite first, then smoke tests:

```bash
# initial runnable chunk (smoke + dosage differential guard)
tests/workflows/run_initial_tests.sh

# rename selected templates from *.nf.test.template to *.nf.test after fixture finalization
# example:
# mv tests/workflows/rare_var_assoc_semantic_dosage.nf.test.template tests/workflows/rare_var_assoc_semantic_dosage.nf.test
nf-test test tests/workflows/*semantic*.nf.test
nf-test test tests/workflows/skip_preparation_skip_reporting.nf.test
```

The exact pattern can be adjusted after test file naming is finalized.

## 9. Interpreting Current Delta Metrics

Current harness prints row-level differences for paired REGENIE outputs:

- `p_true`: `LOG10P` from run with `use_dosage=true`.
- `p_false`: `LOG10P` from run with `use_dosage=false`.
- `delta_p`: absolute difference `|p_true - p_false|`.
- `delta_beta`: absolute difference in effect estimate where `BETA` is available and numeric; missing values (e.g. `NA`) are treated as unavailable, not as failure.

Interpretation guidance for this stage:

- This is currently a **sensitivity guard**, not yet a strict truth check: we verify that changing dosage mode changes downstream association statistics.
- Large `delta_p` means notable evidence-shift in association significance between modes.
- Large `delta_beta` means effect estimate depends materially on dosage handling.

What remains to make this a full semantic oracle (next iterations):

- Introduce engineered fixture truth for specific expected changed rows and expected unchanged rows.
- Replace generic “at least one change” with thresholded acceptance criteria bound to scenario contracts (`S01`, then `S03/S04`).
