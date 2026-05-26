# Test Quality and Cleanup Plan (Design)

Date: 2026-05-21 (updated 2026-05-25)

This is the index. Detailed plans for each work area live in
[test-quality-and-cleanup/](test-quality-and-cleanup/). When working on a task, read this
file plus the sub-doc the task points to — you don't need to load the others.

## Working agreement

- One task at a time. Finish it, update status here and in the sub-doc, then stop for
  user review. Don't chain tasks.
- If a task is too large, propose a split into the sub-doc and hand back before doing
  the work.
- Keep design docs terse: they are loaded into AI context every session. State each rule,
  status, and decision once. Delete superseded resume notes and obsolete detail as soon
  as the task lands. Cite source-of-truth files (paths/line numbers) instead of duplicating
  them.
- **No git usage by the AI**. After completing a task the user will commit the code manually as part of reviewing the task.
- **ASCII only in code**. Do not use non-ASCII characters in any source file (`.nf`, `.py`,
  `.R`, `.config`, `.json`, `.yml`, shell scripts, test files). Non-ASCII is fine in design
  docs (`.md`). If non-ASCII is necessary for some reason, ask the user before using it.

## Goals

1. Stabilize and fix all failing nf-test unit tests.
2. Add missing unit tests for modules with non-trivial logic (Python/R); skip tests for
   dead code and modules that merely wrap an external tool with no custom logic.
3. Design and implement integration tests for subworkflows and the main workflow — the
   most valuable tests in this codebase.
4. Remove dead code (unused modules, benchmark workflow, orphan test files).
5. Update `nextflow_schema.json` for all params defined in `nextflow.config` but absent
   from the schema.
6. Remove nf-core template boilerplate comments throughout the codebase.

**Order of work**: fix failing tests → remove dead-code imports + rename dead dirs → write
missing unit tests for complex logic → write integration tests → update schema → (after all
integration tests pass) delete the renamed `__to_delete` directories → remove comments.
Do not write tests for modules confirmed for deletion.

---

## Testing philosophy

Module-level tests that simply invoke an external tool (plink, bcftools) have low value —
we are not testing whether those tools work. The high-value tests are:
(a) modules containing complex Python or R logic that we wrote, and
(b) subworkflow and workflow integration tests that verify channel wiring and business logic.
Tests for wrapper-only modules should be marked low-priority and run only in the full suite.
All tests must run with `nf-test test --profile podman`. Stub tests are not used — they
verify nothing meaningful.

### Quality bar for assertions

Every test must assert at least one *business-meaningful* property of its output — a
property that would change if the code under test were broken in a realistic way. The
following are **process-success checks in disguise** and do not count as meaningful
assertions on their own:

- "File exists" / "file is non-empty" / "file size > N KB"
- "Output contains the literal string `<html`" or other 4–6 byte markers a templating
  engine emits regardless of content
- "At least one row/variant/sample appears in the output" (trivially true unless the
  process errored, which `process.success` already checks)
- Snapshot-matching the entire output when the snapshot has never been audited for
  correctness

Each test must include at least one assertion that names a concrete invariant: a row count
matched to an input row count, a specific value in a known cell, a column the consumer
relies on, a sample/variant ID that *should* (or *should not*) appear given the test
fixture. When fixing a previously-failing test, audit existing assertions against this bar
before declaring the test "done"; do not leave a green test that only checks
`process.success`. An example of meaningful assertions in a test is the test developed as 
part of T3c task: modules/local/regenie/step1/tests/main.nf.test

### Reporting carve-out

The quality bar above is intentionally relaxed for modules and subworkflows whose only
output is a human-readable artifact (PNG/SVG plots, HTML reports, Sankey diagrams).
Reporting is not a top-priority correctness surface and is easy to verify manually; writing
automated assertions that validate plot contents or HTML structure is disproportionately
costly. For reporting outputs we accept smoke-level checks: "the expected set of output
files is produced" (asserting filenames, not content) and "versions.yml lists the expected
dependencies". Do not invest in parsing HTML, validating PNG byte content, or asserting
plot semantics.

Reporting-tagged tests: §5a `python/eda`, `python/draw_pc_plot`,
`python/generate_tracking_report` (see [unit-tests.md](test-quality-and-cleanup/unit-tests.md));
IT-5 reporting subworkflow; plot-file and HTML assertions inside IT-2 and IT-7 (see
[integration-tests.md](test-quality-and-cleanup/integration-tests.md)).

---

## §0 Conventions

These conventions exist so individual tasks can be terse without losing precision. If a
task contradicts a convention, the task wins.

### 0.1 Test file locations

| Test type | Path pattern |
|---|---|
| Module unit test | `modules/local/<group>/<name>/tests/main.nf.test` |
| Subworkflow integration test | `subworkflows/local/<name>/tests/main.nf.test` |
| Workflow integration test | `workflows/tests/<scenario_name>.nf.test` |
| Snapshot file (nf-test default) | sibling `main.nf.test.snap` next to the test |

Snapshot files (`*.nf.test.snap`) are **gitignored** and not committed. Use explicit
assertions instead of snapshots for all correctness checks. The only acceptable use of
`snapshot()` is for `versions.yml` (a small, deterministic, human-auditable YAML); all
other output properties must be asserted directly.

### 0.2 Required invocation

All tests must be runnable with:

```
nf-test test --profile podman <path/to/main.nf.test>
```

The full suite is `nf-test test --profile podman`. Never use `--shard` (caused prior
failures). Never write stub tests (`process "STUB"`, `workflow "STUB"`).

### 0.3 Tagging convention

Every test declares one of two nf-test tags inside the `test(...)` block:

| Tag | Meaning | What runs in CI |
|---|---|---|
| `tag "ci"` | Module unit tests, plus integration tests with a target wall-clock ≤5 min each | yes |
| `tag "full"` | Slow (>5 min) or low-value — wrapper-only modules, IT-7 | no, manual only |

CI command: `nf-test test --profile podman --tag ci`.
Full nightly: `nf-test test --profile podman`.

Wall-clock budgets are aspirational, not hard gates.

### 0.4 Canonical fixtures

**Unprepared fixtures** — raw 1kGP VCFs: duplicate variant IDs, multiallelic sites, no VEP
annotation. Use these for tests that exercise the preparation steps themselves (bcftools
norm/annotate, VEP) or for integration tests that run the full pipeline.

| Use case | Fixture path |
|---|---|
| Primary VCF input (3 chr inc. X, 3202 samples, ~500 variants) | `assets/three_chr_unprepared/unprepared_rand_500.vcf.gz` |
| Same dataset, larger when Regenie complains about low variance | `assets/three_chr_unprepared/unprepared_rand_{1k,2k,5k,10k}.vcf.gz` |
| Pre-prepared pgen/pvar/psam for downstream subworkflow tests | `assets/three_chr_unprepared/prepared_500/` *(created in T8)* |
| Cases/controls sample lists matching the 3202-sample fixture | `assets/three_chr_unprepared/cases.txt`, `controls.txt` *(created in T8)* |

**Prepared fixtures** — already through bcftools norm + VEP annotation: unique variant IDs,
biallelic only, CSQ-annotated. Use these for tests of downstream modules (plink2, regenie,
etc.) where the module under test should not need to handle messy raw-VCF quirks. Note that
`--write-snplist` and similar plink2 commands **require** unique variant IDs and will refuse
to run on unprepared data.

| Use case | Fixture path |
|---|---|
| Full prepared 3-chr VCF (5000 variants, chr12+chr22+chrX, 3202 samples) | `assets/medium_data/prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz` |
| Small chr12-only subset (100 variants, 3202 samples) for fast module tests | `assets/medium_data/prepared_chr12_100.vcf.gz` |

**Other fixtures**

| Use case | Fixture path |
|---|---|
| Local cache of nf-core sarscov2 test VCF (9 records, no multiallelic sites) | `assets/sarscov2/test.vcf.gz` |
| Local cache of nf-core sarscov2 genome FASTA (chromosome MT192765.1, ~30 kb) | `assets/sarscov2/genome.fasta.gz` + `.fai` + `.gzi` |
| Hand-crafted edge-case inputs for one specific test | `<test-dir>/fixtures/` next to the `.nf.test` file |

Do not create test-specific VCFs unless none of the above can be coerced into the scenario.
Document any new fixture in a `README.md` next to it. **Do not use the legacy top-level
`tests/fixtures/` directory** — it is deleted in T2.5
(see [dead-code.md §1g](test-quality-and-cleanup/dead-code.md)).

### 0.5 What `tracking` actually emits

The five subworkflows (`prepare`, `pca`, `filter_missing_per_pheno`,
`f_coefficient_filtering`, `reporting`) all emit a channel named `tracking` (not
`tracking_out` — that is the per-module output name). When asserting tracking-JSON content,
read the file with `path(workflow.out.tracking[0])`.

### 0.6 Developer tooling

`bcftools` is **not installed** on the development machine and is not going to be. To inspect VCF files locally, use:

```bash
podman run --rm -v $PWD/:/wd/:z quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    sh -c "bcftools <subcommand> /wd/<path>"
```

This is the same container image used by the `BCFTOOLS_*` modules in tests.
The same is true for samtools, plink and other such tools, which are available in respective containers (plink2 in docker.io/psuszynski/plink:2.0-alpha.6.9, samtools in biocontainers/samtools:v1.9-4-deb_cv1).

---

## Task list

Each task links to the sub-doc containing its detail (assertions, files to touch, verify
command, done-when). Conventions in §0 are assumed throughout.

| Task | Status | Detail in |
|---|---|---|
| T1 — Confirm dead-code list | ✅ Done 2026-05-24 | [dead-code.md](test-quality-and-cleanup/dead-code.md) |
| T2 — Phase 1 dead-code removal | ✅ Done 2026-05-24 | [dead-code.md](test-quality-and-cleanup/dead-code.md) |
| T2.5 — Delete abandoned `tests/` scaffold | ✅ Done 2026-05-25 (partial revert — see T2.6) | [dead-code.md](test-quality-and-cleanup/dead-code.md) |
| T2.6 — Prune `tests/nextflow.config` contents | ✅ Done 2026-05-25 | [dead-code.md](test-quality-and-cleanup/dead-code.md) |
| T3 — Fix Cat-B failures (wrong process name / path) — split into T3a–T3d on 2026-05-25 | ✅ Done 2026-05-25 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T3a — `bcftools/filter` (process rename + tuple reshape) | ✅ Done 2026-05-25 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T3b — `rscript/buildreports` (strip dead setup, use fixtures) | ✅ Done 2026-05-25 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T3c — `regenie/step1` (rewrite against current production shape) | ✅ Done 2026-05-25 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T3d — `regenie/step2` (rewrite against current production shape) | ✅ Done 2026-05-25 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| T4 — Fix Cat-A failures (signature mismatches) — split into T4a–T4g on 2026-05-25 | pending | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T4a — `cmds/merge_results` (tuple consolidation; assertions already meet bar) | ✅ Done 2026-05-26 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T4b — `bcftools/annotate` (collapse 5 paths into tuple; many test variants) | ✅ Done 2026-05-26 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T4c — `bcftools/norm` (collapse 4 paths into tuple; large test file) | ✅ Done 2026-05-26 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T4d — `bcftools/view` (collapse 7 paths into tuple) | ✅ Done 2026-05-26 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T4e — `plink2/makebed` (8-path tuple + 4 vals; **unblocks T5**) | ✅ Done 2026-05-26 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T4f — `vep/annotate` (4-path tuple; needs `../vep_cachedir` reference) | ✅ Done 2026-05-26 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| &nbsp;&nbsp;T4g — `rscript/manhattan_qq_plots` (11-element tuple; **deferred — needs IT-5 fixtures from T11**) | deferred | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| T5 — Fix Cat-D transitive failure (`plink2/write_snplist`) | ✅ Done 2026-05-26 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| T6 — Re-record Cat-E snapshot (`bcftools/vcf2frq`) | ✅ Done 2026-05-26 | [failing-tests.md](test-quality-and-cleanup/failing-tests.md) |
| T7 — Write high-priority unit tests (§5a) — split into T7a–T7f on 2026-05-26 | pending | [unit-tests.md](test-quality-and-cleanup/unit-tests.md) |
| &nbsp;&nbsp;T7a — `python/calc_f_outliers` (hand-rolled `het` fixture + FID/IID branch) | pending | [unit-tests.md](test-quality-and-cleanup/unit-tests.md) |
| &nbsp;&nbsp;T7b — `python/vcf2aaf` (hand-rolled VCF + AF fallback logic) | pending | [unit-tests.md](test-quality-and-cleanup/unit-tests.md) |
| &nbsp;&nbsp;T7c — `bcftools/assign_annotations` (hand-rolled VCF + masks TSV) | pending | [unit-tests.md](test-quality-and-cleanup/unit-tests.md) |
| &nbsp;&nbsp;T7d — `python/eda` (smoke; reuses `unprepared_rand_500.vcf.gz`) | pending | [unit-tests.md](test-quality-and-cleanup/unit-tests.md) |
| &nbsp;&nbsp;T7e — `python/draw_pc_plot` (smoke; hand-rolled `.sscore` + pheno) | pending | [unit-tests.md](test-quality-and-cleanup/unit-tests.md) |
| &nbsp;&nbsp;T7f — `python/generate_tracking_report` (smoke; hand-rolled tracking JSONs) | pending | [unit-tests.md](test-quality-and-cleanup/unit-tests.md) |
| T8 — Implement IT-1 + IT-1b; record `prepared_500/` fixture | pending | [integration-tests.md](test-quality-and-cleanup/integration-tests.md) |
| T9 — Implement IT-2 (PCA) | pending | [integration-tests.md](test-quality-and-cleanup/integration-tests.md) |
| T10 — Implement IT-3 and IT-4 in parallel | pending | [integration-tests.md](test-quality-and-cleanup/integration-tests.md) |
| T11 — Implement IT-5 (reporting subworkflow) | pending | [integration-tests.md](test-quality-and-cleanup/integration-tests.md) |
| T12 — Implement IT-6 (workflow, fast path) | pending | [integration-tests.md](test-quality-and-cleanup/integration-tests.md) |
| T13 — Implement IT-7 (workflow, reporting path) | pending | [integration-tests.md](test-quality-and-cleanup/integration-tests.md) |
| T14 — Phase 2 dead-code removal | pending | [dead-code.md](test-quality-and-cleanup/dead-code.md) |
| T15 — Update `nextflow_schema.json` | pending | [schema-and-config.md](test-quality-and-cleanup/schema-and-config.md) |
| T16 — Remove nf-core template comments | pending | [dead-code.md](test-quality-and-cleanup/dead-code.md) |

---

## Task dependency graph

```
T1 ──► T2 ──► T2.5
                 ├──► T3 ──┐
                 ├──► T4 ──┼──► T5 ──┐
                 └──► T6   │         │
                           └────►    │
                                T7   │
                                T8 ──┼──► T9
                                     ├──► T10
                                     ├──► T11
                                     └──► T12 ──► T13 ──► T14 ──► T15 ──► T16
```

T3/T4/T5/T6/T7/T8 are independent of each other and may be parallelised. T9..T12 all depend
on T8 (they need `prepared_500/`). T14 is the only task that strictly requires every
preceding test task to be green.

T4 sub-tasks T4a..T4f are independent of each other. T4e (`plink2/makebed`) is a strict
prerequisite for T5. T4g (`rscript/manhattan_qq_plots`) is deferred until T11 produces the
IT-5 reporting fixtures.
