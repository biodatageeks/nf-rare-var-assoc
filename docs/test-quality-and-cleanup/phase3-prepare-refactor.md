# Phase 3 — Delegate VCF preparation to `nf-prepare-vcf` (P3a–P3e)

See [../test-quality-and-cleanup-plan.md](../test-quality-and-cleanup-plan.md) for goals and
conventions. Phase 3 runs **between T12 and T13** and **unblocks T13** (IT-7).

## Why

T13/IT-7 exercises the `skip_preparation=false, skip_reporting=false, use_dosage=false` path.
That path runs `FIX_ZERO_PL`, whose `comp_ds_htslib` Rust binary panics
(`index out of bounds: the len is 2 but the index is 2`, `bioinf_combo:1.1.1`). The fix is
architectural, not a patch to the Rust binary: the sibling pipeline `../nf-prepare-vcf`
already does VCF preparation with a maintained polars-bio implementation
(`CALC_DOSAGE_POLARSBIO`, `python_tools:1.0.11`) and supersedes the duplicated NORM /
ANNOTATE / VEP / dosage logic in this repo. Delegating to it removes `FIX_ZERO_PL` and a
large block of duplicated, drifting prep code in one move.

This was always planned as "phase 3" of the cleanups; the T13 blocker just forces it earlier.

## Composition strategy — three options (the pivotal decision)

The developer's concern is real: with a plain cross-repo `include`, `nf-prepare-vcf` does
**not** get its own params dictionary or its own `projectDir`. Below are the three viable
ways to compose the two pipelines, with the trade-off that actually matters here —
**config/params isolation vs. operability on the PLGrid HPC executor (slurm / HyperQueue)**.

### Option A — cross-repo `include` (in-process subworkflow composition)

```groovy
include { NF_PREPARE_VCF } from '../../nf-prepare-vcf/workflows/nf-prepare-vcf'
```

(Import the **workflow**, not `main.nf`: `main.nf` only defines the entry wrapper
`PSUSZYNS_NF_PREPARE_VCF`; the reusable `NF_PREPARE_VCF` lives in
`workflows/nf-prepare-vcf.nf`. Path is relative to `workflows/rare-var-assoc.nf`.)

One Nextflow DAG, native channels, one `-resume`, one executor submission path.

Friction (this is exactly what worries the developer):
- **No separate params dict.** Nextflow loads only the *running* pipeline's config.
  `NF_PREPARE_VCF` reads `params.calc_dosage_polarsbio_script_path`, `params.calc_ds_min_gq`,
  `params.plink2_makepgen_options`, etc. — undefined here unless we add them to *our* config.
- **`projectDir` re-resolves to this repo**, so `"${projectDir}/.../calc_dosage.py"` and
  `"${projectDir}/assets/rename_chr.txt"` point at the wrong tree. Must override with explicit
  sibling-repo paths.
- **`withName` selectors collide by simple name** (`BCFTOOLS_NORM`, `VEP_ANNOTATE`, …); needs
  fully-qualified `.*:NF_PREPARE_VCF:*` selectors.
- **`NF_PREPARE_VCF` emits only `versions`** — must add a `prepared_vcf` emit.
- Tight version coupling between the two repos.
- **HPC: works.** Single Nextflow run, normal slurm/hq submission. No nesting.

### Option B — `NEXTFLOW_RUN` process (nested run; the nf-cascade mechanism)

This is precisely the "module that calls `nextflow` in its script section" idea. A local
process runs the child pipeline as a subprocess:

```
nextflow run ../nf-prepare-vcf/main.nf -params-file prep.yml -c prep.config \
    --input_vcf <vcf> --outdir $task.workDir/results
```

Downstream channels then pick files out of the published results dir
(`.map { dir -> file(dir.resolve('bcftools_reheader/..._reheader.vcf.gz')) }`).
[nf-cascade](https://github.com/mahesh-panchal/nf-cascade) formalizes this pattern
(`NEXTFLOW_RUN` module, per-pipeline `*_params_file` / `*_opts` / `*_add_config`, runs the
child on the same node via `String.execute()`, dedicated child work dir for resume).

- **Full isolation — solves the developer's concern directly.** The child run loads
  `nf-prepare-vcf`'s own `nextflow.config`, its own `projectDir`, its own params (via
  `-params-file`). **Zero changes to `nf-prepare-vcf`** — no emits, no scoped selectors, no
  param porting. It runs byte-for-byte as it does standalone.
- No process-name collisions; no `projectDir` foot-guns.
- Cost: outputs arrive as a **published directory**, not native channels — we glob
  `prepared_vcf`/`tbi` (and the tracking JSONs) out of `results/`.
- Cost: **two of everything** — work dirs, trace/timeline/report, resume caches. Resume is
  the fragile part (nf-cascade works around it with a fixed child cache dir).
- Maturity: nf-cascade is an explicit **proof-of-concept**, "not tested in cloud", and needs
  `nextflow` available inside the process environment.
- **HPC: the risk.** A child `nextflow run` launched from inside a task will itself submit
  jobs via the configured executor. Nested slurm submission (sbatch-within-sbatch) is
  disallowed or discouraged on many clusters; HyperQueue's server/worker model may tolerate
  it better. **This must be validated on PLGrid before committing** — if nested submission is
  blocked, Option B is dead for production even though it is the cleanest for isolation.

### Option C — vendor `nf-prepare-vcf` modules into a local `subworkflows/local/prepare_v2`

Copy the handful of modules + the workflow body into this repo.
- Full local config control, no runtime coupling, no nested runner, HPC-safe.
- Cost: code duplication and drift from upstream — directly against the "single source of
  truth" goal that motivated using `nf-prepare-vcf` at all.

### Recommendation

Pick based on whether nested Nextflow submission is acceptable on PLGrid:

- If you want the **cleanest isolation and zero edits to `nf-prepare-vcf`**, and a quick
  PLGrid spike shows a child `nextflow run` can submit hq/slurm jobs from within a job →
  **Option B**. It is the option that actually answers "give `nf-prepare-vcf` its own params
  and `projectDir`".
- If nested submission is blocked or risky on PLGrid (likely for plain slurm) →
  **Option A**, accepting the param-porting + scoped-selector work below. It is the most
  HPC-robust and keeps one resumable DAG.
- **Option C** only if both A's coupling and B's nesting prove unworkable.

The task list below is written for **Option A** (the HPC-safe default). If you choose Option
B, P3a/P3b are replaced by a single "build the `NEXTFLOW_RUN` wrapper + per-child params
file" task and P3e shrinks (nothing to port); the rest (P3c wiring, P3d tests) stay similar.
**This is the one decision to lock before P3c** — see D1.

## Current vs target architecture

**Current `skip_preparation=false` path** (the code to remove), spread across two files:

- `subworkflows/local/prepare/main.nf`: REPLACE_SAMPLE_NAMES -> INDEX_1 -> VIEW_1 (subset to
  all_samples) -> NORM -> ANNOTATE -> (FILTER_AND_ENHANCE_VCF if `use_dosage`).
- `workflows/rare-var-assoc.nf`: FIX_ZERO_PL (if `!use_dosage && !skip_reporting`) ->
  BCFTOOLS_FILTER_1 + FILTER_2 (if `!use_dosage`) -> VEP_ANNOTATE + INDEX_3 (if
  `!skip_preparation`) -> VCF2FRQ/VCF2PSAM/MAKEPGEN... and VEP_UPDATECACHE at the top.

**`nf-prepare-vcf` `NF_PREPARE_VCF`** does: VEP_UPDATECACHE -> BCFTOOLS_SORT -> BCFTOOLS_NORM
-> BCFTOOLS_ANNOTATE -> VEP_ANNOTATE -> CALC_DOSAGE_POLARSBIO -> BCFTOOLS_REHEADER (adds the
`DS` FORMAT header) -> INDEX_2, plus downstream COMPUTE_GENE_RANGES / VCF2PSAM / MAKEPGEN /
LD_REPORT that this repo does **not** need (it redoes VCF2PSAM/MAKEPGEN/VCF2FRQ itself).

**Target** — both branches converge on the same tail (`REPLACE_SAMPLE_NAMES -> ... ->
BCFTOOLS_VIEW_AND_FILTER2`), differing only by whether `NF_PREPARE_VCF` runs in the middle:

```
skip_preparation=true  (tuner; prep done upstream out-of-band):
    input -> REPLACE_SAMPLE_NAMES -> INDEX_1 -> VIEW_AND_FILTER2 -> downstream

skip_preparation=false (this repo drives prep):
    input -> REPLACE_SAMPLE_NAMES -> INDEX_1 -> NF_PREPARE_VCF -> VIEW_AND_FILTER2 -> downstream
```

`BCFTOOLS_VIEW_AND_FILTER2` becomes the single point that subsets to the analysis cohort
(`--samples-file`) **and** applies the QUAL/GQ/DP quality filters — replacing
`BCFTOOLS_VIEW_1`, `BCFTOOLS_FILTER_1/2`, `FILTER_AND_ENHANCE_VCF`, and `FIX_ZERO_PL` at once.
After the refactor the prepared VCF always carries `DS`, so `use_dosage` no longer gates a
dosage-compute step — it only controls the downstream PLINK2_EXPORT/IMPORT dosage path and
EDA's dosage plots.

**Modules removed from `workflows/rare-var-assoc.nf`** (per the developer): VEP_UPDATECACHE,
FIX_ZERO_PL, BCFTOOLS_FILTER_1, BCFTOOLS_FILTER_2, VEP_ANNOTATE, BCFTOOLS_INDEX_3.
**Kept in PREPARE**: BCFTOOLS_REPLACE_SAMPLE_NAMES, BCFTOOLS_INDEX_1, BCFTOOLS_VIEW_AND_FILTER2.
**Dropped from PREPARE**: BCFTOOLS_VIEW_1, BCFTOOLS_NORM, BCFTOOLS_ANNOTATE,
FILTER_AND_ENHANCE_VCF (and the local INDEX_2).

## Decisions

- **D1 — Composition strategy.** OPEN, and the one blocker before P3c. Choose Option A
  (include), B (`NEXTFLOW_RUN`/nf-cascade), or C (vendor) from the assessment above. Drives
  whether the task list stays as written (Option A) or is reshaped (B). Pending a PLGrid
  spike on nested Nextflow submission if Option B is attractive.
- **D2 — Prep runs at the front.** CONFIRMED 2026-05-28:
  `NF_PREPARE_VCF -> REPLACE_SAMPLE_NAMES -> INDEX_1 -> VIEW_AND_FILTER2`, so both
  `skip_preparation` values share an identical tail (the tuner case is just a manual
  `nf-prepare-vcf` run done out-of-band).
- **D3 — Quality filtering home.** `NF_PREPARE_VCF` does **no** QUAL/GQ/DP filtering;
  `BCFTOOLS_VIEW_AND_FILTER2` is the single quality gate for both branches, applying the
  seven `filter_and_enhance_vcf_*` thresholds (IT-1b proves this) and replacing
  `BCFTOOLS_FILTER_1/2` + `FILTER_AND_ENHANCE_VCF`. (Confirm the thresholds are the intended
  gate for the full path.)
- **D4 — EDA input.** CONFIRMED 2026-05-28: EDA runs on **pre-quality-filter** data — i.e.
  the `NF_PREPARE_VCF` output (CSQ + DS present), before `VIEW_AND_FILTER2`.
- **D5 — Run all of NF_PREPARE_VCF.** CONFIRMED 2026-05-28: do **not** add a `skip_downstream`
  flag for now; run the whole `nf-prepare-vcf` pipeline (COMPUTE_GENE_RANGES / VCF2PSAM /
  MAKEPGEN / LD_REPORT included) and ignore the outputs this repo doesn't need. Revisit only
  if IT wall-clock becomes a problem.

---

## Tasks

> The breakdown below assumes **Option A (include)** — the HPC-safe default. If D1 selects
> **Option B (`NEXTFLOW_RUN`)**, P3a + P3b collapse into a single "build the `NEXTFLOW_RUN`
> wrapper + per-child params file + glob outputs" task (no emits, no param porting, no scoped
> selectors), P3c wiring stays similar, and P3e shrinks (nothing to port). Lock D1 first.

> **Cross-repo note**: under Option A, P3a edits `../nf-prepare-vcf`, a **separate git repo**.
> Per the working agreement the AI does not run git; the developer commits each repo
> separately. Call out sibling-repo edits explicitly when handing back. (Option B needs no
> edits to `../nf-prepare-vcf` at all.)

### P3a — Add prepared-VCF emits to `NF_PREPARE_VCF`  *(Option A; edits `../nf-prepare-vcf`)*

- **File**: `../nf-prepare-vcf/workflows/nf-prepare-vcf.nf`.
- **Add to `emit:`** (channels already exist in the workflow body):
  - `prepared_vcf      = ch_with_correct_header_vcf`
  - `prepared_vcf_tbi  = ch_with_correct_header_vcf_tbi`
  - `tracking          = ch_tracking`
  - keep `versions`.
- Adding emits is backward-compatible: the standalone `main.nf` entry ignores them.
- Per D5 we run the whole pipeline — **no** `skip_downstream` gating.
- **Done-when**: standalone `nf-prepare-vcf` still runs (its own nf-test, if any, green) and
  the new emits are visible to an importing workflow.

### P3b — Port params + scoped module config into this repo  *(Option A)*

- **File**: `nextflow.config` — add the params `NF_PREPARE_VCF` reads, pointing python-script
  paths at the sibling repo so they resolve regardless of `projectDir`:
  - `calc_dosage_polarsbio_script_path = "${projectDir}/../nf-prepare-vcf/modules/local/python/calc_dosage_polarsbio/assets/calc_dosage.py"`
  - `compute_ranges_py_script_path     = "${projectDir}/../nf-prepare-vcf/modules/local/python/compute_gene_ranges/assets/compute_ranges.py"`
  - `calc_ds_min_gq`, `plink2_makepgen_options`, `plink2_makepgen_vcf_input_options`,
    `plink2_ld_report_options`, `skip_ld_report`, `cpu_support_avx2`.
  - Verify the `vep_*` params already present here match `nf-prepare-vcf`'s (they do as of
    2026-05-28); reuse ours.
- **File**: `conf/modules.config` — add **fully-qualified** selectors so the imported
  processes get `nf-prepare-vcf`'s required `ext.args` without colliding with our own
  same-named processes. Copy the `ext.args` from `../nf-prepare-vcf/conf/modules.config`:
  - `withName: '.*:NF_PREPARE_VCF:BCFTOOLS_SORT'` -> `--output-type z --write-index=tbi`
  - `.*:NF_PREPARE_VCF:BCFTOOLS_NORM` -> `--rm-dup exact -m -any --do-not-normalize --output-type z --write-index=tbi`
  - `.*:NF_PREPARE_VCF:BCFTOOLS_ANNOTATE` -> `--set-id '%CHROM\_%POS\_%REF\_%ALT' --output-type z --write-index=tbi`
  - `.*:NF_PREPARE_VCF:BCFTOOLS_INDEX` -> `-t`
  - `.*:NF_PREPARE_VCF:VEP_ANNOTATE`, `:VEP_UPDATECACHE`, `:CALC_DOSAGE_POLARSBIO`,
    `:BCFTOOLS_REHEADER`, `:BCFTOOLS_VCF2PSAM`, `:PLINK2_MAKEPGEN`, `:PLINK2_LD_REPORT`,
    `:COMPUTE_GENE_RANGES` — publishDir + any args as in the upstream config.
- **Audit** `projectDir`/`baseDir` references inside the imported modules + `NF_PREPARE_VCF`
  body. `"${projectDir}/assets/rename_chr.txt"` resolves to our copy (fine, identical file);
  `"${projectDir}/../vep_cachedir"` resolves to the shared sibling cache (fine). Anything
  else that breaks must be passed in explicitly. **This audit is the make-or-break check for
  Option A** — if it surfaces references that cannot be cleanly overridden, fall back to
  Option B or C (D1).
- **Done-when**: `nextflow config -profile podman` resolves with no undefined-param errors
  and the scoped selectors are in place.

### P3c — Wire `NF_PREPARE_VCF` in; remove the old prep code

- **File**: `workflows/rare-var-assoc.nf`.
  - Add `include { NF_PREPARE_VCF } from '../../nf-prepare-vcf/workflows/nf-prepare-vcf'`
    (Option A) — or the `NEXTFLOW_RUN`-as-`PREPARE_VCF` alias (Option B).
  - **Remove** includes + call sites + guarding conditionals for: `VEP_UPDATECACHE`,
    `FIX_ZERO_PL`, `BCFTOOLS_FILTER_1`, `BCFTOOLS_FILTER_2`, `VEP_ANNOTATE`, `BCFTOOLS_INDEX_3`.
    Also remove the now-dead `trackingFirstOrEmpty`/`trackingLastOrEmpty` uses tied to them
    if they become unused.
  - New flow (per D2 — prep at the front): for `skip_preparation=false`, run
    `NF_PREPARE_VCF(ch_input_vcf)` first, then feed its `prepared_vcf`/`prepared_vcf_tbi` into
    PREPARE (REPLACE_SAMPLE_NAMES -> INDEX_1 -> VIEW_AND_FILTER2). For `skip_preparation=true`,
    feed `ch_input_vcf` straight into PREPARE. Both then converge into the existing
    `ch_vep_vcf_with_index` channel that `VCF2FRQ`/`VCF2PSAM`/`MAKEPGEN_1` consume.
  - Reconcile `use_dosage`: the prepared VCF always has `DS` now; keep the downstream
    `if (params.use_dosage)` EXPORT/IMPORT block, drop the dosage-compute branch in EDA setup.
  - Per D4, point EDA at the `NF_PREPARE_VCF` output (pre-`VIEW_AND_FILTER2`), and re-point
    the `vep_annotated_vcf` test emit (added during T13 drafting) at the prepared VCF.
- **File**: `subworkflows/local/prepare/main.nf`.
  - Keep REPLACE_SAMPLE_NAMES, INDEX_1, VIEW_AND_FILTER2; make them the single shared body for
    both `skip_preparation` values (the branch difference now lives in the outer workflow:
    whether `NF_PREPARE_VCF` ran upstream).
  - Drop VIEW_1 / NORM / ANNOTATE / FILTER_AND_ENHANCE_VCF (moved to `NF_PREPARE_VCF` +
    VIEW_AND_FILTER2). Keep the `NF_PREPARE_VCF` call in the **outer** workflow (top level),
    not inside PREPARE, so the cross-repo include sits at the top and PREPARE stays a thin
    sample-rename/subset/filter step.
- **Verify**: a manual `skip_preparation=false, skip_reporting=false, use_dosage=false` run on
  a small fixture completes through Regenie with no `FIX_ZERO_PL`.
- **Done-when**: the full path runs end to end; `rg -n 'FIX_ZERO_PL|BCFTOOLS_FILTER_[12]\b'
  workflows/` is empty.

### P3d — Rework IT-1 / IT-1b for the new PREPARE shape

- **IT-1** (`subworkflows/local/prepare/tests/main.nf.test`) currently asserts NORM-split /
  unique-ID / chr-rename invariants that now belong to `NF_PREPARE_VCF`. Options:
  - move those invariant assertions to a new `NF_PREPARE_VCF` integration test (in the
    sibling repo, or a thin wrapper test here), and
  - shrink IT-1 to what PREPARE still does (sample-name replace + subset/filter), or fold it
    into IT-1b if the two branches now share the same PREPARE body.
- **IT-1b** (skip_preparation=true) exercises VIEW_AND_FILTER2 and should be largely
  unchanged; re-run to confirm.
- Update the IT-1 / IT-1b descriptions in
  [integration-tests.md](integration-tests.md) and the §1f note in
  [dead-code.md](dead-code.md) (FILTER_AND_ENHANCE_VCF is no longer exercised by IT-1).
- **Done-when**: prepare-subworkflow tests pass against the new shape; no test still asserts
  behavior that moved to `nf-prepare-vcf`.

### P3e — Delete prep modules orphaned by the refactor

After P3c, ref-check and remove (coordinate with T14 — same `git rm` pass):
- `modules/local/python/fix_zero_PL`
- `modules/local/combo/filter_and_enhance_vcf`
- `modules/local/python/filter_and_enhance_vcf_polarsbio` (already a §1f candidate)
- this repo's now-unused prep copies **iff** nothing else references them: `bcftools/norm`,
  `bcftools/annotate`, `vep/annotate`, `vep/updatecache`, `bcftools/filter`. Grep before
  deleting — some (e.g. `bcftools/filter`) may still be used elsewhere or only have tests.
- **Verify**: `nf-test test --profile podman --tag ci` green; `rg` finds no dangling imports.
- **Done-when**: orphaned prep modules gone, CI green. Then proceed to T13.

## Risks / watch-list

- **Container drift**: `CALC_DOSAGE_POLARSBIO` pins `python_tools:1.0.11`; ensure it is
  pullable in CI (podman) like the other images.
- **Tracking JSON continuity**: downstream `GENERATE_TRACKING_REPORT` expects a tracking
  chain. `NF_PREPARE_VCF` emits `tracking` from NORM + MAKEPGEN; make sure the handoff into
  this repo's `ch_tracking` keeps the Sankey report non-empty (IT-6/IT-7 assert it exists).
- **`projectDir` surprises**: the single most likely runtime failure mode is an imported
  module resolving a `${projectDir}` path into this repo. Audit in P3b; if it bites, that is
  the trigger to fall back to vendoring (D1).
- **chrX / `--split-par`**: `nf-prepare-vcf`'s `PLINK2_MAKEPGEN` uses `--split-par b38`; only
  relevant if D5 keeps MAKEPGEN. This repo does its own chrX handling downstream, so prefer
  D5 = gate-off to avoid double PAR handling.
