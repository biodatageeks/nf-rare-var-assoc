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

### Decision (2026-05-28): Option B selected

**D1 is LOCKED to Option B (`NEXTFLOW_RUN` / nested run).** Both `skip_preparation` values are
**real production paths** (neither is dev/CI-only):

- **`skip_preparation=false` — the main production path.** Actual association testing, usually run
  **on clinical-center machines** (patient data must stay within that infrastructure). This is the
  path that runs the nested `nf-prepare-vcf`. It is **not** run on PLGrid HPC.
- **`skip_preparation=true` — the parameter-tuning path on PLGrid HPC.** Tuning assoc params
  (e.g. GQ thresholds) means many runs; the un-tuned transformations were factored into
  `nf-prepare-vcf`, run **once** by a separate upstream driver, then assoc runs many times with
  `skip_preparation=true`.

The classic Option-B risk — nested Nextflow submission (sbatch-within-sbatch) being disallowed on
a cluster — therefore does **not** bite: the nested run only happens on the `skip_preparation=false`
path, which runs on clinical-center machines, **not** under slurm/HyperQueue on PLGrid HPC.

Option B is chosen for its decisive advantage: **full config/params isolation and zero edits
to `nf-prepare-vcf`** — the child loads its own `nextflow.config`, `projectDir`, and params via
`-params-file`, running byte-for-byte as it does standalone. This sidesteps the `projectDir`
foot-guns and scoped-selector work Option A requires (the make-or-break audit in old P3b).

The **Option B task list (PB1–PB4) is below**; the original **Option A breakdown (P3a–P3e) is
retained as a backup plan** further down, in case Option B's isolation cost (two work dirs,
glob-from-results) proves worse than Option A's coupling. Option C stays the last resort if both
fail.

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

- **D1 — Composition strategy.** LOCKED 2026-05-28 -> **Option B (`NEXTFLOW_RUN` / nested run)**.
  Rationale: the nested run only fires on the `skip_preparation=false` (main) path, which runs on
  clinical-center machines, **not** under slurm/HyperQueue on PLGrid HPC (HPC uses the
  `skip_preparation=true` tuning path), so the nested-submission risk is moot; Option B buys full
  param/`projectDir` isolation and zero edits to `nf-prepare-vcf`. Both `skip_preparation` values
  are production. Option A (P3a–P3e) is retained as the backup plan. See the Decision note above
  and the **Tasks (Option B — SELECTED)** list below.
- **D6 — Nested-run details.** ADDED 2026-05-28 (Option B): (a) **Resume for the child run —
  OPEN, decide later.** Options: no `-resume` (fresh each time), or nf-cascade fixed-cache
  workaround for persistent resume. (b) **Forward `cpu_support_avx2`** from this pipeline's
  params into the child ([nextflow.config:150](../../nextflow.config#L150)) — it selects the
  AVX2 vs non-AVX2 plink container in the child's `PLINK2_MAKEPGEN` (which still runs per D5).
- **D2 — Prep runs at the front.** CONFIRMED 2026-05-28:
  `NF_PREPARE_VCF -> REPLACE_SAMPLE_NAMES -> INDEX_1 -> VIEW_AND_FILTER2`, so both
  `skip_preparation` values share an identical tail (the tuner case is just a manual
  `nf-prepare-vcf` run done out-of-band).
- **D3 — Quality filtering home.** `NF_PREPARE_VCF` does **no** QUAL/GQ/DP filtering;
  `BCFTOOLS_VIEW_AND_FILTER2` is the single quality gate for both branches, applying the
  seven `filter_and_enhance_vcf_*` thresholds and replacing `BCFTOOLS_FILTER_1/2` +
  `FILTER_AND_ENHANCE_VCF`. (The threshold wiring was originally proven by the retired IT-1b;
  it is now exercised by the workflow-level IT-6/IT-7.)
- **D4 — EDA input.** CONFIRMED 2026-05-28: EDA runs on **pre-quality-filter** data — i.e.
  the `NF_PREPARE_VCF` output (CSQ + DS present), before `VIEW_AND_FILTER2`.
- **D5 — Run all of NF_PREPARE_VCF.** CONFIRMED 2026-05-28: do **not** add a `skip_downstream`
  flag for now; run the whole `nf-prepare-vcf` pipeline (COMPUTE_GENE_RANGES / VCF2PSAM /
  MAKEPGEN / LD_REPORT included) and ignore the outputs this repo doesn't need. Revisit only
  if IT wall-clock becomes a problem.

---

## Tasks (Option B — SELECTED)

> Composition via a local `NEXTFLOW_RUN`-style process (the nf-cascade pattern) that runs
> `nf-prepare-vcf` as a child `nextflow run`. **No edits to `../nf-prepare-vcf`.** The child
> loads its own config/params/`projectDir`; we hand it `--input_vcf` / `--outdir` and glob the
> published prepared VCF back into our channels. Grounding facts (verified 2026-05-28):
> child entry point is `../nf-prepare-vcf/main.nf`; the DS-carrying, CSQ-annotated output lands
> at `results/bcftools_reheader/*_reheader.vcf.gz` with its index at
> `results/bcftools_index/*_reheader.vcf.gz.tbi`.
>
> PB1 replaces Option-A P3a+P3b (no emits, no param porting, no scoped selectors). PB2/PB3/PB4
> mirror Option-A P3c/P3d/P3e. PB2..PB4 depend on PB1. PB2 must land before PB3/PB4.

### PB1 — Build the `PREPARE_VCF` wrapper process (child `nextflow run` + glob outputs)

**Status: Done **

- **Module**: `modules/local/nextflow_run/prepare_vcf/main.nf` (process `PREPARE_VCF`).
  Input: `tuple val(meta), path(vcf)` — no tbi (child starts with BCFTOOLS_SORT which
  does not need a pre-built index).
- **Test**: `modules/local/nextflow_run/prepare_vcf/tests/main.nf.test` (tag "full").
  Verifies, per child-pipeline step: DS FORMAT + CSQ INFO header lines (CSQ Format
  field order checked, since downstream parses by position); 3202 samples preserved;
  chr-prefix stripped + CHROM in {12,22,X} (ANNOTATE rename_chr); every record biallelic
  (NORM -m -any), with spanning-deletion star alleles (ALT='*') present as expected;
  ID == CHROM_POS_REF_ALT and globally unique (ANNOTATE set-id); non-empty CSQ on every
  non-star record (VEP); and **dosage correctness recomputed from the spec** — DS
  reconstructed from each genotype's PL (10^(-PL/10), normalized, P(het)+2*P(hom_alt)),
  asserting mean |DS-spec| < 1e-3 and >0.01-violation fraction < 5e-4 (only the PL=255
  saturation tail deviates; ~17/1.54M on the fixture). Plus a DP-plausibility check
  (DP>20 -> round(DS)==GT >= 99.9%; intermediate dosages from shallower calls) and a
  record-count cross-check against NORM tracking (inputs 500 / 3202 samples; output count
  == NORM outputs.variants). Requires `tests/fixtures/child_podman.config` as add_config.
- **Per-child params**: `modules/local/nextflow_run/prepare_vcf/assets/prep.yml`:
  - `skip_ld_report: true` — saves wall-clock; LD report unused by this repo.
  - `publish_intermediate: true` — required so BCFTOOLS_NORM and PLINK2_MAKEPGEN publish
    their tracking JSONs to `results/bcftools_norm/` and `results/plink2_makepgen/`.
  - Everything else comes from nf-prepare-vcf's own `nextflow.config`.
- **Runtime**: `nextflow` must be on PATH; no container on the wrapper process itself.
- **Resume strategy** (D6 resolved): no `-resume`; child runs fresh each time the parent cache
  misses. Child work dir lives at `${task.workDir}/child_work` (ephemeral, cleaned with parent).
- **Tracking glob confirmed**: `results/**/*tracking*.json` catches `*_norm_tracking.json` from
  `results/bcftools_norm/` and `*_makepgen_tracking.json` from `results/plink2_makepgen/`.
- **add_config**: optional path input; pass `[]` when no extra config is needed.
  Caller must supply a config enabling the right container manager (e.g. `podman.enabled=true`)
  for the child run; `add_config_arg` is suppressed when `[]`.
- **Done-when**: `PREPARE_VCF(ch_input_vcf)` on `unprepared_rand_500.vcf.gz` emits a DS-carrying,
  CSQ-annotated `prepared_vcf` + `.tbi` and a non-empty tracking JSON, with **zero changes to
  `../nf-prepare-vcf`**.

### PB2 — Wire `PREPARE_VCF` in; remove the old prep code

**Status: Done 2026-05-29.** Implemented with a structural change vs. the original sketch: the
`PREPARE` subworkflow was **eliminated** rather than kept as a thin shared body. Its three
remaining steps (REPLACE_SAMPLE_NAMES -> INDEX -> VIEW_AND_FILTER2) are now **inlined directly**
into `workflows/rare-var-assoc.nf`, which is simpler than a one-call wrapper subworkflow.

As built:
- **File**: `workflows/rare-var-assoc.nf`.
  - `include { PREPARE_VCF } from '../modules/local/nextflow_run/prepare_vcf'`; also inlined
    `BCFTOOLS_REPLACE_SAMPLE_NAMES`, `BCFTOOLS_INDEX`, `BCFTOOLS_VIEW_AND_FILTER2` includes.
  - **Removed** includes + call sites + guards for VEP_UPDATECACHE, FIX_ZERO_PL, BCFTOOLS_FILTER_1,
    BCFTOOLS_FILTER_2, VEP_ANNOTATE, BCFTOOLS_INDEX_2/_3, and the `PREPARE` subworkflow import.
  - Flow (D2 — prep at the front): `if (skip_preparation==false)` runs `PREPARE_VCF(ch_input_vcf,
    ch_prep_params_file)` and feeds its `prepared_vcf` into REPLACE_SAMPLE_NAMES; else feeds
    `ch_input_vcf` straight in. Both converge through INDEX -> VIEW_AND_FILTER2 into
    `ch_filtered_vcf_with_index` (renamed from `ch_vep_vcf_with_index`) consumed by VCF2FRQ /
    VCF2PSAM / MAKEPGEN_1 / VIEW_2 / VCFTOAAF / EXPORT_OTHER.
  - `use_dosage`: prepared VCF always carries `DS` (child computes it; on `skip_preparation=true`
    the manual upstream `nf-prepare-vcf` run already did — D4/item 5). Kept the downstream
    `if (use_dosage)` EXPORT/IMPORT block; dropped the dosage-compute branch.
  - Per D4, EDA runs on `ch_vcf_with_sample_names_corrected` (post-reheader + sample-rename,
    pre-VIEW_AND_FILTER2). `vep_annotated_vcf` test emit re-pointed at that VCF; added a
    `vep_annotated_vcf_tbi` emit so IT-7's naive-LOG10P helper can colocate the index.
- **Child params file**: `conf/nf_prepare_params.yml` (`skip_ld_report: true`,
  `publish_intermediate: false`). `publish_intermediate: false` is safe — the child's
  BCFTOOLS_REHEADER + BCFTOOLS_INDEX publish unconditionally (no `enabled:` guard), so the
  `prepared_vcf`/`tbi` globs always resolve; only the NORM/MAKEPGEN tracking JSONs are suppressed,
  so the wrapper's `tracking` output is now `optional: true`. (The PB1 module test keeps its own
  `assets/prep.yml` at `publish_intermediate: true` so its NORM-tracking assertions still hold.)
- **Tracking continuity**: `PREPARE_VCF.out.tracking` is threaded into `ch_tracking`, but is empty
  under `publish_intermediate: false`; the Sankey stays non-empty from VIEW_AND_FILTER2 onward
  (IT-6/IT-7 assert the report exists, not its node set). Child prep steps therefore do **not**
  appear in the production Sankey — acceptable; set `publish_intermediate: true` in the conf file
  if they are ever wanted.
- **Resourcing**: the deleted dev-only `child_podman.config` (7 GB cap) is replaced by forwarding
  the `low_resources` profile to the child via `workflow.profile`; `nf-test.config` now sets
  `profile "podman,low_resources"`.
- **Verified**: green — `nf-test test workflows/tests/skip_prep_skip_reporting.nf.test` (IT-6),
  `workflows/tests/full_reporting.nf.test` (IT-7), and
  `modules/local/nextflow_run/prepare_vcf/tests/main.nf.test` (PB1).
- **Done-when**: met. `rg -n 'FIX_ZERO_PL|BCFTOOLS_FILTER_[12]\b' workflows/` is empty.

### PB3 — Drop subworkflow-level IT-1 / IT-1b (PREPARE subworkflow eliminated)

**Status: Done 2026-05-29.** Docs updated (integration-tests.md IT-1/IT-1b marked superseded;
dead-code.md §1f updated). The retired test files still live under `prepare__to_delete/` for PB5
migration; to keep them from running and failing against the now-thin PREPARE body, they are
excluded via `ignore "subworkflows/local/prepare__to_delete/**"` in `nf-test.config`
(`nf-test list` confirms only the PB1 wrapper test remains under "prepare").

PB2 removed the `PREPARE` subworkflow, so its subworkflow-level tests (IT-1 full-prep branch,
IT-1b skip-prep branch) no longer have a target. The invariants they checked are now covered
elsewhere:
- NORM-split / unique-ID / chr-rename / DS-dosage invariants -> the **PB1 wrapper test**
  (`modules/local/nextflow_run/prepare_vcf/tests/main.nf.test`) and, ultimately, tests in
  `../nf-prepare-vcf` (PB5).
- sample-name replace + cohort subset + quality-filter (VIEW_AND_FILTER2) -> exercised by the
  **workflow-level** IT-6 (`skip_preparation=true`) and IT-7 (`skip_preparation=false`).
- The old `subworkflows/local/prepare/` tree was renamed to `subworkflows/local/prepare__to_delete/`
  and is deleted in PB5 (its test logic is being migrated upstream, not discarded).

Tasks:
- Update the IT-1 / IT-1b sections in [integration-tests.md](integration-tests.md) to mark them
  superseded (subworkflow gone; coverage moved to PB1 wrapper test + IT-6/IT-7).
- Update the §1f note in [dead-code.md](dead-code.md): FILTER_AND_ENHANCE_VCF is no longer
  referenced by any prepare-path code.
- **Done-when**: no doc still describes IT-1/IT-1b as live subworkflow tests; no test asserts
  behavior that moved to `nf-prepare-vcf`.

### PB5 — Migrate `prepare__to_delete/tests` upstream, then delete the dir

The renamed `subworkflows/local/prepare__to_delete/` dir is kept **only** to preserve its test
fixtures and assertions (`tests/main.nf.test`, `tests/skip_prep_skip_reporting.nf.test`,
`tests/fixtures/`) for reuse — the NORM/ANNOTATE/VEP invariants they encode belong to
`../nf-prepare-vcf`, which currently lacks equivalent integration coverage.
- Port the still-relevant assertions/fixtures into a `nf-prepare-vcf` integration test
  (sibling repo; the AI does not commit there — hand back the edits).
- Then delete `subworkflows/local/prepare__to_delete/` entirely (coordinate with T14).
- **Done-when**: the dir is gone; equivalent coverage lives in `../nf-prepare-vcf`.

### PB4 — Delete prep modules orphaned by the refactor

Same scope as Option-A **P3e** (under Option B nothing was imported from the sibling, so all of
this repo's now-unused prep copies can go once ref-checked). After PB2, `rg`-check and remove,
coordinating with T14:
- `modules/local/python/fix_zero_PL`
- `modules/local/combo/filter_and_enhance_vcf`
- `modules/local/python/filter_and_enhance_vcf_polarsbio` (already a §1f candidate)
- this repo's now-unused prep copies **iff** unreferenced: `bcftools/norm`, `bcftools/annotate`,
  `vep/annotate`, `vep/updatecache`, `bcftools/filter` (grep first — some may have only tests).
- **Verify**: `nf-test test --profile podman --tag ci` green; `rg` finds no dangling imports.
- **Done-when**: orphaned prep modules gone, CI green. Then proceed to T13.

---

## Tasks (Option A — BACKUP PLAN)

> **NOT the selected plan.** Retained per D1 in case Option B's nested-run isolation cost proves
> worse than Option A's coupling. Option A composes via a cross-repo `include` of
> `NF_PREPARE_VCF` (one Nextflow DAG, native channels) and pays for it with param porting,
> scoped `withName` selectors, and a `projectDir` audit.

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
