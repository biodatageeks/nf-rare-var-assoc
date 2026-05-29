# T17 — GitHub Actions CI + README update — ✅ Done 2026-05-29

See [../test-quality-and-cleanup-plan.md](../test-quality-and-cleanup-plan.md) for goals and
conventions. T17 runs after the test work is stable (logically near T16); it is independent of
PB4/PB5 and the schema task and can be done at any point once the `--tag ci` suite is green.

## What was done

- Created `.github/workflows/ci.yml`: triggers on push/PR to master + `workflow_dispatch`;
  installs Nextflow 25.10.2 via `nf-core/setup-nextflow@v2`; caches and installs nf-test
  0.9.3 via `https://get.nf-test.com`; runs `nf-test test --profile podman,low_resources --tag ci`.
- Rewrote `README.md`: real pipeline description, required params table, key toggles, two
  production modes, Testing section with `--tag ci` / full commands. Dropped dead badges
  (Zenodo, Seqera launch, linting). Updated Nextflow badge to >=25.10.2. CI badge now
  points at the real workflow.

## Why

The README already advertises CI badges (`.github/workflows/ci.yml`, `linting.yml`) but **no
such workflows exist** — the badges are dead nf-core boilerplate. There is no automated test
run on push/PR. T17 makes CI real: a GitHub Actions workflow that runs the fast test suite
(`nf-test test --profile podman --tag ci`, per §0.3), and brings the README in line with what
the pipeline actually is and how CI actually runs.

## Part A — GitHub Actions CI workflow

- **File**: `.github/workflows/ci.yml` (new).
- **Triggers**: `push` to `master` and `pull_request` targeting `master` (and the active
  feature branch while developing, optional). Add `workflow_dispatch` for manual runs.
- **Runner**: `ubuntu-latest`. Podman is preinstalled on GitHub's ubuntu runners, matching the
  `--profile podman` convention in §0.2 — no need to switch the suite to Docker.
- **Steps**:
  1. `actions/checkout@v4`.
  2. Install Nextflow (>=24.04.2 per README badge) — `nf-core/setup-nextflow` action or a
     direct `curl -s https://get.nextflow.io | bash` + cache.
  3. Install nf-test (`nf-test`'s install script or a setup action), pinned to a known version.
  4. Run `nf-test test --profile podman --tag ci`.
- **Scope — only `tag "ci"` runs** (§0.3). The `tag "full"` tests are excluded: they are slow
  and, critically, the `PREPARE_VCF` wrapper test (`modules/local/nextflow_run/prepare_vcf`)
  spawns a **nested `nextflow run`** and needs `nextflow` on PATH inside the process — out of
  scope for CI. Keep those manual/nightly. (A nightly `workflow_dispatch`/`schedule` job that
  runs the full suite is optional, not required for T17.)
- **Container pull caveat**: the `--tag ci` tests pull biocontainer/quay/docker.io images
  (plink2 at `docker.io/psuszynski/plink:2.0-alpha.6.9`, etc.). Confirm all are public and
  pullable from a clean runner; CI has no local image cache. Consider caching the podman image
  store between runs to cut wall-clock.
- **Done-when**: the workflow runs on PR/push, executes `--tag ci`, and goes green; the
  README CI badge resolves to a real run (not 404).

> **Linting badge** (`linting.yml`): the README also shows an nf-core *linting* badge. T17 may
> either add a minimal linting workflow or drop that badge (coordinate with T16, which removes
> nf-core boilerplate). Decide when implementing; not required for the CI run itself.

## Part B — README update

The current `README.md` is almost entirely nf-core template boilerplate (samplesheet/fastq
examples that do not apply, `MultiQC` mention, TODO comments, dead Zenodo/Seqera badges). Bring
it in line with the real pipeline. Coordinate with **T16** (remove nf-core template comments) so
the two passes do not fight — T17 owns README *content*, T16 owns stray template comments
repo-wide; doing them together for the README is fine.

- Replace the placeholder Introduction with a real 2-3 sentence description: a rare-variant
  association pipeline ingesting a (multi-sample) VCF + case/control sample lists, running
  preparation (delegated to `nf-prepare-vcf`), PCA, per-phenotype missingness and F-coefficient
  filtering, Regenie step1/step2 association, and a reporting/tracking artifact.
- Fix the Usage section: drop the fastq/samplesheet example (irrelevant); document the real
  required params (`--input_vcf`, `--input_cases`, `--input_controls`, `--outdir`,
  `--project_name`) and the key toggles (`skip_preparation`, `skip_reporting`, `use_dosage`).
- Document the two production modes (see [[project_hpc-execution-model]] in memory):
  `skip_preparation=false` (main path, runs nested `nf-prepare-vcf`) vs `skip_preparation=true`
  (PLGrid HPC param-tuning, prep done out-of-band upstream).
- Add a **Testing** section: `nf-test test --profile podman --tag ci` (fast/CI suite) vs
  `nf-test test --profile podman` (full suite), and note the `ci`/`full` tag convention.
- Remove or correct the dead badges (Zenodo `XXXXXXX`, Seqera launch, conda/singularity if not
  actually supported). Keep the CI badge once Part A makes it real.
- Strip the nf-core `<!-- TODO nf-core: ... -->` blocks (with T16).

- **Done-when**: README describes this pipeline (no fastq/MultiQC/samplesheet boilerplate),
  documents the real run command + test commands, and every remaining badge points at something
  real.
