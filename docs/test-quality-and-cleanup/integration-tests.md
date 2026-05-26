# Integration Tests (T8 — T13)

See [../test-quality-and-cleanup-plan.md](../test-quality-and-cleanup-plan.md) for goals,
conventions, the **quality bar for assertions**, and the **reporting carve-out**.

**Why these matter most**: subworkflow and workflow tests exercise channel wiring, business
logic, and data transformations that are impossible to cover with module-level tests.

All integration tests must:
- Run with `--profile podman` and real (non-stub) execution.
- Use `assets/three_chr_unprepared/unprepared_rand_500.vcf.gz` (or a 500-variant-derived
  asset) as the primary fixture, unless the scenario requires a synthetic edge case.
- Have assertions focused on business outcomes, not just `process.success`.
- Tag `"ci"` unless explicitly marked `"full"` below.

**Subworkflow emit names**: `prepared_vcf`, `prepared_vcf_tbi`, `versions`, `tracking` (PREPARE);
`sscore`, `plot_file`, `king_cutoff_prune_in`, `tracking` (PCA);
`pgen_pvar_psam_out`, `tracking` (FILTER_MISSING_PER_PHENO, F_COEFFICIENT_FILTERING).
Read tracking JSON with `path(workflow.out.tracking[0]).text` and `JsonSlurper`.

---

## IT-1 — `subworkflows/local/prepare` (full preparation branch)

- **Test file**: `subworkflows/local/prepare/tests/main.nf.test`
- **Inputs**: `unprepared_rand_500.vcf.gz` (+ `.tbi`),
  `assets/three_chr_unprepared/cases.txt` (created in T8), VEP cache reuse via existing
  test config.
- **Run with**: `params.skip_preparation=false`, `params.skip_reporting=false`,
  `params.use_dosage=true` (this exercises NORM → ANNOTATE → FILTER_AND_ENHANCE_VCF).
- **Phase-2 caveat**: `FILTER_AND_ENHANCE_VCF` is a phase-2 deletion candidate (§1f). When
  the VCF-prep refactor lands, this test will need updating to either drop the
  `use_dosage=true` exercise from here or point at the replacement upstream module.
- **Assertions** on `workflow.out.prepared_vcf`:
  - `bcftools view -H <vcf> | awk '$5 !~ /,/' | wc -l` equals `bcftools view -H <vcf> | wc -l`
    (no multiallelics after NORM)
  - `bcftools view -H <vcf> | awk '{print $3}' | sort -u | wc -l` equals the variant count
    (unique IDs after ANNOTATE — load-bearing invariant per the comment in
    [prepare/main.nf:76-78](../../subworkflows/local/prepare/main.nf#L76-L78))
  - `bcftools view -h <vcf>` contains `##INFO=<ID=CSQ,` (VEP annotated)
  - Sample count in the output equals 3202 (sample-preservation invariant — replace_sample_names
    must not drop samples)
  - `prepared_vcf_tbi` is emitted
- **Tracking-JSON assertion**: keys for `bcftools_replace_sample_names`, `bcftools_norm`,
  and `bcftools_annotate` exist; at least one has a non-zero `variants` count; the
  `samples` count in the final entry equals 3202.
- **Tag**: `"ci"`.

## IT-1b — `subworkflows/local/prepare` (skip_preparation + skip_reporting branch)

This branch is the main production path used by `nextflow-gene-assoc-tuner` for parameter
optimization (see the comment at
[prepare/main.nf:38](../../subworkflows/local/prepare/main.nf#L38)). It uses a completely
different code path (`BCFTOOLS_VIEW_AND_FILTER2`) and must have its own coverage.

- **Test file**: `subworkflows/local/prepare/tests/skip_prep_skip_reporting.nf.test`
- **Inputs**: a **pre-prepared** VCF (because `skip_preparation=true` means the caller
  guarantees NORM/ANNOTATE/VEP have already run upstream). Use
  `assets/medium_data/prepared_chr12_100.vcf.gz` (3202 samples, 100 variants, chr12); same
  1000G sample IDs as `cases.txt`/`controls.txt`. The joined all.samples fixture lives at
  `subworkflows/local/prepare/tests/fixtures/all_samples.txt`.
- **Run with**: `params.skip_preparation=true`, `params.skip_reporting=true`,
  `params.use_dosage=true`. The seven filter thresholds are overridden in a sibling
  `skip_prep_skip_reporting.config` so each is active against the fixture (variant-level
  filters drop a known subset; sample-level filters with `--set-GTs '.'` rewrite bad
  genotypes to missing).
- **Assertions** on `workflow.out.prepared_vcf` (all derived from one zcat of the output):
  - Output variant count equals the deterministic survivor count from the fixture +
    threshold combination (currently 82 of 100)
  - A specific high-QUAL variant ID is present in the output, a specific low-QUAL variant
    ID is absent — pins identity-level filter wiring, not just count
  - Sample count in the output equals 3202, and a known sample ID (`HG00096`) appears in
    the `#CHROM` header
  - Missing-genotype count (`./.`) in the output body is well above the input baseline,
    proving the sample-level `--set-GTs '.'` filter actually ran
  - `prepared_vcf_tbi` is emitted with a `.tbi` filename suffix
- **Tracking-JSON assertion**: `process_name` references `BCFTOOLS_VIEW_AND_FILTER2`;
  `inputs.variants`/`outputs.variants`/`inputs.samples`/`outputs.samples` cross-check
  against the VCF-derived counts above.
- **Tag**: `"ci"`.

## IT-2 — `subworkflows/local/pca`

- **Test file**: `subworkflows/local/pca/tests/main.nf.test`
- **Inputs**: pgen/pvar/psam recorded from IT-1 output (see T8 follow-up step). Committed
  to `assets/three_chr_unprepared/prepared_500/{prepared_500.pgen,.pvar,.psam}`.
- **Run with**: default `params.plink2_pca_settings`,
  `regenie_step1_kinship_filtering=true`.
- **Assertions**:
  - `workflow.out.sscore` first line contains all `PC1..PCN` columns derived from
    `params.plink2_pca_settings` (assert the exact N — not just "contains PC1")
  - `workflow.out.king_cutoff_prune_in` is a non-empty list of sample IDs
  - **Channel-wiring**: `tracking` JSON contains entries for
    `plink2_indep_pairwise → plink2_king_cutoff → plink2_pca → plink2_projection_score`
    in that order (the projection step must consume the eigenvec_allele output, which is
    the most subtle wiring in this subworkflow — see
    [pca/main.nf:92-100](../../subworkflows/local/pca/main.nf#L92-L100))
  - `tracking` for `plink2_king_cutoff` has `samples_before > samples_after`
- **Reporting carve-out**: `workflow.out.plot_file` is emitted with the expected filename;
  do not assert byte size or plot contents.
- **Tag**: `"ci"`.

## IT-3 — `subworkflows/local/filter_missing_per_pheno`

- **Test file**: `subworkflows/local/filter_missing_per_pheno/tests/main.nf.test`
- **Inputs**: pgen/pvar/psam from `assets/three_chr_unprepared/prepared_500/` + a
  hand-crafted phenotype file with **two phenotypes**, where ~20% of samples are present
  for one phenotype only. Author at
  `subworkflows/local/filter_missing_per_pheno/tests/fixtures/two_pheno.tsv` with a short
  `README.md` documenting the missingness pattern.
- **Assertions**:
  - `workflow.out.pgen_pvar_psam_out` is a list of 2 entries (one per phenotype) — this
    tests the per-phenotype `flatMap` split in
    [filter_missing_per_pheno/main.nf:33-39](../../subworkflows/local/filter_missing_per_pheno/main.nf#L33-L39)
  - Sample counts in the two output `.psam` files differ by the expected ~20% margin
  - Both output entries' `meta.id` matches the source cohort `meta.id` (verifies the
    `orig_id` round-trip in [main.nf:50-53](../../subworkflows/local/filter_missing_per_pheno/main.nf#L50-L53))
  - `tracking` JSON contains per-phenotype `samples_before` / `samples_after` entries
- **Tag**: `"ci"`.

## IT-4 — `subworkflows/local/f_coefficient_filtering`

- **Test file**: `subworkflows/local/f_coefficient_filtering/tests/main.nf.test`
- **Inputs**: pgen/pvar/psam from `assets/three_chr_unprepared/prepared_500/`, with ≥2
  artificial high-F outlier samples. Commit either the spiked pgen/pvar/psam OR a list of
  outlier IIDs the test can re-inject under
  `subworkflows/local/f_coefficient_filtering/tests/fixtures/` with a short `README.md`
  describing the spike procedure.
- **Run with**: `params.inbreeding_outliers_range_stds=3` (tight enough to catch the spike).
- **Assertions**:
  - The exclusion output (`calc_f_outliers` outliers file) contains every spiked sample IID
    (not just any N IIDs — pin to the IIDs you spiked)
  - `workflow.out.pgen_pvar_psam_out` `.psam` has exactly `N_input - N_spiked` samples
- **Tag**: `"ci"`.

## IT-5 — `subworkflows/local/reporting` *(reporting carve-out)*

- **Test file**: `subworkflows/local/reporting/tests/main.nf.test`
- **Inputs**: pre-computed (and committed) Regenie step-2 output for one phenotype
  (`.regenie.gz`) + a tiny annotations file + setlist + phenotype file. Place under
  `subworkflows/local/reporting/tests/fixtures/`.
- **Assertions** (smoke-only per carve-out):
  - Process succeeds.
  - The expected output report file is emitted with the expected filename.
  - Do not parse HTML or validate plot contents.
- **Tag**: `"ci"`.

## IT-6 — `workflows/rare-var-assoc` — fast path

- **Test file**: `workflows/tests/skip_prep_skip_reporting.nf.test`
- **Config**: `skip_preparation=true`, `skip_reporting=true`, `use_dosage=true`.
- **Inputs**: same shape as `conf/test_skip_preparation_and_reporting.config` but pointed
  at the 500-variant fixture set (pre-prepared VCF + psam + phenotype file + masks).
  Commit any missing pieces to `assets/three_chr_unprepared/prepared_500/`.

### File-shape and presence assertions

The fixture has **one phenotype** (`Y1`, derived from cases/controls).

- Regenie step 1 emits exactly one `*_pred.list` (a per-phenotype LOCO manifest) and one
  `*_1.loco` file (one LOCO file per phenotype; chromosomes are stored inside).
- Regenie step 2 emits exactly one **`*_Y1.regenie`** file (uncompressed; covers all
  chromosomes in a single file because `--chr` is not used).
- `MERGE_RESULTS` output: a single `Y1.regenie.gz` (the `.gz` extension belongs to this
  aggregate, not to step 2's output). After `zcat`, the file has a `##MASKS=<...>`
  metadata line followed by a header line with these exact whitespace-separated columns:
  `CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA`.
- MultiQC HTML report emitted under `outdir` (filename check only — reporting carve-out).
- `outdir/pipeline_info/*tracking*.json` present.

### Statistical-soundness assertion (the load-bearing one)

The point of this assertion is to fail loudly if a regression turns Regenie's output into
nonsense — for example, wrong sample-to-phenotype mapping, miswired covariates, dropped
masks, or cases/controls swapped. Per-row exact numerical comparison is impossible because
Regenie's LOCO step-1 corrections + SKAT-style burden aggregation cannot be equated to a
single-variant chi-square. The test therefore uses a **broad tolerance** sanity check: for
the genes Regenie thinks are most significant, a naive chi-square on the raw genotype
counts should at least *also* point at those genes.

**Helper script**: commit
`workflows/tests/helpers/recompute_naive_log10p.py` that, given the merged regenie output,
the input VCF, the per-gene setlist, the variant→mask annotations, and the cases/controls
lists, computes a naive `LOG10P` per gene.

**Naive LOG10P recipe** (deliberately simple; document this in the script's docstring):

1. Parse the merged regenie output. The ID column encodes the gene as the substring
   before the first `.`, e.g. `HOXC4.Mask_Mod.0.05` → gene `HOXC4`. Group rows by gene and
   keep `max(LOG10P)` and the mask category from the winning row.
2. Take the **top 5 genes** by `max(LOG10P)`. Skip rows where `LOG10P` is `NA`.
3. For each top-5 gene:
   a. Look up that gene's variants in the setlist file produced by
      `BCFTOOLS_ASSIGN_ANNOTATIONS`.
   b. Filter to variants whose annotation mask category matches the winning row's mask
      (e.g. `Mask_Mod`).
   c. For each remaining variant, count `0/0`, `0/1`, `1/1` genotypes separately in cases
      and controls (read from the input VCF via `pysam`).
   d. Build a 2×3 contingency table per variant and compute a chi-square
      `-log10(p)` (use `scipy.stats.chi2_contingency`).
   e. Take the maximum `-log10(p)` across the gene's variants → `naive_log10p[gene]`.

**Assertions** on the top-5 genes:

- For every gene in the top 5, `naive_log10p[gene]` is finite and ≥ 0.5. This is a *very*
  weak floor — it just says "for the genes Regenie ranks highest, a naive test also sees a
  signal stronger than coin-flip". Floor failure would mean Regenie is highlighting genes
  whose raw genotype distribution is statistically indistinguishable between cases and
  controls — that is almost certainly a wiring bug.
- The naive `-log10(p)` and Regenie's `LOG10P` agree on order of magnitude within a broad
  tolerance: `|naive_log10p[gene] - regenie_log10p[gene]| ≤ 2.0`. Two orders of magnitude
  is intentionally generous; this is a smoke check, not a numerical equivalence test.
- If either floor is violated for ≥3 of the 5 genes, fail the test with a message listing
  all five (gene, regenie LOG10P, naive LOG10P) triples so the regression is debuggable
  from CI logs alone. Single-gene violations are tolerated to avoid flakiness on a
  500-variant fixture where one top-5 entry may be borderline.

Tolerances above are first-pass; the implementer should tune them after seeing what a
healthy fixture run produces, and document the chosen numbers + reasoning in the helper
script's docstring.

- **Tag**: `"ci"`. This is the main parameter-optimization scenario; keep it fast.

## IT-7 — `workflows/rare-var-assoc` — full reporting path

- **Test file**: `workflows/tests/full_reporting.nf.test`
- **Config**: `skip_preparation=false`, `skip_reporting=false`, `use_dosage=false`.
- **Inputs**: `unprepared_rand_500.vcf.gz` + cases/controls (or phenotype file) committed
  to `assets/three_chr_unprepared/`.
- **Assertions**: all IT-6 assertions (including the naive-LOG10P statistical-soundness
  check — reuses the same `recompute_naive_log10p.py` helper) PLUS
  - `bcftools view -h <prepared_vcf>` contains `##INFO=<ID=CSQ,`
  - The expected EDA plot files are emitted under `outdir` (filename check only)
  - The GWAS HTML report files are emitted under `outdir` (filename check only)
- **Tag**: `"full"`. This is the slowest test; nightly only.

---

## Tasks

### T8 — Implement IT-1 and IT-1b (prepare subworkflow); prepare downstream fixture

Split into T8a and T8b on 2026-05-26. T8a and T8b are independent and may be parallelised;
the prepared-fixture step depends on T8b's IT-1 producing the prepared VCF.

#### T8a — IT-1b (skip-preparation + skip-reporting branch)

- **Test file to create**: `subworkflows/local/prepare/tests/skip_prep_skip_reporting.nf.test`.
- **Asset to commit**: `assets/three_chr_unprepared/cases.txt` and `controls.txt` — full
  3202-sample list split however convenient (e.g. first 2000 controls, last 1202 cases).
  Reused by T8b and IT-6/IT-7.
- **Verify**: `nf-test test --profile podman subworkflows/local/prepare/tests/skip_prep_skip_reporting.nf.test`.
- **Done-when**: IT-1b passes AND `cases.txt` / `controls.txt` are committed.

#### T8b — IT-1 (full preparation branch); record `prepared_500/` fixture

- **Prerequisite**: T8a (reuses the same `cases.txt`).
- **Test file to create**: `subworkflows/local/prepare/tests/main.nf.test` (IT-1).
- **After IT-1 passes**: copy the prepared VCF + run
  `plink2 --vcf <prepared> --make-pgen --out prepared_500` outside the test; commit the
  resulting `assets/three_chr_unprepared/prepared_500/prepared_500.{pgen,pvar,psam}` for
  reuse by IT-2..IT-7. Add a short `README.md` next to the files explaining their origin.
- **Verify**: `nf-test test --profile podman subworkflows/local/prepare/tests/main.nf.test`.
- **Done-when**: IT-1 passes AND `prepared_500.{pgen,pvar,psam}` are committed.

### T9 — Implement IT-2 (PCA)

- **Test file**: `subworkflows/local/pca/tests/main.nf.test` per IT-2.
- **Prerequisite**: T8 (prepared_500 assets committed).
- **Verify**: `nf-test test --profile podman subworkflows/local/pca/tests/main.nf.test`.
- **Done-when**: IT-2 passes.

### T10 — Implement IT-3 and IT-4 in parallel

Both reuse `prepared_500/` from T8.

- **IT-3 test file**: `subworkflows/local/filter_missing_per_pheno/tests/main.nf.test`;
  fixture: `subworkflows/local/filter_missing_per_pheno/tests/fixtures/two_pheno.tsv`.
- **IT-4 test file**: `subworkflows/local/f_coefficient_filtering/tests/main.nf.test`;
  fixture: spiked-pgen variant of `prepared_500` under
  `subworkflows/local/f_coefficient_filtering/tests/fixtures/`, OR an inline IID-injection
  in the test setup.
- **Verify**:
  `nf-test test --profile podman subworkflows/local/filter_missing_per_pheno/tests/main.nf.test subworkflows/local/f_coefficient_filtering/tests/main.nf.test`.
- **Done-when**: both pass.

### T11 — Implement IT-5 (reporting subworkflow)

- **Test file**: `subworkflows/local/reporting/tests/main.nf.test` per IT-5.
- **Fixtures**: a tiny `.regenie.gz`, annotations file, setlist, and phenotype file under
  `subworkflows/local/reporting/tests/fixtures/`.
- **Verify**: `nf-test test --profile podman subworkflows/local/reporting/tests/main.nf.test`.
- **Done-when**: passes (smoke-only per reporting carve-out).

### T12 — Implement IT-6 (full workflow, fast path)

- **Test file**: `workflows/tests/skip_prep_skip_reporting.nf.test` per IT-6.
- **Helper to author**: `workflows/tests/helpers/recompute_naive_log10p.py` implementing
  the naive LOG10P recipe documented in IT-6. Reused by IT-7. Include a docstring with the
  recipe and chosen tolerances.
- **Prerequisite**: T8.
- **Verify**: `nf-test test --profile podman workflows/tests/skip_prep_skip_reporting.nf.test`.
- **Done-when**: passes within the ≤5 min `ci` budget (including the
  naive-LOG10P recomputation step). If slower, profile the helper first; if Regenie
  itself is the bottleneck, downgrade to `tag "full"`.

### T13 — Implement IT-7 (full workflow, reporting path)

- **Test file**: `workflows/tests/full_reporting.nf.test` per IT-7.
- **Tag**: `"full"` (slow).
- **Verify**: `nf-test test --profile podman workflows/tests/full_reporting.nf.test`.
- **Done-when**: passes. Document expected wall-clock at the top of the test file.
