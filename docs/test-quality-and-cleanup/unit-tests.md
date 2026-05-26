# Module Unit Tests (T7)

See [../test-quality-and-cleanup-plan.md](../test-quality-and-cleanup-plan.md) for goals,
conventions, the **quality bar for assertions**, and the **reporting carve-out**.

## §5a — High-priority: modules with complex Python/R logic

Each row specifies the test file to create, the fixture/input to use, and what to assert.
All tests live at `modules/local/<group>/<name>/tests/main.nf.test` and run with
`nf-test test --profile podman <path>`.

### `python/calc_f_outliers`

- **Test file**: `modules/local/python/calc_f_outliers/tests/main.nf.test`
- **Fixture**: a hand-rolled `het` file with ≥20 normal samples (F drawn from a tight
  distribution) and two spiked outliers (F well outside mean ± 2σ). Commit under
  `modules/local/python/calc_f_outliers/tests/fixtures/` with a short `README.md`
  documenting the spike values.
- **Assertions**:
  - `range_stds=2`: output has exactly 2 lines AND both lines list the spiked sample IIDs
    (not just any 2 IIDs — pin to the IIDs you spiked).
  - `range_stds=10`: output is an empty file.
  - **FID/IID column-selection branch** ([main.nf:44-47](../../modules/local/python/calc_f_outliers/main.nf#L44-L47)):
    include one fixture variant whose `het` header omits `FID`; assert the output still
    contains the spiked IID correctly.

### `python/vcf2aaf`

- **Test file**: `modules/local/python/vcf2aaf/tests/main.nf.test`
- **Fixture**: a tiny hand-rolled VCF (~10 variants) committed under
  `modules/local/python/vcf2aaf/tests/fixtures/`. Engineer it so:
  - some variants have only `AF` in INFO,
  - some have only `AF_nfe`,
  - some have neither,
  - at least one variant has a `chr`-prefixed chromosome (to test prefix stripping).
- **Assertions** (run with `tag_name="AF_nfe"`, `default_tag_name="AF"`):
  - Output row count equals input variant count.
  - `pos` column format is `<chr>_<pos>_<ref>_<alt>` with the `chr` prefix stripped — verify
    on the `chr`-prefixed input variant.
  - **Fallback logic** (the script's main custom behaviour, see
    [vcf2aaf.py:120-165](../../modules/local/python/vcf2aaf/assets/vcf2aaf.py#L120-L165)):
    pick three specific input variants and assert each produces the expected AF value:
    a variant with only `AF_nfe` → its `AF_nfe` value;
    a variant with only `AF` → its `AF` value;
    a variant with neither → `"0"`.

### `bcftools/assign_annotations`

- **Test file**: `modules/local/bcftools/assign_annotations/tests/main.nf.test`
- **Fixture**: small VCF (~20 variants spanning ≥2 genes) + a hand-crafted masks TSV with
  two distinct mask categories (e.g. `LoF`, `missense`). Commit under
  `modules/local/bcftools/assign_annotations/tests/fixtures/`. Engineer it so some variants
  land in mask A, some in mask B, and at least one variant is intentionally NOT in any mask.
- **Assertions**:
  - The variant intentionally outside any mask does NOT appear in `.annotations`.
  - For one variant known to be in mask A, the `.annotations` row has the expected
    `(variant_id, gene, mask_category)` triple.
  - `.setlist` groups variants by gene — for a gene with N input variants in the fixture,
    the setlist row for that gene lists exactly those N variant IDs.

### Reporting carve-out — smoke-only tests

The three modules below produce only human-readable plots/HTML. Per the reporting carve-out
in the main plan, do not parse PNG/SVG content or validate HTML structure; assert only that
the expected output files appear with the expected filename patterns.

| Module | Test file | Fixture | Assertions |
|---|---|---|---|
| `python/eda` | `modules/local/python/eda/tests/main.nf.test` | `assets/three_chr_unprepared/unprepared_rand_500.vcf.gz` + `.tbi` + a small 2-column phenotype TSV committed next to the test | (1) process succeeds; (2) at least one `plots/*.png` and one `plots/*.svg` emitted; (3) `versions.yml` contains `EXPLORATORY_DATA_ANALYSIS:`. Run once with `use_dosage=false`. |
| `python/draw_pc_plot` | `modules/local/python/draw_pc_plot/tests/main.nf.test` | Hand-rolled `.sscore` (3 PCs, 10 samples) + matching pheno TSV under `tests/fixtures/` | (1) process succeeds; (2) `plot_file` (PNG) and `plot_file_svg` both emitted with the expected filename pattern. |
| `python/generate_tracking_report` | `modules/local/python/generate_tracking_report/tests/main.nf.test` | Two hand-rolled `*.tracking.json` files under `tests/fixtures/`, where file B's `predecessor` field names file A's `process_name` | (1) process succeeds; (2) both `report_html_file` and `report_txt_file` are emitted; (3) the `.txt` report mentions both step names from the input JSONs (one-line grep). |

## §5b — Low-priority: wrapper-only modules

These modules contain little or no custom logic — they invoke an external tool with
configured arguments. Testing them has lower value (but still has value); if tests exist they should be marked
`tag "full"`.

- `plink2/makepgen`, `plink2/indep_pairwise`, `plink2/king_cutoff`, `plink2/pca`,
  `plink2/projection_score`, `plink2/het`, `plink2/import_dosage`, `plink2/export_other`,
  `plink19/makeset`, `bcftools/replace_sample_names`, `bcftools/vcf2psam`,
  `bcftools/view_and_filter2`, `cmds/check_x_chrom_present`,
  `cmds/extract_phenotypes_and_samples`, `cmds/rename`, `rscript/build_phenotypes`.

## §5c — Skip writing tests

- `python/filter_and_enhance_vcf_polarsbio` — likely removed in phase 2.
- `python/fix_zero_PL` — likely removed in phase 2.
- `combo/filter_and_enhance_vcf` — likely removed in phase 2.
- `python/view_and_filter2_polarsbio` — dead code (replaced by bcftools version, §1c).
- `vep/updatecache` — downloads external data only; tested manually.

---

## T7 — Write high-priority unit tests (§5a)

Split into T7a–T7f on 2026-05-26 — one sub-task per §5a row. The sub-tasks are
independent (separate modules, no shared fixtures); they may be done in any order.

For each sub-task, follow this recipe against the row referenced in §5a above:

1. Create `modules/local/<group>/<name>/tests/` if missing.
2. If the row needs a hand-rolled fixture, create
   `modules/local/<group>/<name>/tests/fixtures/` with the fixture(s) plus a `README.md`
   explaining how each was constructed.
3. Write `main.nf.test` with `tag "ci"`, the listed input, and the listed assertions.
4. Verify: `nf-test test --profile podman modules/local/<group>/<name>/tests/main.nf.test`.

| Sub-task | §5a row | Notes |
|---|---|---|
| T7a | `python/calc_f_outliers` | Three test cases (range_stds=2, range_stds=10, FID-absent branch). Business-meaningful assertions on outlier IIDs. |
| T7b | `python/vcf2aaf` | Hand-rolled VCF spans `AF`-only, `AF_nfe`-only, neither, and a `chr`-prefixed variant. Assert AF fallback for three named variants. |
| T7c | `bcftools/assign_annotations` | Hand-rolled VCF (≥2 genes) + masks TSV (2 categories). Assert one variant absent, one triple, one setlist row. |
| T7d | `python/eda` | Reporting carve-out — filename-only smoke. Reuses `assets/three_chr_unprepared/unprepared_rand_500.vcf.gz`; pheno TSV hand-rolled next to the test. |
| T7e | `python/draw_pc_plot` | Reporting carve-out — filename-only smoke. Hand-rolled `.sscore` (3 PCs, 10 samples) + pheno TSV. |
| T7f | `python/generate_tracking_report` | Reporting carve-out — filename-only smoke. Two hand-rolled `*.tracking.json` chained via `predecessor`. |

**Done-when (T7 overall)**: all six §5a tests pass; total CI runtime added ≤ ~5 min.
**Done-when (per sub-task)**: the named module's `main.nf.test` passes under `nf-test
test --profile podman <path>`; fixtures (if any) committed with a `README.md`; index +
this sub-doc updated with status.
