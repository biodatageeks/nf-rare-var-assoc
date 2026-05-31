# EDA Memory Reduction (T18)

See [../test-quality-and-cleanup-plan.md](../test-quality-and-cleanup-plan.md) for conventions
and the **reporting carve-out** (EDA is a reporting surface: smoke-level tests, no plot-content
assertions).

## Problem

[`EXPLORATORY_DATA_ANALYSIS`](../../modules/local/python/eda/main.nf) is the most
memory-hungry process in the pipeline: `process_9` = 60 GB, retrying to 120 GB
([conf/base.config:61-64](../../conf/base.config#L61-L64)).

Root cause is [`load_vcf`](../../modules/local/python/eda/main.nf#L57-L111):
it builds a dict of Python lists with one column **per sample per stat** (`4 x n_samples + 4`
columns; ~12,800 for 3202 samples), so every cell is a boxed Python object (~36 B) during
accumulation, then `pl.DataFrame(data, schema=dtypes)` holds the dict and the Polars copy
simultaneously. The full variants x samples matrix is materialized even though every plot
consumes only an aggregate. Memory is O(variants x samples) -- ~230 GB as Python lists / ~26 GB
as Float32 for 500k variants x 3202 samples. Compact dtypes alone cannot fix WGS scale; the
matrix has to go.

## Decision (LOCKED)

**Approach: two-pass, chunked pysam scan (v5).** Do not materialize the full
variants x samples matrix in memory. Use pysam to iterate VCF records in chunks
of N variants (or per-contig windows), compute per-variant statistics within
each chunk, and accumulate per-sample running sums/counters. A second pass is
used only for heatmaps to preserve exact `np.histogram2d` binning (global or
per-phenotype max values must be known to set bin edges). Memory target
O(variants + samples + bins), plus one bounded chunk.

**Why two passes are needed:** per-variant percentiles are local to a variant
and can be computed in the first pass, but the heatmap bins depend on global
maxima. Exact binning requires those maxima before computing the histogram.

**Implementation stance:** keep logic aligned with the previous numpy
computations (exact `np.quantile`, exact heatmaps). Avoid polars-bio and SQL
aggregation for v5.

## Plot -> SQL aggregation map

Every plot is expressible as a streaming group-by; Python only receives small aggregated frames.

| Plot (current fn) | Aggregation | Output groups |
|---|---|---|
| `plot_variant_stats` mean | `GROUP BY variant_idx -> AVG(stat)` | O(variants) |
| `plot_variant_stats` percentile | `GROUP BY variant_idx -> approx_percentile_cont(stat, p)` | O(variants) |
| per-phenotype variant stats | join sample->pheno, filter, same group-by | O(variants) |
| `plot_sample_stats`, `plot_boxplots` | `GROUP BY sample_idx -> AVG(stat)` | O(samples) |
| `plot_missingness` | `AVG(gt IS NULL)` grouped by variant_idx and by sample_idx | O(variants)+O(samples) |
| `plot_dp_differences` | `GROUP BY variant_idx -> AVG(DP) FILTER (case/control)` | O(variants) |
| `plot_heterozygosity` | `GROUP BY sample_idx -> AVG(is_het)` | O(samples) |
| `plot_allele_frequency`, `plot_variant_types`, `plot_chrom_density` | variant-level columns / counts, no unnest | O(variants) |
| `plot_stat_vs_stat` | `GROUP BY floor(s1/bw), floor(s2/bw) -> COUNT(*)` (+ per-pheno) | O(bins) |

## Subtasks

Split so the refactor is validated against the old behavior:

- **T18a -- golden baseline + equivalence harness (no behavior change).** Instrument the
  *current* code to emit a structured stats artifact, capture goldens, and write the
  comparison test. Test passes exactly on current code (code vs itself).
  **Status: ✅ Done 2026-05-30.**
- **T18b -- refactor to v5 (pysam two-pass chunking).** Rewrite data production; keep plotting
  and the stats-emission schema identical. The T18a test then validates new ~= old within
  tolerance. Lower the resource label once memory improves. **Status: DONE 2026-05-31. v5
  (two-pass pysam, O(variants + samples + bins + chunk)) is the active script; passes the
  T18a equivalence test and the T7d smoke test against the committed goldens. Resource label
  is `process_2`.**

## v5 outcome (2026-05-31)

v5 ([assets/eda_v5.py](../../modules/local/python/eda/assets/eda_v5.py)) implements Option 2
(two-pass chunked pysam) and is wired into both tests. Verified facts:

- **Equivalence + smoke both pass** with the script path pointing at v5; spot-checked that they
  also pass against v1 (the de-escaped baseline) so the goldens are reproducible from either.
- **Pass 1 streams** per-variant series (O(variants)) and per-sample running sums (O(samples));
  no variants x samples matrix is ever resident. **Pass 2** (only when `use_dosage=true`) refills
  `np.histogram2d` bins from bounded chunks (`--chunk-size`, default 500) using the per-group
  maxima recorded in pass 1, so percentiles and heatmaps stay exact.
- **`mean_horizontal` vs `nanmean`:** missing FORMAT values arrive as pysam `None` -> polars
  `null`, which `mean_horizontal` skips -- equal to v5's `np.nanmean`. Percentiles propagate NaN
  in both. v5 reproduces v1's mean-kept / percentile-dropped asymmetry per variant.
- **Plots 1_/1b_ ("Mean {stat} Across Variants") are emitted unconditionally** (overall mean over
  all samples), matching v1. Earlier v5 dropped them in the <=5 branch; fixed via a dedicated
  `overall_means` accumulator (avoids computing unused `'all'`-group percentiles in the <=5 path).
- **Smoke test now asserts the full expected file set** (28 PNG + 28 SVG stems for the
  use_dosage=false binary scenario), not just `count >= 1`.

### File state
- `assets/eda_v5.py` -- **active** (two-pass pysam). Used by both tests.
- `assets/eda_v1.py` -- de-escaped + argparse'd baseline kept **only** for v1-vs-v5 comparison
  (same CLI as v5). Delete once comparison is no longer needed.
- `assets/eda_v2.py` (polars-bio, slow), `assets/eda_v4.py` (chunked load, full matrix) --
  **superseded by v5; safe to delete.**
- `main.nf` invokes `python3 ${python_script} --vcf .. --phenotype .. --use-dosage ..
  --process-name ..`; input tuple's 4th element is the script path. Caller
  [workflows/rare-var-assoc.nf](../../workflows/rare-var-assoc.nf) must point `eda_script_ch` at
  `eda_v5.py`.

## Historical: v1 outcome + v5 redesign (superseded by the v5 outcome above)

### What happened with v1
The data layer was rewritten from inline `main.nf` Python into an external script
[modules/local/python/eda/assets/eda_v2.py](../../modules/local/python/eda/assets/eda_v2.py)
(decision: external for testability / SQL escaping; matches the `calc_dosage.py` reference).
It uses polars-bio `register_vcf` + per-plot SQL `GROUP BY` reductions. **Its correctness is
UNCONFIRMED** -- the equivalence test was never run to completion against it (individual SQL
queries were spot-checked against goldens during development, but the full per-series diff never
passed). It was abandoned purely on speed: the old code did **one** VCF scan; v1 issues a separate
`pb.sql(...).collect()` for
nearly every plot, and *each `collect()` re-decompresses and re-parses the whole VCF* (DataFusion
does not cache the registered VCF across queries). The binary/`use_dosage=true` path runs **~26
scans**: 3 variant-stats (DP/GQ/DS) + 6 sample-means (recomputed for boxplots) + missingness x2 +
het + dp_diff + variant_level + **12 heatmap cell-pair collects**. Old equivalence run ~112 s;
v1 is substantially slower.

**Hard rule for v3 (user, LOCKED): one VCF scan; two only if it demonstrably helps.** Compute
every aggregate needed by every plot during that pass, then plot from the small results.

### Verified v4/v5 facts (reuse in v5; do not re-derive)
- `eda_v4.py` loads VCF in chunks and extends a Polars DataFrame; the full matrix still exists
  after load, but peak memory during construction is lower and dtypes are narrower.
- Percentiles and heatmaps remain **exact** in v4 because it still uses `np.quantile` and
  `np.histogram2d` on numpy arrays.
- The equivalence test is currently wired to `eda_v4.py` to avoid the approximate percentile
  drift seen in the polars-bio v1 path.
- The `>5`-phenotype branch is not asserted by the equivalence test, but it must remain
  implemented and functional.

### Verified polars-bio / SQL facts (archived; reuse if we revisit that path)
- Container `python_tools:1.0.11` has polars-bio 0.26.0, polars 1.39.3, pysam 0.23.3,
  matplotlib 3.10.8, seaborn 0.13.2.
- `register_vcf(path, name=, format_fields=['GT','DP','GQ','DS'], info_fields=['AF'])`. Schema:
  `chrom, start, end, id, ref, alt, qual, filter, AF (List), genotypes (Struct)`.
- Access struct fields as `genotypes."GT"` etc.; INFO `AF` is a List -> `"AF"[1]` (1-based);
  **`AF` must be double-quoted** in SQL (bare `af` fails).
- `start` is off by one (polars-bio quirk) -> POS = `start + 1`.
- Unnest pattern: `CROSS JOIN generate_series(0, N-1) AS s` (yields column `value`) with
  `arr[CAST(s.value AS BIGINT)+1]`. Gives deterministic `sample_idx` in VCF header order.
- GT arrives as strings: `'0/0'`, `'0|1'`, `'0'`, `'1'`, `'./.'`, `'.'` (haploid present on chrX).
  Encode f(j,k)=k(k+1)/2+j via `split_part`; missing (any `.`) -> NULL; het = two alleles differ
  (haploid -> 0). See `_gt_encode_sql` / `_gt_is_het_sql` in v1 -- carry them over.
- `approx_percentile_cont(col, p)` and `AGG(...) FILTER (WHERE ...)` both work. Per-variant
  percentile NULL-guard (NULL if any group value is NULL) reproduces golden row counts exactly
  (DP p1 pheno0 = 4860, mean = 4955). Old code used exact `np.quantile`; approx lives in the
  **loose** tolerance tier (see tiers table).
- Heatmap binning exactness is fragile: replicating `np.histogram2d` edges with SQL `floor(v/bw)`
  mismatches at exact internal edges (e.g. DS=1.0 -> `1.0/0.01`=99.9999 floors to 99, np puts it
  in bin 100). **numpy binning is exact by construction.** Also note the OLD overall
  heatmap hardcodes `stat2_bins = linspace(0,2,201)` (so for GQ-vs-DP, DP>2 cells are dropped) --
  a quirk that MUST be preserved to match goldens; per-pheno uses `linspace(0, nanmax, 201)`.
- Goldens are **binary-pheno only**; the equivalence test runs only that path with
  `use_dosage=true`.

### Memory reality (drives the design)
Both a **wide** matrix (old) and a **full long frame** (one row per cell) are O(variants x samples):
~250 MB-1 GB at the 5000-variant fixture, but ~25-40 GB at WGS scale (500k variants x 3202). So
"scan once, collect all cells, do everything in pandas/polars" is fast and simple **but fails at
WGS** -- it is NOT an acceptable v3 on its own. v3 memory must be O(variants + samples + bins),
i.e. only aggregated results (and at most one bounded chunk of cells) ever resident.

### v5 design options

**Option 1 -- single pass with fixed bins.**
Use hard-coded heatmap edges and compute all stats in one pass. Fast, but binning is only
approximate to the old code, so exact-tier heatmap comparisons may fail.

**Option 2 -- two-pass chunked pysam (selected).**
Pass 1: iterate records in chunks; compute per-variant arrays (means/percentiles), accumulate
per-sample running sums/counts, missingness, het; track global and per-phenotype maxima for
heatmap edges. Pass 2 (only when `use_dosage=true`): re-scan and fill `np.histogram2d` bins
exactly with the correct edges, accumulating bin counts per chunk. This preserves exact
percentiles and heatmaps while bounding memory by chunk size.

**Option 3 -- full matrix in memory. REJECTED** (O(cells) = 25-40 GB at WGS).

### Prior v3 design options (polars-bio path, archived)

**Option 1 -- consolidated GROUP BY scans (closest to original LOCKED approach; ~3 scans).**
Merge v1's many queries by grouping key:
- Scan A: one `GROUP BY variant_idx` over the unnested stream computing **all** per-variant series
  (DP/GQ/DS mean+p1+p50 per pheno + 'all', missingness, dp_diff) plus per-variant `MAX(stat)` and
  the per-variant scalars (AF, variant_type via `MAX`/`ANY_VALUE`; chrom). Global/per-pheno maxes =
  max over this 5000-row result (no extra scan).
- Scan B: one `GROUP BY sample_idx` computing all per-sample series (DP/GQ/DS means, missingness,
  het). Feeds both sample-stat and boxplot plots (compute once, not 6x).
- Scan C (`use_dosage` only): all 12 heatmaps in one scan by expanding each cell into the relevant
  (pair_id, group_label, bin_x, bin_y) rows (array-of-structs `UNNEST`) then `GROUP BY` -> COUNT;
  edges precomputed from Scan A maxes. State O(bins); cells never materialized.
- Pros: pure polars-bio streaming, lowest memory, scales to WGS, smallest conceptual delta from v1.
- Cons: 3 scans (3x decompress); heatmap row-expansion SQL is intricate; **SQL binning vs
  `np.histogram2d` exactness risk** on the exact tier (may force relaxing the heatmap tolerance).

**Option 2 -- single (or double) chunked pass in Python (pre-v5 recommendation).**
Loop over genomic chunks (per-chromosome, or fixed N-variant windows via region predicate /
`pysam.fetch`). Per chunk build a *small* typed frame (chunk_variants x samples; e.g. 500 variants
-> ~26 MB) and reuse the OLD per-chunk reductions:
- per-variant series: complete within the chunk (a variant never spans chunks) -> append; exact
  `np.quantile` percentiles (back to tight tier).
- per-sample series: accumulate running sum/count (means), null counts (missingness), het counts.
- heatmaps: accumulate `np.histogram2d` counts per chunk and sum (additive) -> **exact**.
Heatmap edges need global/per-pheno maxes first; track them as running maxes in the pass and do
heatmap binning in a **second** chunked pass (so `use_dosage=false` => 1 pass, `true` => 2 passes).
- Pros: effectively one VCF traversal; memory bounded by chunk; **heatmaps and percentiles exact**
  (can tighten test tiers); large code reuse from `eda__old.py`; scales to WGS.
- Cons: manual chunk orchestration; two passes when `use_dosage`.

**Option 3 -- one scan, collect full long frame, aggregate in-memory. REJECTED** (O(cells) =
25-40 GB at WGS). Fine only as a throwaway correctness oracle at fixture scale.

## Equivalence test design (data layer, not pixels)

PNG/SVG comparison is the wrong layer: matplotlib rendering is not reproducible enough and a
failing pixel diff is non-diagnostic. Instead compare the **numeric series/matrices behind each
plot**. This is a complete equivalence proof for the plots **iff the refactor keeps the
plotting functions unchanged and only replaces data production** -- T18b must respect that.

**Stats artifact.** Add a permanent **optional** `stats` emit channel
(`path("eda_stats/*"), emit: stats, optional: true`); not wired into the workflow, not
published by default. Right before each `sns.*`/`plt.*` call, write the exact array/matrix that
feeds it, keyed by plot id, as a **human-readable CSV** (one file per plotted series/matrix).
Every series is O(variants), O(samples), or O(bins) -- never the full matrix. This also upgrades
EDA from smoke-only to genuinely tested at the data layer.

**Implementation constraints for T18a (LOCKED -- decided 2026-05-30).**
- **CSV, not parquet.** Goldens must be human-readable/diffable plain text. The equivalence
  comparison then runs in the nf-test Groovy `then {}` block (readLines + parse floats +
  tolerance check) -- no extra container/process needed.
- **Minimal, surgical edits to the *current* inline Python in `main.nf`.** Do **NOT** refactor
  the inline script into an external `assets/eda.py`/template, and do not restructure existing
  functions. Rationale: (1) we are editing untested code *before* the test exists, so every edit
  is risky -- keep the diff tiny; (2) T18b replaces this whole script anyway, so refactoring it
  now is wasted work. Add only: a small `emit_*` helper + one emit call right before each plot,
  and the `stats` output channel.
- Generate goldens by running the instrumented module via nf-test and copying the emitted
  `eda_stats/` CSVs into `tests/fixtures/golden_*/`.

**Fixtures.** Use `assets/medium_data/prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz`
(5000 variants, chr12+22+X, 3202 samples; has GT/DP/GQ/DS in FORMAT and AF in INFO). Cover all
branches across two runs sharing this VCF: the binary case/control path (existing 2-group pheno)
**and** the `>5`-phenotype path (a >5-group pheno file). Both runs use `use_dosage=true` to
exercise DS and the `plot_stat_vs_stat` heatmaps.

**Tolerance tiers** (the approx-percentile change lives only in the loose tier):

| Tier | Series | Comparison |
|---|---|---|
| Exact | chrom-density counts, variant-type counts, 2-D histogram bin matrices, total variant count | integer-equal |
| Tight `allclose` | per-variant & per-sample means (DP/GQ/DS), missingness rates, het rates, AF, DP case-vs-control diffs | `rtol~1e-4, atol~1e-6` (float order / F32-vs-F64 only) |
| Loose | per-variant percentile series from `approx_percentile_cont` only | bound max deviation AND >=99% within tight rtol (or small total-variation bound on the 50-bin histogram of the series) |

The tight tier also guards NaN/NULL semantics: current code filters `~np.isnan`; SQL `AVG`
skips `NULL` -- a divergence shows up here.

## Watch list

1. **Chunk size.** Pick a default chunk size that keeps peak memory under the target tier
  (start with 500-1000 variants). Make it easy to tune if needed.
2. **Heatmap edges.** Pass 1 must record max values per stat (overall + per-pheno) so pass 2
  uses the exact edges from the old code.
3. **GT / het.** Preserve the old integer encoding and het logic; haploid calls on chrX must
  keep working.
4. **Resource label.** After v5, confirm the lower tier (`process_2`) on the medium fixture.

## Done when

**T18a**
- Current EDA code emits the `stats` artifact (optional output, **CSV**) for every plotted
  series/matrix, via minimal surgical edits to the existing inline Python (no refactor).
- Golden stats files (CSV) committed under `modules/local/python/eda/tests/fixtures/golden_*/`
  with a `README.md`.
- New nf-test compares emitted stats to goldens with the tolerance tiers (comparison in the
  Groovy `then {}` block); passes exactly on current code (`nf-test test --profile podman` +
  `low_resources`). Existing T7d smoke test
  ([modules/local/python/eda/tests/main.nf.test](../../modules/local/python/eda/tests/main.nf.test))
  still passes.

**T18b -- DONE 2026-05-31**
- [x] `EXPLORATORY_DATA_ANALYSIS` rewritten to two-pass chunked pysam scans (v5); no full
  variants x samples matrix materialized in Python; plotting fns and stats-emission schema
  unchanged (plots 1_/1b_ restored to emit unconditionally, matching v1).
- [x] T18a equivalence test passes against the goldens within tolerance.
- [x] Smoke test (T7d) asserts the full expected output-file set, not just `count >= 1`.
- [x] Resource label lowered to `process_2` and exercised on the medium fixture via the
  equivalence test (5000 variants, 3202 samples) under `low_resources`.
