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

**Approach: polars-bio SQL `UNNEST` + standard `GROUP BY` reductions in DataFusion.** Never
materialize the matrix in Python. polars-bio reads the multisample VCF as a nested `genotypes`
struct of per-FORMAT-field lists; `scan_vcf` is itself a DataFusion streaming source, so the
engine streams the scan through the aggregate -- no Python-side batch/record iteration needed.
We `UNNEST` to a long table `(variant_idx, sample_idx, GT, DP, GQ, DS)` and express each plot's
reduction as a streaming hash-aggregate. Python only receives the small aggregated result
frames. Memory target O(variants + samples + bins). Consistent with `nf-prepare-vcf`'s
`calc_dosage_polarsbio` module. Reference usage:
[../nf-prepare-vcf/.../calc_dosage.py](../../../nf-prepare-vcf/modules/local/python/calc_dosage_polarsbio/assets/calc_dosage.py).

**Why standard GROUP BY is fine here (and why the FusedArrayTransform optimizer is irrelevant).**
That optimizer ([datafusion-bio-functions optimizer_rule.rs](https://github.com/biodatageeks/datafusion-bio-functions/tree/feature/optimize_unnest_groupby/datafusion/bio-function-vcftools),
authored for `calc_dosage.py`) fires *only* for the `UNNEST -> transform -> array_agg`-back
pattern that reconstructs per-variant lists -- that round-trip is what blows up memory, hence
the fused operator. Our queries are plain reductions (`AVG`/`COUNT`/`approx_percentile_cont`),
so they ride DataFusion's default streaming hash-aggregate: state is bounded by group count
(O(samples) or O(variants)), input rows flow through in batches and are never all buffered. We
do **not** need and will not trigger FusedArrayTransform.

**Implementation stance:** implement directly against this design; do not benchmark first. If
peak memory turns out bad in practice, the fallback is region-chunked reads (loop genomic
windows via predicate pushdown on the TBI index, reduce each window in polars/numpy, accumulate)
-- but only refactor to that if the straightforward GROUP BY approach proves insufficient.

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
- **T18b -- refactor to polars-bio streaming.** Rewrite data production; keep plotting and the
  stats-emission schema identical. The T18a test then validates new ~= old within tolerance.
  Lower the resource label; update container/env.

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

1. **Container.** Bump from `python_tools:1.0.4` to an image with polars-bio
   (`1.0.11` ships polars-bio 0.26.1). **Verify the image still has matplotlib + seaborn +
   polars + pysam** that EDA needs; update `environment.yml`.
2. **Approx percentiles.** `approx_percentile_cont` replaces exact `np.quantile` -- acceptable
   for EDA histograms; state it as an intentional behavior change.
3. **`start` -1 quirk.** Apply the same `+1` workaround as the reference module if POS is used.
4. **GT / het.** GT arrives as strings (`'0/1'`, `'0|0'`, `'0'`). Replace the
   integer-encode-then-`is_in([0,2,5,9])` het test with a direct two-differing-alleles check.
   The one numeric use of GT (`plot_stat_vs_stat` GT-vs-DS) needs a deliberate encoding choice.
5. **`datafusion.execution.target_partitions`** defaults to 1; bump for parallelism (trades
   memory for speed).
6. **2-D histogram bins** need a global `MAX(stat)` first (one cheap aggregate), as the numpy
   code does today.
7. **Resource label.** After the rewrite, drop `process_9` to a lower tier (target `process_2`
   / 16 GB) and confirm against the medium fixture.

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

**T18b**
- `EXPLORATORY_DATA_ANALYSIS` rewritten to polars-bio streaming aggregation; no full
  variants x samples matrix materialized in Python; plotting fns and stats-emission schema
  unchanged.
- T18a equivalence test passes against the goldens within tolerance.
- Container/env updated (polars-bio present; matplotlib/seaborn/polars/pysam still present).
- Resource label lowered (target `process_2`) and verified on the medium fixture.
