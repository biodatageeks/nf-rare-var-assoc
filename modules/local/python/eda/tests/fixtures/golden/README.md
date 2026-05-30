# EDA Golden Stats

Golden CSV files that the `equivalence.nf.test` compares against.

## Fixture setup (one-time, before the test can pass)

Run the equivalence test once to generate `eda_stats/`:

```bash
nf-test test --profile podman,low_resources \
  modules/local/python/eda/tests/equivalence.nf.test
```

The test will fail at the "goldens populated" assertion — that is expected.
Find the emitted `eda_stats/` directory in `.nf-test/…/work/` and copy:

```bash
cp <work>/eda_stats/*.csv modules/local/python/eda/tests/fixtures/golden/
```

Then re-run — the test should now pass.  Commit all CSV files in this directory.

## File inventory

One CSV per plotted series/matrix.  Produced with:
- VCF: `assets/medium_data/prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz`
  (5000 variants, chr12+chr22+chrX, 3202 samples)
- Phenotype: `assets/three_chr_unprepared/pheno.tsv` (Y1=0: 2000 controls, Y1=1: 1202 cases)
- `use_dosage=true`

| File | Tolerance tier | Description |
|---|---|---|
| `variant_stats_{DP,GQ,DS}_p{1,50}_pheno{0,1}.csv` | loose | per-variant percentile series per phenotype |
| `variant_stats_{DP,GQ,DS}_pmean_pheno{0,1}.csv` | tight | per-variant mean series per phenotype |
| `sample_stats_{DP,GQ,DS}_pheno{0,1}.csv` | tight | per-sample mean per phenotype |
| `boxplot_{DP,GQ,DS}_pheno{0,1}.csv` | tight | per-sample mean per phenotype (boxplot source) |
| `missingness_variants_pheno{0,1}.csv` | tight | per-variant missingness rate per phenotype |
| `missingness_samples_pheno{0,1}.csv` | tight | per-sample missingness rate per phenotype |
| `het_pheno{0,1}.csv` | tight | per-sample heterozygosity rate per phenotype |
| `dp_diff.csv` | tight | per-variant abs DP diff (cases vs controls) |
| `allele_freq.csv` | tight | per-variant AF (NaN-filtered) |
| `variant_types.csv` | exact | [snp_count, indel_count] |
| `chrom_density.csv` | exact | CHROM,len (header row + one row per chromosome) |
| `heatmap_{stat2}_vs_{stat1}.csv` | exact | 2D histogram matrix (all samples) |
| `heatmap_{stat2}_vs_{stat1}_pheno{0,1}.csv` | exact | 2D histogram matrix per phenotype |

stat pairs: (GQ,DS), (DP,DS), (GQ,DP), (GT,DS)
