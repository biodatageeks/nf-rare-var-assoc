# `vcf2aaf` test fixtures

One hand-rolled VCF file for `PYTHON_VCFTOAAF` unit tests.

## `test_vcf2aaf.vcf`

A 10-variant VCF with no sample/genotype columns (polars_bio does not need them).
Both `AF` and `AF_nfe` are declared in the header so that per-variant presence
controls which fallback branch the script takes.

| Variant | CHROM | POS | REF | ALT | INFO | Purpose |
|---|---|---|---|---|---|---|
| var01 | 1 | 1000 | A | T | AF_nfe=0.05;AF=0.08 | Both tags present |
| var02 | 1 | 2000 | G | C | AF_nfe=0.15;AF=0.12 | Both tags present |
| var03 | 2 | 1000 | T | A | AF_nfe=0.25 | **AF_nfe-only** — primary-tag fallback target |
| var04 | 3 | 1000 | C | G | AF=0.10 | **AF-only** — default-tag fallback target |
| var05 | 4 | 1000 | G | A | . | **Neither tag** — "0" fallback target |
| var06 | chr5 | 1000 | A | C | AF_nfe=0.20 | **chr-prefix stripping** — CHROM prefixed with "chr" |
| var07 | 1 | 3000 | T | G | AF=0.03 | Padding |
| var08 | 1 | 4000 | C | T | AF_nfe=0.50 | Padding |
| var09 | 1 | 5000 | A | G | AF=0.08 | Padding |
| var10 | 1 | 6000 | G | T | . | Padding |

## Assertions enabled by this fixture

1. **Row count** — 10 input variants produce exactly 10 output rows.
2. **chr-prefix stripping** — var06's `CHROM=chr5` must appear as `5_...` in the
   `pos` column, never as `chr5_...`
   ([vcf2aaf.py:107](../../assets/vcf2aaf.py#L107)).
3. **AF fallback logic** ([vcf2aaf.py:120-165](../../assets/vcf2aaf.py#L120-L165)):
   - var03 (chr2, AF_nfe-only): output AF = "0.25"
   - var04 (chr3, AF-only): output AF = "0.1"
   - var05 (chr4, neither): output AF = "0"

Variants are placed on distinct chromosomes (2, 3, 4, chr5) so the test can
identify each target row by its chromosome prefix without knowing the exact
coordinate that polars_bio assigns.
