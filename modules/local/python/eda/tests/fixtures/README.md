# Fixture: pheno.tsv

A hand-crafted 2-column phenotype file (IID, Y1) for the `python/eda` smoke test.

## Construction

Ten sample IDs taken from the first 20 samples listed in
`assets/three_chr_unprepared/unprepared_rand_500.vcf.gz`
(retrieved with `bcftools query -l`):

| IID | Y1 |
|---|---|
| HG00096 | 1 |
| HG00097 | 1 |
| HG00099 | 1 |
| HG00100 | 1 |
| HG00101 | 1 |
| HG00102 | 2 |
| HG00103 | 2 |
| HG00105 | 2 |
| HG00106 | 2 |
| HG00107 | 2 |

The remaining 3192 VCF samples receive `Y1 = null` via left-join inside the
EDA script. With exactly 2 distinct phenotype values (1 and 2), the script
takes the `len(phenotypes) <= 5` code path (per-phenotype histograms, DP
difference plots between the two groups).
