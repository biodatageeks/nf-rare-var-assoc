# IT-3 fixtures — filter_missing_per_pheno

## two_pheno.tsv

3-column phenotype file (FID, IID, Y1) for 3202 samples from the `prepared_500` fixture.
Engineered so that variant `12_52171632_T_C` has deliberately different per-group missing
rates, enabling a business-meaningful assertion that the intersection logic actually filters
the right variant.

### Group design

`12_52171632_T_C` has 1127 total missing genotypes across 3202 samples (35.2% overall).
Samples are assigned to groups as follows:

| Group | Y1 | Total | Missing for 12_52171632_T_C | F_MISS |
|---|---|---|---|---|
| Group 1 | 1 | 1202 | 180 | ~15% |
| Group 2 | 2 | 2000 | 947 | ~47% |

Construction:
- Group 1 = first 180 missing samples + first 1022 non-missing samples (by psam order)
- Group 2 = remaining 947 missing + remaining 1053 non-missing

### Threshold choice

With `--geno 0.30`:
- Group 1: 14.97% < 30% → `12_52171632_T_C` PASSES (stays in Group 1 snplist)
- Group 2: 47.35% > 30% → `12_52171632_T_C` FAILS (filtered from Group 2 snplist)
- `--extract-intersect` → `12_52171632_T_C` absent from output pvar

Control: `12_553708_G_A` has 0% missing in both groups, passes all filters, present in output.

### Regeneration

To regenerate this file from scratch (e.g. after the pgen fixture changes):
```bash
# get missing sample IDs for the target variant
podman run --rm -v $PWD:/wd:z quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    sh -c "bcftools view -r 12:52171632 /wd/assets/three_chr_unprepared/prepared_500/prepared_500.vcf.gz | \
           bcftools query -f '[%SAMPLE\t%GT\n]' | awk '\$2==\"./.\"' | cut -f1" \
    > missing_samples.txt
# get non-missing sample IDs
awk 'NR>1 {print $2}' assets/three_chr_unprepared/prepared_500/prepared_500.psam \
    | sort > all_sorted.txt
sort missing_samples.txt > missing_sorted.txt
comm -23 all_sorted.txt missing_sorted.txt > nonmissing_samples.txt
# build phenotype file
printf 'FID\tIID\tY1\n' > two_pheno.tsv
head -180  missing_samples.txt    | awk '{print $1"\t"$1"\t1"}' >> two_pheno.tsv
head -1022 nonmissing_samples.txt | awk '{print $1"\t"$1"\t1"}' >> two_pheno.tsv
tail -n +181  missing_samples.txt    | awk '{print $1"\t"$1"\t2"}' >> two_pheno.tsv
tail -n +1023 nonmissing_samples.txt | awk '{print $1"\t"$1"\t2"}' >> two_pheno.tsv
```
