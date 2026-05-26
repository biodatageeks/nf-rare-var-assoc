# assets/medium_data

Contains medium-sized test fixtures derived from 1000 Genomes Project (1kGP) data.

## Prepared VCFs (unique IDs, biallelic, VEP-annotated)

These files have been through the full preparation pipeline (bcftools norm → bcftools
annotate → VEP) and are suitable for tests of downstream modules (plink2, regenie, etc.)
where the module under test should not need to handle raw-VCF quirks such as duplicate
variant IDs or multiallelic sites.

### prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz

- Chromosomes: chr12, chr22, chrX
- Variants: ~5000 (2000 CSQ-selected + 3000 random)
- Samples: 3202 (1kGP)
- VEP-annotated; unique variant IDs; no multiallelic sites

### prepared_chr12_100.vcf.gz

- First 100 variant records from `prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz`
  (all chr12)
- Samples: 3202 (1kGP)
- Intended for fast module-level tests where the full 5000-variant file is unnecessarily
  large (e.g. plink2/write_snplist, plink2/makebed unit tests)

Created with:
```bash
zcat prepared_chr_12_22_X_csq_filtered_2k_rand_3k.vcf.gz \
    | head -$((3497 + 100)) \
    > prepared_chr12_100.vcf
# 3497 = number of header lines in the source VCF
podman run --rm -v $PWD:/wd/:z quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    sh -c "bcftools view -Oz -o /wd/prepared_chr12_100.vcf.gz /wd/prepared_chr12_100.vcf \
           && bcftools index -t /wd/prepared_chr12_100.vcf.gz"
rm prepared_chr12_100.vcf
```

## Unprepared VCFs

Raw 1kGP VCFs: may contain duplicate variant IDs, multiallelic sites, no VEP annotation.
Use for tests that exercise the preparation steps or integration tests that run the full
pipeline.

### unprepared_rand_{500,1k,2k,5k,10k}.vcf.gz

Random subsets of the 1kGP dataset. The 500-variant file additionally exists as a plain
`.vcf` for tools that do not accept bgzipped input.
