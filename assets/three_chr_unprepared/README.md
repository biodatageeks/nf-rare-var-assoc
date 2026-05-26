# `assets/three_chr_unprepared/`

Canonical 1000 Genomes test fixtures used by the pipeline's test profiles and by integration
tests (IT-1..IT-7 in `docs/test-quality-and-cleanup/integration-tests.md`).

## Sample lists

- `samples.psam` — full 3202-sample plink2 psam (FID/IID/SEX) for the 1000G high-coverage
  release.
- `cases.txt` — 1202 IIDs designated as cases. Taken as the **last 1202 entries** from
  `bcftools query -l unprepared_rand_500.vcf.gz`.
- `controls.txt` — 2000 IIDs designated as controls. Taken as the **first 2000 entries** from
  the same `bcftools query -l` ordering.

The cases/controls split is arbitrary (sample order, not biology). It exists only to provide
a deterministic ~37%/63% split with 3202 total samples for tests that need a binary phenotype.

Regenerate:

```bash
podman run --rm -v $PWD/:/wd/:z quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    sh -c "bcftools query -l /wd/assets/three_chr_unprepared/unprepared_rand_500.vcf.gz" \
    > /tmp/all_samples.txt
head -2000 /tmp/all_samples.txt > assets/three_chr_unprepared/controls.txt
tail -1202 /tmp/all_samples.txt > assets/three_chr_unprepared/cases.txt
```

## VCFs

- `unprepared_rand_500.vcf.gz` (+ `.tbi`) — primary unprepared test VCF: 3 chromosomes
  (incl. X), 3202 samples, ~500 variants. Contains multiallelic sites and duplicate variant
  IDs; **not** VEP-annotated. Use for tests of preparation steps and full-pipeline runs.
- `unprepared_rand_{1k,2k,5k,10k}.vcf.gz` (+ `.tbi`) — same shape, larger variant counts
  for tests where Regenie complains about low variance at 500 variants.
- `unprepared_rand_500.vcf` — uncompressed copy of the 500-variant VCF.
