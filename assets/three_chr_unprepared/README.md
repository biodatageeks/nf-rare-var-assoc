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

## Phenotype and joined-sample fixtures

- `all_samples.txt` — `cases.txt` concatenated with `controls.txt` (3202 IIDs, cases first).
  Mirrors the `all.samples` file `JOIN_CASES_AND_CONTROLS` builds in production; used as the
  `ch_all_samples` input (the `--samples-file` of `BCFTOOLS_VIEW_AND_FILTER2`) in the IT-6
  fast-path workflow test.
- `pheno_binary.tsv` — three-column `FID IID Y1`, controls=0 (2000), cases=1 (1202). This is
  the shape Regenie consumes (and the shape `RSCRIPT_BUILD_PHENOTYPES` emits in production).
  Used by IT-6.
- `pheno.tsv` — two-column `IID Y1`, same 0=control / 1=case coding, used by the PCA
  subworkflow test (IT-2), where the phenotype is only a per-sample group label (plot
  colouring / covariate merge), not a Regenie trait. The 0/1 coding is what the pipeline's
  plink2 steps assume (`--1`: `0=control, 1=case`) and matches Regenie's `--bt`; it differs
  from `pheno_binary.tsv` only in lacking the `FID` column.

Regenerate `all_samples.txt` / `pheno_binary.tsv`:

```bash
cat assets/three_chr_unprepared/cases.txt assets/three_chr_unprepared/controls.txt \
    > assets/three_chr_unprepared/all_samples.txt
{ printf 'FID\tIID\tY1\n';
  awk '{print $1"\t"$1"\t0"}' assets/three_chr_unprepared/controls.txt;
  awk '{print $1"\t"$1"\t1"}' assets/three_chr_unprepared/cases.txt;
} > assets/three_chr_unprepared/pheno_binary.tsv
```

## VCFs

- `unprepared_rand_500.vcf.gz` (+ `.tbi`) — primary unprepared test VCF: 3 chromosomes
  (incl. X), 3202 samples, ~500 variants. Contains multiallelic sites and duplicate variant
  IDs; **not** VEP-annotated. Use for tests of preparation steps and full-pipeline runs.
- `unprepared_rand_{1k,2k,5k,10k}.vcf.gz` (+ `.tbi`) — same shape, larger variant counts
  for tests where Regenie complains about low variance at 500 variants.
- `unprepared_rand_500.vcf` — uncompressed copy of the 500-variant VCF.
