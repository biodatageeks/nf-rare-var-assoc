# `assets/three_chr_unprepared/prepared_500/`

Pre-prepared plink2 pgen/pvar/psam derived from the PREPARE subworkflow's full-prep
branch run on `assets/three_chr_unprepared/unprepared_rand_500.vcf.gz`. Used as the
input fixture for downstream-subworkflow integration tests (IT-2 PCA, IT-3
filter_missing_per_pheno, IT-4 f_coefficient_filtering) and the IT-6 fast-path
workflow test. Avoids re-running NORM/ANNOTATE/FILTER_AND_ENHANCE_VCF in every
downstream test.

## Files

- `prepared_500.{pgen,pvar,psam}` — canonical fixture. 3202 samples, 417 variants
  on `12`, `22`, `X`, `PAR1` (after `--split-par hg38`). Variant IDs are
  `<CHROM>_<POS>_<REF>_<ALT>` (assigned by BCFTOOLS_ANNOTATE).
- `prepared_500.vcf.gz` — the prepared VCF that produced the pgen above. Kept for
  future tests that need a VCF-shaped prepared input rather than pgen.

## Variant-count derivation

| Stage | Variant count |
|---|---|
| Input (`unprepared_rand_500.vcf.gz`) | 500 |
| Post BCFTOOLS_NORM `-m -any` (split 35 multi-ALT sites into 82 records, 0 exact dups) | 547 |
| Post BCFTOOLS_ANNOTATE (no count change; only renames `chr*` and assigns IDs) | 547 |
| Post FILTER_AND_ENHANCE_VCF (default thresholds; drops 130 on QUAL/AVG_GQ/AVG_DP) | 417 |

The 417 figure is asserted exactly in
`subworkflows/local/prepare/tests/main.nf.test` (IT-1).

## Regenerate

```bash
# 1. Run IT-1 to produce the prepared VCF
nf-test test --profile podman subworkflows/local/prepare/tests/main.nf.test
# 2. Copy the FILTER_AND_ENHANCE_VCF output out of the nf-test work directory
cp .nf-test/tests/<hash>/work/<XX>/<YYY...>/test_filterhance.vcf.gz \
    assets/three_chr_unprepared/prepared_500/prepared_500.vcf.gz

# 3. Build an input .psam in VCF sample order, with SEX looked up from
#    assets/three_chr_unprepared/samples.psam (the prepared VCF reorders samples
#    so the BCFTOOLS_VIEW_1 --samples-file order is preserved: cases first, then
#    controls; see subworkflows/local/prepare/tests/fixtures/all_samples.txt).
cd assets/three_chr_unprepared/prepared_500
podman run --rm -v $PWD/..:/wd:z quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    sh -c "bcftools query -l /wd/prepared_500/prepared_500.vcf.gz" > /tmp/vcf_order.txt
python3 -c "
sex = {}
for line in open('../samples.psam'):
    if line.startswith('#'): continue
    fid, iid, s = line.strip().split('\t'); sex[iid] = s
with open('source.psam', 'w') as o:
    o.write('#FID\tIID\tSEX\n')
    for line in open('/tmp/vcf_order.txt'):
        iid = line.strip(); o.write(f'{iid}\t{iid}\t{sex[iid]}\n')
"

# 4. Convert to pgen with PAR splitting (chrX requires sex info; without --split-par
#    plink2 cannot pseudo-autosomal-correct chrX). One sample (HG02300) has
#    ambiguous sex and ends up with SEX=NA in the output psam; this is expected.
podman run --rm -v $PWD:/wd:z -w /wd docker.io/psuszynski/plink:2.0-alpha.6.9 \
    sh -c "plink2 --vcf prepared_500.vcf.gz --psam source.psam --make-pgen \
        --split-par hg38 --out prepared_500"

rm source.psam prepared_500-temporary.* prepared_500.log
```
