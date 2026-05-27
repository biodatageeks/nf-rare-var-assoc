# `f_coefficient_filtering` test fixture

Synthetic fixture for IT-4 (`subworkflows/local/f_coefficient_filtering`).

## Files

- `gen_fixture.py` — generates `spiked_fixture.vcf` (run to regenerate)
- `spiked_fixture.vcf` — source VCF (text; committed for auditability)
- `spiked_fixture.pgen` / `.pvar` / `.psam` — plink2 binary files used by the test

## Fixture design

| Group | Samples | Genotype | Expected F |
|---|---|---|---|
| Normal | NORM_01..NORM_60 (60) | 25 hom-ref + 25 hom-alt + 50 het per sample (random per-sample seed) | ~0 |
| Spike | SPIKE_HIGH1, SPIKE_HIGH2 (2) | 1/1 at all 100 sites | ~1.0 |

100 biallelic SNPs on chr1 (positions 1000..100000).

plink2 `--het` output (all 100 variants survive `--indep-pairwise 50 5 0.2`):
- Normal: F = -0.00713 (O_HOM=50, E_HOM=50.35, OBS_CT=100)
- Spike:  F = 1.0     (O_HOM=100, E_HOM=50.35, OBS_CT=100)

With `inbreeding_outliers_range_stds=3`:
- mean = 0.025, std = 0.178
- upper_bound = 0.025 + 3*0.178 = 0.559
- Spike F=1.0 > 0.559 → caught
- Normal F=-0.007 > lower_bound=-0.509 → not caught

## Why --double-id is needed

The VCF does not carry FID metadata. Without `--double-id`, plink2 creates a
psam with header `#IID SEX` (no FID). `calc_f_outliers` then writes the outliers
file in `IID IID` (duplicate) format, and plink2 `--remove` cannot match those
rows against the no-FID psam. `--double-id` sets FID=IID in the psam so that
`calc_f_outliers` uses the `FID IID` branch and `--remove` works correctly.

## Regeneration recipe

```bash
cd subworkflows/local/f_coefficient_filtering/tests/fixtures
python3 gen_fixture.py > spiked_fixture.vcf
podman run --rm -v $PWD:/wd/:z docker.io/psuszynski/plink:2.0-alpha.6.9 \
    plink2 --vcf /wd/spiked_fixture.vcf --double-id --out /wd/spiked_fixture --make-pgen
```
