# Expected RSCRIPT_BUILDREPORTS outputs

These files are the byte-for-byte expected outputs that `build_reports.R` produces from the
sibling input fixtures. The test asserts equality against them. Update them only when a
deliberate change in `build_reports.R` semantics has been agreed.

## Why each value is what it is

Input fixtures (one mask, three variants, four samples, two cases + two controls):

| Variant         | S001 (case) | S002 (case) | S003 (ctrl) | S004 (ctrl) |
|-----------------|-------------|-------------|-------------|-------------|
| chr22:100:A:T   | het         | hom-ref     | hom-ref     | hom-ref     |
| chr22:200:G:C   | hom-alt     | het         | hom-ref     | hom-ref     |
| chr22:300:C:G   | hom-ref     | hom-ref     | het         | hom-ref     |

Allele frequency = (#het + 2·#hom-alt) / (2 · #genotyped). With 2 cases / 2 controls and no
missing calls, the denominator is 4 in both groups. So:

- `annotated_snps.csv`:
  - `chr22:100:A:T` → cases_ac = 1, cases_af = 0.25; controls_ac = 0.
  - `chr22:200:G:C` → cases_ac = 1 + 2·1 = 3, cases_af = 0.75; controls_ac = 0.
  - `chr22:300:C:G` → controls_ac = 1, controls_af = 0.25; cases_ac = 0.
- `res_log10p_1_annotated.csv`: the regenie fixture has a single row for set
  `GENE1.test.0.01` with LOG10P = 2.5 (> 1), and that mask covers all three variants — so we
  get three rows joining each variant's annotation and per-group AF to the regenie row.
- `annotated_snps_with_sample_ids.csv`: same three rows as `annotated_snps.csv` plus the
  sample-ID columns naming the carriers (e.g. `cases_htz=S002`, `cases_hmz=S001` for the
  pos-200 variant).
