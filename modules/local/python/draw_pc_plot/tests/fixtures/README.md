# draw_pc_plot test fixtures

## samples.sscore

Tab-separated plink2 `.sscore` output (3 PCs, 10 samples).

Columns: `#FID`, `IID`, `ALLELE_CT`, `NAMED_ALLELE_DOSAGE_SUM`, `PC1_AVG`, `PC2_AVG`, `PC3_AVG`.

Samples `HG00096`–`HG00101` form a tight cluster with positive PC1 (approximately +0.012) and
negative PC2 (approximately -0.004), representing one ancestry group.

Samples `HG00102`–`HG00107` form a second cluster with negative PC1 (approximately -0.031) and
positive PC2 (approximately +0.021), representing a second ancestry group.

The two-cluster structure ensures the scatter plot in `draw_pc_plot` exercises the `hue`
colour-coding path correctly.

## pheno.tsv

Two-column TSV (`IID`, `Y1`) matching the 10 samples in `samples.sscore`:
- `Y1=1`: samples HG00096–HG00101 (5 samples, group 1)
- `Y1=2`: samples HG00102–HG00107 (5 samples, group 2)

The phenotype split intentionally aligns with the PC cluster separation to verify that hue
colouring corresponds to the correct phenotype groups.
