# Reporting subworkflow test fixtures

Minimal hand-rolled fixtures for IT-5 (smoke test; reporting carve-out).

## Y1.regenie.gz

Tab-separated gzip file in MERGE_RESULTS output format (no `##MASKS=` header).
Two genes (GeneA, GeneB), two masks (Mask_Mod, Mask_High), two AAF thresholds (0.05, 0.01).
Includes one BURDEN-SKATO-joint row to exercise the joint-test code path in the Rmd.

## masks.txt

Two masks matching the IDs used in Y1.regenie.gz: `Mask_Mod` and `Mask_High`.

## pheno.tsv

Six samples (3 cases Y1=1, 3 controls Y1=2) from 1kGP sample IDs.

## pc_plot.png / eda_plot.png

Minimal valid 1x1 white pixel PNGs used as placeholder inputs for the R Markdown report.
The report embeds them via `knitr::include_graphics`; content is irrelevant for the smoke test.

## test.setlist

Two-gene Regenie setlist: gene, chr, pos, comma-separated variant IDs.
Matches the gene names in Y1.regenie.gz.
