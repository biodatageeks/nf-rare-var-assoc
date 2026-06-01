# bcftools/view test fixtures

Fixture files for `modules/local/bcftools/view/tests/main.nf.test`.

## regions_first_5kb.bed

BED-format (0-based) regions file covering `MT192765.1:1-5000`.  
Used by the regions-filter test to verify that only records within the
first 5 kb of the sarscov2 chromosome pass the filter.

Source: hand-crafted from the known positions in the nf-core sarscov2
`test.vcf.gz` (9 records on MT192765.1).  Two records fall in this window:
position 197 and position 4788.  The remaining 7 records (positions ≥ 8236)
are outside the window and are filtered out.

## samples_test.txt

Plain-text samples file containing one sample name: `test`.  
Used by the samples-filter test to exercise the `--samples-file` code path.  
The sarscov2 `test.vcf.gz` has exactly one sample (`test`), so the output
sample count stays at 1; the test asserts the header names that sample
correctly.
