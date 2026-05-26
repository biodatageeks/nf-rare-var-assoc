# Test fixtures for `bcftools/norm`

Both fixtures use chromosome **MT192765.1** (sarscov2) so they can be used with the nf-core
sarscov2 genome FASTA (`genomics/sarscov2/genome/genome.fasta.gz`) as `--fasta-ref`.
REF alleles were extracted from that FASTA at the respective positions.

| File | Records | Description |
|---|---|---|
| `multiallelic_sarscov2.vcf.gz` | 8 | 5 biallelic + 3 multiallelic sites (pos 200: 2 ALTs, 400: 2 ALTs, 600: 3 ALTs). After `bcftools norm -m -any` → 12 records. |
| `split_sarscov2.vcf.gz` | 12 | All-biallelic, pre-split form of the above. After `bcftools norm -m +any` → 8 records. |
