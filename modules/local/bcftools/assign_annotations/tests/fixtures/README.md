# Fixtures for bcftools/assign_annotations test

## test_variants.vcf.gz (.tbi)

A hand-crafted VCF with 21 variants spanning 6 genes (and 6 distinct VEP consequences) on
chromosome 1. Compressed with bgzip, indexed with tabix.

### Design

The CSQ INFO field uses the same 52-field VEP format as the production prepared VCFs
(Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|...|DISTANCE|...).
Only the fields used by `assign_annotations.py` carry meaningful values (SYMBOL,
Consequence, Feature_type, Feature, DISTANCE); the rest are empty.

Frequencies are deliberately distinct so frequency tie-breaking never matters:

| Gene  | Consequence            | # variants | Mask category | In mask? |
|-------|------------------------|------------|---------------|----------|
| GeneA | stop_gained            | 6          | LoF           | yes      |
| GeneB | missense_variant       | 5          | missense      | yes      |
| GeneC | synonymous_variant     | 4          | (none)        | no       |
| GeneD | intron_variant         | 3          | (none)        | no       |
| GeneE | upstream_gene_variant  | 2          | (none)        | no       |
| GeneF | 5_prime_UTR_variant    | 1          | (none)        | no       |

With `quantile_threshold=0.25` the (R type=7) threshold over the sorted frequencies
`[1, 2, 3, 4, 5, 6]` is `2.25`, so:
- above-quantile: stop_gained (6), missense_variant (5), synonymous_variant (4), intron_variant (3)
- below-quantile: upstream_gene_variant (2), 5_prime_UTR_variant (1)

The two important consequences are both also above the quantile, so the fixture naturally
exercises the overlap branch of `filter_annotations`.

### Why several test scenarios

Each test in `main.nf.test` drives `filter_annotations` into a different branch by varying
`--min-top-annotations` / `--max-annotations` while keeping the same fixture:

| Test | Args | Branch under test |
|------|------|-------------------|
| 1 | `--min-top-annotations 4` | normal range + important/quantile overlap (no padding, no trimming) |
| 2 | `--min-top-annotations 2 --max-annotations 3` | union exceeds max: quantile set trimmed by frequency |
| 3 | `--min-top-annotations 6 --max-annotations 10` | union below min: padding pulls in below-quantile consequences |

## test_masks.tsv

Tab-separated masks file with two categories:
- `LoF`: stop_gained, frameshift_variant, splice_acceptor_variant, splice_donor_variant
- `missense`: missense_variant

The remaining four consequences (synonymous, intron, upstream, 5' UTR) are intentionally
absent from every mask so the mask-membership branch is exercised.
