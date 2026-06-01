#!/usr/bin/env python3
"""
Generate synthetic VCF fixture for IT-4 f_coefficient_filtering test.

Design
------
100 biallelic SNPs on chr1 (positions 1000, 2000, ..., 100000).
62 samples: 60 normals (NORM_01..NORM_60) + 2 high-F spikes (SPIKE_HIGH1, SPIKE_HIGH2).
Minimum 50 samples required by plink2 for --indep-pairwise LD estimation.

Genotype assignment
  Normals: each sample gets exactly 25 homozygous-ref (0/0), 25 homozygous-alt (1/1),
    and 50 heterozygous (0/1) sites, shuffled with a per-sample random seed. This gives
    each normal a reproducible but distinct pattern across variants, minimising inter-variant
    LD so that most variants survive --indep-pairwise 50 5 0.2.
  Spikes: always 1/1 at all 100 sites.

Expected F statistics (plink2 method-of-moments)
  For each variant, approximate allele freq p ~ 0.5 (normals balanced, spikes shift slightly).
  E(HOM) per variant ~ 0.5.
  E(HOM) summed over 100 variants ~ 50.

  Normal sample: O_HOM = 50 -> F ~ (50-50)/(100-50) = 0
  Spike sample:  O_HOM = 100 -> F ~ (100-50)/(100-50) = 1.0

  With inbreeding_outliers_range_stds=3, the distribution of 60 normals (F~0) + 2 spikes
  (F=1.0) has std ~ 0.18. upper_bound = mean + 3*std ~ 0.55. Both spikes (F=1.0) are
  caught; all normals (F~0) are not.

Usage
-----
python3 gen_fixture.py > spiked_fixture.vcf
Then convert to pgen (must be run from this directory or adjust paths).
--double-id sets FID=IID so that calc_f_outliers writes FID+IID columns and
plink2 --remove can match samples by the FID field:
  plink2 --vcf spiked_fixture.vcf --double-id --out spiked_fixture --make-pgen
"""

import random

N_VARIANTS = 100
N_NORMALS = 60
SPIKE_IIDS = ["SPIKE_HIGH1", "SPIKE_HIGH2"]
SEED = 42

samples = [f"NORM_{i:02d}" for i in range(1, N_NORMALS + 1)] + SPIKE_IIDS

# Build per-normal genotype columns: each normal has exactly 25 hom-ref, 25 hom-alt, 50 het.
normal_gts_by_sample = []  # index: [sample_idx][variant_idx] -> gt string
for j in range(N_NORMALS):
    rng = random.Random(SEED + j * 997)
    indices = list(range(N_VARIANTS))
    rng.shuffle(indices)
    gts = ["0/1"] * N_VARIANTS
    for i in indices[:25]:
        gts[i] = "0/0"
    for i in indices[25:50]:
        gts[i] = "1/1"
    normal_gts_by_sample.append(gts)

header_lines = [
    "##fileformat=VCFv4.1",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
]

col_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)

lines = header_lines + [col_header]

for i in range(N_VARIANTS):
    pos = (i + 1) * 1000
    variant_id = f"chr1_{pos}_A_T"
    gts = [normal_gts_by_sample[j][i] for j in range(N_NORMALS)]
    for _ in SPIKE_IIDS:
        gts.append("1/1")
    row = f"chr1\t{pos}\t{variant_id}\tA\tT\t.\tPASS\t.\tGT\t" + "\t".join(gts)
    lines.append(row)

print("\n".join(lines))
