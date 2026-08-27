# Tool comparison

This directory contains everything needed to reproduce the comparison between
nf-rare-var-assoc and other rare-variant analysis software: the scripts that run each
tool, the container definitions for the tools that ship none, and the scoring and
statistical comparison code.

It is not part of the pipeline. Nothing here is needed to run nf-rare-var-assoc; see
[../usage.md](../usage.md) for that.

## How the comparison is organised

No single existing tool covers the same range of steps as nf-rare-var-assoc, which
processes a raw VCF all the way to annotated rare-variant association results. Tools fall
into two groups: those that do thorough quality control but have no rare-variant
association test, and those with a strong association test but little or no data
preparation. Comparing against either group alone would be unfair in one direction or the
other.

The comparison therefore combines two tools into one complete alternative and compares
that against this pipeline. Two such combinations were run, plus one earlier,
narrower comparison:

| Directory | What it compares | Notes |
|---|---|---|
| [chained-benchmark/](chained-benchmark/) | RICOPILI (quality control and principal components) combined with STAARpipeline (association) | STAARpipeline brings its own functional annotations, so this is the only comparison where the annotation step is entirely the other tool's. Run in two versions, one removing related samples and one modelling them. |
| [chained-benchmark-nf-gwas/](chained-benchmark-nf-gwas/) | RICOPILI combined with [nf-gwas](https://github.com/genepi/nf-gwas) | Both halves are complete pipelines a user configures and runs, with no analysis code of ours between them. The gene groupings are still borrowed, because nf-gwas cannot build them. |
| [nf-gwas-benchmark/](nf-gwas-benchmark/) | nf-rare-var-assoc against nf-gwas alone | The earliest comparison. nf-gwas is given this pipeline's gene groupings *and* principal components, so it measures only the quality-control and dosage stages, not the tools as a whole. The two combined comparisons above exist because of this limitation. Kept unchanged so its result remains reproducible. |
| [benchmark-common/](benchmark-common/) | -- | The parts that do not depend on which tool is being compared: running this pipeline as the reference, scoring any method's results, and comparing two methods statistically. |

Each directory has its own README with prerequisites, the exact commands, and the
practical problems worth knowing about in advance. Start with
[benchmark-common/](benchmark-common/), which explains how results are scored and how the
known answers in the simulated datasets must be read.

## Common ground rules

These apply to every comparison here.

- **Each tool does its own work.** A compared method receives only the raw VCF, the
  phenotype files, and the list of causal genes used for scoring. Its own quality control,
  population-structure correction, annotation and statistical test are used, configured as
  well as we could configure them. Where a tool is given something it cannot produce
  itself, this is stated explicitly in that directory's README and counted as a limitation
  of the comparison, not hidden.

- **Everything runs in containers.** Where a tool publishes no usable container image, a
  Dockerfile is included here.

- **Correcting for population structure is a precondition, not a refinement.** Two
  independent tools measured a genomic inflation factor above 8 on this data with no
  structure correction. Scoring an uncorrected run against a corrected one would be
  meaningless, so the inflation factor is measured for every dataset in every method,
  every time.

- **Results are compared dataset by dataset**, using recall-scaled average precision as
  the single measure chosen in advance, with paired t-tests and Wilcoxon tests. AUC-PR and
  AUC-ROC are reported as supporting checks.

## Working data

The genotype data, the simulated datasets, the annotation databases and all run outputs
live outside this repository and are not committed. Every script defaults to the layout of
the machine they were developed on, but all locations are environment variables -- in most
cases setting `DATA` and `RVA_REPO` is enough to run them elsewhere unchanged.
