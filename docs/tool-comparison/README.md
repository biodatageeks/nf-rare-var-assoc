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

The genotype data, the annotation databases and all run outputs live outside this
repository and are not committed. The simulated datasets are: see the next section. Every
script defaults to the layout of the machine they were developed on, but all locations are
environment variables -- in most cases setting `DATA` and `RVA_REPO` is enough to run them
elsewhere unchanged.

## The simulated datasets

The simulated phenotypes and the known answers they are scored against are committed, in
[`assets/tool_comparison_datasets/`](../../assets/tool_comparison_datasets/). They were
produced by nf-gene-sim-assoc from the 1000 Genomes high-coverage exome VCF (chromosomes
12, 22 and X); the genotypes themselves are not committed.

There are 30 datasets, numbered `run_0` to `run_31` with 1 and 6 absent. The number `<N>`
in a directory name is the same number that appears as `dataset_idx_<N>` in every filename
and in every result the comparison produces.

```
assets/tool_comparison_datasets/
├── gen_params_list.json                    the parameter space the datasets were drawn from
└── run_<N>/
    ├── run_random_params.json              the parameters this dataset was drawn with,
    │                                       including its --simu-cc, --simu-hsq and --simu-k
    ├── run_random_params_pdfs.json         the distributions they were drawn from
    ├── gcta_simu/
    │   └── tuner_base_run_<N>_dataset_idx_<N>_gcta_simu.phenotype.txt
    └── select_genes/
        ├── tuner_base_run_<N>_select_genes_snps_dataset_idx_<N>_in_<GENES>.snplist
        └── tuner_base_run_<N>_select_genes_genes_dataset_idx_<N>.txt
```

| File | Contents |
|---|---|
| `*_gcta_simu.phenotype.txt` | The simulated phenotype: `FID`, `IID`, `Y1`, tab separated, with a header. `Y1` is 0 for a control and 1 for a case. Between 1,344 and 2,000 samples per dataset, with 11 to 667 cases. |
| `*_select_genes_snps_dataset_idx_<N>.snplist` | **The causal variants**: one identifier per line, `CHROM_POS_REF_ALT` with the `chr` prefix stripped, matching the identifiers the pipeline assigns. Between 2 and 196 per dataset (median 41). |
| `*_select_genes_genes_dataset_idx_<N>.txt` | `Gene`, `Variant`, tab separated, with a header. **These are not the causal genes.** It lists every gene the causal variants were drawn from and all of that gene's candidate variants -- up to 2,680 genes and 87,000 rows in one dataset. |

### Deriving the causal genes

A gene is causal when one of its variants in the `*_genes_*.txt` file appears in that
dataset's `.snplist`. The gene list on its own says nothing about which genes are causal,
and scoring against it directly would treat thousands of untouched genes as true answers.

```bash
cd assets/tool_comparison_datasets
awk 'NR==FNR {causal[$1]; next} FNR>1 && ($2 in causal) {print $1}' \
    run_0/select_genes/*_snps_dataset_idx_0_in_*.snplist \
    run_0/select_genes/*_genes_dataset_idx_0.txt | sort -u
```

Do not take the gene names from the `.snplist` filename either. It records the genes the
selection started from, which is not always the set that ends up carrying a causal
variant: `run_0`'s file is named `..._in_RBMXL3_PTPRQ.snplist`, while the join above
returns `LRCH2`, `PTPRQ` and `RBMXL3`, because its chromosome X variant falls inside two
overlapping genes.

Done this way there are 2 to 27 causal genes per dataset, median 9, and 319 across all 30
-- the number the scoring reports as `total_causal_genes`. The run scripts hand both file
patterns to the scorer (`CAUSAL_SNPLIST_GLOB` and `CAUSAL_GENES_GLOB`, see
[benchmark-common/](benchmark-common/)), which performs this join itself.

The scripts read the datasets from `DATASETS_DIR`, which defaults to `${DATA}/datasets`,
so point that at this directory (or copy it there) when reproducing a run.
