# benchmark-common -- shared parts of the benchmark

The parts that do not depend on which tool is being compared: running this pipeline as
the reference method, scoring the results of any method, and comparing two methods
statistically.

The scripts in [../nf-gwas-benchmark/](../nf-gwas-benchmark/) are earlier **copies** of
some of these, deliberately left unchanged so that the nf-gwas comparison can still be
reproduced exactly as it was originally run.

## Files

| File | What it does |
|---|---|
| `run_nf_rare_var_assoc.sh` | Runs this pipeline over all datasets in a single `nextflow run`, with `--skip_preparation`. Writes results to `runs/nf_rare_var_assoc/results/...`, then scores them with `run_eval.sh`. |
| `run_eval.sh` | Scores results. Runs `nf-eval-gene-assoc` once, giving it file patterns that match every dataset; it matches results to the correct dataset by the `dataset_idx` in each filename. Its settings are passed as environment variables, documented in the file header. |
| `filter_causal_autosomal.py` | Makes a copy of the known-answer files that keeps only genes on the numbered chromosomes. Needed when a tool cannot test chromosome X. Only the copy is modified; the original files and the pipeline are untouched. |
| `pairwise_compare.py` | Compares two methods dataset by dataset (`--arm-a` / `--arm-b`). Writes `pairwise_per_dataset.csv` and `pairwise_summary.csv`, and, unless you pass `--no-plots`, the five figures listed below. `--missing {drop,zero}` decides what to do with datasets that one of the two methods failed to produce a result for. Requires pandas, scipy and matplotlib. |
| `project_power.py` | Estimates how many datasets a comparison would need before a difference this small could reach statistical significance, starting from an existing `pairwise_per_dataset.csv`. Uses both an analytical calculation and a bootstrap, and writes `power_projection_<metric>.{png,pdf}`. It assumes the true difference equals the one measured so far, and states this on the figure. |

## Figures written by `pairwise_compare.py`

Each figure is written twice, as .png (300 dpi) and .pdf (vector). Colours are a
diverging pair that remains distinguishable with colour-vision deficiency: blue means the
first method scored higher, red means the second did. Differences are always calculated
as first method minus second method.

| Figure name ends with | What it shows |
|---|---|
| `__paired_diff_ap` | The main figure: the difference in average precision for each dataset, one bar per dataset, coloured by which method won, with the mean and its 95% confidence interval |
| `__scatter_ap` | Average precision of one method against the other, one point per dataset, with a diagonal line marking equal performance |
| `__forest_metrics` | Mean difference and 95% confidence interval for average precision, AUC-PR and AUC-ROC |
| `__pvalue_vs_n` | How the p-value develops as datasets are added one at a time, for both a paired t-test and a Wilcoxon test |
| `__pvalue_vs_n_bootstrap` | The same, but repeated over many random dataset orderings (2000 by default, set with `--bootstrap`), showing the range of possible paths to the final p-value |

The p-values are not corrected for multiple testing. Average precision was chosen in
advance as the single main measure; AUC-PR, AUC-ROC and the Wilcoxon test are reported
only as supporting checks, so no correction applies.

The plotting scripts need matplotlib, which the system `python3` may not have. Use a
virtual environment that does -- the run scripts accept one through the `PAIRWISE_PYTHON`
environment variable.

## Scoring when a tool cannot test chromosome X

Some tools have no statistical model for chromosome X. To keep the comparison fair, both
methods are then scored against a version of the known answers that contains only genes
on the numbered chromosomes:

```bash
# 1) build the reduced set of known answers once
python filter_causal_autosomal.py \
    --datasets-dir .../tools_comparison/datasets \
    --out-dir      .../tools_comparison/datasets_autosomal

# 2) re-score existing reference results against it
#    (this only re-scores; the pipeline itself is not run again)
CAUSAL_SNPLIST_GLOB=".../datasets_autosomal/run_*/select_genes/*_dataset_idx_*_in_*.snplist" \
CAUSAL_GENES_GLOB=".../datasets_autosomal/run_*/select_genes/*_genes_dataset_idx_*.txt" \
EVAL_RUN_DIR=".../runs/nf_rare_var_assoc_autosomal_eval" \
  bash run_nf_rare_var_assoc.sh   # or call run_eval.sh directly on results you already have
```

The two resulting score directories are what you pass to `pairwise_compare.py` as
`--arm-a` and `--arm-b`.

## How to read the known-answer files

Each simulated dataset has two files in `datasets/run_<N>/select_genes/`, and it is easy
to mistake one for the other:

- `*_genes_dataset_idx_<N>.txt` lists **every gene that a causal gene could have been
  chosen from**, together with its variants. It is not the list of causal genes.
- `*_in_<GENE>_<GENE>...snplist` lists the **causal variants**.

To get the causal genes, filter the first file using the variants in the second.
`nf-eval-gene-assoc` receives both file patterns and does this itself.
