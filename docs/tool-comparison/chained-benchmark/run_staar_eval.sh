#!/usr/bin/env bash
#
# Score the RICOPILI + STAARpipeline results with nf-eval-gene-assoc, against the
# COMPLETE list of causal genes
# ============================================================================
# Runs benchmark-common/run_eval.sh once per relatedness version (filtered, full),
# each against the FULL causal truth: chrX and non-coding causal genes STAAR
# structurally cannot test therefore count as misses (recall penalty). This is
# This is the honest whole-method number, and the one that shares an axis with
# nf-rare-var-assoc and with the RICOPILI + nf-gwas comparison, which are both
# scored the same way. A number restricted to only the genes STAAR can test
# (which would measure the engine alone) is not produced here.
#
# nf-eval matches the known answers to the result tables by dataset_idx, so pointing
# pointing it at the complete run_*/ known answers but only at this method's result
# tables scores exactly the datasets it actually produced.
#
# proxy_scoring_mode defaults to 'none' (nf-eval-gene-assoc/nextflow.config), so AP
# is gene-level max-LOG10P vs causal genes and INPUT_VCF is not used for proxy
# mapping -- consistent with how the reference and nf-gwas arms were scored.
#
# Env overrides: DATA, EVAL_REPO, COMMON. Defaults are this workstation's layout.
set -euo pipefail

DATA="${DATA:-/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc}"
T="$DATA/tools_comparison"
COMMON="${COMMON:-/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc/docs/tool-comparison/benchmark-common}"
EVAL_REPO="${EVAL_REPO:-/data/git/doktorat_pw/wum_pims/nf-eval-gene-assoc}"
RAW_VCF="$T/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz"

export EVAL_REPO
export INPUT_VCF="$RAW_VCF"
export SKIP_PREP="true"
# The complete list of causal genes.
export CAUSAL_SNPLIST_GLOB="$T/datasets/run_*/select_genes/*_dataset_idx_*_in_*.snplist"
export CAUSAL_GENES_GLOB="$T/datasets/run_*/select_genes/*_genes_dataset_idx_*.txt"

for V in "${VERSIONS:-filtered full}"; do
  for VER in $V; do
    echo "========================================================"
    echo " STAAR eval: version=$VER   $(date -Is)"
    echo "========================================================"
    export EVAL_RUN_DIR="$T/runs/ricopili_staar_qcmatched_${VER}_eval"
    export EVAL_PROJECT="ricopili_staar_qcmatched_${VER}"
    export REGENIE_GLOB="$T/runs/ricopili_staar_qcmatched/regenie_per_dataset/${VER}/ricopili_staar_qcmatched_${VER}_dataset_idx_*_step2_Y1.regenie"
    echo "REGENIE tables:"; ls -1 $REGENIE_GLOB
    bash "$COMMON/run_eval.sh"
    echo " STAAR eval version=$VER DONE  $(date -Is)"
  done
done
echo "ALL STAAR EVALS DONE $(date -Is)"
