#!/usr/bin/env bash
#
# Score the RICOPILI + STAARpipeline results with nf-eval-gene-assoc, against the
# COMPLETE list of causal genes
# ============================================================================
# Runs benchmark-common/run_eval.sh once per ARM, where an arm is a
# "<relatedness version>_<spa mode>" label produced by run_chain_ricopili_staar.sh
# (default: full_nospa and full_spa). Each is scored against the FULL causal truth:
# chrX and non-coding causal genes STAAR structurally cannot test therefore count as
# misses (recall penalty). This is
# This is the honest whole-method number, and the one that shares an axis with
# nf-rare-var-assoc and with the RICOPILI + nf-gwas comparison, which are both
# scored the same way. A number restricted to only the genes STAAR can test
# (which would measure the engine alone) is not produced here.
#
# CHOOSING BETWEEN THE SPA MODES: score both, then pick ONE for the whole arm by
# mean recall-scaled AP over every dataset -- never per dataset, which would be
# cherry-picking. Both are reported; the paper states STAAR was given the better of
# two configurations, which biases against the reference pipeline by construction.
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

# Which run to score. Must match run_chain_ricopili_staar.sh's RUN_DATE/RUN_DIR --
# METHOD is the single name every table, eval project and pairwise arm derives from.
RUN_DATE="${RUN_DATE:-2026_09_12}"
RUN_DIR="${RUN_DIR:-$T/runs/ricopili_staar_${RUN_DATE}}"
METHOD="$(basename "$RUN_DIR")"
# Arm labels, not bare versions: "<relatedness version>_<spa mode>".
ARMS="${ARMS:-full_nospa full_spa}"
read -r -a ARM_ARR <<< "$ARMS"

[[ -d "$RUN_DIR" ]] || { echo "ERROR: no run dir at $RUN_DIR (set RUN_DATE or RUN_DIR)" >&2; exit 1; }

for ARM_LABEL in "${ARM_ARR[@]}"; do
    echo "========================================================"
    echo " STAAR eval: arm=$ARM_LABEL   $(date -Is)"
    echo "========================================================"
    export EVAL_RUN_DIR="$T/runs/${METHOD}_${ARM_LABEL}_eval"
    export EVAL_PROJECT="${METHOD}_${ARM_LABEL}"
    export REGENIE_GLOB="${RUN_DIR}/regenie_per_dataset/${ARM_LABEL}/${METHOD}_${ARM_LABEL}_dataset_idx_*_step2_Y1.regenie"
    # shellcheck disable=SC2086  # deliberate glob expansion
    if ! ls -1 $REGENIE_GLOB >/dev/null 2>&1; then
        echo "  no result tables for arm=${ARM_LABEL} -- skipping" >&2
        continue
    fi
    echo "REGENIE tables:"; ls -1 $REGENIE_GLOB
    bash "$COMMON/run_eval.sh"
    echo " STAAR eval arm=$ARM_LABEL DONE  $(date -Is)"
done
echo "ALL STAAR EVALS DONE $(date -Is)"
