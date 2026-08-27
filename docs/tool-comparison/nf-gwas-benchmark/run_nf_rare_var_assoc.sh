#!/usr/bin/env bash
#
# Benchmark run: nf-rare-var-assoc (this pipeline) -- MULTI-DATASET
# ================================================================
# Runs the pipeline across MANY simulated phenotype datasets in a SINGLE
# `nextflow run`. nf-rare-var-assoc natively fans out over a comma-separated list
# of phenotype files (input_phenotype is tokenized on ',' and Cartesian-producted
# with the VCF -- see subworkflows/local/utils_nfcore_rare-var-assoc_pipeline/main.nf),
# so each `_dataset_idx_N_` draw runs the full QC/PCA/grouping/association under
# meta.id = "${project_name}_dataset_idx_${N}".
#
# PREPARATION IS SKIPPED. We feed the already-prepared (split + CSQ-annotated +
# DS-bearing) VCF and set `--skip_preparation true`, so the heavy VEP/normalisation
# stage does NOT rerun per dataset. This is the param-tuning fast path (mirrors
# conf/test_skip_preparation_and_reporting.config). The filter_vcf_* genotype QC
# still runs -- it lives in the main pipeline, not in prep.
#
# This runs nf-rare-var-assoc as the reference method of the nf-gwas comparison.
# See README.md in this directory.
#
# IMPORTANT: launched with `--publish_intermediate true` so the per-dataset
# REGENIE-format mask files and PC covariate files are on disk for the nf-gwas arm
# to REUSE (one set PER dataset_idx):
#   - <outdir>/bcftools/<id>.annotations            (REGENIE --anno-file)
#   - <outdir>/bcftools/<id>.setlist                (REGENIE --set-list)
#   - <outdir>/merge/<id>_sex_covar.txt             (PC covariate table)
#   - <outdir>/regenie_step2/<id>_step2_Y1.regenie  (association results, per dataset)
# where <id> = tools_comparison_dataset_idx_<N>. The mask definition itself is the
# static assets/default.masks in this repo.
#
# After the pipeline, nf-eval-gene-assoc is run ONCE over ALL datasets: its inputs
# are passed as GLOBS (it parses dataset_idx from filenames and inner-joins
# causal_snplist/causal_genes/regenie_results by meta.id), so a single invocation
# scores every dataset. The inner join is constrained by the regenie glob, which
# only matches datasets this run actually produced.
#
# CACHING/SPACE: this pipeline runs with the `nocache` profile and WITHOUT
# `-resume` (a fresh run every time). After it finishes, the big intermediate
# result dirs (see CLEANUP_SUBDIRS) are deleted to reclaim space; none of them are
# needed by the eval arm or by nf-gwas.
#
# DATASET SELECTION: set DATASET_IDXS (space-separated) to override; default is all
# run_<N> dirs found under $DATASETS_DIR.
#
set -euo pipefail

# ----------------------------------------------------------------------------
# Configuration (edit paths here)
# ----------------------------------------------------------------------------

DATA="${DATA:-/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/tools_comparison}"
RVA_REPO="${RVA_REPO:-/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc}"
EVAL_REPO="${EVAL_REPO:-/data/git/doktorat_pw/wum_pims/nf-eval-gene-assoc}"
DATASETS_DIR="${DATA}/datasets"

# Prepared (split + CSQ + DS) VCF -> input for skip_preparation=true.
PREPARED_VCF="${DATA}/prepared.vcf.gz"
# Original unprepared exome VCF -> input for the eval arm (it runs its own VEP).
INPUT_VCF="${DATA}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz"

PROJECT="tools_comparison"
RUN_DIR="${DATA}/runs/nf_rare_var_assoc"
EVAL_RUN_DIR="${DATA}/runs/nf_rare_var_assoc_eval"
RVA_PROFILE="podman,medium_resources,nocache"   # nocache: fresh run, no work-dir cache reuse
EVAL_PROFILE="podman,medium_resources"

# Big intermediate result dirs to delete after the run (not needed downstream).
CLEANUP_SUBDIRS=(
    bcftools_view_and_filter2
    bcftools_replace_sample_names
    bcftools_view
    check_x_chrom_present
    plink2_export_other
    plink2_makepgen
)
CLEANUP="${CLEANUP:-true}"   # set CLEANUP=false to keep intermediates

# ----------------------------------------------------------------------------
# Dataset selection + per-dataset file resolution
# ----------------------------------------------------------------------------
# Each datasets/run_<N> dir holds (N == dataset_idx):
#   gcta_simu/tuner_base_run_<N>_dataset_idx_<N>_gcta_simu.phenotype.txt
#   select_genes/tuner_base_run_<N>_select_genes_genes_dataset_idx_<N>.txt
#   select_genes/tuner_base_run_<N>_select_genes_snps_dataset_idx_<N>_in_*.snplist
if [[ -z "${DATASET_IDXS:-}" ]]; then
    DATASET_IDXS="$(ls -1d "${DATASETS_DIR}"/run_*/ 2>/dev/null \
        | sed -E 's#.*/run_([0-9]+)/#\1#' | sort -n | tr '\n' ' ')"
fi
read -r -a IDXS <<< "$DATASET_IDXS"
[[ ${#IDXS[@]} -gt 0 ]] || { echo "ERROR: no datasets selected (DATASETS_DIR=$DATASETS_DIR)" >&2; exit 1; }

pheno_path() {
    echo "${DATASETS_DIR}/run_$1/gcta_simu/tuner_base_run_$1_dataset_idx_$1_gcta_simu.phenotype.txt"
}

# Build the comma-separated phenotype list (one entry per selected dataset).
PHENO_LIST=""
for idx in "${IDXS[@]}"; do
    p="$(pheno_path "$idx")"
    [[ -e "$p" ]] || { echo "ERROR: missing phenotype file for run_${idx}: $p" >&2; exit 1; }
    PHENO_LIST="${PHENO_LIST:+${PHENO_LIST},}${p}"
done

# ----------------------------------------------------------------------------
# Sanity checks
# ----------------------------------------------------------------------------
for f in "$PREPARED_VCF" "$INPUT_VCF" "$RVA_REPO/main.nf" "$EVAL_REPO/main.nf"; do
    [[ -e "$f" ]] || { echo "ERROR: missing required path: $f" >&2; exit 1; }
done
mkdir -p "$RUN_DIR" "$EVAL_RUN_DIR"

# ----------------------------------------------------------------------------
# Parameter file
# ----------------------------------------------------------------------------
PARAMS_FILE="${RUN_DIR}/params.nf_rare_var_assoc.yaml"
cat > "$PARAMS_FILE" <<EOF
# Auto-generated by run_nf_rare_var_assoc.sh -- edit the script, not this file.
input_vcf: "${PREPARED_VCF}"
input_phenotype: "${PHENO_LIST}"
project_name: "${PROJECT}"
outdir: "${RUN_DIR}/results"

# Prepared VCF in, preparation skipped (genotype filter_vcf_* QC still runs).
skip_preparation: true

# Expose intermediates (mask files + covar) for the nf-gwas arm to reuse.
publish_intermediate: true

use_dosage: true

filter_vcf_qual_min: 23
filter_vcf_avg_gq_min: 23
filter_vcf_avg_dp_min: 23
filter_vcf_avg_dp_max: 107
filter_vcf_sample_gq_min: 11
filter_vcf_sample_dp_min: 22
filter_vcf_sample_dp_max: 339

inbreeding_outliers_range_stds: 6

plink2_makepgen_3_options: "--geno 0.250 --hwe 1e-9 0.01 --mac 16 --maf 0.045000"
plink2_makepgen_4_options: "--geno 0.100"
plink2_makepgen_5_options: "--mind 0.200"
plink2_write_snplist_qc_options: "--mind 0.150"
plink2_indep_pairwise_options: "--mind 0.300"
plink2_indep_pairwise_window: "800 80 0.2"
plink2_missing_per_pheno_options: "--geno 0.350"
plink2_indep_pairwise_window_pca: "800 80 0.2"
plink2_king_cutoff_threshold_pca: 0.19
plink2_write_snplist_step2_options: "--mind 0.200"

regenie_step1_options: "--bt --bsize 400 --covarColList PC1_AVG,PC2_AVG,PC3_AVG,PC4_AVG"
regenie_step2_options: "--bt --minMAC 7 --ref-first --firth --approx --bsize 200 --aaf-bins 0.2 --vc-tests skato --covarColList PC1_AVG,PC2_AVG,PC3_AVG,PC4_AVG"
EOF

# Extra config: make podman robust on this host (verified: rootless podman needs
# keep-id userns so containers can write to host-owned work dirs).
EXTRA_CFG="${RUN_DIR}/podman_userns.config"
cat > "$EXTRA_CFG" <<'EOF'
podman.runOptions = '--userns=keep-id'
EOF

# ----------------------------------------------------------------------------
# Launch
# ----------------------------------------------------------------------------
echo "=================================================================="
echo " nf-rare-var-assoc benchmark run (multi-dataset, skip_preparation)"
echo "   project  : ${PROJECT}"
echo "   input    : ${PREPARED_VCF} (prepared; --skip_preparation)"
echo "   datasets : ${IDXS[*]}  (${#IDXS[@]} total)"
echo "   outdir   : ${RUN_DIR}/results"
echo "   params   : ${PARAMS_FILE}"
echo "   started  : $(date -Is)"
echo "=================================================================="
echo ""

cd "$RVA_REPO"
nextflow run "$RVA_REPO/main.nf" \
    -profile "$RVA_PROFILE" \
    -c "$EXTRA_CFG" \
    -params-file "$PARAMS_FILE" \
    -work-dir "${RUN_DIR}/work" \
    -with-trace "${RUN_DIR}/trace.txt" \
    -with-report "${RUN_DIR}/report.html" \
    -with-timeline "${RUN_DIR}/timeline.html"

echo ""
echo "Finished: $(date -Is)"
echo ""
echo "Per-dataset files the nf-gwas arm will reuse (one set per dataset_idx):"
ls -1 "${RUN_DIR}/results/bcftools/"*_dataset_idx_*.annotations 2>/dev/null || echo "  MISSING: *.annotations"
ls -1 "${RUN_DIR}/results/bcftools/"*_dataset_idx_*.setlist     2>/dev/null || echo "  MISSING: *.setlist"
ls -1 "${RUN_DIR}/results/merge/"*_dataset_idx_*_sex_covar.txt   2>/dev/null || echo "  MISSING: *_sex_covar.txt"

# ----------------------------------------------------------------------------
# Cleanup: drop big intermediate result dirs (not needed by eval or nf-gwas).
# ----------------------------------------------------------------------------
if [[ "$CLEANUP" == "true" ]]; then
    echo ""
    echo "[cleanup] removing big intermediate result dirs under ${RUN_DIR}/results/ ..."
    for sub in "${CLEANUP_SUBDIRS[@]}"; do
        d="${RUN_DIR}/results/${sub}"
        if [[ -d "$d" ]]; then
            echo "  rm -rf ${d}"
            rm -rf "$d"
        fi
    done
fi

# ----------------------------------------------------------------------------
# Evaluation: ONE nf-eval-gene-assoc run over ALL datasets (glob inputs).
# ----------------------------------------------------------------------------
# nf-eval-gene-assoc parses dataset_idx from filenames and inner-joins
# causal_snplist/causal_genes/regenie_results by meta.id, so globs fan out per
# dataset in a single invocation. The regenie glob constrains the join to the
# datasets this run produced (broad causal globs over run_*/ are inner-joined away).
EXTRA_CFG_EVAL="${EVAL_RUN_DIR}/podman_userns.config"
cat > "$EXTRA_CFG_EVAL" <<'EOF'
podman.runOptions = '--userns=keep-id'
EOF

CAUSAL_SNPLIST_GLOB="${DATASETS_DIR}/run_*/select_genes/*_dataset_idx_*_in_*.snplist"
CAUSAL_GENES_GLOB="${DATASETS_DIR}/run_*/select_genes/*_genes_dataset_idx_*.txt"
REGENIE_GLOB="${RUN_DIR}/results/regenie_step2/*_dataset_idx_*_step2_Y1.regenie"

echo ""
echo "=================================================================="
echo " nf-eval-gene-assoc -- scoring ALL datasets in one fan-out run"
echo "   regenie : ${REGENIE_GLOB}"
echo "   outdir  : ${EVAL_RUN_DIR}/results"
echo "=================================================================="

cd "$EVAL_REPO"
nextflow run "$EVAL_REPO/main.nf" \
    -profile "$EVAL_PROFILE" \
    -c "$EXTRA_CFG_EVAL" \
    -work-dir "${EVAL_RUN_DIR}/work" \
    --project_name "${PROJECT}_eval" \
    --input_vcf "${PREPARED_VCF}" \
    --causal_snplist "${CAUSAL_SNPLIST_GLOB}" \
    --causal_genes "${CAUSAL_GENES_GLOB}" \
    --regenie_results "${REGENIE_GLOB}" \
    --outdir "${EVAL_RUN_DIR}/results" \
    --skip_preparation true

echo ""
echo "nf-rare-var-assoc evaluation Finished: $(date -Is)"
