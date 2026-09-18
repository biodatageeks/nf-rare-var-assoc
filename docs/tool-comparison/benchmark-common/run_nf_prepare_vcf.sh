#!/usr/bin/env bash
#
# Prepares the raw exome VCF once with nf-prepare-vcf, using the same parameters file as
# this pipeline's own preparation step. The output,
# ${OUTDIR}/bcftools_reheader/prepare_vcf_reheader.vcf.gz, is the PREPARED_VCF input of
# run_nf_rare_var_assoc.sh.
#
set -euo pipefail

DATA="${DATA:-/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/tools_comparison}"
PV_REPO="${PV_REPO:-/data/git/doktorat_pw/wum_pims/nf-prepare-vcf}"
RVA_REPO="${RVA_REPO:-/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc}"
INPUT_VCF_RAW="${INPUT_VCF_RAW:-${DATA}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz}"
RUN_DATE="${RUN_DATE:-2026_09_17}"
OUTDIR="${OUTDIR:-${DATA}/nf_prepare_vcf_results_${RUN_DATE}}"
WORK_DIR="${WORK_DIR:-${DATA}/nf_prepare_vcf_work_${RUN_DATE}}"

for f in "$INPUT_VCF_RAW" "$PV_REPO/main.nf" "$RVA_REPO/conf/nf_prepare_params.yml"; do
    [[ -e "$f" ]] || { echo "ERROR: missing required path: $f" >&2; exit 1; }
done

# Run from the output directory so .nextflow.log lands beside the results.
mkdir -p "$OUTDIR"
cd "$OUTDIR"

nextflow run "${PV_REPO}/main.nf" \
    -profile podman,medium_resources,nocache \
    -params-file "${RVA_REPO}/conf/nf_prepare_params.yml" \
    --input_vcf "${INPUT_VCF_RAW}" \
    --outdir "${OUTDIR}" \
    --calc_ds_min_gq 3 \
    --project_name "prepare_vcf" \
    -work-dir "${WORK_DIR}" \
    -ansi-log false
