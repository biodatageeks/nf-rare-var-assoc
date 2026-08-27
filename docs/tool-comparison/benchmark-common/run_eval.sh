#!/usr/bin/env bash
#
# Shared evaluation fan-out: nf-eval-gene-assoc over ALL datasets in ONE run
# ==========================================================================
# Tool-agnostic. Hoisted from the identical eval block that used to be inlined
# at the tail of BOTH nf-gwas-benchmark/run_nf_rare_var_assoc.sh and
# nf-gwas-benchmark/run_nf_gwas.sh. Any comparison calls this once, after it
# has produced a per-dataset REGENIE-format results glob.
#
# nf-eval-gene-assoc parses dataset_idx from filenames and inner-joins
# causal_snplist / causal_genes / regenie_results by meta.id, so a single
# invocation with GLOB inputs scores every dataset. The inner join is
# constrained by the regenie file pattern (only datasets that actually produced a result).
#
# Interface is ENV VARS (so the caller can export them and just `bash run_eval.sh`):
#   EVAL_REPO            nf-eval-gene-assoc checkout (has main.nf)          [required]
#   EVAL_RUN_DIR        output and working directory for this scoring run  [required]
#   EVAL_PROJECT        project_name (a "_eval" suffix is appended)        [required]
#   INPUT_VCF            VCF handed to the scoring pipeline                 [required]
#   REGENIE_GLOB         glob of *_dataset_idx_<N>_..._Y1.regenie          [required]
#   CAUSAL_SNPLIST_GLOB  glob of causal *_in_*.snplist                     [required]
#   CAUSAL_GENES_GLOB    glob of causal *_genes_dataset_idx_*.txt          [required]
#   INPUT_VCF_TBI        VCF index (only if VEP has to run here)           [optional]
#   SKIP_PREP            "true" -> pass --skip_preparation true            [optional]
#   EVAL_PROFILE         nextflow profile (default: podman,medium_resources)
#
# NOTE on the two historical callers (kept identical here):
#   - nf-rare-var-assoc : SKIP_PREP=true, INPUT_VCF=prepared.vcf.gz, no index
#   - nf-gwas           : SKIP_PREP unset, INPUT_VCF=raw exome VCF, INPUT_VCF_TBI set
# The RICOPILI + STAARpipeline comparison uses the first set of settings but is
# pointed at a set of known answers restricted to the numbered chromosomes (see
# filter_causal_autosomal.py).
#
set -euo pipefail

: "${EVAL_REPO:?set EVAL_REPO}"
: "${EVAL_RUN_DIR:?set EVAL_RUN_DIR}"
: "${EVAL_PROJECT:?set EVAL_PROJECT}"
: "${INPUT_VCF:?set INPUT_VCF}"
: "${REGENIE_GLOB:?set REGENIE_GLOB}"
: "${CAUSAL_SNPLIST_GLOB:?set CAUSAL_SNPLIST_GLOB}"
: "${CAUSAL_GENES_GLOB:?set CAUSAL_GENES_GLOB}"
EVAL_PROFILE="${EVAL_PROFILE:-podman,medium_resources}"

[[ -e "$EVAL_REPO/main.nf" ]] || { echo "ERROR: no main.nf under EVAL_REPO=$EVAL_REPO" >&2; exit 1; }
[[ -e "$INPUT_VCF" ]] || { echo "ERROR: INPUT_VCF not found: $INPUT_VCF" >&2; exit 1; }
mkdir -p "$EVAL_RUN_DIR"

# Rootless podman on this host needs keep-id userns to write host-owned work dirs.
EXTRA_CFG_EVAL="${EVAL_RUN_DIR}/podman_userns.config"
cat > "$EXTRA_CFG_EVAL" <<'EOF'
podman.runOptions = '--userns=keep-id'
EOF

# Assemble optional args.
opt_args=()
[[ -n "${INPUT_VCF_TBI:-}" ]] && opt_args+=(--input_vcf_tbi "$INPUT_VCF_TBI")
[[ "${SKIP_PREP:-}" == "true" ]] && opt_args+=(--skip_preparation true)

echo "=================================================================="
echo " nf-eval-gene-assoc -- scoring ALL datasets in one fan-out run"
echo "   method  : ${EVAL_PROJECT}"
echo "   vcf     : ${INPUT_VCF}${INPUT_VCF_TBI:+ (+tbi)}${SKIP_PREP:+  [skip_preparation]}"
echo "   regenie : ${REGENIE_GLOB}"
echo "   causal  : ${CAUSAL_GENES_GLOB}"
echo "   outdir  : ${EVAL_RUN_DIR}/results"
echo "=================================================================="

cd "$EVAL_REPO"
nextflow run "$EVAL_REPO/main.nf" \
    -profile "$EVAL_PROFILE" \
    -c "$EXTRA_CFG_EVAL" \
    -work-dir "${EVAL_RUN_DIR}/work" \
    --project_name "${EVAL_PROJECT}_eval" \
    --input_vcf "${INPUT_VCF}" \
    --causal_snplist "${CAUSAL_SNPLIST_GLOB}" \
    --causal_genes "${CAUSAL_GENES_GLOB}" \
    --regenie_results "${REGENIE_GLOB}" \
    --outdir "${EVAL_RUN_DIR}/results" \
    "${opt_args[@]}"

echo ""
echo "eval (${EVAL_PROJECT}) finished: $(date -Is)"
