#!/usr/bin/env bash
#
# Benchmark run: genepi/nf-gwas (comparator) -- MULTI-DATASET
# ==========================================================
# nf-gwas has NO multi-phenotype-file input (phenotypes_columns are columns of a
# SINGLE file) and emits a single Y1.regenie, so -- unlike nf-rare-var-assoc, which
# fans out over a comma-separated phenotype list in one run -- nf-gwas must be
# LOOPED once per simulated dataset. Each iteration:
#   PHASE 1 (podman bcftools/plink2): produce the inputs nf-gwas cannot make itself
#           for THIS dataset -- a sample-filtered, multiallelic-split, ID-assigned
#           VCF (association input) + a PLINK1 bed/bim/fam (prediction input).
#   PHASE 2 (nextflow): run nf-gwas in gene-based (RVAS) mode, REUSING the
#           per-dataset mask files and PC covariates nf-rare-var-assoc produced.
#   Then the single Y1.regenie is renamed to embed _dataset_idx_<N>_ and copied to
#   a shared dir so a SINGLE nf-eval-gene-assoc fan-out run can score all datasets.
#
# WHY nf-gwas READS THE RAW VCF (not the prepared one nf-rare-var-assoc now uses):
#   The comparison's premise is that nf-gwas lacks normalisation + PL-dosage. We
#   give it only the minimal split + ID assignment it cannot do itself (needed for
#   mask reuse), starting from the RAW VCF, which has no DS field -- so plink2
#   `dosage=DS` silently falls back to HARD CALLS (verified). Feeding nf-gwas the
#   prepared (DS-bearing, CSQ-filtered) VCF would hand it nf-rare-var-assoc's
#   dosage + normalisation for free and erase the documented intrinsic difference.
#
# WHY THE SPLIT + ID ASSIGNMENT IS MANDATORY:
#   We reuse nf-rare-var-assoc's gene mask files, whose variant keys are
#   CHROM_POS_REF_ALT (chr-stripped) built AFTER a multiallelic split. plink2
#   --make-pgen preserves the VCF ID column, so we set those IDs here. Skipping the
#   split / ID assignment => REGENIE matches zero variants. SPLIT_MULTIALLELIC must
#   stay true; REMOVE_DUPS is a toggle.
#
# PREREQUISITE: run_nf_rare_var_assoc.sh has completed (it produces the per-dataset
#   mask files and covar tables, named tools_comparison_dataset_idx_<N>.*).
#
# SPACE: each dataset's prep/ and work/ dirs are deleted after its REGENIE result
#   is extracted (set CLEANUP=false to keep them).
#
# DATASET SELECTION: set DATASET_IDXS (space-separated) to override; default is all
#   run_<N> dirs found under $DATASETS_DIR.
#
set -euo pipefail

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

DATA="${DATA:-/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/tools_comparison}"
RVA_REPO="${RVA_REPO:-/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc}"
EVAL_REPO="${EVAL_REPO:-/data/git/doktorat_pw/wum_pims/nf-eval-gene-assoc}"
NFGWAS_REPO="${NFGWAS_REPO:-/data/git/doktorat_pw/nf_rare_var_assoc__tools_comparison/genepi/nf-gwas}"
DATASETS_DIR="${DATA}/datasets"

INPUT_VCF="${DATA}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz"

RVA_PROJECT="tools_comparison"          # project_name used by run_nf_rare_var_assoc.sh
PROJECT="tools_comparison_nfgwas"
RUN_DIR="${DATA}/runs/nf_gwas"
EVAL_RUN_DIR="${DATA}/runs/nf_gwas_eval"
REGENIE_OUT_DIR="${RUN_DIR}/regenie_per_dataset"   # renamed per-dataset REGENIE results
PROFILE="podman,medium_resources"

# Reused outputs from the nf-rare-var-assoc arm:
RVA_RESULTS="${DATA}/runs/nf_rare_var_assoc/results"
GENE_MASKS="${RVA_REPO}/assets/default.masks"      # static mask definition

# Prep behaviour:
SPLIT_MULTIALLELIC=true   # MUST be true for mask reuse. Splits -m -any.
REMOVE_DUPS=true          # bcftools norm --rm-dup exact (matches nf-rare-var-assoc).
THREADS=8
CLEANUP="${CLEANUP:-true}"   # set CLEANUP=false to keep per-dataset prep/work dirs

# Container images (podman).
BCFTOOLS_IMG="community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f"
PLINK2_IMG="quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
# nf-gwas container: code at $NFGWAS_REPO is v1.0.11; only v1.0.9 is pulled locally.
NFGWAS_CONTAINER="quay.io/genepi/nf-gwas:v1.0.9"

# ----------------------------------------------------------------------------
# Dataset selection
# ----------------------------------------------------------------------------
if [[ -z "${DATASET_IDXS:-}" ]]; then
    DATASET_IDXS="$(ls -1d "${DATASETS_DIR}"/run_*/ 2>/dev/null \
        | sed -E 's#.*/run_([0-9]+)/#\1#' | sort -n | tr '\n' ' ')"
fi
read -r -a IDXS <<< "$DATASET_IDXS"
[[ ${#IDXS[@]} -gt 0 ]] || { echo "ERROR: no datasets selected (DATASETS_DIR=$DATASETS_DIR)" >&2; exit 1; }

pheno_path() {
    echo "${DATASETS_DIR}/run_$1/gcta_simu/tuner_base_run_$1_dataset_idx_$1_gcta_simu.phenotype.txt"
}

# ----------------------------------------------------------------------------
# Sanity checks
# ----------------------------------------------------------------------------
for f in "$INPUT_VCF" "$NFGWAS_REPO/main.nf" "$EVAL_REPO/main.nf" "$GENE_MASKS"; do
    [[ -e "$f" ]] || { echo "ERROR: missing required path: $f" >&2; exit 1; }
done
[[ "$SPLIT_MULTIALLELIC" == "true" ]] || {
    echo "ERROR: SPLIT_MULTIALLELIC=false desyncs IDs from the reused masks -> 0 variants." >&2; exit 1; }
mkdir -p "$RUN_DIR" "$EVAL_RUN_DIR" "$REGENIE_OUT_DIR"

VCFBASE="$(basename "$INPUT_VCF")"

# ----------------------------------------------------------------------------
# Shared prep (dataset-independent): VCF sample list + chr rename map. Done once.
# ----------------------------------------------------------------------------
SHARED="${RUN_DIR}/prep_shared"
mkdir -p "$SHARED"
bcft_shared() { podman run --rm --userns=keep-id -v "$DATA":/indata:z -v "$SHARED":/w:z -w /w "$BCFTOOLS_IMG" "$@"; }

echo "[shared] VCF sample list + chr rename map ..."
bcft_shared bcftools query -l "/indata/${VCFBASE}" | sort -u > "${SHARED}/vcf_ids.txt"
bcft_shared bcftools view -h "/indata/${VCFBASE}" \
    | grep '^##contig' \
    | sed -E 's/.*ID=([^,>]+).*/\1/' \
    | awk '{n=$1; sub(/^chr/,"",n); print $1"\t"n}' > "${SHARED}/chr_map.txt"

NORM_OPTS="-m -any"
[[ "$REMOVE_DUPS" == "true" ]] && NORM_OPTS="${NORM_OPTS} --rm-dup exact"

# ============================================================================
# PER-DATASET LOOP
# ============================================================================
for idx in "${IDXS[@]}"; do
    PHENO="$(pheno_path "$idx")"
    [[ -e "$PHENO" ]] || { echo "ERROR: missing phenotype for run_${idx}: $PHENO" >&2; exit 1; }

    RVA_ID="${RVA_PROJECT}_dataset_idx_${idx}"
    ANNO_FILE="${RVA_RESULTS}/bcftools/${RVA_ID}.annotations"
    SETLIST_FILE="${RVA_RESULTS}/bcftools/${RVA_ID}.setlist"
    COVAR_SRC="${RVA_RESULTS}/merge/${RVA_ID}_sex_covar.txt"
    for v in "$ANNO_FILE:.annotations" "$SETLIST_FILE:.setlist" "$COVAR_SRC:_sex_covar.txt"; do
        p="${v%%:*}"; name="${v##*:}"
        [[ -e "$p" ]] || { echo "ERROR: reused ${name} for run_${idx} not found: $p" >&2
                           echo "       Run run_nf_rare_var_assoc.sh first (project ${RVA_PROJECT})." >&2; exit 1; }
    done

    PREP="${RUN_DIR}/prep/run_${idx}"
    RES="${RUN_DIR}/results/run_${idx}"
    WORK="${RUN_DIR}/work/run_${idx}"
    mkdir -p "$PREP"
    PREPPED_VCF="prepped.vcf.gz"   # inside $PREP

    plink2c() { podman run --rm --userns=keep-id -v "$PREP":/w:z "$PLINK2_IMG" "$@"; }

    echo "=================================================================="
    echo " nf-gwas run_${idx} -- PHASE 1: preparation"
    echo "   pheno : ${PHENO}"
    echo "   started : $(date -Is)"
    echo "=================================================================="

    # 1) keep-list = phenotype IIDs present in the VCF
    awk 'NR>1{print $2}' "$PHENO" | sort -u > "${PREP}/pheno_ids.txt"
    comm -12 "${PREP}/pheno_ids.txt" "${SHARED}/vcf_ids.txt" > "${PREP}/keep_samples.txt"
    echo "       kept samples: $(wc -l < "${PREP}/keep_samples.txt")"

    # 2) sample-filter -> split multiallelics [-> rm-dup exact] -> rename chrs ->
    #    assign CHROM_POS_REF_ALT IDs -> bgzip + index. (Raw VCF has no DS field.)
    #    chr_map lives in $SHARED; copy into $PREP so the one-container pipeline sees it.
    cp "${SHARED}/chr_map.txt" "${PREP}/chr_map.txt"
    echo "[prep] writing ${PREPPED_VCF} (bcftools norm ${NORM_OPTS}) ..."
    podman run --rm --userns=keep-id -v "$DATA":/indata:z -v "$PREP":/w:z -w /w "$BCFTOOLS_IMG" \
        sh -euc "
            bcftools view -S /w/keep_samples.txt --force-samples /indata/${VCFBASE} -Ou \
            | bcftools norm ${NORM_OPTS} --threads ${THREADS} -Ou \
            | bcftools annotate --rename-chrs /w/chr_map.txt -Ou \
            | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \
                  --threads ${THREADS} -Oz -o /w/${PREPPED_VCF}
            bcftools index -t /w/${PREPPED_VCF}
        "

    # 3) prediction bed/bim/fam. --update-sex needed because chrX present; we supply
    #    the sex nf-rare-var-assoc imputed (from this dataset's _sex_covar.txt).
    awk -F'\t' 'NR==1{sub(/^#/,"",$1); for(i=1;i<=NF;i++) col[$i]=i; next}
                {print $col["FID"]"\t"$col["IID"]"\t"$col["SEX"]}' \
        "$COVAR_SRC" > "${PREP}/sex_update.txt"
    echo "[prep] exporting prediction bed/bim/fam ..."
    plink2c plink2 \
        --vcf "/w/${PREPPED_VCF}" \
        --double-id --split-par b38 \
        --update-sex "/w/sex_update.txt" \
        --make-bed \
        --out /w/prediction \
        --threads "${THREADS}"

    # 4) PC covariate file (FID, IID, PC1_AVG..PC4_AVG; no SEX -- covarColList is PCs).
    COVAR_OUT="${PREP}/covar_pc.txt"
    awk -F'\t' '
        NR==1 { sub(/^#/,"",$1); for(i=1;i<=NF;i++) col[$i]=i;
                print "FID\tIID\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG"; next }
        { print $col["FID"]"\t"$col["IID"]"\t"$col["PC1_AVG"]"\t"$col["PC2_AVG"]"\t"$col["PC3_AVG"]"\t"$col["PC4_AVG"] }
    ' "$COVAR_SRC" > "$COVAR_OUT"

    # 5) nf-gwas params for this dataset
    PARAMS_FILE="${PREP}/params.nf_gwas.yaml"
    cat > "$PARAMS_FILE" <<EOF
project: "${PROJECT}_dataset_idx_${idx}"
outdir: "${RES}"

genotypes_prediction: "${PREP}/prediction.{bed,bim,fam}"
genotypes_association: "${PREP}/${PREPPED_VCF}"
genotypes_association_format: "vcf"
genotypes_build: "hg38"
association_build: "hg38"

phenotypes_filename: "${PHENO}"
phenotypes_columns: "Y1"
phenotypes_binary_trait: true

covariates_filename: "${COVAR_OUT}"
covariates_columns: "PC1_AVG,PC2_AVG,PC3_AVG,PC4_AVG"

# Single QC call (mapped from nf-rare-var-assoc makepgen_3 + snplist_qc mind)
qc_maf: 0.045
qc_mac: 16
qc_geno: 0.25
qc_hwe: "1e-9"
qc_mind: 0.15

# REGENIE
regenie_bsize_step1: 400
regenie_bsize_step2: 200
regenie_ref_first: true
regenie_firth: true
regenie_firth_approx: true
regenie_min_mac: 7

# Gene-based (RVAS) -- reuse nf-rare-var-assoc's per-dataset mask files
regenie_run_gene_based_tests: true
regenie_gene_anno: "${ANNO_FILE}"
regenie_gene_setlist: "${SETLIST_FILE}"
regenie_gene_masks: "${GENE_MASKS}"
regenie_gene_aaf: 0.2
regenie_gene_test: "skato"
EOF

    EXTRA_CFG="${PREP}/nfgwas_podman.config"
    cat > "$EXTRA_CFG" <<EOF
podman.enabled    = true
docker.enabled    = false
podman.runOptions = '--userns=keep-id'
podman.mountFlags = 'z'
process.container = '${NFGWAS_CONTAINER}'
EOF

    echo "=================================================================="
    echo " nf-gwas run_${idx} -- PHASE 2: nextflow"
    echo "=================================================================="
    cd "$NFGWAS_REPO"
    nextflow run "$NFGWAS_REPO/main.nf" \
        -c "$EXTRA_CFG" \
        -params-file "$PARAMS_FILE" \
        -work-dir "${WORK}" \
        -with-trace "${RES}/trace.txt" \
        -with-report "${RES}/report.html" \
        -with-timeline "${RES}/timeline.html"

    # 6) Extract + rename REGENIE result so dataset_idx is in the filename (the
    #    eval fan-out parses it). nf-eval reads sep=' ', comment='#'.
    OUT_REGENIE="${REGENIE_OUT_DIR}/${PROJECT}_dataset_idx_${idx}_step2_Y1.regenie"
    zcat "${RES}/results/Y1.regenie.gz" > "${REGENIE_OUT_DIR}/.tmp_Y1.regenie"
    echo "## ignored line" > "$OUT_REGENIE"
    awk '{ gsub(/\t/, " "); print }' "${REGENIE_OUT_DIR}/.tmp_Y1.regenie" >> "$OUT_REGENIE"
    rm -f "${REGENIE_OUT_DIR}/.tmp_Y1.regenie"
    echo "[run_${idx}] REGENIE result -> ${OUT_REGENIE}"

    # 7) Per-dataset cleanup (keep the extracted REGENIE; drop bulky prep/work).
    if [[ "$CLEANUP" == "true" ]]; then
        echo "[cleanup] rm -rf ${PREP} ${WORK}"
        rm -rf "$PREP" "$WORK"
    fi

    echo "[run_${idx}] done: $(date -Is)"
done

# ----------------------------------------------------------------------------
# Evaluation: ONE nf-eval-gene-assoc run over ALL datasets (glob inputs).
# ----------------------------------------------------------------------------
EXTRA_CFG_EVAL="${EVAL_RUN_DIR}/podman_userns.config"
cat > "$EXTRA_CFG_EVAL" <<'EOF'
podman.runOptions = '--userns=keep-id'
EOF

CAUSAL_SNPLIST_GLOB="${DATASETS_DIR}/run_*/select_genes/*_dataset_idx_*_in_*.snplist"
CAUSAL_GENES_GLOB="${DATASETS_DIR}/run_*/select_genes/*_genes_dataset_idx_*.txt"
REGENIE_GLOB="${REGENIE_OUT_DIR}/*_dataset_idx_*_step2_Y1.regenie"

echo ""
echo "=================================================================="
echo " nf-eval-gene-assoc -- scoring ALL nf-gwas datasets in one fan-out run"
echo "   regenie : ${REGENIE_GLOB}"
echo "   outdir  : ${EVAL_RUN_DIR}/results"
echo "=================================================================="

cd "$EVAL_REPO"
nextflow run "$EVAL_REPO/main.nf" \
    -profile "$PROFILE" \
    -c "$EXTRA_CFG_EVAL" \
    -work-dir "${EVAL_RUN_DIR}/work" \
    --project_name "${PROJECT}_eval" \
    --input_vcf "${INPUT_VCF}" \
    --input_vcf_tbi "${INPUT_VCF}.tbi" \
    --causal_snplist "${CAUSAL_SNPLIST_GLOB}" \
    --causal_genes "${CAUSAL_GENES_GLOB}" \
    --regenie_results "${REGENIE_GLOB}" \
    --outdir "${EVAL_RUN_DIR}/results"

echo ""
echo "nf-gwas evaluation Finished: $(date -Is)"
