#!/usr/bin/env bash
#
# Runs RICOPILI (quality control + principal components) combined with STAARpipeline
# (association) as one method, for comparison against nf-rare-var-assoc.
# ================================================================================
# This method receives only the raw VCF, the phenotype files, and the causal genes
# used for scoring. Quality control, population structure, annotation and the
# statistical test are all done by RICOPILI and STAARpipeline themselves.
#
# See README.md in this directory for prerequisites, the one-time setup (containers
# and the FAVOR annotation database), and the practical problems worth knowing about.
#
# Runs per dataset rather than once: RICOPILI's quality control depends on the
# case/control split, and every dataset has its own sample subset.
#
# TWO INDEPENDENT AXES produce the ARMS this script runs. An arm is labelled
# "<version>_<spa mode>" and gets its own result table and its own eval.
#
#   VERSIONS -- how relatedness is handled:
#     full     = all samples kept. Relatedness is modelled in the null model with a
#                GENESIS PC-Relate matrix, and principal components come from GENESIS
#                PC-AiR. This is how STAARpipeline is designed to be used, and it
#                matches nf-rare-var-assoc's own handling of relatedness. DEFAULT.
#     filtered = RICOPILI as shipped. pcaer removes related samples; its smartpca
#                principal components carry into the association test. Kept runnable
#                (VERSIONS="filtered full") but no longer run by default: it scored
#                worse, for the marker-starvation reason documented in README.md.
#
#   SPA_MODES -- which test STAAR runs. These are DIFFERENT TESTS, not a calibration
#                tweak, and both are run so the better one can be chosen:
#     nospa    = use_SPA=FALSE -> STAAR-O, the full omnibus (SKAT + Burden + ACAT-V).
#                The like-for-like analogue of REGENIE's SKAT-O.
#     spa      = use_SPA=TRUE  -> STAAR-B, a BURDEN-ONLY omnibus. STAARpipeline skips
#                the SKAT and ACAT-V branches entirely under SPA (coding.R:502-516).
#                The saddlepoint approximation is what these datasets' extreme
#                case/control imbalance (as few as 11 cases) calls for.
#     The p-value column differs per mode and is handed to staar_to_eval.py --p-col;
#     getting that wrong is a hard failure at stage G, not a silent one.
#
# STAGES PER DATASET:
#   A raw VCF -> sample subset -> genotypes only -> split multi-allelic -> PLINK bed
#   B preimp_dir  (sample and variant quality control)
#   C pcaer       (relatedness removal + smartpca components; Filtered version)
#   D GENESIS PC-AiR + PC-Relate matrix (Full version)  genesis_pcair_pcrelate.R
#   E per chromosome: bed -> VCF -> GDS -> annotated GDS (FAVORannotator)
#   F association test per arm and chromosome           staar_gene_centric_coding.R
#   G staar_to_eval.py -> one REGENIE-shaped result table per arm
# Stages A-E are shared by every arm; only F and G repeat per arm.
# With SCORE=true, benchmark-common/run_eval.sh then scores each arm.
#
# RESUME: a dataset whose result tables all exist already is skipped untouched, and
#   one dataset's failure is recorded and skipped rather than aborting the run. A
#   re-launch after an interruption picks up where it stopped.
#
# THIS SCRIPT DELETES NOTHING by default. It refuses to start a dataset whose work
#   dir is dirty (RICOPILI will not re-run in one) and tells you the path to remove.
#   Set CLEANUP=true to have it drop each dataset's work dir once that dataset's
#   tables are safely written.
#
# QUALITY-CONTROL THRESHOLDS are set explicitly rather than left at RICOPILI's
#   defaults, which were designed for genotyping arrays and remove so much from exome
#   data that some datasets produce no result at all. The four PREIMP_* values below
#   are preimp_dir's own command-line options -- nothing is patched, and every other
#   threshold is left as shipped. Set PREIMP_GENO/MIND/PRE_GENO/MIDI to
#   0.02/0.02/0.05/0.02 to reproduce a run with RICOPILI's defaults. The reasoning,
#   including why raising --geno alone has no effect, is in README.md.
#
# CHROMOSOMES: CHRS defaults to "12 22". The FAVOR database for each chromosome must
#   be extracted beforehand (about 196 GB for chromosome 12, 53 GB for chromosome 22);
#   the check below fails loudly if one is missing. Set CHRS="22" for chromosome 22 only.
#
# DATASETS: set DATASET_IDXS (space-separated); the default is every run_<N> found.
#
# QUICK CHECK: MAX_GENES=N limits how many genes are tested per chromosome, turning a
#   ~45 minute run into about a minute while still exercising every stage. The
#   resulting tables CANNOT be scored -- leave MAX_GENES unset for any real run.
#
set -euo pipefail

# ----------------------------------------------------------------------------
# Configuration
#
# The location variables default to one particular machine's layout but can all be
# set from the environment, so the same script runs elsewhere without editing.
# Setting DATA and RVA_REPO is usually enough; ARM, COMMON and REFDIR only need
# setting if the repository or the reference genome is somewhere unusual. All paths
# are passed into containers, so they must be absolute.
# ----------------------------------------------------------------------------
DATA="${DATA:-/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/tools_comparison}"
RVA_REPO="${RVA_REPO:-/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc}"
DATASETS_DIR="${DATASETS_DIR:-${DATA}/datasets}"
ARM="${ARM:-${RVA_REPO}/docs/tool-comparison/chained-benchmark-staar}"
COMMON="${COMMON:-${RVA_REPO}/docs/tool-comparison/benchmark-common}"

INPUT_VCF="${INPUT_VCF:-${DATA}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz}"
PED="${PED:-${RVA_REPO}/assets/integrated_call_samples_v3.20250704.ALL.ped}"   # 1000G sex source
REFDIR="${REFDIR:-/data/doktorat/biodatageeks/genome_in_a_bottle/reference}"
REF="${REF:-GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz}"  # chr-prefixed
FAVOR_DB="${FAVOR_DB:-${DATA}/favor/db}"          # chr<N>_*.csv(.idx), extracted per chr
GENES_INFO="${GENES_INFO:-${ARM}/genes_info_hgnc.tsv}"

# Run directory. Dated, so every re-run lands beside its predecessors instead of
# overwriting them -- nothing on disk is ever replaced by a later run. METHOD is
# derived from the directory name and is what every result table, eval project and
# pairwise arm is named after, so pointing RUN_DIR somewhere else renames all of them
# consistently and cannot produce a mislabelled table.
RUN_DATE="${RUN_DATE:-2026_09_12}"
RUN_DIR="${RUN_DIR:-${DATA}/runs/ricopili_staar_${RUN_DATE}}"
METHOD="$(basename "$RUN_DIR")"
REGENIE_OUT_DIR="${RUN_DIR}/regenie_per_dataset"  # per-arm subdirs (<version>_<spa mode>)

# RICOPILI quality-control thresholds. All four are preimp_dir's own command-line
# options -- nothing is patched -- and they are the only thresholds moved off
# RICOPILI's defaults. Hardy-Weinberg (1e-6 controls / 1e-10 cases), the inbreeding
# threshold (0.2) and the frequency floor (0) are left as shipped.
#   geno 0.10 : default 0.02. Matches nf-rare-var-assoc's --geno 0.100 on its
#               association set.
#   mind 0.20 : default 0.02. Matches nf-rare-var-assoc's --mind 0.200.
#   pre_geno  : default 0.05. RAISING --geno ALONE HAS NO EFFECT. The pre-filter list
#     0.35      (variants whose missingness exceeds pre_geno, measured before any
#               sample filtering) is excluded from every later step and never returns,
#               so the effective cut is min(pre_geno, geno). Left at 0.05 it would
#               silently bind and undo --geno 0.10. RICOPILI's own defaults keep
#               pre_geno looser than geno for exactly this reason: the pre-filter
#               guards against catastrophically bad variants inflating sample
#               missingness, it is not meant to be the binding cut.
#   midi 0.35 : default 0.02, and not an exact equivalent of anything in
#               nf-rare-var-assoc. RICOPILI limits the case/control missingness
#               DIFFERENCE; nf-rare-var-assoc limits each group's ABSOLUTE rate
#               (--geno 0.350 within each phenotype group, intersected). 0.35 is the
#               largest difference that filter can let through, so it is the loosest
#               setting still implied by this pipeline rather than an invented one.
#               The default 0.02 is severe here: with about 11 cases in the most
#               imbalanced datasets, one missing genotype is a 9% difference.
# pre_geno and midi both land on 0.35 because one nf-rare-var-assoc step does both
# jobs at once: a loose absolute pre-pass and an implied cap on the difference.
# INVARIANT: keep PREIMP_PRE_GENO >= PREIMP_GENO, or --geno stops being the real cut.
PREIMP_GENO="${PREIMP_GENO:-0.10}"         # variant missingness (preimp_dir --geno)
PREIMP_MIND="${PREIMP_MIND:-0.20}"         # sample missingness  (preimp_dir --mind)
PREIMP_PRE_GENO="${PREIMP_PRE_GENO:-0.35}" # pre-filter variant missingness (--pre_geno)
PREIMP_MIDI="${PREIMP_MIDI:-0.35}"         # case/control differential missingness (--midi)
PREIMP_EXTRA_ARGS="${PREIMP_EXTRA_ARGS:-}" # anything else to hand preimp_dir verbatim
PREIMP_QC_ARGS="--geno ${PREIMP_GENO} --mind ${PREIMP_MIND} --pre_geno ${PREIMP_PRE_GENO} --midi ${PREIMP_MIDI} ${PREIMP_EXTRA_ARGS}"
awk -v p="$PREIMP_PRE_GENO" -v g="$PREIMP_GENO" 'BEGIN{exit !(p+0 >= g+0)}' || {
    echo "ERROR: PREIMP_PRE_GENO (${PREIMP_PRE_GENO}) < PREIMP_GENO (${PREIMP_GENO}):" >&2
    echo "       the pre-filter would bind instead of --geno. Raise PREIMP_PRE_GENO." >&2
    exit 1; }

# Run knobs.
CHRS="${CHRS:-12 22}"              # chromosomes to test; each needs its FAVOR database extracted
VERSIONS="${VERSIONS:-full}"       # relatedness axis; "filtered full" to run both
SPA_MODES="${SPA_MODES:-nospa spa}"  # test axis; see the header. "nospa" or "spa" alone works
NPCS_COVAR="${NPCS_COVAR:-4}"      # principal components used as covariates
# STAAR's rare_maf_cutoff: the frequency ceiling below which a variant enters the
# gene aggregate. HARMONIZED WITH THE REFERENCE, which is the entire point of this
# value -- nf-rare-var-assoc aggregates at regenie's --aaf-bins 0.1 and the nf-gwas
# arm at regenie_gene_aaf 0.1, so anything else here means the arms are not testing
# the same variant sets. (Every run before 2026-09-11 silently used STAAR's own
# default of 0.01, because this variable existed but was never passed.)
# CAVEAT: regenie's bin is on ALT allele frequency and STAAR's cutoff is on MINOR
# allele frequency. With --ref-first and the REF restoration stage E performs these
# agree for everything rare; they part company only where ALT is the major allele,
# which at a 0.1 ceiling is a rounding-error set.
RARE_MAF="${RARE_MAF:-0.1}"        # STAAR rare_maf_cutoff (--rare-maf)
THREADS="${THREADS:-4}"
MAX_GENES="${MAX_GENES:-}"         # smoke-test cap on genes tested per chr; empty = all genes
# This script deletes nothing unless you ask it to. With CLEANUP=false (the default)
# it prints the size of each dataset's work dir and the exact command to remove it;
# the result tables, the retention TSVs and the logs all live outside it.
CLEANUP="${CLEANUP:-false}"        # true -> drop each dataset's work dir once its tables exist
SCORE="${SCORE:-false}"            # true -> run run_eval.sh per arm at the tail

# Container images.
BCFTOOLS_IMG="docker.io/psuszynski/bioinf_combo:1.5.1"     # bcftools + pandas (staar_to_eval)
RICOPILI_IMG="localhost/ricopili:2025_Feb_20.001"
FAVOR_IMG="localhost/favorannotator:1.0.0"
STAAR_IMG="docker.io/zilinli/staarpipeline:0.9.7"
GENESIS_IMG="${GENESIS_IMG:-docker.io/uwgac/topmed-roybranch:latest}"  # SNPRelate+GWASTools+GENESIS

# ----------------------------------------------------------------------------
# Dataset selection
# ----------------------------------------------------------------------------
if [[ -z "${DATASET_IDXS:-}" ]]; then
    DATASET_IDXS="$(ls -1d "${DATASETS_DIR}"/run_*/ 2>/dev/null \
        | sed -E 's#.*/run_([0-9]+)/#\1#' | sort -n | tr '\n' ' ')"
fi
read -r -a IDXS <<< "$DATASET_IDXS"
[[ ${#IDXS[@]} -gt 0 ]] || { echo "ERROR: no datasets selected (DATASETS_DIR=$DATASETS_DIR)" >&2; exit 1; }
read -r -a CHR_ARR <<< "$CHRS"
read -r -a VER_ARR <<< "$VERSIONS"
read -r -a SPA_ARR <<< "$SPA_MODES"

# The arms: the cross product of the relatedness axis and the test axis. An arm label
# is "<version>_<spa mode>" and is the name of its result subdir, its result tables,
# its eval project and its pairwise arm -- one string, so those can never disagree.
# Each mode also fixes which p-value column staar_to_eval.py must read: STAAR omits
# STAAR-O entirely under SPA and STAAR-B entirely without it, so reading the wrong one
# is a hard failure at stage G rather than a quietly empty table.
ARM_LABELS=()
declare -A ARM_VER ARM_SPA_FLAG ARM_P_COL
for v in "${VER_ARR[@]}"; do
    for s in "${SPA_ARR[@]}"; do
        case "$s" in
            nospa) flag=FALSE; pcol="STAAR-O" ;;
            spa)   flag=TRUE;  pcol="STAAR-B" ;;
            *) echo "ERROR: SPA_MODES must contain only 'nospa' and/or 'spa', got '${s}'" >&2
               exit 1 ;;
        esac
        label="${v}_${s}"
        ARM_LABELS+=("$label")
        ARM_VER["$label"]="$v"
        ARM_SPA_FLAG["$label"]="$flag"
        ARM_P_COL["$label"]="$pcol"
    done
done

# Smoke-test cap: --max-genes tests only the first N genes of a chromosome, which
# turns a ~45 min STAAR call into ~1 min. NEVER set it for a scored run.
MAXG_ARG=()
if [[ -n "$MAX_GENES" ]]; then
    MAXG_ARG=(--max-genes "$MAX_GENES")
    echo "WARNING: MAX_GENES=${MAX_GENES} -- SMOKE TEST ONLY, results are not scoreable." >&2
fi

pheno_path() {
    echo "${DATASETS_DIR}/run_$1/gcta_simu/tuner_base_run_$1_dataset_idx_$1_gcta_simu.phenotype.txt"
}

# ----------------------------------------------------------------------------
# Sanity checks
# ----------------------------------------------------------------------------
for f in "$INPUT_VCF" "$PED" "${REFDIR}/${REF}" "$GENES_INFO" \
         "${ARM}/staar_gene_centric_coding.R" "${ARM}/staar_to_eval.py" \
         "${ARM}/convert_vcf_to_gds.R" "${ARM}/favorannotator_csv_essential.R" \
         "${ARM}/genesis_pcair_pcrelate.R"; do
    [[ -e "$f" ]] || { echo "ERROR: missing required path: $f" >&2; exit 1; }
done
for chr in "${CHR_ARR[@]}"; do
    ls "${FAVOR_DB}/chr${chr}_"*.csv >/dev/null 2>&1 || {
        echo "ERROR: FAVOR DB for chr${chr} not extracted under ${FAVOR_DB}" >&2
        echo "       (chr12 needs ~196 GB; see README next-action #7). Set CHRS to omit it," >&2
        echo "       or extract it first." >&2; exit 1; }
done
mkdir -p "$RUN_DIR" "$REGENIE_OUT_DIR"
for a in "${ARM_LABELS[@]}"; do mkdir -p "${REGENIE_OUT_DIR}/${a}"; done

# The result table one arm produces for one dataset. Single definition, used by the
# resume check, by stage G and by the scoring tail, so those three cannot drift apart.
arm_table() {   # $1 = arm label, $2 = dataset idx
    echo "${REGENIE_OUT_DIR}/${1}/${METHOD}_${1}_dataset_idx_${2}_step2_Y1.regenie"
}

VCFBASE="$(basename "$INPUT_VCF")"

# Container helpers. Each mounts what its stage needs; :z relabels, --userns=keep-id
# lets rootless podman write host-owned dirs.
bcft()  { podman run --rm --userns=keep-id "$@"; }
# The GENESIS image ships USER=topmed (uid 2049); under --userns=keep-id that lands on
# a subuid with no write access to our host-owned work dir, so its first output file
# fails with "Permission denied". Run it as our own uid instead.
genesis() { podman run --rm --userns=keep-id --user "$(id -u):$(id -g)" "$@"; }

# ============================================================================
# PER-DATASET LOOP
#
# Each dataset runs in its own subshell so ONE dataset's failure does not abort the
# whole benchmark -- at ~2 h per dataset a single bad draw must not cost the night.
# The subshell is a standalone command (not inside an `if`/`||`), so `set -e` stays
# ACTIVE within it and the dataset still stops at its first real error; only the outer
# `set +e` keeps that from killing the loop. Failures are collected and reported.
#
# RESUME: a dataset whose tables for EVERY selected arm already exist is skipped
# untouched, so a re-launch after an interruption picks up where it stopped without
# redoing finished datasets. A dataset with only some arms done is redone in full --
# stages A-E are shared, and re-deriving them is what makes the arms comparable.
# ============================================================================
FAILED_IDXS=()
for idx in "${IDXS[@]}"; do
    # ---- resume: every selected arm already has its table for this dataset? -------
    N_HAVE=0
    for a in "${ARM_LABELS[@]}"; do
        [[ -e "$(arm_table "$a" "$idx")" ]] && N_HAVE=$((N_HAVE + 1))
    done
    if [[ "$N_HAVE" -eq "${#ARM_LABELS[@]}" ]]; then
        echo "[run_${idx}] all ${#ARM_LABELS[@]} arm table(s) present -- skipping"
        continue
    fi
    if [[ "$N_HAVE" -gt 0 ]]; then
        echo "[run_${idx}] ${N_HAVE}/${#ARM_LABELS[@]} arm tables present -- redoing the dataset in full" >&2
    fi

    set +e
    ( set -e
    PHENO="$(pheno_path "$idx")"
    [[ -e "$PHENO" ]] || { echo "ERROR: missing phenotype for run_${idx}: $PHENO" >&2; exit 1; }
    STUDY="run$(printf '%02d' "$idx")"          # 5-char RICOPILI study name (idx<=99)

    # RICOPILI records progress and refuses to repeat a step that "has been done
    # repeatedly without any progress", so a retry has to start from a clean
    # directory. This script will NOT delete one for you -- it stops and says so.
    WD="${RUN_DIR}/work/run_${idx}"
    if [[ -d "$WD" ]] && [[ -n "$(ls -A "$WD" 2>/dev/null)" ]]; then
        echo "ERROR: work dir for run_${idx} already exists and is not empty:" >&2
        echo "         $WD" >&2
        echo "       RICOPILI will not re-run in a dirty directory. Remove it yourself:" >&2
        echo "         rm -rf '$WD'" >&2
        exit 1
    fi
    mkdir -p "$WD"
    RET="${RUN_DIR}/retention/run_${idx}.tsv"; mkdir -p "$(dirname "$RET")"
    : > "$RET"
    # Provenance: the settings that produced these counts, so a table on disk says
    # which configuration it came from rather than relying on the directory name.
    printf 'preimp_qc_args\t%s\n' "$PREIMP_QC_ARGS" >> "$RET"
    printf 'rare_maf\t%s\n'       "$RARE_MAF"       >> "$RET"
    printf 'variant_type\t%s\n'   "SNV"             >> "$RET"
    printf 'arms\t%s\n'           "${ARM_LABELS[*]}" >> "$RET"
    printf 'npcs_covar\t%s\n'     "$NPCS_COVAR"     >> "$RET"

    echo "=================================================================="
    echo " chain run_${idx}  (study ${STUDY}, chrs ${CHRS})"
    echo "   pheno   : ${PHENO}"
    echo "   QC      : ${PREIMP_QC_ARGS}"
    echo "   arms    : ${ARM_LABELS[*]}"
    echo "   rare_maf: ${RARE_MAF}"
    echo "   out     : ${RUN_DIR}"
    echo "   started : $(date -Is)"
    echo "=================================================================="

    # ------------------------------------------------------------------ Stage A
    # keep-list, sex (1000G ped), plink pheno; raw VCF -> subset -> GT-only ->
    # split -> bed. Two plink passes (plink2 converts, plink1.9 re-sorts the
    # --split-par output that plink1.9 otherwise refuses -- README gotcha).
    echo "[A] RICOPILI input prep ..."
    awk 'NR>1{print $2}' "$PHENO" | sort -u > "${WD}/keep_iids.txt"
    awk -F'\t' 'NR>1{print $2"\t"$5}' "$PED" > "${WD}/sex_all.tsv"
    { printf '#FID\tIID\tSEX\n'
      awk 'NR==FNR{s[$1]=$2; next} ($1 in s){print $1"\t"$1"\t"s[$1]}' \
          "${WD}/sex_all.tsv" "${WD}/keep_iids.txt"; } > "${WD}/sex.txt"
    { printf '#FID\tIID\tY1\n'
      awk 'NR>1{print $1"\t"$2"\t"$3}' "$PHENO"; } > "${WD}/pheno_plink.txt"

    bcft -v "$DATA":/d:z,ro -v "$WD":/w:z "$BCFTOOLS_IMG" bash -lc "
        bcftools view -S /w/keep_iids.txt --force-samples -Ou /d/${VCFBASE} \
        | bcftools annotate -x ^FORMAT/GT -Ou \
        | bcftools norm -m -any --threads ${THREADS} -Oz -o /w/rico.vcf.gz
        bcftools index -t /w/rico.vcf.gz"

    # VARID template passed via env so the in-container bash never expands \$r/\$a.
    # --keep-allele-order on the plink 1.9 re-sort: PLINK 1.9 otherwise makes A1 the
    # minor allele, which loses which allele was REF (stage E has to restore it).
    bcft -v "$WD":/w:z -e "VARID=@:#:\$r:\$a" "$RICOPILI_IMG" bash -lc '
        /opt/rp_dep/plink2/plink2 --vcf /w/rico.vcf.gz --double-id \
          --set-all-var-ids "$VARID" --new-id-max-allele-len 1000 \
          --update-sex /w/sex.txt --split-par hg38 --output-chr 26 \
          --pheno /w/pheno_plink.txt --pheno-name Y1 --1 \
          --make-bed --out /w/conv
        /opt/rp_dep/plink/plink --bfile /w/conv --make-bed --allow-no-sex \
          --keep-allele-order --out /w/sorted'
    [[ "$(awk '{print $2}' "${WD}/sorted.bim" | sort | uniq -d | wc -l)" -eq 0 ]] \
        || { echo "ERROR: duplicate variant IDs in sorted.bim (run_${idx})" >&2; exit 1; }
    printf 'n_variants_converted\t%s\n' "$(wc -l < "${WD}/sorted.bim")" >> "$RET"
    printf 'n_samples_converted\t%s\n'  "$(wc -l < "${WD}/sorted.fam")" >> "$RET"

    # ------------------------------------------------------------------ Stage B
    # preimp_dir runs TWICE: the first call writes a template name file and exits,
    # we set the 5-char study name, the second call does the QC (README gotcha).
    # PREIMP_QC_ARGS goes to BOTH calls: preimp_dir re-reads its whole command line on
    # the second pass, so omitting them there would silently run the QC at defaults.
    echo "[B] preimp_dir (sample + variant QC) ..."
    echo "    QC thresholds: ${PREIMP_QC_ARGS}"
    mkdir -p "${WD}/preimp"
    for e in bed bim fam; do cp "${WD}/sorted.${e}" "${WD}/preimp/${STUDY}_raw.${e}"; done
    bcft -v "$WD":/w:z "$RICOPILI_IMG" bash -lc \
        "cd /w/preimp && preimp_dir --disease sim --outname ${STUDY} --popname mix --serial ${PREIMP_QC_ARGS}" || true
    # first column of the matching line (placeholder 'sim1') -> the study name
    sed -i "s/^sim1\t${STUDY}_raw/${STUDY}\t${STUDY}_raw/" "${WD}/preimp/sim.names"
    bcft -v "$WD":/w:z "$RICOPILI_IMG" bash -lc \
        "cd /w/preimp && preimp_dir --disease sim --outname ${STUDY} --popname mix --serial ${PREIMP_QC_ARGS}" \
        2>&1 | tee "${WD}/preimp/preimp.log"

    QCPRE="${WD}/preimp/qc/sim_${STUDY}_mix_rp-qc1"
    for e in bed bim fam; do
        [[ -e "${QCPRE}.${e}" ]] || { echo "ERROR: preimp_dir produced no ${QCPRE}.${e}" >&2; exit 1; }
    done
    # Per-step retention counts, taken from preimp_dir's own .meta file.
    if [[ -e "${QCPRE}.meta" ]]; then
        grep -E '^(nsnpex_(mono|prefilter|hwe-co|miss|midi|hwe-ca|prekno)|nidex_(fhet|miss|sexcheck_ex))' \
            "${QCPRE}.meta" | tr -s ' ' '\t' >> "$RET" || true
    fi
    N_VAR_QC="$(wc -l < "${QCPRE}.bim")"
    printf 'n_variants_after_qc\t%s\n' "$N_VAR_QC" >> "$RET"
    printf 'n_samples_after_qc\t%s\n'  "$(wc -l < "${QCPRE}.fam")" >> "$RET"
    # Running out of variants is the failure mode this method hits silently: too few
    # survive, STAAR finds nothing testable, and the dataset returns an empty table
    # several stages later. Say so here instead.
    [[ "$N_VAR_QC" -ge 10000 ]] || \
        echo "WARNING: run_${idx} has only ${N_VAR_QC} variants after QC -- expect few or no testable genes." >&2

    # --------------------------------------------------------- Stage C (Filtered)
    # pcaer: IBD relatedness removal + smartpca PCs on the 454 unrelated survivors.
    if [[ " ${VERSIONS} " == *" filtered "* ]]; then
        echo "[C] pcaer (Filtered: relatedness removal + smartpca PCs) ..."
        mkdir -p "${WD}/pca"
        for e in bed bim fam; do cp "${QCPRE}.${e}" "${WD}/pca/"; done
        bcft -v "$WD":/w:z "$RICOPILI_IMG" bash -lc \
            "cd /w/pca && pcaer --out ${STUDY}pca --serial sim_${STUDY}_mix_rp-qc1" \
            2>&1 | tee "${WD}/pca/pcaer.log"
        MDS="${WD}/pca/pcaer_${STUDY}pca/${STUDY}pca.menv.mds"
        [[ -e "$MDS" ]] || { echo "ERROR: pcaer produced no ${MDS}" >&2; exit 1; }
        printf 'n_samples_pca_unrelated\t%s\n' "$(($(wc -l < "$MDS") - 1))" >> "$RET"
    fi

    # ------------------------------------------------------------- Stage D (Full)
    # GENESIS PC-AiR + PC-Relate over the whole quality-controlled bed (chr12+22+X).
    if [[ " ${VERSIONS} " == *" full "* ]]; then
        echo "[D] GENESIS PC-AiR + PC-Relate (Full: all samples kept) ..."
        mkdir -p "${WD}/genesis"
        for e in bed bim fam; do cp "${QCPRE}.${e}" "${WD}/genesis/qc.${e}"; done
        genesis -v "$WD":/w:z -v "$ARM":/arm:z,ro "$GENESIS_IMG" \
            Rscript /arm/genesis_pcair_pcrelate.R \
              --bed /w/genesis/qc --out-dir /w/genesis --npcs 20
        [[ -e "${WD}/genesis/pcair.tsv" && -e "${WD}/genesis/grm.rds" ]] \
            || { echo "ERROR: GENESIS step produced no pcair.tsv / grm.rds (run_${idx})" >&2; exit 1; }
        printf 'n_samples_full\t%s\n' "$(($(wc -l < "${WD}/genesis/pcair.tsv") - 1))" >> "$RET"
    fi

    # ------------------------------------------------------------------ Stage E
    # Per chr: QC'd bed -> chr-prefixed VCF -> FAVOR preprocess (left-normalize) ->
    # GDS -> aGDS. aGDS is per-dataset (the QC'd variant set is dataset-specific).
    for chr in "${CHR_ARR[@]}"; do
        echo "[E] chr${chr}: QC'd bed -> aGDS ..."
        AG="${WD}/c${chr}.agds"
        # Restore the REF allele before exporting. A PLINK bed does not record which
        # allele was REF, and preimp_dir's own PLINK 1.9 passes make A1 the minor
        # allele, so ~1.4% of sites come out of QC with REF and ALT exchanged and
        # bcftools norm -f then aborts on them. The variant IDs carry the answer --
        # stage A set them to <chr>:<pos>:<REF>:<ALT> -- so field 3 of the ID is the
        # true REF and --ref-allele restores it exactly, with no reference-genome
        # guessing. (Field 4 can be '*': a spanning deletion from splitting
        # multiallelics. Those parse the same way.)
        awk -F'\t' -v chr="$chr" '$1==chr {n=split($2,a,":"); if (n==4 && a[3]!="") print $2"\t"a[3]}' \
            "${QCPRE}.bim" > "${WD}/ref_alleles_chr${chr}.txt"
        # Export just this chromosome, chr-prefixed (FAVOR expects chr22, not 22).
        # id-paste=iid is REQUIRED: plink2 names VCF samples FID_IID by default, and
        # preimp_dir rewrites FID to 'con_sim_<study>_mix_rp_*<IID>', so the default
        # would give sample names no phenotype or PC file can join to.
        bcft -v "$WD":/w:z "$RICOPILI_IMG" bash -lc "
            /opt/rp_dep/plink2/plink2 --bfile /w/preimp/qc/sim_${STUDY}_mix_rp-qc1 \
              --chr ${chr} --ref-allele force /w/ref_alleles_chr${chr}.txt 2 1 \
              --output-chr chrM --export vcf bgz id-paste=iid \
              --out /w/c${chr}.export"
        # FAVORannotator's own preprocessing recipe, minus its GT-only step: a
        # plink-exported VCF already carries GT and nothing else, and bcftools fails
        # ("No matching tag in -x ^FORMAT/GT") when the strip would remove nothing.
        # The trailing filter step marks every record PASS. It is REQUIRED, not
        # cosmetic: STAAR reads QC_label ("annotation/filter") and keeps only
        # variants whose FILTER is exactly PASS ([R/coding.R] SNVlist <- filter ==
        # "PASS" & isSNV). PLINK writes FILTER as '.', so without this every gene
        # tests zero variants and the whole run comes back empty. PASS is also the
        # honest label here -- these are precisely the variants that survived
        # RICOPILI's QC, which is what the column means in this chain.
        bcft -v "$WD":/w:z -v "$REFDIR":/ref:z,ro "$BCFTOOLS_IMG" bash -lc "
            set -e
            bcftools norm -m -any --threads ${THREADS} -Ou /w/c${chr}.export.vcf.gz \
            | bcftools norm -f /ref/${REF} --threads ${THREADS} -Ou \
            | bcftools filter -e 'N_ALT<0' --threads ${THREADS} -Oz -o /w/c${chr}.pp.vcf.gz
            bcftools index -t /w/c${chr}.pp.vcf.gz"
        # VCF -> GDS (in place: the annotator modifies this file to become the aGDS)
        bcft -v "$WD":/w:z -v "$ARM":/arm:z,ro "$FAVOR_IMG" \
            Rscript /arm/convert_vcf_to_gds.R "/w/c${chr}.pp.vcf.gz" "/w/c${chr}.agds" "${THREADS}"
        # GDS -> aGDS (annotate against the chr's FAVOR DB)
        bcft -v "$WD":/w:z -v "$FAVOR_DB":/favordb:z,ro -v "$ARM":/arm:z,ro "$FAVOR_IMG" \
            Rscript /arm/favorannotator_csv_essential.R \
              "/w/c${chr}.agds" "${chr}" /favordb "/w/anno_out/"
        [[ -e "$AG" ]] || { echo "ERROR: no aGDS for chr${chr} (run_${idx})" >&2; exit 1; }
    done

    # ---------------------------------------------------------- Stages F + G
    # The augmented phenotype (Y1 + PCs) depends only on the RELATEDNESS version, so
    # it is built once per version and shared by that version's SPA modes -- the two
    # modes must see byte-identical covariates or they are not comparable.
    # Then, per arm: STAAR per chromosome, merge the category CSVs into one eval
    # table, and measure the inflation factor.
    for ver in "${VER_ARR[@]}"; do
        echo "[F/G] version=${ver}: augmented phenotype ..."
        AUG="${WD}/pheno_${ver}.tsv"
        GRM_ARG=()
        if [[ "$ver" == "filtered" ]]; then
            # menv.mds: FID IID SOL C1..C20  -> join Y1 on IID, take C1..C<NPCS_COVAR>.
            MDS="${WD}/pca/pcaer_${STUDY}pca/${STUDY}pca.menv.mds"
            awk -v n="$NPCS_COVAR" '
                NR==FNR { if (FNR>1) y[$2]=$3; next }               # pheno: IID->Y1
                FNR==1  { printf "FID\tIID\tY1"; for(i=1;i<=n;i++) printf "\tPC%d", i; printf "\n"; next }
                ($2 in y) { printf "%s\t%s\t%s", $2, $2, y[$2]
                            for(i=1;i<=n;i++) printf "\t%s", $(3+i); printf "\n" }
            ' "$PHENO" "$MDS" > "$AUG"
        else
            # pcair.tsv: IID PC1..PC20 -> join Y1 on IID, take PC1..PC<NPCS_COVAR>.
            PCAIR="${WD}/genesis/pcair.tsv"
            awk -v n="$NPCS_COVAR" '
                NR==FNR { if (FNR>1) y[$2]=$3; next }               # pheno: IID->Y1
                FNR==1  { printf "FID\tIID\tY1"; for(i=1;i<=n;i++) printf "\tPC%d", i; printf "\n"; next }
                ($1 in y) { printf "%s\t%s\t%s", $1, $1, y[$1]
                            for(i=1;i<=n;i++) printf "\t%s", $(1+i); printf "\n" }
            ' "$PHENO" "$PCAIR" > "$AUG"
            GRM_ARG=(--grm /w/genesis/grm.rds)
        fi
        printf 'n_samples_staar_%s\t%s\n' "$ver" "$(($(wc -l < "$AUG") - 1))" >> "$RET"

        for spa in "${SPA_ARR[@]}"; do
            label="${ver}_${spa}"
            SPA_FLAG="${ARM_SPA_FLAG[$label]}"
            P_COL="${ARM_P_COL[$label]}"
            echo "[F/G] arm=${label}: STAAR (use_SPA=${SPA_FLAG}, p-col=${P_COL}, rare_maf=${RARE_MAF}) ..."

            CSVS=()
            for chr in "${CHR_ARR[@]}"; do
                OUTSUB="staar_${label}_chr${chr}"
                bcft -v "$WD":/w:z -v "$ARM":/arm:z,ro "$STAAR_IMG" \
                    Rscript /arm/staar_gene_centric_coding.R \
                      --agds "/w/c${chr}.agds" --pheno "/w/pheno_${ver}.tsv" --chr "${chr}" \
                      --out-dir "/w/${OUTSUB}" --pheno-col Y1 --id-col IID \
                      --use-spa "${SPA_FLAG}" --rare-maf "${RARE_MAF}" \
                      "${GRM_ARG[@]}" "${MAXG_ARG[@]}" \
                    2>&1 | tee "${WD}/${OUTSUB}.log"
                for c in plof plof_ds missense disruptive_missense synonymous; do
                    [[ -e "${WD}/${OUTSUB}/${c}.csv" ]] && CSVS+=("/w/${OUTSUB}/${c}.csv")
                done
            done
            [[ ${#CSVS[@]} -gt 0 ]] || { echo "ERROR: arm=${label} produced no STAAR CSVs" >&2; exit 1; }

            # --p-col is NOT optional. STAAR emits STAAR-O only without SPA and
            # STAAR-B only with it, and staar_to_eval.py's default candidate list
            # names STAAR-O alone -- so an SPA run without this dies here rather than
            # writing a wrong table. That is the intended behaviour: loud, not silent.
            OUT_TABLE="$(arm_table "$label" "$idx")"
            TABLE_BASE="$(basename "$OUT_TABLE")"
            bcft -v "$WD":/w:z -v "$ARM":/arm:z,ro -v "$REGENIE_OUT_DIR":/out:z "$BCFTOOLS_IMG" \
                python3 /arm/staar_to_eval.py \
                  --staar-results "${CSVS[@]}" \
                  --genes-info /arm/genes_info_hgnc.tsv \
                  --p-col "${P_COL}" \
                  --out "/out/${label}/${TABLE_BASE}"

            # Inflation factor. The README states it is measured for every dataset in
            # every method, every time -- this arm is the one that was not doing it.
            # LOG10P is -log10(p) of the arm's omnibus column, so lambda is
            # qchisq(median p, 1, upper) / qchisq(0.5, 1). With no structure
            # correction both tools measured >8 on this data; near 1.0 is the evidence
            # that the principal components did their job.
            LAMBDA="$(bcft -v "$REGENIE_OUT_DIR":/out:z,ro "$STAAR_IMG" Rscript -e \
                'a <- commandArgs(TRUE); d <- read.table(a[1], header = TRUE, comment.char = "#"); p <- 10^(-as.numeric(d$LOG10P)); p <- p[is.finite(p) & p > 0 & p <= 1]; cat(if (length(p) == 0) "NA" else sprintf("%.4f", qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1)))' \
                "/out/${label}/${TABLE_BASE}" 2>/dev/null | tail -1)"
            N_ROWS="$(awk 'NR>2' "$OUT_TABLE" | wc -l)"
            N_GENES="$(awk 'NR>2{split($3,a,"."); print a[1]}' "$OUT_TABLE" | sort -u | wc -l)"
            printf 'use_spa_%s\t%s\n'        "$label" "$SPA_FLAG" >> "$RET"
            printf 'p_col_%s\t%s\n'          "$label" "$P_COL"    >> "$RET"
            printf 'n_result_rows_%s\t%s\n'  "$label" "$N_ROWS"   >> "$RET"
            printf 'n_genes_tested_%s\t%s\n' "$label" "$N_GENES"  >> "$RET"
            printf 'lambda_%s\t%s\n'         "$label" "${LAMBDA:-NA}" >> "$RET"
            echo "[run_${idx}] ${label} eval table -> ${OUT_TABLE}"
            echo "             ${N_ROWS} rows, ${N_GENES} genes, lambda ${LAMBDA:-NA}"
        done
    done

    # This script deletes nothing unless CLEANUP=true. The work dir is the only bulky
    # thing here; the tables, the retention TSV and the logs live outside it.
    if [[ "$CLEANUP" == "true" ]]; then
        echo "[cleanup] CLEANUP=true -- removing ${WD}"
        rm -rf "$WD"
    else
        echo "[run_${idx}] disposable now that the tables above are written:"
        du -sh "$WD" 2>/dev/null | sed 's/^/    /'
        echo "    remove by hand with:  rm -rf '${WD}'"
    fi
    echo "[run_${idx}] done: $(date -Is)"
    )                                # end of the per-dataset subshell
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
        echo "[run_${idx}] FAILED (exit ${rc}) -- recorded and skipped; continuing." >&2
        echo "             (leftover work under ${RUN_DIR}/work/run_${idx} may need removing" >&2
        echo "              by hand before this dataset can be retried.)" >&2
        FAILED_IDXS+=("$idx")
    fi
done

if [[ ${#FAILED_IDXS[@]} -gt 0 ]]; then
    echo ""
    echo "NOTE: ${#FAILED_IDXS[@]} dataset(s) did not produce a full set of arm tables: ${FAILED_IDXS[*]}"
    echo "      Scoring below simply omits them (the regenie glob only matches tables that exist)."
fi

# ----------------------------------------------------------------------------
# Scoring (optional). One run_eval.sh call per version, over that version's result
# tables. CAUSAL_*_GLOB default to the known answers restricted to the numbered
# chromosomes (the output of filter_causal_autosomal.py), so causal genes on
# chromosome X count as misses. Point them elsewhere to score differently.
# Then compare with benchmark-common/pairwise_compare.py, passing
# --arm-b <METHOD> <METHOD>_<arm>_eval  and --missing zero.
#
# The eval dir is a SIBLING of RUN_DIR under runs/, not a child: pairwise_compare.py
# and three_method_ap_bar.py both resolve eval subdirs directly under runs/.
# ----------------------------------------------------------------------------
if [[ "$SCORE" == "true" ]]; then
    : "${CAUSAL_SNPLIST_GLOB:?set CAUSAL_SNPLIST_GLOB (autosomal truth) to score}"
    : "${CAUSAL_GENES_GLOB:?set CAUSAL_GENES_GLOB (autosomal truth) to score}"
    RUNS_DIR="$(dirname "$RUN_DIR")"
    for label in "${ARM_LABELS[@]}"; do
        echo ""
        echo "=== scoring arm=${label} ==="
        export EVAL_REPO="${EVAL_REPO:-/data/git/doktorat_pw/wum_pims/nf-eval-gene-assoc}"
        export EVAL_RUN_DIR="${RUNS_DIR}/${METHOD}_${label}_eval"
        export EVAL_PROJECT="${METHOD}_${label}"
        export EVAL_PROFILE="${EVAL_PROFILE:-podman,medium_resources}"
        export INPUT_VCF SKIP_PREP="true"
        export REGENIE_GLOB="${REGENIE_OUT_DIR}/${label}/${METHOD}_${label}_dataset_idx_*_step2_Y1.regenie"
        export CAUSAL_SNPLIST_GLOB CAUSAL_GENES_GLOB
        bash "${COMMON}/run_eval.sh"
    done
fi

echo ""
echo "RICOPILI + STAARpipeline run finished: $(date -Is)"
