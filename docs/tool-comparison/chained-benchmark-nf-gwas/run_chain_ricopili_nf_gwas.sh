#!/usr/bin/env bash
#
# Runs RICOPILI (quality control + principal components) combined with nf-gwas
# (association) as one method, for comparison against nf-rare-var-assoc.
# ========================================================================================
# Both halves are complete pipelines a user configures and runs, with no analysis code
# of ours between them. That is why this uses RICOPILI's own menv.trans.mds principal
# components rather than the GENESIS components the STAARpipeline comparison needs.
#
# See README.md in this directory for prerequisites, the conversion recipe and its
# pitfalls, and the practical problems worth knowing about.
#
# The only thing borrowed from nf-rare-var-assoc is the gene groupings (.annotations
# and .setlist), which nf-gwas cannot build. Sample subsetting, quality control,
# principal components and the association engine are all the other tools'.
#
# ONE VERSION, unlike the STAARpipeline comparison. RICOPILI's quality-controlled
# sample set is kept WITH the relatives, because REGENIE models relatedness in its
# whole-genome step. Relatedness shapes only the principal-component axes, which is
# also how nf-rare-var-assoc handles it.
#
# STAGES PER DATASET:
#   A raw VCF -> sample subset -> genotypes only -> split multi-allelic -> PLINK bed
#   B preimp_dir  (sample and variant quality control)
#   C pcaer       -> menv.trans.mds (unrelated samples plus projected relatives)
#   D bed -> association VCF, prediction bed and covariates, then the identifier
#     check -- the one place this can fail silently
#   E nf-gwas (REGENIE steps 1 and 2, gene-based, on the borrowed groupings)
#   F rename the result table; compute the inflation factor and diagnostics
# With SCORE=true, benchmark-common/run_eval.sh then scores the results.
#
# STAGES A-C ARE COPIED VERBATIM from run_chain_ricopili_staar.sh and deliberately
# not factored into benchmark-common/: the repository should hold the exact code that
# produced each comparison's numbers. That includes the quality-control settings
# below -- both comparisons must run RICOPILI at the same thresholds, or their
# quality-control halves stop being the same pipeline.
#
# QUALITY-CONTROL THRESHOLDS are set explicitly rather than left at RICOPILI's
#   array-era defaults, which remove so much from exome data that the smaller datasets
#   run out of variants. The PREIMP_* values below are preimp_dir's own command-line
#   options; nothing is patched and every other threshold stays as shipped. Set them
#   to 0.02/0.02/0.05/0.02 to reproduce a run with RICOPILI's defaults. The reasoning
#   is in ../chained-benchmark/README.md.
#
# QUICK CHECK: EXPORT_ONLY=true runs stages A-D and the identifier check, then stops
#   (about 5 minutes per dataset instead of a full Nextflow run). This is the cheap
#   way to confirm the identifier convention before spending time on association.
#
# DATASETS: set DATASET_IDXS (space-separated); the default is every run_<N> found.
#
set -euo pipefail

# ----------------------------------------------------------------------------
# Configuration. Same convention as the STAARpipeline script: every location
# defaults to one machine's layout but can be set from the environment, so this
# runs elsewhere without editing. Absolute paths only, since they are passed into
# containers.
# ----------------------------------------------------------------------------
DATA="${DATA:-/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/tools_comparison}"
RVA_REPO="${RVA_REPO:-/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc}"
NFGWAS_REPO="${NFGWAS_REPO:-/data/git/doktorat_pw/nf_rare_var_assoc__tools_comparison/genepi/nf-gwas}"
DATASETS_DIR="${DATASETS_DIR:-${DATA}/datasets}"
COMMON="${COMMON:-${RVA_REPO}/docs/tool-comparison/benchmark-common}"

INPUT_VCF="${INPUT_VCF:-${DATA}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz}"
PED="${PED:-${RVA_REPO}/assets/integrated_call_samples_v3.20250704.ALL.ped}"   # 1000G sex source

# Reference for the stage D normalisation. This is the FASTA nf-prepare-vcf uses, and
# the association VCF's identifiers only match the borrowed masks when both are
# normalised against the same sequence. Chr-prefixed, so stage D exports chr-prefixed
# names and strips the prefix afterwards. The .fai must sit next to it.
REFDIR="${REFDIR:-/data/doktorat/biodatageeks/1000g/GRCh38_reference_genome}"
REF="${REF:-GRCh38_full_analysis_set_plus_decoy_hla.fa}"

# The one borrowed input: the per-dataset gene groupings published by a run of
# nf-rare-var-assoc with --publish_intermediate true.
RVA_RESULTS="${RVA_RESULTS:-${DATA}/runs/nf_rare_var_assoc/results}"
RVA_PROJECT="${RVA_PROJECT:-tools_comparison}"
GENE_MASKS="${GENE_MASKS:-${RVA_REPO}/assets/default.masks}"   # static mask definition

# NOTE: the default moved to .../ricopili_nf_gwas_qcmatched with the QC harmonization,
# so the original defaults run under .../ricopili_nf_gwas stays on disk for comparison.
# Set RUN_DIR=${DATA}/runs/ricopili_nf_gwas to overwrite it instead.
RUN_DIR="${RUN_DIR:-${DATA}/runs/ricopili_nf_gwas_qcmatched}"
REGENIE_OUT_DIR="${RUN_DIR}/regenie_per_dataset"

# RICOPILI QC thresholds. All four are preimp_dir's own command-line options -- no
# source patch -- and are the ONLY thresholds we move off RICOPILI's defaults; HWE
# (1e-6 co / 1e-10 ca), F-het (0.2) and maf (0) stay as shipped.
# KEEP THESE IN SYNC WITH run_chain_ricopili_staar.sh.
#   geno 0.10 : preimp_dir default 0.02. Matches plink2_makepgen_4_options --geno 0.100,
#               the variant missingness nf-rare-var-assoc applies to its association set.
#   mind 0.20 : preimp_dir default 0.02. Matches plink2_makepgen_5_options --mind 0.200 /
#               plink2_write_snplist_step2_options --mind 0.200.
#   pre_geno  : preimp_dir default 0.05, and RAISING --geno ALONE IS NOT ENOUGH. In
#     0.35      rep_qc2_14 the pre-filter list (variants with missingness > pre_geno,
#               measured before any sample filtering) is --exclude'd from every later
#               step and never comes back, so the effective variant cut is
#               min(pre_geno, geno). Left at 0.05 it would silently bind at 0.05 and
#               undo --geno 0.10. RICOPILI's own defaults keep pre_geno LOOSER than geno
#               (0.05 vs 0.02) precisely because the pre-filter is a guard against
#               catastrophically bad variants inflating sample missingness, not the
#               binding cut. 0.35 restores that ordering and matches
#               nf-rare-var-assoc's own loose pre-pass, --geno 0.350.
#   midi 0.35 : preimp_dir default 0.02, and NOT a like-for-like port -- RICOPILI bounds
#               the case/control missingness DIFFERENCE, while nf-rare-var-assoc's
#               FILTER_MISSING_PER_PHENO bounds each group's ABSOLUTE rate separately
#               (--geno 0.350 per phenotype value, intersected). 0.35 is the largest
#               differential our own filter can let through, i.e. the loosest cut that is
#               still implied by this pipeline rather than invented. The default 0.02
#               is punishing here: with ~11 cases in the most imbalanced datasets, a
#               single missing genotype is a 9% differential.
# pre_geno and midi both come out at 0.35 because ONE nf-rare-var-assoc step
# (FILTER_MISSING_PER_PHENO, --geno 0.350 per phenotype value, intersected) does both
# jobs at once: a loose absolute pre-pass and an implied cap on the differential.
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
NPCS_COVAR="${NPCS_COVAR:-4}"        # principal components used as covariates
THREADS="${THREADS:-8}"
MIN_ID_MATCH="${MIN_ID_MATCH:-0.5}"  # least fraction of exported variants that must be
                                     # named in .annotations, or the run stops
EXPORT_ONLY="${EXPORT_ONLY:-false}"  # true -> stop after stage D and the identifier check
# nf-gwas asks for 8 CPUs / 16 GB for its REGENIE processes (its conf/base.config) and
# Nextflow aborts on a host with less. Set either of these on a smaller machine, e.g.
# NFGWAS_MAX_CPUS=4 NFGWAS_MAX_MEM=12.GB (Nextflow memory literal). Empty = no clamp.
NFGWAS_MAX_CPUS="${NFGWAS_MAX_CPUS:-}"
NFGWAS_MAX_MEM="${NFGWAS_MAX_MEM:-}"
SCORE="${SCORE:-false}"              # true -> run run_eval.sh at the tail
# NOTE: this script deletes nothing. It prints the size of the two disposable
# per-dataset directories at the end of each dataset; remove them by hand.

# Container images.
BCFTOOLS_IMG="${BCFTOOLS_IMG:-docker.io/psuszynski/bioinf_combo:1.5.1}"
RICOPILI_IMG="${RICOPILI_IMG:-localhost/ricopili:2025_Feb_20.001}"
NFGWAS_CONTAINER="${NFGWAS_CONTAINER:-quay.io/genepi/nf-gwas:v1.0.9}"

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
for f in "$INPUT_VCF" "$PED" "$GENE_MASKS" "${NFGWAS_REPO}/main.nf" \
         "${REFDIR}/${REF}" "${REFDIR}/${REF}.fai"; do
    [[ -e "$f" ]] || { echo "ERROR: missing required path: $f" >&2; exit 1; }
done
mkdir -p "$RUN_DIR" "$REGENIE_OUT_DIR"

VCFBASE="$(basename "$INPUT_VCF")"

# Container helpers (as in the STAARpipeline script): :z relabels, --userns=keep-id lets
# rootless podman write host-owned dirs.
bcft() { podman run --rm --userns=keep-id "$@"; }

# median of a numeric stream (empty stream -> NA)
median() { sort -g | awk '{v[n++]=$1} END{ if(n==0){print "NA"} else if(n%2){print v[(n-1)/2]} else {print (v[n/2-1]+v[n/2])/2} }'; }

# ============================================================================
# PER-DATASET LOOP
#
# Each dataset runs in its own subshell so that ONE dataset's failure does not
# abort the whole benchmark. Some simulated phenotypes are legitimately
# unrunnable (e.g. REGENIE step 1 aborts with "all phenotypes have less than 10
# cases" when a draw has too few cases) -- those are recorded and skipped, not
# fatal. The subshell is a standalone command (not in an `if`/`||`), so `set -e`
# stays ACTIVE inside it and the dataset still stops at its first real error;
# only the outer `set +e` keeps that from killing the loop.
#
# RESUME: a dataset whose final eval table already exists is skipped untouched,
# so a re-launch after an interruption picks up where it stopped without redoing
# finished datasets or deleting anything.
# ============================================================================
FAILED_IDXS=()
for idx in "${IDXS[@]}"; do
    OUT_TABLE="${REGENIE_OUT_DIR}/ricopili_nf_gwas_qcmatched_dataset_idx_${idx}_step2_Y1.regenie"
    if [[ -e "$OUT_TABLE" ]]; then
        echo "[run_${idx}] eval table already present -- skipping (${OUT_TABLE})"
        continue
    fi

    set +e
    ( set -e
    PHENO="$(pheno_path "$idx")"
    [[ -e "$PHENO" ]] || { echo "ERROR: missing phenotype for run_${idx}: $PHENO" >&2; exit 1; }
    STUDY="run$(printf '%02d' "$idx")"          # 5-char RICOPILI study name (idx<=99)

    # The borrowed input: the per-dataset gene groupings.
    RVA_ID="${RVA_PROJECT}_dataset_idx_${idx}"
    ANNO_FILE="${RVA_RESULTS}/bcftools/${RVA_ID}.annotations"
    SETLIST_FILE="${RVA_RESULTS}/bcftools/${RVA_ID}.setlist"
    for p in "$ANNO_FILE" "$SETLIST_FILE"; do
        [[ -e "$p" ]] || { echo "ERROR: reused mask file for run_${idx} not found: $p" >&2
                           echo "       These come from an nf-rare-var-assoc run (project ${RVA_PROJECT})." >&2
                           exit 1; }
    done

    WD="${RUN_DIR}/work/run_${idx}"
    # RICOPILI records progress and refuses to repeat a step that "has been done
    # repeatedly without any progress", so a retry has to start from a clean
    # directory. This script will not delete one for you -- it stops and says so.
    if [[ -d "$WD" ]] && [[ -n "$(ls -A "$WD" 2>/dev/null)" ]]; then
        echo "ERROR: work dir for run_${idx} already exists and is not empty:" >&2
        echo "         $WD" >&2
        echo "       RICOPILI will not re-run in a dirty directory. Remove it yourself" >&2
        echo "       (rm -rf that path) and start again." >&2
        exit 1
    fi
    mkdir -p "$WD"
    RES="${RUN_DIR}/results/run_${idx}"
    NFWORK="${RUN_DIR}/nfwork/run_${idx}"
    RET="${RUN_DIR}/retention/run_${idx}.tsv"; mkdir -p "$(dirname "$RET")"
    : > "$RET"
    # Provenance: which QC thresholds produced these counts (defaults vs harmonized).
    printf 'preimp_qc_args\t%s\n' "$PREIMP_QC_ARGS" >> "$RET"

    echo "=================================================================="
    echo " chain run_${idx}  (study ${STUDY}) -- RICOPILI -> nf-gwas"
    echo "   pheno   : ${PHENO}"
    echo "   masks   : ${SETLIST_FILE}"
    echo "   QC      : ${PREIMP_QC_ARGS}"
    echo "   out     : ${RUN_DIR}"
    echo "   started : $(date -Is)"
    echo "=================================================================="

    # ------------------------------------------------------------------ Stage A
    # COPIED VERBATIM from the STAARpipeline script. Keep-list, sex (1000G pedigree),
    # plink phenotype;
    # raw VCF -> subset -> GT-only -> split -> bed. Two plink passes (plink2
    # converts, plink1.9 re-sorts the --split-par output that plink1.9 otherwise
    # refuses -- see ../chained-benchmark/README.md).
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
    # minor allele, which loses which allele was REF (stage D has to restore it).
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
    # COPIED VERBATIM from the STAARpipeline script. preimp_dir runs TWICE: the first call
    # writes a template name file and exits, we set the 5-char study name, the
    # second call does the QC.
    # PREIMP_QC_ARGS goes to BOTH calls: preimp_dir re-reads its whole command line on
    # the second pass, so omitting them there would silently run the QC at defaults.
    echo "[B] preimp_dir (sample + variant QC) ..."
    echo "    QC thresholds: ${PREIMP_QC_ARGS}"
    mkdir -p "${WD}/preimp"
    for e in bed bim fam; do cp "${WD}/sorted.${e}" "${WD}/preimp/${STUDY}_raw.${e}"; done
    bcft -v "$WD":/w:z "$RICOPILI_IMG" bash -lc \
        "cd /w/preimp && preimp_dir --disease sim --outname ${STUDY} --popname mix --serial ${PREIMP_QC_ARGS}" || true
    sed -i "s/^sim1\t${STUDY}_raw/${STUDY}\t${STUDY}_raw/" "${WD}/preimp/sim.names"
    bcft -v "$WD":/w:z "$RICOPILI_IMG" bash -lc \
        "cd /w/preimp && preimp_dir --disease sim --outname ${STUDY} --popname mix --serial ${PREIMP_QC_ARGS}" \
        2>&1 | tee "${WD}/preimp/preimp.log"

    QCPRE="${WD}/preimp/qc/sim_${STUDY}_mix_rp-qc1"
    for e in bed bim fam; do
        [[ -e "${QCPRE}.${e}" ]] || { echo "ERROR: preimp_dir produced no ${QCPRE}.${e}" >&2; exit 1; }
    done
    if [[ -e "${QCPRE}.meta" ]]; then
        grep -E '^(nsnpex_(mono|prefilter|hwe-co|miss|midi|hwe-ca|prekno)|nidex_(fhet|miss|sexcheck_ex))' \
            "${QCPRE}.meta" | tr -s ' ' '\t' >> "$RET" || true
    fi
    N_VAR_QC="$(wc -l < "${QCPRE}.bim")"
    printf 'n_variants_after_qc\t%s\n' "$N_VAR_QC" >> "$RET"
    printf 'n_samples_after_qc\t%s\n'  "$(wc -l < "${QCPRE}.fam")" >> "$RET"
    # Running out of variants is the failure mode this hits silently: too few variants
    # survive quality control, the identifier check or REGENIE finds nothing testable, and the dataset
    # returns an empty table stages later. Say so here instead.
    [[ "$N_VAR_QC" -ge 10000 ]] || \
        echo "WARNING: run_${idx} has only ${N_VAR_QC} variants after QC -- expect few or no testable genes." >&2

    # ------------------------------------------------------------------ Stage C
    # pcaer. We take menv.trans.mds -- the unrelated set PLUS the relatives
    # projected onto the same axes -- because this keeps the relatives and
    # lets REGENIE model them. smartpca still drops ancestry outliers, so this
    # file covers fewer samples than the quality-controlled fam. That gap is a result,
    # not a defect: smartpca drops some samples as ancestry outliers.
    echo "[C] pcaer (relatedness-aware smartpca PCs, projected onto all samples) ..."
    mkdir -p "${WD}/pca"
    for e in bed bim fam; do cp "${QCPRE}.${e}" "${WD}/pca/"; done
    bcft -v "$WD":/w:z "$RICOPILI_IMG" bash -lc \
        "cd /w/pca && pcaer --out ${STUDY}pca --serial sim_${STUDY}_mix_rp-qc1" \
        2>&1 | tee "${WD}/pca/pcaer.log"
    MDS="${WD}/pca/pcaer_${STUDY}pca/${STUDY}pca.menv.trans.mds"
    MDS_UNREL="${WD}/pca/pcaer_${STUDY}pca/${STUDY}pca.menv.mds"
    [[ -e "$MDS" ]] || { echo "ERROR: pcaer produced no ${MDS}" >&2; exit 1; }
    printf 'n_samples_pca_projected\t%s\n' "$(($(wc -l < "$MDS") - 1))" >> "$RET"
    [[ -e "$MDS_UNREL" ]] && printf 'n_samples_pca_unrelated\t%s\n' "$(($(wc -l < "$MDS_UNREL") - 1))" >> "$RET"

    # ------------------------------------------------------------------ Stage D
    # QC'd bed -> the two genotype inputs nf-gwas needs, with MASK-COMPATIBLE
    # variant identifiers, plus the covariate file. Then the identifier check.
    echo "[D] export mask-compatible genotypes + covariates ..."

    # Restore the REF allele before exporting. A PLINK bed does
    # not record which allele was REF, and preimp_dir's PLINK 1.9 passes make A1 the
    # minor allele, so a few percent of sites come out of QC with REF and ALT
    # exchanged. Our own variant IDs carry the answer -- stage A set them to
    # <chr>:<pos>:<REF>:<ALT> -- so field 3 is the true REF and --ref-allele restores
    # it exactly, with no reference-genome guessing. Without this the identifiers
    # bcftools re-derives below would be CHROM_POS_ALT_REF and would not match the masks,
    # and the normalisation below would abort on the exchanged sites.
    awk -F'\t' '{n=split($2,a,":"); if (n==4 && a[3]!="") print $2"\t"a[3]}' \
        "${QCPRE}.bim" > "${WD}/ref_alleles.txt"
    printf 'n_variants_ref_allele_forced\t%s\n' "$(wc -l < "${WD}/ref_alleles.txt")" >> "$RET"

    # The association VCF, in two plink2 passes. Both passes are forced on us; the
    # single-pass version is wrong in two different ways (measured, see the README):
    #
    #  pass 1  --merge-x folds the pseudo-autosomal variants (chromosome 25 after
    #          stage A's --split-par, exported as "XY") back into X, because the masks
    #          call them plain X. NOT --merge-par: our PAR sits under the XY code, and
    #          plink2 answers --merge-par with "had no effect (no PAR1/PAR2 chromosome
    #          codes present)" -- silently leaving 934 variants unmatchable. --merge-x
    #          needs --sort-vars, and --sort-vars needs --make-pgen, hence a second pass.
    #  pass 2  exports the VCF. --update-sex marks every sample female, which makes
    #          chromosome X export DIPLOID. That is deliberate: plink2 refuses to export
    #          a VCF straight after --merge-x precisely because the male calls would come
    #          out haploid, while nf-rare-var-assoc reads chromosome X
    #          from the input VCF, where 1000G codes males diploid. Exporting haploid
    #          would change the chromosome-X burden dosages relative to both, mixing an
    #          encoding artefact into the comparison. The true sex is untouched
    #          everywhere it is used -- RICOPILI's QC ran on it, and it stays in the
    #          prediction genotypes below.
    #
    # --output-chr chrM writes chr12 / chr22 / chrX. The prefix is what the reference
    # FASTA uses, so it has to be present for the normalisation below; it is stripped
    # straight afterwards, leaving the 12 / 22 / X the gene groupings use.
    # id-paste=iid is REQUIRED: plink2 names VCF samples FID_IID and preimp_dir
    # rewrites FID to 'con_sim_<study>_mix_rp_*<IID>'.
    { printf '#FID\tIID\tSEX\n'; awk '{print $1"\t"$2"\t2"}' "${QCPRE}.fam"; } > "${WD}/sex_diploid.txt"
    bcft -v "$WD":/w:z "$RICOPILI_IMG" bash -lc "
        set -e
        /opt/rp_dep/plink2/plink2 --bfile /w/preimp/qc/sim_${STUDY}_mix_rp-qc1 \
          --merge-x --sort-vars --ref-allele force /w/ref_alleles.txt 2 1 \
          --make-pgen --out /w/xmerged
        /opt/rp_dep/plink2/plink2 --pfile /w/xmerged --update-sex /w/sex_diploid.txt \
          --output-chr chrM --export vcf bgz id-paste=iid \
          --out /w/assoc.export"

    # The rename map for the step below: every contig of the export, chr prefix stripped.
    bcft -v "$WD":/w:z "$BCFTOOLS_IMG" bash -lc 'bcftools view -h /w/assoc.export.vcf.gz' \
        | grep '^##contig' | sed -E 's/.*ID=([^,>]+).*/\1/' \
        | awk '{n=$1; sub(/^chr/,"",n); print $1"\t"n}' > "${WD}/chr_map.txt"

    # The association VCF, prepared the way nf-prepare-vcf prepares the reference arm's
    # input, in that order: split multi-allelic sites and normalise against the reference
    # FASTA, strip the chr prefix, then set CHROM_POS_REF_ALT. The identifiers therefore
    # carry left-aligned indel positions, which is what makes the borrowed masks match.
    # bcftools norm exits on a REF/ALT mismatch -- the --ref-allele restoration above is
    # what keeps it from doing so.
    bcft -v "$WD":/w:z -v "$REFDIR":/ref:z,ro "$BCFTOOLS_IMG" bash -lc "
        set -e -o pipefail
        bcftools norm --fasta-ref /ref/${REF} -m -any --rm-dup exact \
            --threads ${THREADS} -Ou /w/assoc.export.vcf.gz 2> /w/norm.log \
        | bcftools annotate --rename-chrs /w/chr_map.txt -Ou \
        | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \
              --threads ${THREADS} -Oz -o /w/assoc.vcf.gz
        bcftools index -t /w/assoc.vcf.gz"
    sed 's/^/    [norm] /' "${WD}/norm.log"

    # Prediction genotypes for REGENIE step 1: the same QC'd bed, with FID rewritten
    # back to the sample id so it joins to the phenotype and covariate files. TRUE sex
    # rides along in the .fam (stage A set it from the 1000G pedigree), and the
    # pseudo-autosomal region stays on its own chromosome code -- the same structure
    # the reference run's prediction bed had, since it imported with --split-par b38.
    awk '{print $1"\t"$2"\t"$2"\t"$2}' "${QCPRE}.fam" > "${WD}/update_ids.txt"
    bcft -v "$WD":/w:z "$RICOPILI_IMG" bash -lc "
        set -e
        /opt/rp_dep/plink2/plink2 --bfile /w/preimp/qc/sim_${STUDY}_mix_rp-qc1 \
          --update-ids /w/update_ids.txt --output-chr MT \
          --make-bed --out /w/prediction"

    # PC covariates from RICOPILI's own projected MDS. Layout (whitespace-split):
    # FID IID SOL C1..C20 -> $2=IID, $4=C1. Samples smartpca dropped as ancestry
    # outliers are simply absent, and REGENIE will not analyse them.
    COVAR_OUT="${WD}/covar_pc.txt"
    awk -v n="$NPCS_COVAR" '
        FNR==1 { printf "FID\tIID"; for(i=1;i<=n;i++) printf "\tPC%d", i; printf "\n"; next }
        { printf "%s\t%s", $2, $2; for(i=1;i<=n;i++) printf "\t%s", $(3+i); printf "\n" }
    ' "$MDS" > "$COVAR_OUT"
    COVAR_COLS="$(awk -v n="$NPCS_COVAR" 'BEGIN{for(i=1;i<=n;i++) printf "%sPC%d", (i>1?",":""), i}')"
    printf 'n_samples_covar\t%s\n' "$(($(wc -l < "$COVAR_OUT") - 1))" >> "$RET"

    # ---- CHECK: the identifiers must match the borrowed gene groupings ----------
    bcft -v "$WD":/w:z "$BCFTOOLS_IMG" bash -lc \
        "bcftools query -f '%CHROM\t%ID\n' /w/assoc.vcf.gz > /w/exported_ids.tsv"
    cut -f1 "${WD}/exported_ids.tsv" | sort -u > "${WD}/exported_chroms.txt"
    if grep -qE '^(chr|XY|25|23)' "${WD}/exported_chroms.txt"; then
        echo "ERROR: exported chromosome names are not mask-compatible:" >&2
        tr '\n' ' ' < "${WD}/exported_chroms.txt" >&2; echo >&2
        echo "       masks use 12 / 22 / X with no chr prefix and no PAR code." >&2; exit 1
    fi
    cut -f2 "${WD}/exported_ids.tsv" > "${WD}/exported_ids.txt"
    N_EXPORTED="$(wc -l < "${WD}/exported_ids.txt")"
    N_MATCHED="$(awk -F'\t' 'NR==FNR{a[$1]=1; next} ($1 in a){m++} END{print m+0}' \
                     "$ANNO_FILE" "${WD}/exported_ids.txt")"
    FRAC_MATCHED="$(awk -v m="$N_MATCHED" -v n="$N_EXPORTED" 'BEGIN{printf "%.4f", (n?m/n:0)}')"
    printf 'n_variants_exported\t%s\n'        "$N_EXPORTED"   >> "$RET"
    printf 'n_exported_in_masks\t%s\n'        "$N_MATCHED"    >> "$RET"
    printf 'frac_exported_in_masks\t%s\n'     "$FRAC_MATCHED" >> "$RET"
    echo "    exported variants : ${N_EXPORTED}"
    echo "    named in masks    : ${N_MATCHED} (${FRAC_MATCHED})"
    awk -v f="$FRAC_MATCHED" -v t="$MIN_ID_MATCH" 'BEGIN{exit !(f<t)}' && {
        echo "ERROR: only ${FRAC_MATCHED} of exported variants are named in the masks" >&2
        echo "       (threshold MIN_ID_MATCH=${MIN_ID_MATCH}). The identifier convention" >&2
        echo "       is wrong -- REGENIE would silently test almost nothing. Compare:" >&2
        head -3 "${WD}/exported_ids.txt" >&2; head -3 "$ANNO_FILE" >&2; exit 1; }

    # Per chromosome as well. The overall fraction is dominated by the autosomes and
    # would stay comfortably above the threshold even if chromosome X were named wrongly
    # (X carries only 2k of the ~128k mask variants), so require every exported
    # chromosome to match SOMETHING. The X fraction itself is legitimately low: the
    # masks are coding-only and cover far less of X than of 12 or 22.
    awk -F'\t' 'NR==FNR{a[$1]=1; next} {t[$1]++; if($2 in a) m[$1]++}
                END{for(c in t) printf "chr_%s_exported_in_masks\t%d/%d\t%.4f\n", c, m[c]+0, t[c], (m[c]+0)/t[c]}' \
        "$ANNO_FILE" "${WD}/exported_ids.tsv" | sort >> "$RET"
    BAD_CHR="$(awk -F'\t' 'NR==FNR{a[$1]=1; next} {t[$1]++; if($2 in a) m[$1]++}
                           END{for(c in t) if((m[c]+0)==0) printf "%s ", c}' \
                   "$ANNO_FILE" "${WD}/exported_ids.tsv")"
    [[ -z "$BAD_CHR" ]] || { echo "ERROR: no mask match at all on chromosome(s): ${BAD_CHR}" >&2
                             echo "       the identifier convention is wrong for those." >&2; exit 1; }

    # ---- How many grouping variants survived RICOPILI's quality control ---------
    # The set-lists still name every variant nf-rare-var-assoc saw, so genes lose
    # variants here and some fall below testability. REGENIE skips absent variants;
    # what matters is that we MEASURE the loss.
    awk -F'\t' 'NR==FNR{keep[$1]=1; next}
                { n=split($4,v,","); c=0; for(i=1;i<=n;i++) if (v[i] in keep) c++;
                  print $1"\t"n"\t"c }' \
        "${WD}/exported_ids.txt" "$SETLIST_FILE" > "${WD}/gene_survival.tsv"
    printf 'n_genes_setlist\t%s\n'          "$(wc -l < "${WD}/gene_survival.tsv")" >> "$RET"
    printf 'n_genes_ge1_variant\t%s\n'      "$(awk -F'\t' '$3>=1' "${WD}/gene_survival.tsv" | wc -l)" >> "$RET"
    printf 'n_genes_ge2_variants\t%s\n'     "$(awk -F'\t' '$3>=2' "${WD}/gene_survival.tsv" | wc -l)" >> "$RET"
    printf 'median_variants_per_gene_before\t%s\n' "$(cut -f2 "${WD}/gene_survival.tsv" | median)" >> "$RET"
    printf 'median_variants_per_gene_after\t%s\n'  "$(cut -f3 "${WD}/gene_survival.tsv" | median)" >> "$RET"
    printf 'frac_setlist_variants_surviving\t%s\n' \
        "$(awk -F'\t' '{b+=$2; a+=$3} END{printf "%.4f", (b?a/b:0)}' "${WD}/gene_survival.tsv")" >> "$RET"

    if [[ "$EXPORT_ONLY" == "true" ]]; then
        echo "[run_${idx}] EXPORT_ONLY=true -- stopping after the identifier check. Retention: ${RET}"
        exit 0                       # ends this dataset's subshell, not the loop
    fi

    # ------------------------------------------------------------------ Stage E
    # nf-gwas. Parameters are copied from nf-rare-var-assoc's own tuned settings, so
    # the step-1 quality-control layer is identical across the compared methods; they differ only
    # in what they are fed. nf-gwas ALWAYS runs its own QC pass on the prediction
    # genotypes (there is no off switch), so this carries RICOPILI's quality control and
    # nf-gwas's QC in sequence -- which is exactly what chaining the two gives a user.
    mkdir -p "$RES"
    PARAMS_FILE="${WD}/params.nf_gwas.yaml"
    cat > "$PARAMS_FILE" <<EOF
project: "ricopili_nf_gwas_qcmatched_dataset_idx_${idx}"
outdir: "${RES}"

genotypes_prediction: "${WD}/prediction.{bed,bim,fam}"
genotypes_association: "${WD}/assoc.vcf.gz"
genotypes_association_format: "vcf"
genotypes_build: "hg38"
association_build: "hg38"

phenotypes_filename: "${PHENO}"
phenotypes_columns: "Y1"
phenotypes_binary_trait: true

covariates_filename: "${COVAR_OUT}"
covariates_columns: "${COVAR_COLS}"

# nf-gwas's own quality-control call
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

# Gene-based (RVAS) -- THE LOAN: nf-rare-var-assoc's per-dataset mask files
regenie_run_gene_based_tests: true
regenie_gene_anno: "${ANNO_FILE}"
regenie_gene_setlist: "${SETLIST_FILE}"
regenie_gene_masks: "${GENE_MASKS}"
regenie_gene_aaf: 0.2
regenie_gene_test: "skato"
EOF

    EXTRA_CFG="${WD}/nfgwas_podman.config"
    cat > "$EXTRA_CFG" <<EOF
podman.enabled    = true
docker.enabled    = false
podman.runOptions = '--userns=keep-id'
podman.mountFlags = 'z'
process.container = '${NFGWAS_CONTAINER}'
EOF
    if [[ -n "$NFGWAS_MAX_CPUS" || -n "$NFGWAS_MAX_MEM" ]]; then
        LIMITS=""
        [[ -n "$NFGWAS_MAX_CPUS" ]] && LIMITS="cpus: ${NFGWAS_MAX_CPUS}"
        [[ -n "$NFGWAS_MAX_MEM"  ]] && LIMITS="${LIMITS:+${LIMITS}, }memory: ${NFGWAS_MAX_MEM}"
        echo "process.resourceLimits = [ ${LIMITS} ]" >> "$EXTRA_CFG"
        echo "    [resource clamp] ${LIMITS}"
    fi

    echo "[E] nf-gwas (REGENIE step 1 + step 2) ..."
    ( cd "$NFGWAS_REPO" && nextflow run "$NFGWAS_REPO/main.nf" \
        -c "$EXTRA_CFG" \
        -params-file "$PARAMS_FILE" \
        -work-dir "$NFWORK" \
        -with-trace "${RES}/trace.txt" \
        -with-report "${RES}/report.html" \
        -with-timeline "${RES}/timeline.html" )

    # ------------------------------------------------------------------ Stage F
    # Rename the REGENIE table into the shared dir so ONE eval fan-out scores every
    # dataset (nf-eval parses dataset_idx from the filename, reads sep=' ' comment='#').
    OUT_TABLE="${REGENIE_OUT_DIR}/ricopili_nf_gwas_qcmatched_dataset_idx_${idx}_step2_Y1.regenie"
    [[ -e "${RES}/results/Y1.regenie.gz" ]] \
        || { echo "ERROR: nf-gwas produced no ${RES}/results/Y1.regenie.gz" >&2; exit 1; }
    zcat "${RES}/results/Y1.regenie.gz" > "${WD}/.tmp_Y1.regenie"
    echo "## ignored line" > "$OUT_TABLE"
    awk '{ gsub(/\t/, " "); print }' "${WD}/.tmp_Y1.regenie" >> "$OUT_TABLE"

    # Non-degeneracy + inflation. Both tools measured lambda > 8 on this fixture with
    # no structure correction, so a value near 1.0 here is the evidence
    # that RICOPILI's principal components did their job. CHISQ is regenie's own 1-df
    # statistic, present for the burden and SKAT-O rows alike.
    N_ROWS="$(awk 'NR>2' "$OUT_TABLE" | wc -l)"
    N_GENES="$(awk 'NR>2{split($3,a,"."); print a[1]}' "$OUT_TABLE" | sort -u | wc -l)"
    N_DISTINCT_P="$(awk 'NR>2{print $12}' "$OUT_TABLE" | sort -u | wc -l)"
    printf 'n_result_rows\t%s\n'     "$N_ROWS"        >> "$RET"
    printf 'n_genes_tested\t%s\n'    "$N_GENES"       >> "$RET"
    printf 'n_distinct_log10p\t%s\n' "$N_DISTINCT_P"  >> "$RET"
    for t in ADD ADD-SKAT ADD-SKATO; do
        LAM="$(awk -v t="$t" 'NR>2 && $8==t && $11!="NA"{print $11}' "$OUT_TABLE" \
               | median | awk '{ if($1=="NA"){print "NA"} else {printf "%.4f", $1/0.4549364} }')"
        printf 'lambda_%s\t%s\n' "$t" "$LAM" >> "$RET"
        echo "    lambda ${t}: ${LAM}"
    done
    [[ "$N_ROWS" -gt 0 && "$N_DISTINCT_P" -gt 1 ]] \
        || { echo "ERROR: degenerate REGENIE output for run_${idx} (${N_ROWS} rows, ${N_DISTINCT_P} distinct p)" >&2
             exit 1; }
    echo "[run_${idx}] eval table -> ${OUT_TABLE}  (${N_ROWS} rows, ${N_GENES} genes)"

    # This script NEVER deletes anything. The two bulky per-dataset directories below
    # are safe to remove by hand once the eval table exists -- the table, the retention
    # TSV and the nf-gwas reports live outside them. Sizes are printed so the decision
    # is informed.
    echo "[run_${idx}] disposable once the table above is safe (delete by hand):"
    du -sh "$WD" "$NFWORK" 2>/dev/null | sed 's/^/    /'
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
    echo "NOTE: ${#FAILED_IDXS[@]} dataset(s) did not produce an eval table: ${FAILED_IDXS[*]}"
    echo "      Scoring below simply omits them (the regenie glob only matches tables that exist)."
fi

# ----------------------------------------------------------------------------
# Scoring (optional). Unlike the STAARpipeline comparison there is no reduced set of
# known answers:
# REGENIE tests chromosome X and the masks cover all three chromosomes, so this
# this scores against the COMPLETE list of causal genes. That
# is why CAUSAL_*_GLOB default to the unfiltered per-dataset truth files.
#
# Two steps: (1) run_eval.sh fans nf-eval-gene-assoc over every dataset table;
# (2) pairwise_compare.py compares this against nf-rare-var-assoc and against the
# nf-gwas-alone comparison, writing the paired statistics AND the figures (headline
# per-dataset difference, paired scatter, metrics forest -- .png + .pdf each).
# Difference is reported as reference - comparator, so a positive value is the
# reference pipeline scoring higher.
# ----------------------------------------------------------------------------
if [[ "$SCORE" == "true" ]]; then
    export EVAL_REPO="${EVAL_REPO:-/data/git/doktorat_pw/wum_pims/nf-eval-gene-assoc}"
    # Sibling of runs/ricopili_nf_gwas_qcmatched, not a child: pairwise_compare.py takes eval
    # subdirs directly under runs/ (--arm-b ricopili_nf_gwas_qcmatched ricopili_nf_gwas_qcmatched_eval).
    export EVAL_RUN_DIR="${EVAL_RUN_DIR:-${DATA}/runs/ricopili_nf_gwas_qcmatched_eval}"
    export EVAL_PROJECT="${EVAL_PROJECT:-ricopili_nf_gwas_qcmatched}"
    export EVAL_PROFILE="${EVAL_PROFILE:-podman,medium_resources}"
    export INPUT_VCF SKIP_PREP="true"
    export REGENIE_GLOB="${REGENIE_OUT_DIR}/ricopili_nf_gwas_qcmatched_dataset_idx_*_step2_Y1.regenie"
    export CAUSAL_SNPLIST_GLOB="${CAUSAL_SNPLIST_GLOB:-${DATASETS_DIR}/run_*/select_genes/*_dataset_idx_*_in_*.snplist}"
    export CAUSAL_GENES_GLOB="${CAUSAL_GENES_GLOB:-${DATASETS_DIR}/run_*/select_genes/*_genes_dataset_idx_*.txt}"
    echo ""
    echo "=== scoring (full causal truth) ==="
    bash "${COMMON}/run_eval.sh"

    # ---- Pairwise comparison + figures --------------------------------------
    # pairwise_compare.py needs pandas + scipy + matplotlib. PAIRWISE_PYTHON points
    # at an interpreter that has them; fall back to python3 if the preferred venv is
    # not present on this host.
    PAIRWISE_PYTHON="${PAIRWISE_PYTHON:-/data/git/playground_all/python/.venv/bin/python}"
    [[ -x "$PAIRWISE_PYTHON" ]] || PAIRWISE_PYTHON="python3"
    RUNS_DIR="$(dirname "$EVAL_RUN_DIR")"                 # .../runs
    EVAL_SUBDIR="$(basename "$EVAL_RUN_DIR")"             # ricopili_nf_gwas_qcmatched_eval
    PW_OUT="${PW_OUT:-${RUNS_DIR}/pairwise_ricopili_nf_gwas_qcmatched}"

    # --missing zero: a dataset this could not analyse (e.g. run_27, whose quality control left
    # too few cases for REGENIE) scores 0 rather than being dropped from the pairing --
    # excluding it would reward the method for failing to run on a dataset the reference
    # pipeline handled. See pairwise_compare.py's --missing docs.
    MISSING="${PAIRWISE_MISSING:-zero}"
    echo ""
    echo "=== pairwise comparison + figures (missing=${MISSING}) -> ${PW_OUT} ==="
    # (a) headline: nf-rare-var-assoc against this combined method.
    "$PAIRWISE_PYTHON" "${COMMON}/pairwise_compare.py" \
        --runs "$RUNS_DIR" --missing "$MISSING" \
        --arm-a nf_rare_var_assoc nf_rare_var_assoc_eval \
        --arm-b ricopili_nf_gwas_qcmatched "$EVAL_SUBDIR" \
        --out "${PW_OUT}/vs_reference"
    # (b) cross-check: nf-gwas alone (given our quality control and components)
    #     against this one (given RICOPILI's).
    if [[ -d "${RUNS_DIR}/nf_gwas_eval/results/compute_score" ]]; then
        "$PAIRWISE_PYTHON" "${COMMON}/pairwise_compare.py" \
            --runs "$RUNS_DIR" --missing "$MISSING" \
            --arm-a nf_gwas nf_gwas_eval \
            --arm-b ricopili_nf_gwas_qcmatched "$EVAL_SUBDIR" \
            --out "${PW_OUT}/vs_nf_gwas_E1"
    else
        echo "  (skipping the E1 nf-gwas cross-check: no ${RUNS_DIR}/nf_gwas_eval on this host)"
    fi
fi

echo ""
echo "chained (RICOPILI->nf-gwas) arm finished: $(date -Is)"
