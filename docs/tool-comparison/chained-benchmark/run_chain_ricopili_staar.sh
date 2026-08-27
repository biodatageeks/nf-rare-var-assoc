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
# TWO VERSIONS PER DATASET, differing only in how relatedness is handled:
#   Filtered = RICOPILI as shipped. pcaer removes related samples; its smartpca
#              principal components carry into the association test.
#   Full     = all samples kept. Relatedness is modelled in the null model with a
#              GENESIS PC-Relate matrix, and principal components come from GENESIS
#              PC-AiR. This is how STAARpipeline is designed to be used, and it
#              matches nf-rare-var-assoc's own handling of relatedness.
#
# STAGES PER DATASET:
#   A raw VCF -> sample subset -> genotypes only -> split multi-allelic -> PLINK bed
#   B preimp_dir  (sample and variant quality control)
#   C pcaer       (relatedness removal + smartpca components; Filtered version)
#   D GENESIS PC-AiR + PC-Relate matrix (Full version)  genesis_pcair_pcrelate.R
#   E per chromosome: bed -> VCF -> GDS -> annotated GDS (FAVORannotator)
#   F association test per version and chromosome       staar_gene_centric_coding.R
#   G staar_to_eval.py -> one REGENIE-shaped result table per version
# With SCORE=true, benchmark-common/run_eval.sh then scores each version.
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
ARM="${ARM:-${RVA_REPO}/docs/tool-comparison/chained-benchmark}"
COMMON="${COMMON:-${RVA_REPO}/docs/tool-comparison/benchmark-common}"

INPUT_VCF="${INPUT_VCF:-${DATA}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz}"
PED="${PED:-${RVA_REPO}/assets/integrated_call_samples_v3.20250704.ALL.ped}"   # 1000G sex source
REFDIR="${REFDIR:-/data/doktorat/biodatageeks/genome_in_a_bottle/reference}"
REF="${REF:-GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz}"  # chr-prefixed
FAVOR_DB="${FAVOR_DB:-${DATA}/favor/db}"          # chr<N>_*.csv(.idx), extracted per chr
GENES_INFO="${GENES_INFO:-${ARM}/genes_info_hgnc.tsv}"

# NOTE: the default moved to .../ricopili_staar_qcmatched with the QC harmonization, so
# the original defaults run under .../ricopili_staar stays on disk for comparison.
# Set RUN_DIR=${DATA}/runs/ricopili_staar to overwrite it instead.
RUN_DIR="${RUN_DIR:-${DATA}/runs/ricopili_staar_qcmatched}"
REGENIE_OUT_DIR="${RUN_DIR}/regenie_per_dataset"  # per-version subdirs (filtered/full)

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
VERSIONS="${VERSIONS:-filtered full}"
NPCS_COVAR="${NPCS_COVAR:-4}"      # principal components used as covariates
USE_SPA="${USE_SPA:-FALSE}"        # STAAR-O (FALSE) or STAAR-B (TRUE)
RARE_MAF="${RARE_MAF:-0.01}"       # STAAR rare_maf_cutoff
THREADS="${THREADS:-4}"
MAX_GENES="${MAX_GENES:-}"         # smoke-test cap on genes tested per chr; empty = all genes
CLEANUP="${CLEANUP:-true}"         # drop bulky per-dataset work after the eval table is made
SCORE="${SCORE:-false}"            # true -> run run_eval.sh per version at the tail

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
for v in "${VER_ARR[@]}"; do mkdir -p "${REGENIE_OUT_DIR}/${v}"; done

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
# ============================================================================
for idx in "${IDXS[@]}"; do
    PHENO="$(pheno_path "$idx")"
    [[ -e "$PHENO" ]] || { echo "ERROR: missing phenotype for run_${idx}: $PHENO" >&2; exit 1; }
    STUDY="run$(printf '%02d' "$idx")"          # 5-char RICOPILI study name (idx<=99)

    WD="${RUN_DIR}/work/run_${idx}"
    rm -rf "$WD"; mkdir -p "$WD"                 # RICOPILI refuses to re-run in a dirty dir
    RET="${RUN_DIR}/retention/run_${idx}.tsv"; mkdir -p "$(dirname "$RET")"
    : > "$RET"
    # Provenance: which QC thresholds produced these counts (defaults vs harmonized).
    printf 'preimp_qc_args\t%s\n' "$PREIMP_QC_ARGS" >> "$RET"

    echo "=================================================================="
    echo " chain run_${idx}  (study ${STUDY}, chrs ${CHRS}, versions ${VERSIONS})"
    echo "   pheno   : ${PHENO}"
    echo "   QC      : ${PREIMP_QC_ARGS}"
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
    # For each version, build an augmented phenotype (Y1 + PCs), run STAAR-O per chr
    # (Full also passes --grm), then merge all category CSVs into one eval table.
    for ver in "${VER_ARR[@]}"; do
        echo "[F/G] version=${ver}: augmented phenotype + STAAR-O + eval table ..."
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

        CSVS=()
        for chr in "${CHR_ARR[@]}"; do
            OUTSUB="staar_${ver}_chr${chr}"
            bcft -v "$WD":/w:z -v "$ARM":/arm:z,ro "$STAAR_IMG" \
                Rscript /arm/staar_gene_centric_coding.R \
                  --agds "/w/c${chr}.agds" --pheno "/w/pheno_${ver}.tsv" --chr "${chr}" \
                  --out-dir "/w/${OUTSUB}" --pheno-col Y1 --id-col IID \
                  --use-spa "${USE_SPA}" "${GRM_ARG[@]}" "${MAXG_ARG[@]}" \
                2>&1 | tee "${WD}/${OUTSUB}.log"
            for c in plof plof_ds missense disruptive_missense synonymous; do
                [[ -e "${WD}/${OUTSUB}/${c}.csv" ]] && CSVS+=("/w/${OUTSUB}/${c}.csv")
            done
        done
        [[ ${#CSVS[@]} -gt 0 ]] || { echo "ERROR: version=${ver} produced no STAAR CSVs" >&2; exit 1; }

        OUT_TABLE="${REGENIE_OUT_DIR}/${ver}/ricopili_staar_qcmatched_${ver}_dataset_idx_${idx}_step2_Y1.regenie"
        bcft -v "$WD":/w:z -v "$ARM":/arm:z,ro -v "$REGENIE_OUT_DIR":/out:z "$BCFTOOLS_IMG" \
            python3 /arm/staar_to_eval.py \
              --staar-results "${CSVS[@]}" \
              --genes-info /arm/genes_info_hgnc.tsv \
              --out "/out/${ver}/ricopili_staar_qcmatched_${ver}_dataset_idx_${idx}_step2_Y1.regenie"
        echo "[run_${idx}] ${ver} eval table -> ${OUT_TABLE}"
    done

    if [[ "$CLEANUP" == "true" ]]; then
        echo "[cleanup] rm -rf ${WD}  (retention + eval tables are kept elsewhere)"
        rm -rf "$WD"
    fi
    echo "[run_${idx}] done: $(date -Is)"
done

# ----------------------------------------------------------------------------
# Scoring (optional). One run_eval.sh call per version, over that version's result
# tables. CAUSAL_*_GLOB default to the known answers restricted to the numbered
# chromosomes (the output of filter_causal_autosomal.py), so causal genes on
# chromosome X count as misses. Point them elsewhere to score differently.
# Then compare with benchmark-common/pairwise_compare.py, passing
# --arm-b ricopili_staar_qcmatched_<version>.
# ----------------------------------------------------------------------------
if [[ "$SCORE" == "true" ]]; then
    : "${CAUSAL_SNPLIST_GLOB:?set CAUSAL_SNPLIST_GLOB (autosomal truth) to score}"
    : "${CAUSAL_GENES_GLOB:?set CAUSAL_GENES_GLOB (autosomal truth) to score}"
    for ver in "${VER_ARR[@]}"; do
        echo ""
        echo "=== scoring version=${ver} ==="
        export EVAL_REPO="/data/git/doktorat_pw/wum_pims/nf-eval-gene-assoc"
        export EVAL_RUN_DIR="${RUN_DIR}/eval_${ver}"
        export EVAL_PROJECT="ricopili_staar_qcmatched_${ver}"
        export EVAL_PROFILE="podman,medium_resources"
        export INPUT_VCF SKIP_PREP="true"
        export REGENIE_GLOB="${REGENIE_OUT_DIR}/${ver}/ricopili_staar_qcmatched_${ver}_dataset_idx_*_step2_Y1.regenie"
        export CAUSAL_SNPLIST_GLOB CAUSAL_GENES_GLOB
        bash "${COMMON}/run_eval.sh"
    done
fi

echo ""
echo "RICOPILI + STAARpipeline run finished: $(date -Is)"
