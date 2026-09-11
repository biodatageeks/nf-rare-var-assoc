#!/usr/bin/env Rscript
#
# STAARpipeline gene-centric coding analysis over a FAVORannotator aGDS.
#
# The association step of the RICOPILI + STAARpipeline comparison: STAAR's own
# engine, over STAAR's own
# annotation, with nothing borrowed from nf-rare-var-assoc except the raw
# genotypes, the phenotype, and (later) the QC/PC covariates from the QC half.
#
# Runs inside docker.io/zilinli/staarpipeline:0.9.7
# (STAAR/STAARpipeline 0.9.7.2, SeqArray 1.26.2, GMMAT 1.3.2, GENESIS 2.16.1).
#
# WHAT IT DOES
#   1. fits the STAAR null model on the phenotype (no GRM by default -- see
#      --grm; the chain's step-3 design substitutes a GENESIS PC-Relate GRM)
#   2. runs Gene_Centric_Coding(category="all_categories") gene by gene over
#      every gene STAAR ships for the chromosome (STAARpipeline:::genes_info)
#   3. writes one CSV per coding category, shaped for staar_to_eval.py
#
# TWO ENGINE BEHAVIOURS THIS SCRIPT PINS DOWN:
#
#   * use_SPA=FALSE emits columns  ... ACAT-O, STAAR-O   (full omnibus: burden+SKAT)
#     use_SPA=TRUE  emits column   ... STAAR-B           (burden-only omnibus)
#     They are DIFFERENT TESTS, not a calibration tweak. STAAR-O is the
#     like-for-like analogue of REGENIE's SKAT-O; STAAR-B drops the SKAT
#     component. --use-spa selects; the output column is reported either way.
#
#   * Use_annotation_weights applies ONLY when variant_type=="SNV"
#     (STAARpipeline R/coding.R:90-92). With variant_type="variant" the
#     annotation weights are silently NOT applied to ANY variant. So including
#     indels costs the annotation weighting entirely -- there is no
#     "weighted SNVs + unweighted indels" mode. Default here is SNV, which is
#     canonical STAARpipeline usage and keeps aPC weighting live.
#
# USAGE
#   Rscript staar_gene_centric_coding.R \
#     --agds /work/c22.agds --pheno /work/pheno.tsv --chr 22 \
#     --out-dir /work/staar_out [--pheno-col Y1] [--id-col IID] \
#     [--use-spa TRUE|FALSE] [--variant-type SNV|Indel|variant] \
#     [--rare-maf 0.01] [--max-genes N] [--grm <rds>] [--kins-cutoff 0.022]
#
# Covariates: every phenotype column that is not the id column, "FID", or the
# phenotype column is used as a fixed-effect covariate. With none present the
# model is <pheno> ~ 1.

suppressPackageStartupMessages({
    library(gdsfmt)
    library(SeqArray)
    library(STAAR)
    library(STAARpipeline)
})

## ---------------------------------------------------------------------------
## args
## ---------------------------------------------------------------------------
parse_args <- function(argv) {
    opt <- list(agds = NULL, pheno = NULL, chr = NULL, out_dir = NULL,
                pheno_col = NULL, id_col = "IID", use_spa = FALSE,
                variant_type = "SNV", rare_maf = 0.01, max_genes = NA_integer_,
                grm = NULL, kins_cutoff = 0.022)
    i <- 1
    while (i <= length(argv)) {
        k <- argv[i]
        v <- if (i < length(argv)) argv[i + 1] else NA_character_
        switch(k,
            "--agds"         = { opt$agds <- v; i <- i + 2 },
            "--pheno"        = { opt$pheno <- v; i <- i + 2 },
            "--chr"          = { opt$chr <- as.integer(v); i <- i + 2 },
            "--out-dir"      = { opt$out_dir <- v; i <- i + 2 },
            "--pheno-col"    = { opt$pheno_col <- v; i <- i + 2 },
            "--id-col"       = { opt$id_col <- v; i <- i + 2 },
            "--use-spa"      = { opt$use_spa <- toupper(v) %in% c("TRUE","T","YES","1"); i <- i + 2 },
            "--variant-type" = { opt$variant_type <- v; i <- i + 2 },
            "--rare-maf"     = { opt$rare_maf <- as.numeric(v); i <- i + 2 },
            "--max-genes"    = { opt$max_genes <- as.integer(v); i <- i + 2 },
            "--grm"          = { opt$grm <- v; i <- i + 2 },
            "--kins-cutoff"  = { opt$kins_cutoff <- as.numeric(v); i <- i + 2 },
            stop("unknown argument: ", k)
        )
    }
    for (req in c("agds", "pheno", "chr", "out_dir")) {
        if (is.null(opt[[req]]) || is.na(opt[[req]])) {
            stop("missing required argument: --", gsub("_", "-", req))
        }
    }
    opt
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== STAAR gene-centric coding ===\n")
cat("aGDS        :", opt$agds, "\n")
cat("phenotype   :", opt$pheno, "\n")
cat("chromosome  :", opt$chr, "\n")
cat("use_SPA     :", opt$use_spa, "\n")
cat("variant_type:", opt$variant_type, "\n")
cat("rare_maf    :", opt$rare_maf, "\n")

if (opt$variant_type != "SNV") {
    cat("\nWARNING: variant_type != SNV -> STAAR applies NO annotation weights\n",
        "         to any variant (coding.R:90-92). The aPC weighting that is\n",
        "         STAAR's whole premise is inactive in this configuration.\n\n", sep = "")
}

## ---------------------------------------------------------------------------
## annotation catalog: STAAR's expected names -> FAVORannotator channel nodes
##
## STAAR reads seqGetData(genofile, paste0(Annotation_dir, dir)), with
## Annotation_dir = "annotation/info/FunctionalAnnotation", so every dir needs a
## leading slash. The dir values below are the actual child nodes written by
## favorannotator_csv_essential.R (verified against the built aGDS).
##
## "CADD" and "aPC.LocalDiversity" are special-cased by name inside STAAR
## (coding.R:101,106 -- NA->0, and a reversed companion channel) so those two
## strings must be spelled exactly as here.
## ---------------------------------------------------------------------------
Annotation_name_catalog <- data.frame(
    name = c("rs_num", "GENCODE.Category", "GENCODE.Info",
             "GENCODE.EXONIC.Category", "GENCODE.EXONIC.Info", "MetaSVM",
             "GeneHancer", "CAGE", "DHS", "CADD", "LINSIGHT", "FATHMM.XF",
             "aPC.EpigeneticsActive", "aPC.EpigeneticsRepressed",
             "aPC.EpigeneticsTranscription", "aPC.Conservation",
             "aPC.LocalDiversity", "aPC.Mappability",
             "aPC.TranscriptionFactor", "aPC.ProteinFunction"),
    dir = c("/rsid", "/genecode_comprehensive_category",
            "/genecode_comprehensive_info",
            "/genecode_comprehensive_exonic_category",
            "/genecode_comprehensive_exonic_info", "/metasvm_pred",
            "/genehancer", "/cage_tc", "/rdhs", "/cadd_phred", "/linsight",
            "/fathmm_xf", "/apc_epigenetics_active",
            "/apc_epigenetics_repressed", "/apc_epigenetics_transcription",
            "/apc_conservation", "/apc_local_nucleotide_diversity",
            "/apc_mappability", "/apc_transcription_factor",
            "/apc_protein_function"),
    stringsAsFactors = FALSE
)

## the continuous channels used as STAAR weights (aPC.LocalDiversity implies a
## second, reversed channel that STAAR generates itself)
Annotation_name <- c("CADD", "LINSIGHT", "FATHMM.XF",
                     "aPC.EpigeneticsActive", "aPC.EpigeneticsRepressed",
                     "aPC.EpigeneticsTranscription", "aPC.Conservation",
                     "aPC.LocalDiversity", "aPC.Mappability",
                     "aPC.TranscriptionFactor", "aPC.ProteinFunction")

## ---------------------------------------------------------------------------
## open aGDS, verify the catalog resolves against it
## ---------------------------------------------------------------------------
genofile <- seqOpen(opt$agds)
on.exit(try(seqClose(genofile), silent = TRUE), add = TRUE)

anno_root <- "annotation/info/FunctionalAnnotation"
present <- ls.gdsn(index.gdsn(genofile, anno_root))
wanted <- sub("^/", "", Annotation_name_catalog$dir)
absent <- setdiff(wanted, present)
if (length(absent) > 0) {
    stop("aGDS is missing annotation channels the catalog names: ",
         paste(absent, collapse = ", "))
}
cat("annotation channels resolved:", length(wanted), "of", length(wanted), "\n")

gds_samples <- seqGetData(genofile, "sample.id")
cat("aGDS samples:", length(gds_samples),
    " variants:", length(seqGetData(genofile, "variant.id")), "\n")

## ---------------------------------------------------------------------------
## phenotype
## ---------------------------------------------------------------------------
pheno <- read.table(opt$pheno, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE)
if (!opt$id_col %in% names(pheno)) {
    stop("phenotype lacks id column '", opt$id_col, "'; have: ",
         paste(names(pheno), collapse = ", "))
}
if (is.null(opt$pheno_col)) {
    cand <- setdiff(names(pheno), c("FID", opt$id_col))
    if (length(cand) < 1) stop("cannot infer phenotype column")
    opt$pheno_col <- cand[1]
    cat("inferred phenotype column:", opt$pheno_col, "\n")
}

pheno[[opt$id_col]] <- as.character(pheno[[opt$id_col]])
pheno <- pheno[pheno[[opt$id_col]] %in% gds_samples, , drop = FALSE]
pheno <- pheno[!is.na(pheno[[opt$pheno_col]]), , drop = FALSE]
if (nrow(pheno) == 0) stop("no phenotyped samples overlap the aGDS")

covars <- setdiff(names(pheno), c("FID", opt$id_col, opt$pheno_col))
fml <- if (length(covars) > 0) {
    as.formula(paste(opt$pheno_col, "~", paste(covars, collapse = " + ")))
} else {
    as.formula(paste(opt$pheno_col, "~ 1"))
}

yv <- pheno[[opt$pheno_col]]
is_binary <- all(yv %in% c(0, 1))
fam <- if (is_binary) binomial(link = "logit") else gaussian(link = "identity")

cat("samples analysed:", nrow(pheno), "\n")
cat("model           :", deparse(fml), "\n")
cat("family          :", if (is_binary) "binomial" else "gaussian", "\n")
if (is_binary) {
    cat("cases/controls  :", sum(yv == 1), "/", sum(yv == 0),
        sprintf(" (%.1f%% cases)\n", 100 * mean(yv == 1)))
}

## ---------------------------------------------------------------------------
## null model
## ---------------------------------------------------------------------------
kins <- NULL
if (!is.null(opt$grm)) {
    cat("loading GRM:", opt$grm, "\n")
    kins <- readRDS(opt$grm)
    ids <- pheno[[opt$id_col]]
    kins <- kins[ids, ids, drop = FALSE]
}

cat("\nfitting null model", if (is.null(kins)) "(no GRM -> glm)" else "(GRM, sparse)", "...\n")
t0 <- Sys.time()
## use_sparse = TRUE is REQUIRED, not a tuning choice. fit_nullmodel defaults it
## to NULL, which sends a dense GRM down glmmkin's dense path; STAAR then calls
## STAAR_O_SMMAT with a Matrix S4 projection matrix and its C++ side rejects it
## ("Not a matrix."). Every gene comes back empty, and because coding() wraps the
## call in try() the failure is SILENT -- a whole run of NULL results with zero
## reported errors. TRUE takes glmmkin's sparse path (kinship thresholded at
## kins_cutoff), which is also what STAARpipeline's own tutorial does for related
## samples. Ignored when kins is NULL.
obj_nullmodel <- fit_nullmodel(
    fml, data = pheno, kins = kins, id = opt$id_col,
    use_sparse = TRUE, kins_cutoff = opt$kins_cutoff,
    use_SPA = opt$use_spa, family = fam, verbose = FALSE
)
cat("null model fitted in",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")

## ---------------------------------------------------------------------------
## gene loop
## ---------------------------------------------------------------------------
genes_info <- STAARpipeline:::genes_info
genes_chr <- genes_info[genes_info[, 2] == opt$chr, ]
gene_names <- as.character(genes_chr[, 1])
if (!is.na(opt$max_genes) && opt$max_genes < length(gene_names)) {
    gene_names <- gene_names[seq_len(opt$max_genes)]
}
cat("\ngenes on chr", opt$chr, ":", length(gene_names), "\n", sep = "")

## STAAR's coding() builds its result tables with rbind() over mixed-type
## vectors, so each is a MATRIX OF MODE LIST (character gene name alongside
## numeric p-values). as.data.frame() on that yields list columns, which
## write.csv cannot encode. Unlist column by column instead: that keeps each
## column's native type (numeric stays numeric, so no precision is lost
## round-tripping p-values through character).
mat_to_df <- function(m) {
    if (is.data.frame(m)) return(m)
    m <- as.matrix(m)
    cols <- lapply(seq_len(ncol(m)), function(j) {
        v <- m[, j]
        if (is.list(v)) v <- unlist(v, use.names = FALSE)
        v
    })
    names(cols) <- colnames(m)
    as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
}

categories <- c("plof", "plof_ds", "missense", "disruptive_missense", "synonymous")
acc <- setNames(vector("list", length(categories)), categories)
n_err <- 0L
errors <- character(0)

t0 <- Sys.time()
for (gi in seq_along(gene_names)) {
    g <- gene_names[gi]
    res <- tryCatch(
        Gene_Centric_Coding(
            chr = opt$chr, gene_name = g, category = "all_categories",
            genofile = genofile, obj_nullmodel = obj_nullmodel,
            rare_maf_cutoff = opt$rare_maf,
            QC_label = "annotation/filter",
            variant_type = opt$variant_type,
            geno_missing_imputation = "mean",
            Annotation_dir = anno_root,
            Annotation_name_catalog = Annotation_name_catalog,
            Use_annotation_weights = TRUE,
            Annotation_name = Annotation_name,
            silent = TRUE
        ),
        error = function(e) {
            n_err <<- n_err + 1L
            if (length(errors) < 10) errors <<- c(errors, paste0(g, ": ", conditionMessage(e)))
            NULL
        }
    )
    if (!is.null(res)) {
        for (cat_nm in categories) {
            r <- res[[cat_nm]]
            if (!is.null(r) && NROW(r) > 0) {
                acc[[cat_nm]][[length(acc[[cat_nm]]) + 1L]] <- mat_to_df(r)
            }
        }
    }
    if (gi %% 50 == 0 || gi == length(gene_names)) {
        cat(sprintf("  %d/%d genes (%.0fs elapsed, %d errors)\n", gi, length(gene_names),
                    as.numeric(difftime(Sys.time(), t0, units = "secs")), n_err))
    }
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("gene loop finished in", round(elapsed, 1), "s;", n_err, "errors\n")
if (length(errors) > 0) {
    cat("first errors:\n"); for (e in errors) cat("  ", e, "\n")
}

## ---------------------------------------------------------------------------
## write results + diagnostics
## ---------------------------------------------------------------------------
p_col <- if (opt$use_spa) "STAAR-B" else "STAAR-O"
cat("\n=== results (p-value column: ", p_col, ") ===\n", sep = "")

all_p <- numeric(0)
n_rows_total <- 0L
for (cat_nm in categories) {
    parts <- acc[[cat_nm]]
    out_path <- file.path(opt$out_dir, paste0(cat_nm, ".csv"))
    if (length(parts) == 0) {
        cat(sprintf("  %-20s 0 genes (no testable genes)\n", cat_nm))
        next
    }
    df <- do.call(rbind, parts)
    write.csv(df, out_path, row.names = FALSE, quote = TRUE)
    n_rows_total <- n_rows_total + nrow(df)
    if (p_col %in% names(df)) {
        pv <- suppressWarnings(as.numeric(df[[p_col]]))
        all_p <- c(all_p, pv[!is.na(pv)])
        cat(sprintf("  %-20s %4d genes  p: min=%.3g med=%.3g max=%.3g  NA=%d\n",
                    cat_nm, nrow(df), min(pv, na.rm = TRUE),
                    median(pv, na.rm = TRUE), max(pv, na.rm = TRUE), sum(is.na(pv))))
    } else {
        cat(sprintf("  %-20s %4d genes  (no '%s' column; have: %s)\n",
                    cat_nm, nrow(df), p_col, paste(names(df), collapse = ",")))
    }
}

cat("\n=== check that the p-values are not degenerate ===\n")
cat("total result rows :", n_rows_total, "\n")
cat("non-NA p-values   :", length(all_p), "\n")
if (length(all_p) > 0) {
    cat("distinct p-values :", length(unique(all_p)), "\n")
    cat("range             :", sprintf("%.3g .. %.3g", min(all_p), max(all_p)), "\n")
    cat("quantiles         :",
        paste(sprintf("%.3g", quantile(all_p, c(0, .25, .5, .75, 1))), collapse = "  "), "\n")
    cat("p < 0.05          :", sum(all_p < 0.05), "\n")
    cat("p < 0.001         :", sum(all_p < 0.001), "\n")
    degenerate <- length(unique(all_p)) < 10 || all(all_p == 1) || all(all_p == 0)
    cat("VERDICT           :", if (degenerate) "DEGENERATE" else "NON-DEGENERATE", "\n")
} else {
    cat("VERDICT           : DEGENERATE (no p-values produced)\n")
}
cat("\nwrote per-category CSVs to", opt$out_dir, "\n")
