#!/usr/bin/env Rscript
#
# Relatedness-robust principal components (PC-AiR) and a relatedness matrix
# (PC-Relate) from a PLINK bed file, using GENESIS.
#
# Used by the "Full" version of the RICOPILI + STAARpipeline comparison, which keeps
# every sample and models relatedness rather than removing it: the PC-AiR components
# become fixed-effect covariates and the PC-Relate matrix becomes the null model's
# random effect. STAARpipeline needs a relatedness matrix and RICOPILI has no
# equivalent concept, so this is the one substitution that comparison makes; GENESIS
# is what STAARpipeline's own documentation prescribes for it.
#
# Computed on all available markers. The components and the matrix are genome-wide
# covariates, independent of which chromosomes the association test covers.
#
# Runs inside a container with SNPRelate, GWASTools and GENESIS. The default is
# uwgac/topmed-roybranch; the staarpipeline image carries GENESIS but not
# necessarily SNPRelate and GWASTools.
#
# OUTPUTS (into --out-dir):
#   pcair.tsv   IID PC1 PC2 ... PC<npcs>, tab-separated, all samples
#   grm.rds     the PC-Relate matrix, as a SPARSE Matrix indexed by sample id and
#               thresholded at --grm-thresh. It must be sparse: a dense matrix sends
#               STAAR down a code path that silently returns nothing for every gene.
#               The dimnames must be the sample ids, because
#               staar_gene_centric_coding.R subsets the matrix by them.
#
# USAGE
#   Rscript genesis_pcair_pcrelate.R \
#     --bed <prefix> --out-dir <dir> \
#     [--npcs 20] [--pcrelate-npcs 3] [--kin-thresh 0.04419417] [--grm-thresh 0.1767767] \
#     [--ld-threshold 0.31622777] [--maf 0.01] [--missing 0.01] [--seed 100]
#
# --bed is the PLINK prefix (expects <prefix>.bed/.bim/.fam).

suppressPackageStartupMessages({
    library(SNPRelate)
    library(GWASTools)
    library(GENESIS)
})

## ---------------------------------------------------------------------------
## args
## ---------------------------------------------------------------------------
parse_args <- function(argv) {
    opt <- list(bed = NULL, out_dir = NULL, npcs = 20L, pcrelate_npcs = 3L,
                grm_thresh = 2^(-5 / 2),          # ~0.1768 on the 2x-scaled GRM
                                                  #   = 2nd-degree kinship; see README.md
                kin_thresh = 2^(-9 / 2),          # ~0.0442: GENESIS 3rd-degree cut
                ld_threshold = sqrt(0.1),         # ~0.3162, the usual GENESIS default
                maf = 0.01, missing = 0.01, seed = 100L)
    i <- 1
    while (i <= length(argv)) {
        k <- argv[i]
        v <- if (i < length(argv)) argv[i + 1] else NA_character_
        switch(k,
            "--bed"           = { opt$bed <- v; i <- i + 2 },
            "--out-dir"       = { opt$out_dir <- v; i <- i + 2 },
            "--npcs"          = { opt$npcs <- as.integer(v); i <- i + 2 },
            "--pcrelate-npcs" = { opt$pcrelate_npcs <- as.integer(v); i <- i + 2 },
            "--kin-thresh"    = { opt$kin_thresh <- as.numeric(v); i <- i + 2 },
            "--grm-thresh"    = { opt$grm_thresh <- as.numeric(v); i <- i + 2 },
            "--ld-threshold"  = { opt$ld_threshold <- as.numeric(v); i <- i + 2 },
            "--maf"           = { opt$maf <- as.numeric(v); i <- i + 2 },
            "--missing"       = { opt$missing <- as.numeric(v); i <- i + 2 },
            "--seed"          = { opt$seed <- as.integer(v); i <- i + 2 },
            stop("unknown argument: ", k)
        )
    }
    for (req in c("bed", "out_dir")) {
        if (is.null(opt[[req]]) || is.na(opt[[req]])) {
            stop("missing required argument: --", gsub("_", "-", req))
        }
    }
    opt
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

bed <- paste0(opt$bed, ".bed")
bim <- paste0(opt$bed, ".bim")
fam <- paste0(opt$bed, ".fam")
for (f in c(bed, bim, fam)) if (!file.exists(f)) stop("missing bed component: ", f)

cat("=== GENESIS PC-AiR + PC-Relate ===\n")
cat("bed prefix   :", opt$bed, "\n")
cat("out dir      :", opt$out_dir, "\n")
cat("npcs (out)   :", opt$npcs, "\n")
cat("pcrelate pcs :", opt$pcrelate_npcs, "\n")
cat("kin.thresh   :", opt$kin_thresh, "\n")

## ---------------------------------------------------------------------------
## 1. bed -> SNP GDS
## ---------------------------------------------------------------------------
gdsfile <- file.path(opt$out_dir, "geno.gds")
cat("\n[1/5] BED -> GDS ...\n")
snpgdsBED2GDS(bed, fam, bim, gdsfile, verbose = FALSE)
gds <- snpgdsOpen(gdsfile)
on.exit(try(snpgdsClose(gds), silent = TRUE), add = TRUE)

## ---------------------------------------------------------------------------
## 2. LD pruning for the structure markers (KING + PC-AiR run on the pruned set)
## ---------------------------------------------------------------------------
cat("[2/5] LD pruning (maf=", opt$maf, " missing<=", opt$missing,
    " ld<=", round(opt$ld_threshold, 4), ") ...\n", sep = "")
set.seed(opt$seed)
snpset <- snpgdsLDpruning(gds, method = "corr", maf = opt$maf,
                          missing.rate = opt$missing,
                          ld.threshold = opt$ld_threshold,
                          slide.max.bp = 10e6, verbose = FALSE)
pruned <- unlist(snpset, use.names = FALSE)
cat("      pruned markers:", length(pruned), "\n")
if (length(pruned) < 2) stop("LD pruning left < 2 markers -- cannot compute structure")

## ---------------------------------------------------------------------------
## 3. KING-robust kinship (seeds PC-AiR's relatedness + divergence)
## ---------------------------------------------------------------------------
cat("[3/5] KING-robust kinship ...\n")
king <- snpgdsIBDKING(gds, snp.id = pruned, verbose = FALSE)
kingMat <- king$kinship
sample_ids <- king$sample.id
dimnames(kingMat) <- list(sample_ids, sample_ids)
snpgdsClose(gds)
on.exit()   # clear the snpgdsClose handler; gds is now closed

## ---------------------------------------------------------------------------
## 4. PC-AiR (relatedness-robust PCs over all samples)
## ---------------------------------------------------------------------------
cat("[4/5] PC-AiR ...\n")
geno <- GdsGenotypeReader(gdsfile)
genoData <- GenotypeData(geno)
on.exit(try(close(genoData), silent = TRUE), add = TRUE)

pcair_res <- pcair(genoData, kinobj = kingMat, divobj = kingMat,
                   kin.thresh = opt$kin_thresh, div.thresh = -opt$kin_thresh,
                   snp.include = pruned)
n_pcs <- min(opt$npcs, ncol(pcair_res$vectors))
pcs <- pcair_res$vectors[, seq_len(n_pcs), drop = FALSE]
colnames(pcs) <- paste0("PC", seq_len(n_pcs))
cat("      samples:", nrow(pcs), " unrelated set:", length(pcair_res$unrels),
    " related set:", length(pcair_res$rels), "\n")

pcair_df <- data.frame(IID = rownames(pcs), pcs, check.names = FALSE,
                       stringsAsFactors = FALSE)
pcair_path <- file.path(opt$out_dir, "pcair.tsv")
write.table(pcair_df, pcair_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("      wrote", pcair_path, "\n")

## ---------------------------------------------------------------------------
## 5. PC-Relate GRM (ancestry-adjusted kinship, all samples)
##    Adjusted for the top pcrelate-npcs PC-AiR PCs; trained on the unrelated set.
## ---------------------------------------------------------------------------
cat("[5/5] PC-Relate ...\n")
n_pcr <- min(opt$pcrelate_npcs, n_pcs)
genoIterator <- GenotypeBlockIterator(genoData, snpInclude = pruned)
pcrel <- pcrelate(genoIterator,
                  pcs = pcs[, seq_len(n_pcr), drop = FALSE],
                  training.set = pcair_res$unrels)
## thresh zeroes kinship below the cutoff, which is what makes the GRM SPARSE --
## and sparse is not an optimisation here, it is the only working path: STAAR's
## fit_nullmodel sends a dense GRM to glmmkin's dense branch, whose Matrix-class
## projection matrix STAAR's C++ side then rejects, silently returning no results
## for every gene.
## NOTE ON SCALE: pcrelateToMatrix applies thresh to the values it returns, i.e.
## AFTER scaleKin -- so these numbers are 2 x kinship. 2nd-degree kinship (0.0884)
## is thresh 0.1768 = 2^-2.5; the conventional 4th-degree cut (0.0221) would be
## 2^-4.5. Getting this backwards is easy and fails loudly (the run below).
## We cut at 2nd degree rather than the conventional 4th because on an exome the
## conventional value does not sparsify anything. Measured on run_18 (3,149
## LD-pruned markers): at 4th degree the relatedness graph is fully connected --
## one block of all 1,356 samples, i.e. still dense -- and at 3rd degree one block
## of 1,346. Only at 2nd degree does it break up (275 clusters, largest 12, 534
## unrelated singletons), and that structure is the fixture's real one: the pairs
## kept are essentially just its ~600 trios (571 pairs at 2nd degree vs 544 at
## 1st). Too few independent exome markers make weak kinship estimates too noisy
## to threshold at the usual place -- the same marker starvation that collapses
## RICOPILI's relatedness filter, but far less damaging, since every sample is
## kept either way.
grm <- pcrelateToMatrix(pcrel, thresh = opt$grm_thresh, scaleKin = 2)  # 2*kinship, diag ~1
# order rows/cols to the PC-AiR sample order for a stable, id-named matrix
grm <- grm[rownames(pcs), rownames(pcs)]
if (!inherits(grm, "sparseMatrix")) {
    stop("GRM is not sparse at thresh=", opt$grm_thresh,
         " -- STAAR's dense path does not work (see the comment above); ",
         "raise --grm-thresh")
}
grm_path <- file.path(opt$out_dir, "grm.rds")
saveRDS(grm, grm_path)
n_rel <- length(grm@x) - nrow(grm)   # dsCMatrix stores one triangle + diagonal
cat("      GRM:", nrow(grm), "x", ncol(grm), " class ", class(grm)[1],
    ", related pairs kept: ", n_rel, " -> ", grm_path, "\n", sep = "")

cat("\nGENESIS structure step done.\n")
