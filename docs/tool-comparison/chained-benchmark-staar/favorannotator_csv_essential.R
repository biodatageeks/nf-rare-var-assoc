#############################################################################
# favorannotator_csv_essential.R -- GDS -> aGDS using the FAVOR Essential
# Database (CSV variant).
#
# Derived from zhouhufeng/FAVORannotator @ d65f565
# Scripts/CSV/FAVORannotatorCSVEssentialDB.R. The annotation logic is upstream's
# and is kept unchanged: the same three stages (VarInfo split by DB boundary ->
# xsv left-join per split -> cat rows -> select the STAARpipeline column subset
# -> add.gdsn as annotation/info/FunctionalAnnotation), the same column subset
# `anno_colnum`, the same LZMA_ra compression, and the same three renamed
# columns. What we changed, and why:
#
#   * DB_path, output_path, and the xsv binary are arguments/env, not constants
#     pointing at the authors' cluster ("n/holystore01/LABS/xlin/...").
#   * Upstream's step 0 (wget the chromosome tarball from Dataverse into the CWD
#     and untar it, on every run) is REMOVED. We pre-download and verify the
#     tarball by md5 once, out of band, and mount the extracted CSVs read-only.
#     Same bytes, but the run is offline, re-runnable, and does not re-fetch
#     5-18 GiB per invocation.
#   * FAVORdatabase_chrsplit.csv is read from a local path (vendored into the
#     image) rather than raw.githubusercontent.com at runtime.
#   * Added preflight checks: every expected chr<N>_<k>.csv split must exist
#     before any work starts, and the join output is checked for the annotation
#     channels STAAR weights by.
#
# Usage:
#   Rscript favorannotator_csv_essential.R <in.gds> <chr> <db_path> [out_path]
# The input GDS is modified IN PLACE into an aGDS (upstream behaviour: it opens
# readonly=FALSE and adds a node). Work on a copy.
#############################################################################
suppressPackageStartupMessages({
    library(gdsfmt)
    library(SeqArray)
    library(readr)
})

args <- commandArgs(TRUE)
if (length(args) < 3) {
    stop("usage: favorannotator_csv_essential.R <in.gds> <chr> <db_path> [out_path]")
}
gds.file <- args[1]
chr <- as.numeric(args[2])
DB_path <- sub("/+$", "", args[3])
output_path <- if (length(args) >= 4) args[4] else "./"
if (!grepl("/$", output_path)) output_path <- paste0(output_path, "/")

xsv <- Sys.getenv("XSV", unset = "/usr/local/bin/xsv")
chrsplit_file <- Sys.getenv("FAVOR_CHRSPLIT",
                            unset = "/opt/favorannotator/FAVORdatabase_chrsplit.csv")

cat("gds.file   :", gds.file, "\n")
cat("chr        :", chr, "\n")
cat("DB_path    :", DB_path, "\n")
cat("output_path:", output_path, "\n")
cat("xsv        :", xsv, "\n")

stopifnot(file.exists(gds.file), file.exists(chrsplit_file), file.exists(xsv))

### annotation file naming (upstream)
anno_file_name_1 <- "Anno_chr"
anno_file_name_2 <- "_STAARpipeline.csv"

## read DB split info
DB_info <- read.csv(chrsplit_file, header = TRUE)
DB_info_chr <- DB_info[DB_info$Chr == chr, ]
chr_splitnum <- sum(DB_info$Chr == chr)
if (chr_splitnum == 0) stop("no DB split rows for chr ", chr)
cat("DB splits  :", chr_splitnum, "\n")

## preflight: every split CSV must be present (upstream would fail mid-join)
split_files <- file.path(DB_path, paste0("chr", chr, "_", seq_len(chr_splitnum), ".csv"))
missing <- split_files[!file.exists(split_files)]
if (length(missing) > 0) {
    stop("missing FAVOR DB split file(s):\n  ", paste(missing, collapse = "\n  "))
}

start_time <- Sys.time()
dir.create(paste0(output_path, "chr", chr), showWarnings = FALSE, recursive = TRUE)

##########################################################################
### Step 1 (Varinfo_gds)
##########################################################################
genofile <- seqOpen(gds.file, readonly = FALSE)
on.exit(try(seqClose(genofile), silent = TRUE), add = TRUE)

CHR <- as.numeric(seqGetData(genofile, "chromosome"))
position <- as.integer(seqGetData(genofile, "position"))
REF <- as.character(seqGetData(genofile, "$ref"))
ALT <- as.character(seqGetData(genofile, "$alt"))

VarInfo_genome <- paste0(CHR, "-", position, "-", REF, "-", ALT)
cat("variants in GDS:", length(VarInfo_genome), "\n")

for (kk in seq_len(nrow(DB_info_chr))) {
    VarInfo <- VarInfo_genome[(position >= DB_info_chr$Start_Pos[kk]) &
                              (position <= DB_info_chr$End_Pos[kk])]
    VarInfo <- data.frame(VarInfo)
    cat("  split", kk, "variants:", nrow(VarInfo), "\n")
    write.csv(VarInfo,
              paste0(output_path, "chr", chr, "/VarInfo_chr", chr, "_", kk, ".csv"),
              quote = FALSE, row.names = FALSE)
}

##########################################################################
### Step 2 (Annotate) -- xsv left-join against each DB split, then merge
##########################################################################
### anno channel (subset) -- upstream's column selection, unchanged
anno_colnum <- c(1, 8:12, 15, 16, 19, 23, 25:36)

run <- function(cmd) {
    cat("+ ", cmd, "\n", sep = "")
    st <- system(cmd)
    if (st != 0) stop("command failed (status ", st, "): ", cmd)
}

for (kk in seq_len(chr_splitnum)) {
    run(paste0(xsv, " join --left VarInfo ",
               output_path, "chr", chr, "/VarInfo_chr", chr, "_", kk, ".csv",
               " variant_vcf ", DB_path, "/chr", chr, "_", kk, ".csv > ",
               output_path, "chr", chr, "/Anno_chr", chr, "_", kk, ".csv"))
}

## merge info
Anno <- paste0(output_path, "chr", chr, "/Anno_chr", chr, "_", seq_len(chr_splitnum), ".csv ")
merge_command <- paste0(xsv, " cat rows ", paste(Anno, collapse = ""))
merge_command <- paste0(merge_command, "> ", output_path, "chr", chr, "/Anno_chr", chr, ".csv")
run(merge_command)

## subset to the STAARpipeline columns
anno_colnum_xsv <- paste(anno_colnum, collapse = ",")
run(paste0(xsv, " select ", anno_colnum_xsv, " ",
           output_path, "chr", chr, "/Anno_chr", chr, ".csv > ",
           output_path, "chr", chr, "/", anno_file_name_1, chr, anno_file_name_2))

##########################################################################
### Step 3 (gds2agds)
##########################################################################
FunctionalAnnotation <- read_csv(
    paste0(output_path, "chr", chr, "/", anno_file_name_1, chr, anno_file_name_2),
    col_types = list(col_character(), col_double(), col_double(), col_double(), col_double(),
                     col_double(), col_double(), col_double(), col_double(), col_double(),
                     col_character(), col_character(), col_character(), col_double(), col_character(),
                     col_character(), col_character(), col_character(), col_character(), col_double(),
                     col_double(), col_character()))

cat("annotation table dim:", dim(FunctionalAnnotation), "\n")

## rename colnames (upstream)
colnames(FunctionalAnnotation)[2] <- "apc_conservation"
colnames(FunctionalAnnotation)[7] <- "apc_local_nucleotide_diversity"
colnames(FunctionalAnnotation)[9] <- "apc_protein_function"

cat("annotation channels:\n")
print(colnames(FunctionalAnnotation))

## check that the channels STAAR's Use_annotation_weights=TRUE reads are present
required <- c("apc_conservation", "apc_local_nucleotide_diversity", "apc_protein_function")
absent <- required[!required %in% colnames(FunctionalAnnotation)]
if (length(absent) > 0) {
    warning("aGDS is MISSING STAAR weighting channels: ", paste(absent, collapse = ", "))
}
nonmissing <- vapply(FunctionalAnnotation, function(x) sum(!is.na(x)), numeric(1))
cat("non-missing values per channel:\n")
print(nonmissing)
if (any(nonmissing == 0)) {
    warning("all-NA annotation channel(s): ",
            paste(names(nonmissing)[nonmissing == 0], collapse = ", "))
}

Anno.folder <- index.gdsn(genofile, "annotation/info")
add.gdsn(Anno.folder, "FunctionalAnnotation", val = FunctionalAnnotation,
         compress = "LZMA_ra", closezip = TRUE)

print(genofile)
seqClose(genofile)

cat("aGDS written:", gds.file, "\n")
cat("elapsed:", format(Sys.time() - start_time), "\n")
