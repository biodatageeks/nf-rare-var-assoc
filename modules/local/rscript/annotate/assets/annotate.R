library("data.table")
library("dplyr")
library("optparse")

# Define command line options
option_list <- list(
  make_option(c("--fam-path"), type="character", default=NULL,
              help="Path to input fam file", metavar="character"),
  make_option(c("--controls-path"), type="character", default=NULL,
              help="Path to input controls file", metavar="character"),
  make_option(c("--cases-paths"), type="character", default=NULL, action="append",
              help="Path to a cases file (can be specified multiple times)", metavar="character"),
  make_option(c("--vcf-path"), type="character", default=NULL,
              help="Path to input VCF file", metavar="character"),
  make_option(c("--sample-path"), type="character", default=NULL,
              help="Path to input sample file", metavar="character"),
  make_option(c("--masks-path"), type="character", default=NULL,
              help="Path to input masks file", metavar="character"),
  
  make_option(c("--out-fam-path"), type="character", default=NULL,
              help="Path to output fam file", metavar="character"),
  make_option(c("--out-sample-path"), type="character", default=NULL,
              help="Path to output sample file", metavar="character"),
  make_option(c("--out-pheno-path"), type="character", default=NULL,
              help="Path to output phenotype file", metavar="character"),
  make_option(c("--out-anno-path"), type="character", default=NULL,
              help="Path to output annotations file", metavar="character"),
  make_option(c("--out-setlist-path"), type="character", default=NULL,
              help="Path to output setlist file", metavar="character"),
  make_option(c("--min_top_annotations"), type="integer", default=50,
              help="Filter annotations threshold - keep at least that many [default=50]", metavar="integer"),
  make_option(c("--quantile_threshold"), type="double", default=0.25,
              help="Filter annotations threshold - keep based on the quantile of the frequency table [default=0.25]", metavar="double"),
  make_option(c("--include-intergenic"), type="logical", default=FALSE,
              help="Include intergenic variants [default=FALSE]", metavar="logical")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check for required arguments and assign to variables
if (is.null(opt$`fam-path`)) stop("Input fam file path is required")
if (is.null(opt$`controls-path`)) stop("Input controls file path is required")
if (is.null(opt$`cases-paths`)) stop("At least one cases file is required")
if (is.null(opt$`vcf-path`)) stop("Input VCF file path is required")
if (is.null(opt$`sample-path`)) stop("Input sample file path is required")
if (is.null(opt$`masks-path`)) stop("Input masks file path is required")

if (is.null(opt$`out-fam-path`)) stop("Output fam file path is required")
if (is.null(opt$`out-sample-path`)) stop("Output sample file path is required")
if (is.null(opt$`out-pheno-path`)) stop("Output phenotype file path is required")
if (is.null(opt$`out-anno-path`)) stop("Output annotations file path is required")
if (is.null(opt$`out-setlist-path`)) stop("Output setlist file path is required")

r_in_fam_path <- opt$`fam-path`
r_in_controls_path <- opt$`controls-path`
r_in_cases_paths <- opt$`cases-paths`    # List of cases file paths
r_in_vcf_path <- opt$`vcf-path`
r_in_sample_path <- opt$`sample-path`
r_in_masks_path <- opt$`masks-path`

r_out_fam_path <- opt$`out-fam-path`
r_out_sample_path <- opt$`out-sample-path`
r_out_phenotype_path <- opt$`out-pheno-path`
r_out_annotations_path <- opt$`out-anno-path`
r_out_setlist_path <- opt$`out-setlist-path`

include_intergenic <- opt$`include-intergenic`
quantile_threshold <- opt$`quantile_threshold`  # Keep consequences above Nth percentile (adjustable)
min_top_annotations <- opt$`min_top_annotations`   # Ensure at least N top annotations are kept, if available

cat(paste0("quantile_threshold = ", quantile_threshold, "\n"))
cat(paste0("min_top_annotations = ", min_top_annotations, "\n"))
cat(paste0("include_intergenic = ", include_intergenic, "\n"))
cat("Number of cases files:", length(r_in_cases_paths), "\n")


# supporting multiple phenotypes, work in progress, temporarily turned off
#if (FALSE) {
# Process fam and phenotype files
fam <- fread(r_in_fam_path, header=F)
controls <- fread(r_in_controls_path, header=F)

# Clean IDs by removing ".hg38" suffix if present
fam_ids <- gsub("\\.hg38", "", fam$V1)

# Initialize phenotype table with controls (default NA for phenotypes)
pheno <- fam[, .(V1, V2)]
setnames(pheno, c("FID", "IID"))
pheno[, (paste0("Y", seq_along(r_in_cases_paths))) := NA_real_]  # Initialize with NA

# Process each phenotype file
for (i in seq_along(r_in_cases_paths)) {
    pheno_file <- fread(r_in_cases_paths[[i]], header=FALSE)
    pheno_name <- paste0("Y", i)

    setnames(pheno_file, "V1", "ID")
    if (ncol(pheno_file) > 1) setnames(pheno_file, "V2", "VALUE")

    # Check if second column exists and is numeric
    if (ncol(pheno_file) == 1 || !all(grepl("^[0-9.-]+$", pheno_file$VALUE[!is.na(pheno_file$VALUE)]))) {
        # Treat as binary phenotype (1 for presence, 0 for controls, NA otherwise)
        cat(sprintf("Treating %s as binary phenotype (no numeric second column)\n", r_in_cases_paths[[i]]))
        pheno[FID %in% pheno_file$ID, (pheno_name) := 1]
    } else {
        # Treat as numeric phenotype (continuous or categorical)
        cat(sprintf("Treating %s as numeric phenotype\n", r_in_cases_paths[[i]]))
        setkey(pheno, FID)
        setkey(pheno_file, ID)
        pheno[pheno_file, (pheno_name) := as.numeric(i.VALUE)]
    }
}

# Set controls to 0 for all phenotypes where not specified
control_ids <- gsub("\\.hg38", "", controls$V1)
for (pheno_col in paste0("Y", seq_along(r_in_cases_paths))) {
    pheno[FID %in% control_ids & is.na(get(pheno_col)), (pheno_col) := 0]
}

# Update fam file (use first phenotype for V6, as fam has only one phenotype column)
fam <- merge(fam[, .(V1, V2, V3, V4, V5, V6)], pheno[, .(FID, Y1)], by.x="V1", by.y="FID", all.x=TRUE)
fam[, V6 := Y1]
fam[, Y1 := NULL]
fwrite(fam, r_out_fam_path, sep="\t", quote=F, col.names=F)

# Write phenotype file with all phenotypes
fwrite(pheno, r_out_phenotype_path, sep="\t", quote=F, col.names=T)


# Process VCF
dd <- fread(r_in_vcf_path, skip="#CHROM")
dd[, c("A", "CSQ") := tstrsplit(INFO, "CSQ=", fixed=TRUE)]
dd2 <- dd[,"CSQ",with=F]

get_csq_format <- function(vcf_path) {
    header <- fread(vcf_path, skip=0, nrows=4000, header=F, sep="\n", colClasses="character")[[1]]
    csq_line <- header[grep("##INFO=<ID=CSQ", header, fixed=TRUE)]
    if (length(csq_line) == 0) stop("No CSQ format found in VCF header")
    
    # Extract field names
    csq_fields <- strsplit(sub(".*Format: ", "", sub("\">$", "", csq_line)), "\\|")[[1]]
    
    required_fields <- c("Consequence", "SYMBOL", "Feature_type", "Feature", "DISTANCE")
    return(list(
        fields = csq_fields,
        indices = sapply(required_fields, function(f) which(csq_fields == f))
    ))
}

# Function to create variant identifier based on available information
create_variant_id <- function(row, feature_type) {
    if (!is.na(row$SYMBOL) && nzchar(row$SYMBOL)) {
        return(row$SYMBOL)
    } else if (!is.na(row$Feature) && feature_type == "RegulatoryFeature") {
        return(paste0("REG_", row$Feature))
    } else if (row$Consequence == "intergenic_variant" && include_intergenic) {
        # For intergenic variants, include distance to nearest gene if available
        distance_info <- if (!is.na(row$DISTANCE)) paste0("_d", row$DISTANCE) else ""
        return(paste0("INT_", row$CHROM, "_", row$POS, distance_info))
    } else {
        return(NA)
    }
}

csq_info <- get_csq_format(r_in_vcf_path)
processed_csq <- sapply(dd2$CSQ, function(csq_entry) {
    # Split the CSQ entry by pipe character
    fields <- strsplit(csq_entry, "\\|")[[1]]
    expected_fields <- length(csq_info$fields)
    
    # Pad or trim to match expected length
    if (length(fields) < expected_fields) {
        fields <- c(fields, rep("", expected_fields - length(fields)))
    } else if (length(fields) > expected_fields) {
        fields <- fields[1:expected_fields]
    }
    # Select only the required fields using indices and combine with tab separator
    paste(fields[csq_info$indices], collapse="\t")
})

out_csq_file <- "./csq_split.tsv"
writeLines(processed_csq, out_csq_file)
dd3 <- fread(out_csq_file, fill=length(csq_info$indices), header=F)
field_names <- names(csq_info$indices)
setnames(dd3, field_names)

# Update main data table with annotations
dd[, (field_names) := dd3[, field_names, with = FALSE]]

cn <- colnames(dd)
cn[1] <- "CHROM"
setnames(dd, cn)

# Create comprehensive variant identifiers
dd[, Symbol := sapply(1:nrow(dd), function(i) {
    create_variant_id(dd[i,], dd$Feature_type[i])
})]

dd[, key := paste0(CHROM, '_', POS, '_', REF, '_', ALT)]

# Create annotations with only required columns
anno <- dd[!is.na(Symbol), c("key", "Symbol", "Consequence"), with=F]

# Filter out chrM and chrY, replace : with _
anno <- anno[!grepl("^chrM|^chrY", key)]
anno$key <- gsub(":", "_", anno$key)

# remove multiallelic
anno <- anno[!grepl(",", key)]

# Setlist creation with grouping by Symbol
unique_features <- unique(anno$Symbol)
setlist <- rbindlist(lapply(unique_features, function(feature) {
    feature_rows <- dd[Symbol == feature]
    data.table(
        symbol = feature,
        chrom = feature_rows$CHROM[1],
        pos = feature_rows$POS[1],
        variants = paste(feature_rows$key, collapse=",")
    )
}))
# Filter and process setlist (same chrM/chrY filtering as annotations)
setlist <- setlist[!grepl("^chrM|^chrY", chrom)]
setlist$variants <- gsub(":", "_", setlist$variants)


# Read and parse the mask file
mask_data <- read.table(r_in_masks_path, sep="\t", header=FALSE, col.names=c("Mask", "Annotations"), 
                        stringsAsFactors=FALSE)

# Extract all unique annotations, treating '&' combinations as single units
important_annotations <- unique(unlist(strsplit(mask_data$Annotations, ",")))

# Compute and display Consequence value counts
cat("\n")
cat("Consequence value counts:\n")
consequence_counts <- sort(table(anno$Consequence), decreasing=TRUE)
for (conseq in names(consequence_counts)) {
    cat(sprintf("%s: %d\n", conseq, consequence_counts[conseq]))
}
cat("\n")

# Calculate frequency table
freq_table <- table(anno$Consequence)
sorted_freqs <- sort(freq_table, decreasing=TRUE)

# Determine which consequences to keep
# a. Keep all biologically important annotations
keep_important <- names(freq_table) %in% important_annotations

# b. Apply quantile threshold
threshold_freq <- quantile(freq_table, probs=quantile_threshold, names=FALSE)
keep_quantile <- freq_table >= threshold_freq

# c. Ensure at least 50 top annotations (if available)
n_consequences <- length(freq_table)
n_to_keep <- min(min_top_annotations, n_consequences)  # Don't exceed total unique consequences
keep_top <- names(sorted_freqs)[1:n_to_keep]  # Top N by frequency
keep_top_logical <- names(freq_table) %in% keep_top

# Combine filters: keep if important OR above quantile OR in top 50
keep <- keep_important | keep_quantile | keep_top_logical

# Replace consequences not meeting criteria with "NULL"
anno$Consequence[!keep[match(anno$Consequence, names(freq_table))]] <- "NULL"

# Report kept consequences
kept_consequences <- unique(anno$Consequence[anno$Consequence != "NULL"])
cat("Kept consequences (after filtering):\n")
cat(paste(kept_consequences, collapse=", "), "\n")


# Write final annotations file and setlist file
fwrite(anno, r_out_annotations_path, sep="\t", col.names=F, quote=F)
fwrite(setlist, r_out_setlist_path, sep="\t", col.names=F, quote=F)

# need to set missing to 0 in *.sample file
dd <- fread(r_in_sample_path)
dd$missing <- 0
fwrite(dd, r_out_sample_path, sep="\t", quote=F)