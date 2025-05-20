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
  
  make_option(c("--out-fam-path"), type="character", default=NULL,
              help="Path to output fam file", metavar="character"),
  make_option(c("--out-pheno-path"), type="character", default=NULL,
              help="Path to output phenotype file", metavar="character")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check for required arguments and assign to variables
if (is.null(opt$`fam-path`)) stop("Input fam file path is required")
if (is.null(opt$`controls-path`)) stop("Input controls file path is required")
if (is.null(opt$`cases-paths`)) stop("At least one cases file is required")

if (is.null(opt$`out-fam-path`)) stop("Output fam file path is required")
if (is.null(opt$`out-pheno-path`)) stop("Output phenotype file path is required")

r_in_fam_path <- opt$`fam-path`
r_in_controls_path <- opt$`controls-path`
r_in_cases_paths <- opt$`cases-paths`    # List of cases file paths

r_out_fam_path <- opt$`out-fam-path`
r_out_phenotype_path <- opt$`out-pheno-path`

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
