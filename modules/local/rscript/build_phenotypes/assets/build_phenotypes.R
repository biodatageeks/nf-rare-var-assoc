library("data.table")
library("dplyr")
library("optparse")

# Define command line options
option_list <- list(
  make_option(c("--controls-path"), type="character", default=NULL,
              help="Path to input controls file", metavar="character"),
  make_option(c("--cases-paths"), type="character", default=NULL, action="append",
              help="Path to a cases file (can be specified multiple times)", metavar="character"),
  make_option(c("--out-pheno-path"), type="character", default=NULL,
              help="Path to output phenotype file", metavar="character")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opt$`controls-path`)) stop("Input controls file path is required")
if (is.null(opt$`cases-paths`)) stop("At least one cases file is required")
if (is.null(opt$`out-pheno-path`)) stop("Output phenotype file path is required")

r_in_controls_path <- opt$`controls-path`
r_in_cases_paths <- opt$`cases-paths`
r_out_phenotype_path <- opt$`out-pheno-path`

cat("Number of cases files:", length(r_in_cases_paths), "\n")

# Read controls and collect all sample IDs
controls <- fread(r_in_controls_path, header=F)
all_ids <- gsub("\\.hg38", "", controls$V1)

# Collect IDs from all case files
for (case_file in r_in_cases_paths) {
  case_data <- fread(case_file, header=FALSE)
  case_ids <- gsub("\\.hg38", "", case_data$V1)
  all_ids <- unique(c(all_ids, case_ids))
}

# Initialize phenotype table with all unique IDs
pheno <- data.table(FID = all_ids, IID = all_ids)
pheno[, (paste0("Y", seq_along(r_in_cases_paths))) := NA_real_]  # Initialize with NA

# Process each phenotype file
for (i in seq_along(r_in_cases_paths)) {
  pheno_file <- fread(r_in_cases_paths[[i]], header=FALSE)
  pheno_name <- paste0("Y", i)

  setnames(pheno_file, "V1", "ID")
  if (ncol(pheno_file) > 1) setnames(pheno_file, "V2", "VALUE")

  # Check if second column exists and is numeric
  if (ncol(pheno_file) == 1 || !all(grepl("^[0-9.-]+$", pheno_file$VALUE[!is.na(pheno_file$VALUE)]))) {
    cat(sprintf("Treating %s as binary phenotype (no numeric second column)\n", r_in_cases_paths[[i]]))
    pheno[FID %in% pheno_file$ID, (pheno_name) := 1]
  } else {
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

# Write phenotype file
fwrite(pheno, r_out_phenotype_path, sep="\t", quote=F, col.names=T)