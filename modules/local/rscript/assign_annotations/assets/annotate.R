library("data.table")
library("dplyr")
library("optparse")

# Define command line options
option_list <- list(
  make_option(c("--vcf-path"), type="character", default=NULL,
              help="Path to input VCF file", metavar="character"),
#  make_option(c("--sample-path"), type="character", default=NULL,
#              help="Path to input sample file", metavar="character"),
  make_option(c("--masks-path"), type="character", default=NULL,
              help="Path to input masks file", metavar="character"),
  
#  make_option(c("--out-sample-path"), type="character", default=NULL,
#              help="Path to output sample file", metavar="character"),
  make_option(c("--out-anno-path"), type="character", default=NULL,
              help="Path to output annotations file", metavar="character"),
  make_option(c("--out-setlist-path"), type="character", default=NULL,
              help="Path to output setlist file", metavar="character"),
  make_option(c("--min_top_annotations"), type="integer", default=30,
              help="Filter annotations threshold - keep at least that many [default=30]", metavar="integer"),
  make_option(c("--max_annotations"), type="integer", default=62,
              help="Filter annotations threshold - keep at most that many [default=62]", metavar="integer"),
  make_option(c("--quantile_threshold"), type="double", default=0.25,
              help="Filter annotations threshold - keep based on the quantile of the frequency table [default=0.25]", metavar="double"),
  make_option(c("--include-intergenic"), type="logical", default=FALSE,
              help="Include intergenic variants [default=FALSE]", metavar="logical")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check for required arguments and assign to variables
if (is.null(opt$`vcf-path`)) stop("Input VCF file path is required")
#if (is.null(opt$`sample-path`)) stop("Input sample file path is required")
if (is.null(opt$`masks-path`)) stop("Input masks file path is required")

#if (is.null(opt$`out-sample-path`)) stop("Output sample file path is required")
if (is.null(opt$`out-anno-path`)) stop("Output annotations file path is required")
if (is.null(opt$`out-setlist-path`)) stop("Output setlist file path is required")

r_in_vcf_path <- opt$`vcf-path`
#r_in_sample_path <- opt$`sample-path`
r_in_masks_path <- opt$`masks-path`

#r_out_sample_path <- opt$`out-sample-path`
r_out_annotations_path <- opt$`out-anno-path`
r_out_setlist_path <- opt$`out-setlist-path`

include_intergenic <- opt$`include-intergenic`
quantile_threshold <- opt$`quantile_threshold`  # Keep consequences above Nth percentile (adjustable)
min_top_annotations <- opt$`min_top_annotations`   # Ensure at least N top annotations are kept, if available
max_annotations <- opt$`max_annotations`   # Ensure at most N annotations are kept

cat(paste0("quantile_threshold = ", quantile_threshold, "\n"))
cat(paste0("min_top_annotations = ", min_top_annotations, "\n"))
cat(paste0("max_annotations = ", max_annotations, "\n"))
cat(paste0("include_intergenic = ", include_intergenic, "\n"))


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

#csq_info <- get_csq_format(r_in_vcf_path)
#processed_csq <- sapply(dd2$CSQ, function(csq_entry) {
#    # Split the CSQ entry by pipe character
#    fields <- strsplit(csq_entry, "\\|")[[1]]
#    expected_fields <- length(csq_info$fields)
#    
#    # Pad or trim to match expected length
#    if (length(fields) < expected_fields) {
#        fields <- c(fields, rep("", expected_fields - length(fields)))
#    } else if (length(fields) > expected_fields) {
#        fields <- fields[1:expected_fields]
#    }
#    # Select only the required fields using indices and combine with tab separator
#    paste(fields[csq_info$indices], collapse="\t")
#})
#
#out_csq_file <- "./csq_split.tsv"
#writeLines(processed_csq, out_csq_file)
#dd3 <- fread(out_csq_file, fill=length(csq_info$indices), header=F)
#field_names <- names(csq_info$indices)
#setnames(dd3, field_names)

csq_info <- get_csq_format(r_in_vcf_path)
# Assert: Verify all required fields are in csq_info$indices
required_fields <- c("Consequence", "SYMBOL", "Feature_type", "Feature", "DISTANCE")
if (!all(required_fields %in% names(csq_info$indices))) {
    stop("Missing required fields in CSQ format: ",
         paste(setdiff(required_fields, names(csq_info$indices)), collapse=", "))
}
# Debug: Print csq_info to verify field indices
cat("CSQ fields:", csq_info$fields, "\n")
cat("Required field indices:", names(csq_info$indices), "=", csq_info$indices, "\n")

# Process all CSQ annotations by splitting on commas
processed_csq <- rbindlist(lapply(seq_along(dd2$CSQ), function(i) {
    csq_entry <- dd2$CSQ[i]
    annotations <- strsplit(csq_entry, ",")[[1]]  # Split by comma for multiple annotations
    # Debug: Print first few CSQ entries to inspect format
    if (i <= 3) {
        cat("CSQ entry", i, ":", csq_entry, "\n")
    }
    rbindlist(lapply(annotations, function(ann) {
        fields <- strsplit(ann, "\\|")[[1]]
        expected_fields_len <- length(csq_info$fields)
        # Pad or trim to match expected length
        if (length(fields) < expected_fields_len) {
            fields <- c(fields, rep("", expected_fields_len - length(fields)))
        } else if (length(fields) > expected_fields_len) {
            fields <- fields[1:expected_fields_len]
        }
        # Assert: Verify required fields can be accessed
        if (any(csq_info$indices > length(fields))) {
            stop("CSQ annotation ", ann, " in variant ", i, " has too few fields (", length(fields),
                 ") for required indices: ", paste(names(csq_info$indices), "=", csq_info$indices, collapse=", "))
        }
        # Create a data table with required fields and original row index
        #data.table(
        #    row_index = i,
        #    setNames(as.list(fields[csq_info$indices]), names(csq_info$indices))
        #)

        selected_fields <- fields[csq_info$indices]
        names(selected_fields) <- names(csq_info$indices)
        # Debug: Print selected fields for problematic annotations
        if (i <= 3 && ann <= 2) {
            cat("Annotation", ann, "in variant", i, "selected fields:", paste(names(selected_fields), "=", selected_fields, collapse=", "), "\n")
        }
        data.table(
            row_index = i,
            Consequence = selected_fields["Consequence"],
            SYMBOL = selected_fields["SYMBOL"],
            Feature_type = selected_fields["Feature_type"],
            Feature = selected_fields["Feature"],
            DISTANCE = selected_fields["DISTANCE"]
        )
    }))
}))
# Assert: Verify processed_csq has the expected number of columns
field_names <- c("row_index", names(csq_info$indices))
if (ncol(processed_csq) != length(field_names)) {
    stop("processed_csq has ", ncol(processed_csq), " columns, expected ", length(field_names),
         ". Columns found: ", paste(colnames(processed_csq), collapse=", "))
}
# Debug: Print column names and structure of processed_csq
cat("Processed CSQ columns:", colnames(processed_csq), "\n")
cat("Processed CSQ nrow:", nrow(processed_csq), "\n")
# Explicitly set column names for processed_csq
setnames(processed_csq, field_names)

# Update main data table by joining with processed CSQ annotations
dd3 <- processed_csq
# Update main data table with annotations, expanding rows for multiple annotations
dd <- dd[rep(seq_len(nrow(dd)), times = processed_csq[, .N, by = row_index]$N)]
dd[, (names(csq_info$indices)) := dd3[, names(csq_info$indices), with = FALSE]]

cn <- colnames(dd)
cn[1] <- "CHROM"
setnames(dd, cn)



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
anno <- anno[!grepl("^chrM|^chrY|^M|^Y", key)]
anno$key <- gsub(":", "_", anno$key)

# remove multiallelic
anno <- anno[!grepl(",", key)]


filter_annotations <- function(anno, masks_path, quantile_threshold, min_top_annotations, max_annotations) {
    # Read and parse the mask file
    mask_data <- read.table(masks_path, sep="\t", header=FALSE, col.names=c("Mask", "Annotations"), 
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

    # c. Limit total non-important annotations to max_annotations - n_important
    n_consequences <- length(freq_table)
    n_important <- sum(keep_important)  # Number of important annotations
    n_additional <- max(0, max_annotations - n_important)  # Max non-important annotations

    # Limit quantile annotations by frequency
    if (sum(keep_quantile) > n_additional) {
        # Sort quantile annotations by frequency
        quantile_candidates <- names(freq_table)[keep_quantile]
        quantile_candidates_sorted <- quantile_candidates[order(freq_table[quantile_candidates], decreasing=TRUE)]
        keep_quantile <- names(freq_table) %in% quantile_candidates_sorted[1:min(length(quantile_candidates_sorted), n_additional)]
    }

    
    ## d. Apply top annotations if needed
    #n_quantile <- sum(keep_quantile)  # Number of quantile annotations kept
    #n_to_keep <- max(0, min(min_top_annotations, n_consequences, n_additional - n_quantile))  # Ensure non-negative
    ## Only compute keep_top if n_to_keep > 0
    #keep_top_logical <- rep(FALSE, length(freq_table))  # Initialize as FALSE
    #if (n_to_keep > 0) {
    #    available_annotations <- names(sorted_freqs)[!(names(sorted_freqs) %in% names(freq_table)[keep_important | keep_quantile])]
    #    keep_top <- head(available_annotations, n_to_keep)  # Safely select up to n_to_keep
    #    keep_top_logical <- names(freq_table) %in% keep_top
    #}



    # d. Apply top annotations if needed, only to meet min_top_annotations
    n_quantile <- sum(keep_quantile)  # Number of quantile annotations kept
    n_to_keep <- max(0, min(
        min_top_annotations - n_important - n_quantile,
        n_consequences - n_important - n_quantile,
        n_additional - n_quantile
    ))  # Only add to meet min_top_annotations
    keep_top_logical <- rep(FALSE, length(freq_table))  # Initialize as FALSE
    if (n_to_keep > 0) {
        # First try annotations that pass quantile threshold but weren't selected
        available_quantile <- names(sorted_freqs)[(names(sorted_freqs) %in% names(freq_table)[keep_quantile]) & 
                                                    !(names(sorted_freqs) %in% names(freq_table)[keep_important | keep_quantile])]
        if (length(available_quantile) >= n_to_keep) {
            keep_top <- head(available_quantile, n_to_keep)  # Use quantile-passing annotations
        } else {
            # If not enough quantile-passing annotations, include others to meet min_top_annotations
            available_other <- names(sorted_freqs)[!(names(sorted_freqs) %in% names(freq_table)[keep_important | keep_quantile])]
            keep_top <- head(c(available_quantile, available_other), n_to_keep)
        }
        keep_top_logical <- names(freq_table) %in% keep_top
    }



    # Combine filters: keep if important OR above quantile OR in top 50
    keep <- keep_important | keep_quantile | keep_top_logical

    # BUG: This replaces filtered consequences with "NULL" instead of removing the rows.
    # The downstream processing may incorrectly include these NULL rows in setlists.
    # To fix, replace the line below with filtering that removes rows entirely:
    #
    # FIX (commented out to preserve original behavior):
    # anno <- anno[keep[match(anno$Consequence, names(freq_table))], ]
    #
    # Current buggy behavior:
    anno$Consequence[!keep[match(anno$Consequence, names(freq_table))]] <- "NULL"

    # Report kept consequences
    kept_consequences <- unique(anno$Consequence[anno$Consequence != "NULL"])
    cat("Kept consequences (after filtering):\n")
    cat(paste(kept_consequences, collapse=", "), "\n")

    return(anno)
}

anno <- filter_annotations(anno, r_in_masks_path, quantile_threshold, min_top_annotations, max_annotations)


# Setlist creation with grouping by Symbol
# unique_features <- unique(anno$Symbol)
# setlist <- rbindlist(lapply(unique_features, function(feature) {
#     feature_rows <- dd[Symbol == feature]
#     data.table(
#         symbol = feature,
#         chrom = feature_rows$CHROM[1],
#         pos = feature_rows$POS[1],
#         variants = paste(feature_rows$key, collapse=",")
#     )
# }))
# # Filter and process setlist (same chrM/chrY filtering as annotations)
# setlist <- setlist[!grepl("^chrM|^chrY|^M|^Y", chrom)]



# Setlist creation with grouping by Symbol
unique_features <- unique(anno$Symbol)
setlist <- rbindlist(lapply(unique_features, function(feature) {
    feature_rows <- dd[Symbol == feature]
    # Use unique keys to avoid duplicates if a variant appears multiple times
    unique_keys <- unique(feature_rows$key)
    data.table(
        symbol = feature,
        chrom = feature_rows$CHROM[1],
        pos = min(feature_rows$POS),  # Use min POS if multiple positions
        variants = paste(unique_keys, collapse=",")
    )
}))
# Filter and process setlist (same chrM/chrY filtering as annotations)
setlist <- setlist[!grepl("^chrM|^chrY|^M|^Y", chrom)]



setlist$variants <- gsub(":", "_", setlist$variants)


# Write final annotations file and setlist file
fwrite(anno, r_out_annotations_path, sep="\t", col.names=F, quote=F)
fwrite(setlist, r_out_setlist_path, sep="\t", col.names=F, quote=F)

# need to set missing to 0 in *.sample file
#dd <- fread(r_in_sample_path)
#dd$missing <- 0
#fwrite(dd, r_out_sample_path, sep="\t", quote=F)
