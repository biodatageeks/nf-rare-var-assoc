#!/bin/bash
#
# assign_annotations.sh - Extract and filter VEP annotations from VCF files
#
# This script replaces the slower R-based annotate.R script.
# Uses bcftools +split-vep for efficient VEP annotation parsing.
#
# Usage:
#   bash assign_annotations.sh \
#       --vcf-path input.vcf.gz \
#       --masks-path masks.txt \
#       --out-anno-path output.annotations \
#       --out-setlist-path output.setlist \
#       [--min-top-annotations 30] \
#       [--max-annotations 62] \
#       [--quantile-threshold 0.25] \
#       [--include-intergenic]

set -euo pipefail

# Default values
MIN_TOP_ANNOTATIONS=30
MAX_ANNOTATIONS=62
QUANTILE_THRESHOLD=0.25
INCLUDE_INTERGENIC=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf-path)
            VCF_PATH="$2"
            shift 2
            ;;
        --masks-path)
            MASKS_PATH="$2"
            shift 2
            ;;
        --out-anno-path)
            OUT_ANNO_PATH="$2"
            shift 2
            ;;
        --out-setlist-path)
            OUT_SETLIST_PATH="$2"
            shift 2
            ;;
        --min-top-annotations|--min_top_annotations)
            MIN_TOP_ANNOTATIONS="$2"
            shift 2
            ;;
        --max-annotations|--max_annotations)
            MAX_ANNOTATIONS="$2"
            shift 2
            ;;
        --quantile-threshold|--quantile_threshold)
            QUANTILE_THRESHOLD="$2"
            shift 2
            ;;
        --include-intergenic)
            # Accept TRUE/FALSE value (case-insensitive)
            val=$(echo "$2" | tr '[:lower:]' '[:upper:]')
            if [[ "$val" == "TRUE" || "$val" == "T" || "$val" == "YES" || "$val" == "Y" || "$val" == "1" ]]; then
                INCLUDE_INTERGENIC=true
            else
                INCLUDE_INTERGENIC=false
            fi
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "${VCF_PATH:-}" ]]; then
    echo "Error: --vcf-path is required" >&2
    exit 1
fi
if [[ -z "${MASKS_PATH:-}" ]]; then
    echo "Error: --masks-path is required" >&2
    exit 1
fi
if [[ -z "${OUT_ANNO_PATH:-}" ]]; then
    echo "Error: --out-anno-path is required" >&2
    exit 1
fi
if [[ -z "${OUT_SETLIST_PATH:-}" ]]; then
    echo "Error: --out-setlist-path is required" >&2
    exit 1
fi

echo "Parameters:" >&2
echo "  quantile_threshold = $QUANTILE_THRESHOLD" >&2
echo "  min_top_annotations = $MIN_TOP_ANNOTATIONS" >&2
echo "  max_annotations = $MAX_ANNOTATIONS" >&2
echo "  include_intergenic = $INCLUDE_INTERGENIC" >&2

# Create temporary files
TMP_DIR=$(mktemp -d)
trap "rm -rf $TMP_DIR" EXIT

TMP_VEP="$TMP_DIR/vep_annotations.tsv"
TMP_PROCESSED="$TMP_DIR/processed.tsv"
TMP_FREQ="$TMP_DIR/consequence_freq.tsv"
TMP_KEEP="$TMP_DIR/keep_consequences.txt"

# Step 1: Extract VEP annotations using bcftools +split-vep
echo "Extracting VEP annotations..." >&2
bcftools +split-vep \
    -d \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%SYMBOL\t%Feature_type\t%Feature\t%DISTANCE\n' \
    -A tab \
    "$VCF_PATH" > "$TMP_VEP"

echo "Extracted $(wc -l < "$TMP_VEP") annotation records" >&2

# Step 2: Process annotations - create key and Symbol
echo "Processing annotations..." >&2
awk -F'\t' -v include_intergenic="$INCLUDE_INTERGENIC" '
BEGIN { OFS="\t" }
{
    chrom = $1
    pos = $2
    ref = $3
    alt = $4
    consequence = $5
    symbol = $6
    feature_type = $7
    feature = $8
    distance = $9
    
    # Remove chr prefix
    gsub(/^[Cc][Hh][Rr]/, "", chrom)
    
    # Skip chrM and chrY
    if (chrom ~ /^M$|^Y$/) next
    
    # Skip multiallelic (contains comma in alt)
    if (alt ~ /,/) next
    
    # Create key: CHROM_POS_REF_ALT
    key = chrom "_" pos "_" ref "_" alt
    gsub(/:/, "_", key)
    
    # Create Symbol (variant identifier)
    variant_symbol = ""
    if (symbol != "" && symbol != ".") {
        variant_symbol = symbol
    } else if (feature_type == "RegulatoryFeature" && feature != "" && feature != ".") {
        variant_symbol = "REG_" feature
    } else if (consequence == "intergenic_variant" && include_intergenic == "true") {
        distance_info = ""
        if (distance != "" && distance != ".") {
            distance_info = "_d" distance
        }
        variant_symbol = "INT_" chrom "_" pos distance_info
    }
    
    # Skip if no valid symbol
    if (variant_symbol == "") next
    
    # Output: key, symbol, consequence, chrom, pos
    print key, variant_symbol, consequence, chrom, pos
}' "$TMP_VEP" > "$TMP_PROCESSED"

echo "Processed $(wc -l < "$TMP_PROCESSED") valid records" >&2

# Step 3: Calculate consequence frequencies
echo "Calculating consequence frequencies..." >&2
cut -f3 "$TMP_PROCESSED" | sort | uniq -c | sort -rn | \
    awk '{ print $2 "\t" $1 }' > "$TMP_FREQ"

echo "" >&2
echo "Consequence value counts:" >&2
cat "$TMP_FREQ" >&2
echo "" >&2

# Step 4: Determine which consequences to keep
echo "Filtering consequences..." >&2

# Extract important annotations from masks file
IMPORTANT_ANNOTATIONS=$(cut -f2 "$MASKS_PATH" | tr ',' '\n' | sort -u)

# Calculate threshold using awk (matching R's exact logic)
awk -F'\t' -v quantile="$QUANTILE_THRESHOLD" -v min_top="$MIN_TOP_ANNOTATIONS" \
    -v max_anno="$MAX_ANNOTATIONS" -v important="$IMPORTANT_ANNOTATIONS" '
BEGIN {
    # Load important annotations into array
    n_imp = split(important, imp_arr, "\n")
    for (i = 1; i <= n_imp; i++) {
        if (imp_arr[i] != "") {
            important_set[imp_arr[i]] = 1
        }
    }
}
{
    conseq[NR] = $1
    freq[NR] = $2
    n++
}
END {
    # Step 1: Identify important annotations
    n_important = 0
    for (i = 1; i <= n; i++) {
        if (conseq[i] in important_set) {
            keep_important[conseq[i]] = 1
            n_important++
        }
    }
    
    n_additional = max_anno - n_important
    if (n_additional < 0) n_additional = 0
    
    # Step 2: Calculate quantile threshold using R-compatible method (type=7)
    # Create array of sorted frequencies (ascending order for quantile)
    for (i = 1; i <= n; i++) {
        asc_freq[i] = freq[n - i + 1]  # reverse to ascending
    }
    
    # R quantile type=7: index = 1 + (n-1)*p, then interpolate
    idx_float = 1 + (n - 1) * quantile
    idx_lo = int(idx_float)
    idx_hi = idx_lo + 1
    if (idx_hi > n) idx_hi = n
    frac = idx_float - idx_lo
    
    # Linear interpolation
    threshold_freq = asc_freq[idx_lo] + frac * (asc_freq[idx_hi] - asc_freq[idx_lo])
    
    # Step 3: Identify all consequences that pass quantile threshold (including important)
    n_quantile = 0
    for (i = 1; i <= n; i++) {
        if (freq[i] >= threshold_freq) {
            keep_quantile[conseq[i]] = 1
            quantile_order[n_quantile++] = conseq[i]
        }
    }
    
    # Step 4: Limit keep_quantile to top n_additional by frequency (R behavior)
    # Note: This includes important annotations in the count, matching R exactly
    if (n_quantile > n_additional) {
        # Already sorted by frequency (descending), take first n_additional
        delete keep_quantile
        for (i = 0; i < n_additional; i++) {
            keep_quantile[quantile_order[i]] = 1
        }
    }
    
    # Step 5: Combine - keep if important OR in quantile set
    for (c in keep_important) {
        keep[c] = 1
    }
    for (c in keep_quantile) {
        keep[c] = 1
    }
    
    # Step 6: Ensure minimum annotations (R exact logic)
    # R calculates n_to_keep using a formula that prevents adding more than available:
    # n_to_keep = max(0, min(
    #     min_top_annotations - n_important - n_quantile,  # How many more needed
    #     n_consequences - n_important - n_quantile,       # How many available (may be negative!)
    #     n_additional - n_quantile                        # Space in max limit
    # ))
    # Recount n_quantile after limiting
    actual_n_quantile = 0
    for (c in keep_quantile) actual_n_quantile++
    
    n_to_keep_1 = min_top - n_important - actual_n_quantile
    n_to_keep_2 = n - n_important - actual_n_quantile
    n_to_keep_3 = n_additional - actual_n_quantile
    
    n_to_keep = n_to_keep_1
    if (n_to_keep_2 < n_to_keep) n_to_keep = n_to_keep_2
    if (n_to_keep_3 < n_to_keep) n_to_keep = n_to_keep_3
    if (n_to_keep < 0) n_to_keep = 0
    
    added = 0
    for (i = 1; i <= n && added < n_to_keep; i++) {
        if (!(conseq[i] in keep)) {
            keep[conseq[i]] = 1
            added++
        }
    }
    
    # Output kept consequences
    for (c in keep) {
        print c
    }
}' "$TMP_FREQ" > "$TMP_KEEP"

echo "Keeping $(wc -l < "$TMP_KEEP") consequence types" >&2

# Step 5: Create annotations output (key, Symbol, Consequence)
# SKIP rows with consequences not in keep set (not replace with NULL)
echo "Creating annotations file..." >&2
awk -F'\t' '
BEGIN { OFS="\t" }
NR == FNR { keep[$1] = 1; next }
{
    key = $1
    symbol = $2
    consequence = $3
    
    # Skip row if consequence not in keep set
    if (!(consequence in keep)) {
        next
    }
    
    # Output unique key-symbol-consequence combinations
    combo = key "\t" symbol "\t" consequence
    if (!(combo in seen)) {
        seen[combo] = 1
        print key, symbol, consequence
    }
}' "$TMP_KEEP" "$TMP_PROCESSED" > "$OUT_ANNO_PATH"

echo "Written $(wc -l < "$OUT_ANNO_PATH") annotations" >&2

# Step 6: Create setlist output (symbol, chrom, pos, variants)
echo "Creating setlist file..." >&2
awk -F'\t' '
BEGIN { OFS="\t" }
{
    key = $1
    symbol = $2
    chrom = $4
    pos = $5
    
    # Track variants per symbol
    if (!(symbol in first_chrom)) {
        first_chrom[symbol] = chrom
        min_pos[symbol] = pos
    }
    if (pos < min_pos[symbol]) {
        min_pos[symbol] = pos
    }
    
    # Collect unique variants per symbol
    if (!(symbol SUBSEP key in seen_var)) {
        seen_var[symbol SUBSEP key] = 1
        if (symbol in variants) {
            variants[symbol] = variants[symbol] "," key
        } else {
            variants[symbol] = key
        }
    }
}
END {
    for (symbol in variants) {
        # Skip if chrom is M or Y
        if (first_chrom[symbol] ~ /^M$|^Y$/) continue
        
        print symbol, first_chrom[symbol], min_pos[symbol], variants[symbol]
    }
}' "$TMP_PROCESSED" | sort -k2,2 -k3,3n > "$OUT_SETLIST_PATH"

echo "Written $(wc -l < "$OUT_SETLIST_PATH") gene sets" >&2
echo "Done!" >&2
