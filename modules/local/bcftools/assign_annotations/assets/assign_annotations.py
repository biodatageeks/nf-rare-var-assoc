#!/usr/bin/env python3
"""
assign_annotations.py - Extract and filter VEP annotations from VCF files

This script replaces the slower R-based annotate.R script.
Uses bcftools +split-vep for efficient VEP annotation parsing,
then processes with polars for fast data manipulation.

Usage:
    python assign_annotations.py \
        --vcf-path input.vcf.gz \
        --masks-path masks.txt \
        --out-anno-path output.annotations \
        --out-setlist-path output.setlist \
        [--min-top-annotations 30] \
        [--max-annotations 62] \
        [--quantile-threshold 0.25] \
        [--include-intergenic]
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import polars as pl


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract and filter VEP annotations from VCF files"
    )
    parser.add_argument("--vcf-path", required=True, help="Path to input VCF file")
    parser.add_argument("--masks-path", required=True, help="Path to input masks file")
    parser.add_argument(
        "--out-anno-path", required=True, help="Path to output annotations file"
    )
    parser.add_argument(
        "--out-setlist-path", required=True, help="Path to output setlist file"
    )
    parser.add_argument(
        "--min-top-annotations", "--min_top_annotations",
        type=int,
        default=30,
        help="Keep at least this many annotations [default: 30]",
    )
    parser.add_argument(
        "--max-annotations", "--max_annotations",
        type=int,
        default=62,
        help="Keep at most this many annotations [default: 62]",
    )
    parser.add_argument(
        "--quantile-threshold", "--quantile_threshold",
        type=float,
        default=0.25,
        help="Keep annotations above this quantile threshold [default: 0.25]",
    )
    parser.add_argument(
        "--include-intergenic",
        type=str,
        default="FALSE",
        help="Include intergenic variants (TRUE/FALSE) [default: FALSE]",
    )
    
    args = parser.parse_args()
    # Convert include-intergenic string to boolean
    args.include_intergenic = args.include_intergenic.upper() in ("TRUE", "T", "YES", "Y", "1")
    return args


def extract_vep_annotations(vcf_path: str, output_tsv: str) -> None:
    """
    Use bcftools +split-vep to extract VEP annotations efficiently.
    
    bcftools +split-vep is a specialized plugin for parsing VEP CSQ annotations.
    It handles:
    - Multiple transcripts per variant (comma-separated in CSQ)
    - Proper field extraction based on CSQ header format
    - Streaming processing (low memory)
    """
    # bcftools +split-vep extracts VEP annotations and outputs as TSV
    # -f specifies output format, -d splits duplicate variants (one line per transcript)
    # -A tab outputs tab-delimited
    cmd = [
        "bcftools", "+split-vep",
        "-d",  # Duplicate: output one line per transcript/consequence
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%SYMBOL\t%Feature_type\t%Feature\t%DISTANCE\n",
        "-A", "tab",  # Tab-delimited output
        vcf_path
    ]
    
    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    
    with open(output_tsv, 'w') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        print(f"bcftools error: {result.stderr}", file=sys.stderr)
        raise RuntimeError(f"bcftools +split-vep failed with return code {result.returncode}")


def load_important_annotations(masks_path: str) -> set:
    """Load biologically important annotations from masks file."""
    important = set()
    with open(masks_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                # Split annotations by comma, keep '&' combinations as single units
                annotations = parts[1].split(',')
                important.update(annotations)
    return important


def create_variant_id(row: dict, include_intergenic: bool) -> str | None:
    """Create variant identifier based on available information."""
    symbol = row.get("SYMBOL", "")
    feature = row.get("Feature", "")
    feature_type = row.get("Feature_type", "")
    consequence = row.get("Consequence", "")
    chrom = row.get("CHROM", "")
    pos = row.get("POS", "")
    distance = row.get("DISTANCE", "")
    
    if symbol and symbol != ".":
        return symbol
    elif feature and feature != "." and feature_type == "RegulatoryFeature":
        return f"REG_{feature}"
    elif consequence == "intergenic_variant" and include_intergenic:
        distance_info = f"_d{distance}" if distance and distance != "." else ""
        return f"INT_{chrom}_{pos}{distance_info}"
    else:
        return None


def filter_annotations(
    df: pl.DataFrame,
    important_annotations: set,
    quantile_threshold: float,
    min_top_annotations: int,
    max_annotations: int,
) -> pl.DataFrame:
    """
    Filter annotations based on frequency and importance.
    
    Logic:
    1. Always keep biologically important annotations (from masks file)
    2. Keep annotations above quantile threshold
    3. Ensure at least min_top_annotations are kept
    4. Limit total to max_annotations
    """
    # Calculate consequence frequency table
    consequence_counts = (
        df.group_by("Consequence")
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
    )
    
    # Print consequence value counts
    print("\nConsequence value counts:", file=sys.stderr)
    for row in consequence_counts.iter_rows(named=True):
        print(f"  {row['Consequence']}: {row['count']}", file=sys.stderr)
    print("", file=sys.stderr)
    
    # Get frequency values
    freq_values = consequence_counts["count"].to_list()
    consequence_names = consequence_counts["Consequence"].to_list()
    
    # Calculate quantile threshold value using R-compatible method (type=7)
    # R quantile type=7: index = 1 + (n-1)*p, then linear interpolation
    if freq_values:
        sorted_freqs = sorted(freq_values)
        n = len(sorted_freqs)
        # R uses 1-based indexing, we use 0-based
        # R: index = 1 + (n-1)*p  ->  Python: index = (n-1)*p
        idx_float = (n - 1) * quantile_threshold
        idx_lo = int(idx_float)
        idx_hi = min(idx_lo + 1, n - 1)
        frac = idx_float - idx_lo
        # Linear interpolation
        threshold_freq = sorted_freqs[idx_lo] + frac * (sorted_freqs[idx_hi] - sorted_freqs[idx_lo])
    else:
        threshold_freq = 0
    
    # Build frequency lookup for limiting step
    freq_lookup = dict(zip(consequence_names, freq_values))
    
    # Determine which consequences to keep (matching R's exact logic)
    
    # 1. Identify important annotations
    keep_important = set()
    for conseq in consequence_names:
        if conseq in important_annotations:
            keep_important.add(conseq)
    
    n_important = len(keep_important)
    n_additional = max(0, max_annotations - n_important)
    
    # 2. Identify quantile-passing consequences, excluding ones already kept as important.
    # Keeping the two sets disjoint makes the size accounting in step 5 trivially correct
    # when an important consequence is also above the quantile threshold.
    keep_quantile = set()
    for conseq, freq in zip(consequence_names, freq_values):
        if freq >= threshold_freq and conseq not in keep_important:
            keep_quantile.add(conseq)

    # 3. Trim keep_quantile by frequency to fit alongside the important set under max_annotations.
    if len(keep_quantile) > n_additional:
        # Sort by frequency descending
        quantile_sorted = sorted(keep_quantile, key=lambda c: freq_lookup[c], reverse=True)
        keep_quantile = set(quantile_sorted[:n_additional])

    # 4. Combine: keep if important OR in quantile set
    keep_consequences = keep_important | keep_quantile

    # 5. Ensure at least min_top_annotations are kept.
    # min() terms: more-needed-for-floor, available-remaining, room-before-max.
    n_consequences = len(consequence_names)
    n_quantile = len(keep_quantile)
    n_to_keep = max(0, min(
        min_top_annotations - n_important - n_quantile,
        n_consequences - n_important - n_quantile,
        n_additional - n_quantile
    ))
    
    if n_to_keep > 0:
        added = 0
        for conseq in consequence_names:  # Already sorted by frequency descending
            if conseq not in keep_consequences:
                keep_consequences.add(conseq)
                added += 1
                if added >= n_to_keep:
                    break
    
    # Filter: REMOVE rows with consequences not in keep set (not replace with NULL)
    print(f"Kept consequences ({len(keep_consequences)}):", file=sys.stderr)
    print(f"  {', '.join(sorted(keep_consequences))}", file=sys.stderr)
    
    df = df.filter(pl.col("Consequence").is_in(list(keep_consequences)))
    
    return df


def main():
    args = parse_args()
    
    print(f"quantile_threshold = {args.quantile_threshold}", file=sys.stderr)
    print(f"min_top_annotations = {args.min_top_annotations}", file=sys.stderr)
    print(f"max_annotations = {args.max_annotations}", file=sys.stderr)
    print(f"include_intergenic = {args.include_intergenic}", file=sys.stderr)
    
    # Step 1: Extract VEP annotations using bcftools +split-vep
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp:
        tmp_tsv = tmp.name
    
    try:
        extract_vep_annotations(args.vcf_path, tmp_tsv)
        
        # Step 2: Load extracted annotations with polars
        # Force all columns to be read as strings to avoid type inference issues
        df = pl.read_csv(
            tmp_tsv,
            separator='\t',
            has_header=False,
            new_columns=["CHROM", "POS", "REF", "ALT", "Consequence", "SYMBOL", "Feature_type", "Feature", "DISTANCE"],
            schema={
                "CHROM": pl.Utf8,
                "POS": pl.Utf8,
                "REF": pl.Utf8,
                "ALT": pl.Utf8,
                "Consequence": pl.Utf8,
                "SYMBOL": pl.Utf8,
                "Feature_type": pl.Utf8,
                "Feature": pl.Utf8,
                "DISTANCE": pl.Utf8,
            },
        )
        
        print(f"Loaded {len(df)} annotation records", file=sys.stderr)
        
    finally:
        Path(tmp_tsv).unlink(missing_ok=True)
    
    # Step 3: Create variant key
    df = df.with_columns([
        # Remove 'chr' prefix
        pl.col("CHROM").str.replace(r"^[Cc][Hh][Rr]", "").alias("chrom_clean"),
    ]).with_columns([
        # Create key: CHROM_POS_REF_ALT
        pl.concat_str([
            pl.col("chrom_clean"),
            pl.lit("_"),
            pl.col("POS"),
            pl.lit("_"),
            pl.col("REF"),
            pl.lit("_"),
            pl.col("ALT"),
        ]).alias("key"),
    ])
    
    # Step 4: Create Symbol column (variant identifier)
    # Using when/then/otherwise for vectorized operation
    df = df.with_columns([
        pl.when(
            (pl.col("SYMBOL").is_not_null()) & 
            (pl.col("SYMBOL") != "") & 
            (pl.col("SYMBOL") != ".")
        ).then(pl.col("SYMBOL"))
        .when(
            (pl.col("Feature_type") == "RegulatoryFeature") &
            (pl.col("Feature").is_not_null()) &
            (pl.col("Feature") != ".")
        ).then(pl.concat_str([pl.lit("REG_"), pl.col("Feature")]))
        .when(
            pl.lit(args.include_intergenic) &
            (pl.col("Consequence") == "intergenic_variant")
        ).then(
            pl.concat_str([
                pl.lit("INT_"),
                pl.col("chrom_clean"),
                pl.lit("_"),
                pl.col("POS"),
                pl.when(
                    (pl.col("DISTANCE").is_not_null()) & 
                    (pl.col("DISTANCE") != "") &
                    (pl.col("DISTANCE") != ".")
                ).then(pl.concat_str([pl.lit("_d"), pl.col("DISTANCE")]))
                .otherwise(pl.lit(""))
            ])
        )
        .otherwise(pl.lit(None))
        .alias("Symbol")
    ])
    
    # Step 5: Filter out rows with null Symbol, chrM, chrY, and multiallelic
    df = df.filter(
        pl.col("Symbol").is_not_null() &
        ~pl.col("key").str.contains(r"^M_|^Y_|^chrM|^chrY") &
        ~pl.col("key").str.contains(",")
    )
    
    # Replace colons with underscores in key
    df = df.with_columns(
        pl.col("key").str.replace_all(":", "_").alias("key")
    )
    
    # Step 6: Filter annotations based on frequency and importance
    important_annotations = load_important_annotations(args.masks_path)
    df = filter_annotations(
        df,
        important_annotations,
        args.quantile_threshold,
        args.min_top_annotations,
        args.max_annotations,
    )
    
    # Step 7: Create annotations output (key, Symbol, Consequence)
    # Keep all rows (like R does) - don't deduplicate
    # The order is determined by bcftools +split-vep which processes transcripts
    # in VCF order (more severe consequences first within each variant)
    anno = df.select(["key", "Symbol", "Consequence"])
    
    # Step 8: Create setlist output (symbol, chrom, pos, variants)
    setlist = (
        df.group_by("Symbol")
        .agg([
            pl.col("chrom_clean").first().alias("chrom"),
            pl.col("POS").min().alias("pos"),
            pl.col("key").unique().alias("variants_list"),
        ])
        .with_columns(
            pl.col("variants_list").list.join(",").alias("variants")
        )
        .select(["Symbol", "chrom", "pos", "variants"])
        .rename({"Symbol": "symbol"})
        .filter(~pl.col("chrom").str.contains(r"^M$|^Y$"))
    )
    
    # Replace colons with underscores in setlist variants
    setlist = setlist.with_columns(
        pl.col("variants").str.replace_all(":", "_").alias("variants")
    )
    
    # Step 9: Write outputs (no header, tab-separated)
    anno.write_csv(args.out_anno_path, separator='\t', include_header=False)
    setlist.write_csv(args.out_setlist_path, separator='\t', include_header=False)
    
    print(f"Written {len(anno)} annotations to {args.out_anno_path}", file=sys.stderr)
    print(f"Written {len(setlist)} gene sets to {args.out_setlist_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
