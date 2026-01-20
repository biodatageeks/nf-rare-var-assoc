#!/usr/bin/env python3
"""
vcf2aaf.py - Extract allele frequencies from VCF files using polars-bio

This script replaces the slower R-based vcf2aaf.R script.
Uses polars-bio for efficient VCF parsing with lazy evaluation.

Usage:
    python vcf2aaf.py <input.vcf> <output_aaf.tsv> <tag_name> <default_tag_name>

Example:
    python vcf2aaf.py gnomad.vcf output_aaf.tsv AF_nfe AF
"""

import gzip
import re
import sys
import polars as pl

try:
    import polars_bio as pb
except ImportError:
    print("Error: polars-bio is required. Install with: pip install polars-bio", file=sys.stderr)
    sys.exit(1)


def get_vcf_info_tags(vcf_path: str) -> set:
    """
    Parse VCF header to get available INFO tags.
    
    Args:
        vcf_path: Path to VCF file (can be gzipped)
        
    Returns:
        Set of INFO tag names available in the VCF
    """
    info_tags = set()
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                break  # End of header
            if line.startswith('##INFO='):
                # Parse INFO line: ##INFO=<ID=tag_name,...>
                match = re.search(r'ID=([^,>]+)', line)
                if match:
                    info_tags.add(match.group(1))
    
    return info_tags


def extract_aaf(vcf_path: str, output_path: str, tag_name: str, default_tag_name: str) -> None:
    """
    Extract allele frequencies from a VCF file.

    Args:
        vcf_path: Path to the input VCF file
        output_path: Path to the output AAF TSV file
        tag_name: Primary INFO tag to extract (e.g., "AF_nfe")
        default_tag_name: Fallback INFO tag when primary is missing (e.g., "AF")
    """
    # First, check which INFO tags actually exist in the VCF
    available_tags = get_vcf_info_tags(vcf_path)
    print(f"Available INFO tags in VCF: {len(available_tags)} tags", file=sys.stderr)
    
    has_primary = tag_name in available_tags
    has_default = default_tag_name in available_tags
    
    print(f"Primary tag '{tag_name}': {'found' if has_primary else 'NOT FOUND'}", file=sys.stderr)
    print(f"Default tag '{default_tag_name}': {'found' if has_default else 'NOT FOUND'}", file=sys.stderr)
    
    # Build list of fields to request - only include tags that exist
    info_fields = []
    if has_primary:
        info_fields.append(tag_name)
    if has_default:
        info_fields.append(default_tag_name)
    
    if not info_fields:
        print(f"Warning: Neither '{tag_name}' nor '{default_tag_name}' found in VCF. Using 0 for all AF values.", file=sys.stderr)

    # Use scan_vcf for lazy evaluation - more memory efficient for large files
    if info_fields:
        lf = pb.scan_vcf(vcf_path, info_fields=info_fields)
    else:
        lf = pb.scan_vcf(vcf_path, info_fields=[])
    
    # Build select columns dynamically based on available tags
    select_cols = [
        "chrom",
        pl.col("start").alias("pos"),
        "ref",
        "alt",
    ]
    
    if has_primary:
        select_cols.append(pl.col(tag_name).alias("primary_af"))
    if has_default:
        select_cols.append(pl.col(default_tag_name).alias("default_af"))
    
    df = lf.select(select_cols).collect()

    # Process the data using polars expressions
    result = df.with_columns([
        # Remove 'chr' prefix from chromosome (case-insensitive)
        pl.col("chrom").str.replace(r"^[Cc][Hh][Rr]", "").alias("chrom_clean"),
    ]).with_columns([
        # Create position identifier: chr_pos_ref_alt
        pl.concat_str([
            pl.col("chrom_clean"),
            pl.lit("_"),
            pl.col("pos").cast(pl.Utf8),
            pl.lit("_"),
            pl.col("ref"),
            pl.lit("_"),
            pl.col("alt"),
        ]).alias("pos_id"),
    ])
    
    # Build AF expression based on available tags
    if has_primary and has_default:
        # Both tags available - use coalesce logic
        af_expr = (
            pl.when(
                pl.col("primary_af").is_not_null() &
                (pl.col("primary_af").list.first().is_not_null())
            ).then(
                pl.col("primary_af").list.first().cast(pl.Utf8)
            ).when(
                pl.col("default_af").is_not_null() &
                (pl.col("default_af").list.first().is_not_null())
            ).then(
                pl.col("default_af").list.first().cast(pl.Utf8)
            ).otherwise(
                pl.lit("0")
            ).alias("af")
        )
    elif has_primary:
        # Only primary tag available
        af_expr = (
            pl.when(
                pl.col("primary_af").is_not_null() &
                (pl.col("primary_af").list.first().is_not_null())
            ).then(
                pl.col("primary_af").list.first().cast(pl.Utf8)
            ).otherwise(
                pl.lit("0")
            ).alias("af")
        )
    elif has_default:
        # Only default tag available
        af_expr = (
            pl.when(
                pl.col("default_af").is_not_null() &
                (pl.col("default_af").list.first().is_not_null())
            ).then(
                pl.col("default_af").list.first().cast(pl.Utf8)
            ).otherwise(
                pl.lit("0")
            ).alias("af")
        )
    else:
        # No tags available - use 0 for all
        af_expr = pl.lit("0").alias("af")
    
    result = result.with_columns([af_expr]).select([
        pl.col("pos_id").alias("pos"),
        "af",
    ])

    # Write output TSV without header (matching original R script behavior)
    result.write_csv(output_path, separator="\t", include_header=False)

    print(f"Processed {len(result)} variants from {vcf_path}", file=sys.stderr)


def main():
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <input.vcf> <output_aaf.tsv> <tag_name> <default_tag_name>")
        print(f"Example: {sys.argv[0]} gnomad.vcf output_aaf.tsv AF_nfe AF")
        sys.exit(1)

    vcf_path = sys.argv[1]
    output_path = sys.argv[2]
    tag_name = sys.argv[3]
    default_tag_name = sys.argv[4]

    extract_aaf(vcf_path, output_path, tag_name, default_tag_name)


if __name__ == "__main__":
    main()
