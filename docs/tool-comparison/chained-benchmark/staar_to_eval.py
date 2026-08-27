#!/usr/bin/env python3
"""
Convert STAARpipeline gene-centric coding output into the space-separated,
REGENIE-shaped table that nf-eval-gene-assoc's compute_score.py consumes.

WHY SO FEW COLUMNS ARE NEEDED:
compute_score.py reads its results table with `sep=' ', comment='#'` and
REQUIRES columns `ID` and `LOG10P`; it also reads `CHROM, GENPOS, A1FREQ, BETA`
into a groupby.agg. Of these, only LOG10P is load-bearing for scoring: it derives
`gene = ID.split('.')[0]` and takes MAX LOG10P per gene (the exact analogue of
REGENIE's max-over-masks). A1FREQ and GENPOS are aggregated but never read again;
BETA feeds only a volcano scatter that is guarded for emptiness. STAAR-O is an
OMNIBUS test with no effect size or single allele frequency, so we write BETA and
A1FREQ as NA rather than inventing them (na_values=['NA'] on the reader side).

STAAR's `coding()` emits one data frame per functional category with columns:
  "Gene name", "Chr", "Category", "#SNV", "cMAC", ..., "ACAT-O", "STAAR-O"
It OVERWRITES its own position columns with Category/#SNV, so the output carries no
coordinates -- GENPOS is therefore joined from STAAR's `genes_info` (exported from
STAARpipeline/R/sysdata.rda; cols hgnc_symbol|chromosome_name|start_position|
end_position). That export is produced once inside the STAAR container in Step 1;
pass it via --genes-info.

Mapping (one output row per gene x STAAR coding category):
  ID      = "<Gene name>.<Category>"       (gene = ID.split('.')[0] downstream)
  LOG10P  = -log10(STAAR-O p)              (p floored at 1e-300; NA if p missing)
  CHROM   = genes_info.chromosome_name     (fallback: STAAR "Chr" column)
  GENPOS  = genes_info.start_position      (NA if the symbol is absent)
  A1FREQ  = NA
  BETA    = NA

Usage:
  python staar_to_eval.py --staar-results plof.csv missense.csv ... \
      --genes-info genes_info_hgnc.tsv --out tools_comparison_dataset_idx_<N>_step2_Y1.regenie
"""
from __future__ import annotations

import argparse
import math
import sys

import pandas as pd

P_FLOOR = 1e-300  # avoid -log10(0)=inf when STAAR-O underflows

# Candidate column names (STAAR ships spaced/hyphenated names; R writers may mangle).
GENE_COLS = ["Gene name", "Gene.name", "Gene_name", "gene"]
CAT_COLS = ["Category", "category"]
P_COLS = ["STAAR-O", "STAAR.O", "STAAR_O", "STAAR-O.pvalue"]
CHR_COLS = ["Chr", "chr", "chromosome"]


def resolve(df: pd.DataFrame, candidates: list[str], what: str, required: bool = True):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        sys.exit(f"ERROR: no {what} column found; tried {candidates}; "
                 f"have {list(df.columns)}")
    return None


def to_log10p(p) -> str:
    try:
        p = float(p)
    except (TypeError, ValueError):
        return "NA"
    if not (p == p) or p < 0:      # NaN or nonsensical
        return "NA"
    if p <= 0:
        p = P_FLOOR
    return f"{-math.log10(p):.6g}"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--staar-results", nargs="+", required=True,
                    help="one or more STAAR gene-centric coding result CSVs")
    ap.add_argument("--genes-info", required=True,
                    help="genes_info TSV (hgnc_symbol chromosome_name start_position end_position)")
    ap.add_argument("--out", required=True, help="output results table (space-separated)")
    ap.add_argument("--gene-col", default=None, help="override gene-name column")
    ap.add_argument("--cat-col", default=None, help="override category column")
    ap.add_argument("--p-col", default=None, help="override STAAR-O p-value column")
    args = ap.parse_args()

    gi = pd.read_csv(args.genes_info, sep="\t", dtype=str)
    if "hgnc_symbol" not in gi.columns:
        sys.exit(f"ERROR: --genes-info lacks hgnc_symbol; have {list(gi.columns)}")
    gi = gi.drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")
    gi_chrom = gi["chromosome_name"].to_dict()
    gi_start = gi["start_position"].to_dict()

    frames = []
    for path in args.staar_results:
        df = pd.read_csv(path)
        if df.empty:
            continue
        gcol = args.gene_col or resolve(df, GENE_COLS, "gene-name")
        ccol = args.cat_col or resolve(df, CAT_COLS, "category")
        pcol = args.p_col or resolve(df, P_COLS, "STAAR-O p-value")
        chrcol = resolve(df, CHR_COLS, "chromosome", required=False)
        sub = pd.DataFrame({
            "gene": df[gcol].astype(str).str.strip(),
            "category": df[ccol].astype(str).str.strip(),
            "p": df[pcol],
            "chr_staar": df[chrcol].astype(str).str.strip() if chrcol else "NA",
        })
        frames.append(sub)

    if not frames:
        sys.exit("ERROR: all STAAR result files were empty")
    allr = pd.concat(frames, ignore_index=True)

    missing = sorted(set(allr["gene"]) - set(gi.index))
    if missing:
        print(f"WARNING: {len(missing)} STAAR gene symbol(s) absent from genes_info "
              f"(GENPOS/CHROM fall back): {missing[:10]}{' ...' if len(missing) > 10 else ''}",
              file=sys.stderr)

    out_rows = []
    for _, r in allr.iterrows():
        gene = r["gene"]
        chrom = gi_chrom.get(gene, r["chr_staar"] if r["chr_staar"] else "NA")
        genpos = gi_start.get(gene, "NA")
        out_rows.append({
            "CHROM": chrom if chrom else "NA",
            "GENPOS": genpos if genpos else "NA",
            "ID": f"{gene}.{r['category']}",
            "A1FREQ": "NA",
            "BETA": "NA",
            "LOG10P": to_log10p(r["p"]),
        })

    out = pd.DataFrame(out_rows, columns=["CHROM", "GENPOS", "ID", "A1FREQ", "BETA", "LOG10P"])
    with open(args.out, "w") as fh:
        fh.write("## staar_to_eval: STAAR-O gene-centric -> eval table (ID=<SYMBOL>.<category>, LOG10P=-log10 STAAR-O p)\n")
        out.to_csv(fh, sep=" ", index=False)

    print(f"wrote {len(out)} rows ({out['ID'].str.split('.').str[0].nunique()} genes) -> {args.out}")


if __name__ == "__main__":
    main()
