#!/usr/bin/env python3
"""
Restrict the known-answer files to the numbered chromosomes.

WHY: STAAR's shipped gene model (`genes_info`
in STAARpipeline/R/sysdata.rda) is AUTOSOMES-ONLY -- 0 chrX genes. Our fixture is
chr12/22/X and the causal truth DOES include chrX genes (e.g. run_18 names SLC9A7,
REPS2, GAGE12B). Scoring the chain against chrX genes it structurally cannot test
would penalise STAAR for an artifact, not a result. So BOTH arms of the chained
comparison must be scored on autosomal (chr12 + chr22) causal genes only; "no chrX
gene-centric test" is recorded in the functionality matrix instead.

This produces an autosome-only COPY of each dataset's causal files, mirroring the
`run_<N>/select_genes/` layout, so the same globs used by run_eval.sh work against
the filtered tree. It touches nothing else; the original truth files are unchanged
and the pipeline is NOT re-run (only the cheap eval re-scores).

Variant IDs are CHROM_POS_REF_ALT (chr-stripped), e.g. `12_9068058_T_A`, `X_...`;
the chromosome is the token before the first '_'. A variant is kept iff that token
is a plain integer in the autosome set (default 1..22 -> for this fixture: 12, 22).

  causal genes file : TSV `Gene<TAB>Variant` WITH a header line (kept verbatim);
                      each data row filtered by its Variant's chromosome.
  causal snplist    : one variant ID per line, NO header; filtered line-by-line.

Usage:
  python filter_causal_autosomal.py --datasets-dir <dir> --out-dir <dir>
  # optional: --chroms 1-22   (comma list and ranges, e.g. "1-22" or "12,22")
"""
from __future__ import annotations

import argparse
import glob
import os
import sys


def parse_chroms(spec: str) -> set[int]:
    out: set[int] = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-", 1)
            out.update(range(int(lo), int(hi) + 1))
        else:
            out.add(int(part))
    return out


def chrom_of(variant_id: str) -> str:
    """Chromosome token of a CHROM_POS_REF_ALT id (before the first '_')."""
    return variant_id.split("_", 1)[0]


def is_autosomal(variant_id: str, autosomes: set[int]) -> bool:
    tok = chrom_of(variant_id)
    return tok.isdigit() and int(tok) in autosomes


def looks_like_header(line: str) -> bool:
    """A genes-file header row: two tab fields, neither a CHROM_POS_REF_ALT id."""
    cols = line.rstrip("\n").split("\t")
    return len(cols) >= 2 and "_" not in cols[1]


def filter_genes_file(src: str, dst: str, autosomes: set[int]) -> tuple[int, int]:
    kept = total = 0
    with open(src) as fh:
        lines = fh.readlines()
    out = []
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        if i == 0 and looks_like_header(line):
            out.append(line if line.endswith("\n") else line + "\n")
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 2:
            continue
        total += 1
        if is_autosomal(cols[1], autosomes):
            out.append(line if line.endswith("\n") else line + "\n")
            kept += 1
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as fh:
        fh.writelines(out)
    return kept, total


def filter_snplist_file(src: str, dst: str, autosomes: set[int]) -> tuple[int, int]:
    kept = total = 0
    out = []
    with open(src) as fh:
        for line in fh:
            vid = line.strip()
            if not vid:
                continue
            total += 1
            if is_autosomal(vid, autosomes):
                out.append(vid + "\n")
                kept += 1
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as fh:
        fh.writelines(out)
    return kept, total


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datasets-dir", required=True,
                    help="root holding run_<N>/select_genes/ (originals; unchanged)")
    ap.add_argument("--out-dir", required=True,
                    help="root for the autosome-only mirror")
    ap.add_argument("--chroms", default="1-22",
                    help="autosome set to keep (default 1-22)")
    ap.add_argument("--genes-glob",
                    default="run_*/select_genes/*_genes_dataset_idx_*.txt")
    ap.add_argument("--snplist-glob",
                    default="run_*/select_genes/*_dataset_idx_*_in_*.snplist")
    args = ap.parse_args()

    autosomes = parse_chroms(args.chroms)
    datasets_dir = os.path.abspath(args.datasets_dir)
    out_dir = os.path.abspath(args.out_dir)
    if out_dir == datasets_dir:
        sys.exit("ERROR: --out-dir must differ from --datasets-dir (won't overwrite originals)")

    genes = sorted(glob.glob(os.path.join(datasets_dir, args.genes_glob)))
    snps = sorted(glob.glob(os.path.join(datasets_dir, args.snplist_glob)))
    if not genes and not snps:
        sys.exit(f"ERROR: no causal files matched under {datasets_dir}")

    print(f"autosomes kept: {sorted(autosomes)}")
    print(f"{'file':<70} {'kept':>7} {'total':>7} {'dropped':>8}")
    g_kept = g_tot = s_kept = s_tot = 0
    for src in genes:
        rel = os.path.relpath(src, datasets_dir)
        k, t = filter_genes_file(src, os.path.join(out_dir, rel), autosomes)
        g_kept += k; g_tot += t
        print(f"{rel:<70} {k:>7} {t:>7} {t-k:>8}")
    for src in snps:
        rel = os.path.relpath(src, datasets_dir)
        k, t = filter_snplist_file(src, os.path.join(out_dir, rel), autosomes)
        s_kept += k; s_tot += t
        print(f"{rel:<70} {k:>7} {t:>7} {t-k:>8}")

    print("-" * 96)
    print(f"causal gene-variant rows : kept {g_kept}/{g_tot} (dropped {g_tot-g_kept} chrX/non-autosomal)")
    print(f"causal snplist variants  : kept {s_kept}/{s_tot} (dropped {s_tot-s_kept})")
    print(f"filtered truth written under: {out_dir}")


if __name__ == "__main__":
    main()
