#!/usr/bin/env python3
"""
recompute_naive_log10p.py - independent sanity check for the IT-6/IT-7 workflow tests.

Regenie's gene-burden output (LOG10P per gene/mask) is the product of a long pipeline:
LOCO step-1 corrections, dosage import, AAF binning, SKAT-style aggregation. There is no
cheap way to reproduce it exactly. Instead this script recomputes a *deliberately naive*
single-variant chi-square signal for the genes Regenie ranks highest, and checks that the
two broadly agree. The point is to fail loudly on wiring regressions (cases/controls
swapped, phenotype-to-sample mismap, dropped masks/covariates) that would make Regenie's
top hits statistically meaningless, while tolerating the large numerical gap that legitimately
exists between a burden test and a naive per-variant test.

Naive LOG10P recipe (used by the IT-6 integration test):
  1. Parse the merged regenie output. The ID column encodes the gene as the substring
     before the first '.', e.g. 'HOXC4.Mask_Mod.0.05' -> gene 'HOXC4', mask 'Mask_Mod'.
     Group rows by gene; keep max(LOG10P) and the mask of the winning row. Skip NA LOG10P.
  2. Take the top-N genes (default 5) by max(LOG10P).
  3. For each top gene:
     a. Look up the gene's variants in the setlist file.
     b. Keep variants whose annotation consequence (in the annotations file) falls in the
        winning mask's consequence set (from the masks file).
     c. Count 0/0, 0/1, 1/1 (by alt-allele dosage 0/1/2) in cases and controls, read from
        the input VCF via pysam.
     d. Per variant build a 2x3 contingency table and compute chi-square -log10(p)
        (scipy.stats.chi2_contingency, zero-sum columns dropped first).
     e. naive_log10p[gene] = max -log10(p) across the gene's qualifying variants.

Pass/fail (tunable via flags; defaults chosen after observing a healthy fixture run --
see the IT-6 test file for the calibrated values actually used):
  A gene "violates" if its naive LOG10P is non-finite, below --floor, or differs from
  Regenie's LOG10P by more than --tol. The check fails if >= --max-violations of the
  top-N genes violate. Isolated single-gene violations are tolerated because on a
  no-real-signal fixture one borderline top hit is expected noise.

Exit status: 0 = pass, 1 = check failed, 2 = usage / data error.
"""

import argparse
import gzip
import math
import sys

import numpy as np
import pysam
from scipy.stats import chi2_contingency

LN10 = math.log(10.0)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--regenie", required=True, help="regenie step-2 / merged output")
    p.add_argument("--vcf", required=True, help="CSQ-annotated input VCF (indexed)")
    p.add_argument("--vcf-tbi", required=False, default=None, help="index for the input VCF")
    p.add_argument("--setlist", required=True, help="gene\\tchrom\\tpos\\tvar1,var2,...")
    p.add_argument("--annotations", required=True, help="variant_id\\tgene\\tconsequence")
    p.add_argument("--masks", required=True, help="mask_name\\tconseq1,conseq2,...")
    p.add_argument("--cases", required=True, help="case sample IDs, one per line")
    p.add_argument("--controls", required=True, help="control sample IDs, one per line")
    p.add_argument("--top-n", type=int, default=5)
    p.add_argument("--floor", type=float, default=0.5,
                   help="minimum acceptable naive -log10(p) for a top gene")
    p.add_argument("--tol", type=float, default=2.0,
                   help="max allowed |naive - regenie| LOG10P gap per gene")
    p.add_argument("--max-violations", type=int, default=3,
                   help="fail if at least this many of the top-N genes violate")
    p.add_argument("--report-only", action="store_true",
                   help="print the table but always exit 0 (for tolerance calibration)")
    return p.parse_args()


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_id_list(path):
    with open(path) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def parse_regenie(path):
    """Return {gene: (max_log10p, winning_mask)} over set-based rows."""
    best = {}
    with open_text(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            fields = line.split()  # regenie output is whitespace-separated
            if header is None:
                header = fields
                try:
                    id_i = header.index("ID")
                    lp_i = header.index("LOG10P")
                except ValueError:
                    sys.exit("ERROR: regenie output missing ID/LOG10P columns: "
                             f"{header}")
                continue
            if len(fields) <= max(id_i, lp_i):
                continue
            vid = fields[id_i]
            raw = fields[lp_i]
            if raw in ("NA", "nan", "", "."):
                continue
            try:
                lp = float(raw)
            except ValueError:
                continue
            if "." not in vid:
                continue  # not a gene.mask.aaf set row
            parts = vid.split(".")
            gene = parts[0]
            mask = parts[1] if len(parts) > 1 else None
            if gene not in best or lp > best[gene][0]:
                best[gene] = (lp, mask)
    return best


def load_setlist(path):
    """Return {gene: [variant_id, ...]}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            out[f[0]] = [v for v in f[3].split(",") if v]
    return out


def load_annotations(path):
    """Return {(variant_id, gene): set(consequences)}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            out.setdefault((f[0], f[1]), set()).add(f[2])
    return out


def load_masks(path):
    """Return {mask_name: set(consequences)}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            out[f[0]] = set(f[1].split(","))
    return out


def variant_region(vid):
    """'12_553708_G_A' -> (chrom, pos, ref, alt)."""
    chrom, rest = vid.split("_", 1)
    pos_s, ref, alt = rest.split("_", 2)
    return chrom, int(pos_s), ref, alt


def count_genotypes(vcf, vid, case_set, ctrl_set):
    """2x3 contingency [[case 0,1,2 alt],[ctrl 0,1,2 alt]] from the matching record."""
    chrom, pos, ref, alt = variant_region(vid)
    rec = None
    for r in vcf.fetch(chrom, pos - 1, pos):
        if r.pos == pos and r.ref == ref and list(r.alts or []) == [alt]:
            rec = r
            break
    if rec is None:
        return None
    table = np.zeros((2, 3), dtype=float)
    for sample, data in rec.samples.items():
        if sample in case_set:
            row = 0
        elif sample in ctrl_set:
            row = 1
        else:
            continue
        gt = data.get("GT")
        if gt is None or any(a is None for a in gt):
            continue
        n_alt = sum(1 for a in gt if a == 1)
        if n_alt > 2:
            n_alt = 2
        table[row, n_alt] += 1
    return table


def naive_log10p_for_table(table):
    cols = table[:, table.sum(axis=0) > 0]
    if cols.shape[1] < 2 or (cols.sum(axis=1) == 0).any():
        return 0.0
    chi2, p, _, _ = chi2_contingency(cols, correction=False)
    if p <= 0:
        return float("inf")
    return -math.log10(p)


def main():
    a = parse_args()
    cases = load_id_list(a.cases)
    controls = load_id_list(a.controls)
    regenie = parse_regenie(a.regenie)
    if not regenie:
        sys.exit("ERROR: no usable gene-burden rows parsed from regenie output")
    setlist = load_setlist(a.setlist)
    annos = load_annotations(a.annotations)
    masks = load_masks(a.masks)

    top = sorted(regenie.items(), key=lambda kv: kv[1][0], reverse=True)[:a.top_n]

    if a.vcf_tbi is not None and a.vcf_tbi != f"{a.vcf}.tbi":
        import shutil
        shutil.copy(a.vcf_tbi, f"{a.vcf}.tbi")
    
    vcf = pysam.VariantFile(a.vcf)

    rows = []
    for gene, (regenie_lp, mask) in top:
        mask_conseq = masks.get(mask, set())
        variants = setlist.get(gene, [])
        qualifying = [v for v in variants
                      if annos.get((v, gene), set()) & mask_conseq]
        naive = float("nan")
        n_used = 0
        for v in qualifying:
            tbl = count_genotypes(vcf, v, cases, controls)
            if tbl is None:
                continue
            n_used += 1
            lp = naive_log10p_for_table(tbl)
            if math.isnan(naive) or lp > naive:
                naive = lp
        rows.append((gene, mask, regenie_lp, naive, len(qualifying), n_used))

    print(f"{'GENE':<16}{'MASK':<12}{'REGENIE':>10}{'NAIVE':>10}"
          f"{'NVAR':>6}{'NUSED':>6}  VERDICT")
    violations = 0
    for gene, mask, rlp, nlp, nvar, nused in rows:
        bad = (not math.isfinite(nlp)) or (nlp < a.floor) or (abs(nlp - rlp) > a.tol)
        if bad:
            violations += 1
        nlp_s = f"{nlp:.3f}" if math.isfinite(nlp) else "nan"
        print(f"{gene:<16}{str(mask):<12}{rlp:>10.3f}{nlp_s:>10}"
              f"{nvar:>6}{nused:>6}  {'VIOLATION' if bad else 'ok'}")
    print(f"\nfloor={a.floor} tol={a.tol} max_violations={a.max_violations} "
          f"-> {violations}/{len(rows)} genes violated")

    if a.report_only:
        return 0
    if violations >= a.max_violations:
        print("FAIL: too many top genes lack a corroborating naive signal "
              "(possible wiring regression)")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
