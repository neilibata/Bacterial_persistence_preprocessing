#!/usr/bin/env python3
"""Relate alignment rate to depth and detection, per sample and per study.

A low overall alignment rate on its own says nothing. A shallow sample and a
contaminated one both detect few taxa, and comparing one sample against another
picked by hand is worthless -- the first attempt at this compared a 4.59%-aligning
sample against a 121 MB one and concluded, meaninglessly, that both detected six
SGBs.

So this reports four numbers together for every sample:

    read pairs      was there enough input to detect anything
    alignment %     did the reads match the catalogue
    SGBs > 10%      how many taxa cleared the main-text detection threshold
    summed breadth  total genome fraction covered across the catalogue

and groups them by study, because the question is usually whether a whole cohort
sits below the rest (the catalogue covers that population poorly) or whether a
few samples within it are individually broken (contamination, isolates, wrong
body site). Those have opposite consequences and the same symptom.

    python workflow/scripts/alignment_qc.py
    python workflow/scripts/alignment_qc.py --ref lev --out qc_lev.csv
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import statistics as st
import sys


def collect(results_dir: str, logs_dir: str, ref: str):
    rows = []
    pat = os.path.join(logs_dir, ref, "03_align", "*.log")
    for f in sorted(glob.glob(pat)):
        sample = os.path.basename(f)[:-4]
        try:
            txt = open(f, errors="replace").read()
        except OSError:
            continue
        m = re.search(r"([\d]+) reads; of these", txt)
        a = re.search(r"([\d.]+)% overall alignment rate", txt)
        if not (m and a):
            continue
        agg = os.path.join(results_dir, "aggregated", ref,
                           f"{sample}_depth_and_coverage_stats.csv")
        n10 = 0
        summed = 0.0
        if os.path.exists(agg):
            with open(agg) as fh:
                for r in csv.DictReader(fh):
                    try:
                        cov = float(r["Weighted_Mean_Coverage"])
                    except (KeyError, ValueError):
                        continue
                    summed += cov
                    if cov > 10:
                        n10 += 1
        rows.append({
            "sample": sample,
            "read_pairs": int(m.group(1)),
            "align_pct": float(a.group(1)),
            "sgbs_over_10pct": n10,
            "summed_breadth": summed,
        })
    return rows


def study_of(sample: str) -> str:
    """Coarse study key from the accession prefix -- enough to group by."""
    if sample.startswith("SAMN0233"):
        return "PRJNA217052 (Brito Fiji)"
    for pre in ("SAMN", "SAMEA", "ERR", "SRR", "ERS", "SRS"):
        if sample.startswith(pre):
            return f"{pre}* other"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--logs-dir", default="logs")
    ap.add_argument("--ref", default="UHGG")
    ap.add_argument("--low", type=float, default=50.0,
                    help="alignment %% below which a sample is called low")
    ap.add_argument("--out", default=None, help="write the full table here")
    ap.add_argument("--exclude-out", default=None,
                    help="write the low samples as a one-per-line exclusion list")
    args = ap.parse_args()

    rows = collect(args.results_dir, args.logs_dir, args.ref)
    if not rows:
        sys.exit(f"no alignment logs under {args.logs_dir}/{args.ref}/03_align/")
    rows.sort(key=lambda r: r["align_pct"])

    al = [r["align_pct"] for r in rows]
    n10 = [r["sgbs_over_10pct"] for r in rows]
    br = [r["summed_breadth"] for r in rows]
    rp = [r["read_pairs"] for r in rows]
    print(f"  reference {args.ref}: {len(rows)} samples\n")
    print(f"  {'metric':<18}{'min':>12}{'median':>12}{'max':>12}")
    for name, v in (("alignment %", al), ("read pairs", rp),
                    ("SGBs > 10%", n10), ("summed breadth", br)):
        print(f"  {name:<18}{min(v):>12,.0f}{st.median(v):>12,.0f}{max(v):>12,.0f}")

    hdr = f"\n  {'sample':<16}{'read pairs':>13}{'align%':>8}{'SGB>10%':>9}{'breadth':>10}"
    print(hdr)
    print(f"  --- 12 lowest alignment ---")
    for r in rows[:12]:
        print(f"  {r['sample']:<16}{r['read_pairs']:>13,}{r['align_pct']:>8.2f}"
              f"{r['sgbs_over_10pct']:>9}{r['summed_breadth']:>10.0f}")
    mid = len(rows) // 2
    print(f"  --- 6 at the median, for comparison ---")
    for r in rows[mid:mid + 6]:
        print(f"  {r['sample']:<16}{r['read_pairs']:>13,}{r['align_pct']:>8.2f}"
              f"{r['sgbs_over_10pct']:>9}{r['summed_breadth']:>10.0f}")

    # The decisive comparison: is a whole study low, or only part of one?
    print(f"\n  === by study: is the whole cohort low, or only some of it? ===")
    print(f"  {'study':<28}{'n':>5}{'n_low':>7}{'align med':>11}"
          f"{'SGB>10% med':>13}{'breadth med':>13}")
    groups = {}
    for r in rows:
        groups.setdefault(study_of(r["sample"]), []).append(r)
    for key, g in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        low = [r for r in g if r["align_pct"] < args.low]
        print(f"  {key:<28}{len(g):>5}{len(low):>7}"
              f"{st.median([r['align_pct'] for r in g]):>11.1f}"
              f"{st.median([r['sgbs_over_10pct'] for r in g]):>13.0f}"
              f"{st.median([r['summed_breadth'] for r in g]):>13.0f}")

    # Within the study that carries the low tail, split it and compare.
    brito = [r for r in rows if study_of(r["sample"]).startswith("PRJNA217052")]
    if brito:
        lo = [r for r in brito if r["align_pct"] < args.low]
        hi = [r for r in brito if r["align_pct"] >= args.low]
        print(f"\n  === Brito Fiji, split at {args.low:.0f}% alignment ===")
        for label, g in (("low", lo), ("rest", hi)):
            if not g:
                continue
            print(f"    {label:<6}n={len(g):<4} "
                  f"reads med {st.median([r['read_pairs'] for r in g]):>12,.0f}   "
                  f"align med {st.median([r['align_pct'] for r in g]):>6.1f}%   "
                  f"SGB>10% med {st.median([r['sgbs_over_10pct'] for r in g]):>5.0f}   "
                  f"breadth med {st.median([r['summed_breadth'] for r in g]):>7.0f}")
        print("\n    If 'rest' matches the cohort medians above, the low ones are")
        print("    individually defective and can be dropped. If both sit below the")
        print("    cohort, the catalogue covers this population poorly and the whole")
        print("    study would join with depressed richness.")

    if args.out:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"\n  wrote {args.out} ({len(rows)} rows)")

    if args.exclude_out:
        low = [r for r in rows if r["align_pct"] < args.low]
        with open(args.exclude_out, "w") as fh:
            fh.write(f"# samples aligning below {args.low:.0f}% against {args.ref}\n")
            fh.write("# sample\tread_pairs\talign_pct\tsgbs_over_10pct\n")
            for r in low:
                fh.write(f"{r['sample']}\t{r['read_pairs']}\t{r['align_pct']}\t"
                         f"{r['sgbs_over_10pct']}\n")
        print(f"  wrote {args.exclude_out} ({len(low)} samples below {args.low:.0f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
