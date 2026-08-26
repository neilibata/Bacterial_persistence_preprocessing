"""
aggregate_coverage.py
------------------------------------------------------------------------
Aggregates per-contig samtools coverage output to per-species (SGB) level
using contig-length weighting.

For each species computes:
  - Weighted mean depth  + weighted std deviation
  - Weighted mean coverage (breadth) + weighted std deviation

Contigs in the mapping file but absent from the coverage file are counted
as zeros (no reads mapped), which ensures the denominator (total species
genome length) is always correct.

Usage:
    python aggregate_coverage.py <sample_id> <input_dir> <output_dir> <contig_lengths_file>

Arguments:
    sample_id           e.g. "SRR12345678"
    input_dir           directory containing {sample_id}_high_mapq_coverage.txt
    output_dir          where to write {sample_id}_depth_and_coverage_stats.csv
    contig_lengths_file tab-separated: species_id  contig_name  length
------------------------------------------------------------------------
"""

import os
import sys
import pandas as pd
from collections import defaultdict



# ══════════════════════════════════════════════════════════════════════
#  Amplicon gate
# ══════════════════════════════════════════════════════════════════════
# Random shotgun ties breadth to depth: at genome-wide depth d the expected
# breadth is 1 - e^-d, so at the depths seen here breadth ~= depth and the depth
# over the COVERED bases sits near 1x. Amplicon reads break that. They pile onto
# one conserved window of every reference they match, so the covered fraction
# stays tiny while the depth on it climbs without limit.
#
# This was not hypothetical. 16S libraries registered as WGS in MehtaRS 2018
# reached the analysis undetected: 379 bp covered per reference at up to 45,262x
# on that window, against 1-2x for genuine shotgun. They passed the >= 30 Mb
# mapped-bases filter because their mapped bases are real -- they are simply all
# in one place -- so no depth-based check could have caught them.
#
# Thresholds are empirical, from 24,173 samples in the 2026-08-20 cohort:
#
#   depth per covered base   cohort median 1.77, p99 4.29, and FLAT across every
#                            breadth band (1.68-1.99), exactly as the physics
#                            predicts. Confirmed amplicon samples sat at 4.9.
#   median breadth           the elevated-dpcb subpopulation exists ONLY below
#                            0.02%. Above it the p90 pins at ~2.25 in every band;
#                            below it the p90 is 6.6-10.7.
#
# Both conditions are required. Low breadth alone is just a shallow sample, and
# a high ratio alone occurs in deep well-structured samples.
AMPLICON_MAX_BREADTH_PCT = 0.02   # median breadth across covered species
AMPLICON_MIN_DPCB        = 2.5    # depth per covered base; cohort median is 1.77
QC_FLAG_SUFFIX           = "_AMPLICON_SUSPECT.tsv"


def amplicon_signature(results):
    """Breadth/depth consistency for one sample. None if nothing is covered."""
    import statistics as _st

    breadths, dpcb = [], []
    for r in results:
        cov = r["Weighted_Mean_Coverage"]
        if cov <= 0:
            continue
        breadths.append(cov)
        # depth per covered base = genome-wide depth / covered fraction
        dpcb.append(r["Weighted_Mean_Depth"] / (cov / 100.0))
    if not breadths:
        return None

    mb = _st.median(breadths)
    md = _st.median(dpcb)
    suspect = (mb < AMPLICON_MAX_BREADTH_PCT) and (md > AMPLICON_MIN_DPCB)
    return {
        "median_breadth_pct": mb,
        "dpcb_median": md,
        "n_species": len(breadths),
        "suspect": suspect,
        "verdict": "amplicon-like" if suspect else "shotgun-like",
    }


def main():
    if len(sys.argv) != 5:
        print("Usage: python aggregate_coverage.py <sample_id> <input_dir> <output_dir> <contig_lengths_file>")
        sys.exit(1)

    sample_id           = sys.argv[1]
    input_dir           = sys.argv[2]
    output_dir          = sys.argv[3]
    contig_lengths_file = sys.argv[4]

    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Read contig → species mapping and contig lengths
    # ------------------------------------------------------------------
    contig_to_species = {}
    contig_lengths = {}

    with open(contig_lengths_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()  # fallback to whitespace
            if len(parts) != 3:
                continue
            species_id, contig_name, length = parts
            contig_to_species[contig_name] = species_id
            contig_lengths[contig_name] = int(length)

    print(f"Loaded {len(contig_lengths)} contigs across "
          f"{len(set(contig_to_species.values()))} species")

    # ------------------------------------------------------------------
    # 2. Read samtools coverage output
    # ------------------------------------------------------------------
    input_file = os.path.join(input_dir, f"{sample_id}_high_mapq_coverage.txt")
    if not os.path.exists(input_file):
        print(f"No coverage file found for {sample_id}.")
        sys.exit(0)

    df = pd.read_csv(input_file, sep="\t")
    covered_contigs = set(df["#rname"])

    # ------------------------------------------------------------------
    # 3. First pass: accumulate length-weighted sums
    # ------------------------------------------------------------------
    species_data = defaultdict(lambda: {
        "total_meandepth": 0.0,
        "total_coverage": 0.0,
        "total_length": 0,
        "sq_diff_depth": 0.0,
        "sq_diff_coverage": 0.0,
    })

    for _, row in df.iterrows():
        ctg = row["#rname"]
        sp  = contig_to_species.get(ctg, "UNKNOWN")
        cl  = contig_lengths.get(ctg, 0)

        if pd.notna(row["coverage"]) and row["coverage"] > 0:
            species_data[sp]["total_meandepth"] += row["meandepth"] * cl
            species_data[sp]["total_coverage"]  += row["coverage"]  * cl

        species_data[sp]["total_length"] += cl

    # Add zero-coverage contigs (present in reference but not in coverage file)
    for ctg, length in contig_lengths.items():
        if ctg not in covered_contigs:
            sp = contig_to_species.get(ctg, "UNKNOWN")
            species_data[sp]["total_length"] += length

    # ------------------------------------------------------------------
    # 4. Compute weighted means
    # ------------------------------------------------------------------
    for d in species_data.values():
        tl = d["total_length"]
        d["mean_depth"]    = d["total_meandepth"] / tl if tl > 0 else 0.0
        d["mean_coverage"] = d["total_coverage"]  / tl if tl > 0 else 0.0

    # ------------------------------------------------------------------
    # 5. Second pass: weighted variance (for standard deviation)
    # ------------------------------------------------------------------
    for _, row in df.iterrows():
        ctg = row["#rname"]
        sp  = contig_to_species.get(ctg, "UNKNOWN")
        cl  = contig_lengths.get(ctg, 0)
        d   = species_data[sp]

        d["sq_diff_depth"]    += (row["meandepth"] - d["mean_depth"])**2    * cl
        d["sq_diff_coverage"] += (row["coverage"]  - d["mean_coverage"])**2 * cl

    # ------------------------------------------------------------------
    # 6. Build results table
    # ------------------------------------------------------------------
    results = []
    for sp, d in species_data.items():
        tl = d["total_length"]
        results.append({
            "Species_ID":               sp,
            "Weighted_Mean_Depth":      round(d["mean_depth"], 6),
            "Weighted_Std_Dev_Depth":   round((d["sq_diff_depth"]   / tl)**0.5, 6) if tl > 0 else 0.0,
            "Weighted_Mean_Coverage":   round(d["mean_coverage"], 6),
            "Weighted_Std_Dev_Coverage":round((d["sq_diff_coverage"]/ tl)**0.5, 6) if tl > 0 else 0.0,
        })

    # ------------------------------------------------------------------
    # 7. Amplicon gate
    # ------------------------------------------------------------------
    flags = amplicon_signature(results)
    if flags:
        print(f"  [QC] {sample_id}: {flags['verdict']} — "
              f"median breadth {flags['median_breadth_pct']:.4f}%, "
              f"depth/covered-base {flags['dpcb_median']:.2f}x "
              f"over {flags['n_species']} species")
        if flags["suspect"]:
            print(f"  [QC] {sample_id}: SUSPECTED AMPLICON — excluded from the "
                  f"analysis-ready set. See {QC_FLAG_SUFFIX}.")
            with open(os.path.join(output_dir, f"{sample_id}{QC_FLAG_SUFFIX}"), "w") as fh:
                fh.write("# sample\tmedian_breadth_pct\tdpcb_median\tn_species\tverdict\n")
                fh.write(f"{sample_id}\t{flags['median_breadth_pct']:.6f}\t"
                         f"{flags['dpcb_median']:.4f}\t{flags['n_species']}\t"
                         f"{flags['verdict']}\n")

    # ------------------------------------------------------------------
    # 8. Save
    # ------------------------------------------------------------------
    outfile = os.path.join(output_dir, f"{sample_id}_depth_and_coverage_stats.csv")
    pd.DataFrame(results).to_csv(outfile, index=False)
    print(f"Done: {sample_id} → {outfile} ({len(results)} species)")


if __name__ == "__main__":
    main()
