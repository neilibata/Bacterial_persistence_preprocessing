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
    # 7. Save
    # ------------------------------------------------------------------
    outfile = os.path.join(output_dir, f"{sample_id}_depth_and_coverage_stats.csv")
    pd.DataFrame(results).to_csv(outfile, index=False)
    print(f"Done: {sample_id} → {outfile} ({len(results)} species)")


if __name__ == "__main__":
    main()
