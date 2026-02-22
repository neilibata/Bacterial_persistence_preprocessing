"""
combine_aggregated_stats.py
------------------------------------------------------------------------
Combines all per-sample aggregated coverage CSVs into two sample × species
matrices:

  1. depth_matrix.csv      — Weighted_Mean_Depth ± Weighted_Std_Dev_Depth
  2. coverage_matrix.csv   — Weighted_Mean_Coverage ± Weighted_Std_Dev_Coverage

Each matrix has:
  - Rows    = samples (one per input CSV)
  - Columns = Sample_ID, {species}_mean, {species}_std (in consistent order)

Species order is derived from the contig_to_species mapping file to ensure
all matrices are column-identical regardless of which species appear in
individual samples. Missing species default to 0.0.

Usage:
    python combine_aggregated_stats.py <aggregated_dir> <contig_to_species> <output_dir>

Arguments:
    aggregated_dir      directory containing {sample}_depth_and_coverage_stats.csv
    contig_to_species   tab-separated: species_id  contig_name  length
    output_dir          where to write the two output matrices
------------------------------------------------------------------------
"""

import os
import sys
import csv
from collections import OrderedDict


def load_species_order(contig_to_species_path):
    """Extract unique species IDs in stable order from the mapping file."""
    seen = set()
    species_order = []
    with open(contig_to_species_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) >= 1:
                sp = parts[0]
                if sp not in seen:
                    species_order.append(sp)
                    seen.add(sp)
    return species_order


def parse_sample_csv(filepath, species_set):
    """Read one per-sample aggregated CSV and return dicts of values."""
    depth_mean = {}
    depth_std = {}
    cov_mean = {}
    cov_std = {}

    with open(filepath, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sp = row["Species_ID"]
            if sp not in species_set:
                continue
            depth_mean[sp] = row.get("Weighted_Mean_Depth", "0.0")
            depth_std[sp] = row.get("Weighted_Std_Dev_Depth", "0.0")
            cov_mean[sp] = row.get("Weighted_Mean_Coverage", "0.0")
            cov_std[sp] = row.get("Weighted_Std_Dev_Coverage", "0.0")

    return depth_mean, depth_std, cov_mean, cov_std


def main():
    if len(sys.argv) != 4:
        print("Usage: python combine_aggregated_stats.py <aggregated_dir> <contig_to_species> <output_dir>")
        sys.exit(1)

    aggregated_dir = sys.argv[1]
    contig_to_species_path = sys.argv[2]
    output_dir = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load species order from mapping file
    # ------------------------------------------------------------------
    species_order = load_species_order(contig_to_species_path)
    species_set = set(species_order)
    print(f"Species order: {len(species_order)} species from {contig_to_species_path}")

    # ------------------------------------------------------------------
    # 2. Discover input files
    # ------------------------------------------------------------------
    input_files = sorted(
        f for f in os.listdir(aggregated_dir)
        if f.endswith("_depth_and_coverage_stats.csv") or f.endswith("_aggregated_stats.csv")
    )
    print(f"Found {len(input_files)} sample files in {aggregated_dir}")

    if not input_files:
        print("ERROR: No input files found.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 3. Build column headers
    # ------------------------------------------------------------------
    mean_cols = [f"{sp}_mean" for sp in species_order]
    std_cols = [f"{sp}_std" for sp in species_order]
    header = ["Sample_ID"] + mean_cols + std_cols

    # ------------------------------------------------------------------
    # 4. Process all samples, writing both matrices simultaneously
    # ------------------------------------------------------------------
    depth_out = os.path.join(output_dir, "depth_matrix.csv")
    cov_out = os.path.join(output_dir, "coverage_matrix.csv")

    with open(depth_out, "w", newline="") as df, open(cov_out, "w", newline="") as cf:
        dw = csv.writer(df)
        cw = csv.writer(cf)
        dw.writerow(header)
        cw.writerow(header)

        for i, filename in enumerate(input_files, 1):
            filepath = os.path.join(aggregated_dir, filename)
            sample_id = (
                filename
                .replace("_depth_and_coverage_stats.csv", "")
                .replace("_aggregated_stats.csv", "")
            )

            try:
                d_mean, d_std, c_mean, c_std = parse_sample_csv(filepath, species_set)

                # Build rows in consistent species order (0.0 for missing)
                depth_row = [sample_id] + \
                    [d_mean.get(sp, "0.0") for sp in species_order] + \
                    [d_std.get(sp, "0.0") for sp in species_order]

                cov_row = [sample_id] + \
                    [c_mean.get(sp, "0.0") for sp in species_order] + \
                    [c_std.get(sp, "0.0") for sp in species_order]

                dw.writerow(depth_row)
                cw.writerow(cov_row)

                if i % 50 == 0 or i == len(input_files):
                    print(f"  [{i}/{len(input_files)}] processed", flush=True)

            except Exception as e:
                print(f"  WARNING: error reading {filename}: {e}", flush=True)
                continue

    print(f"\nDone:")
    print(f"  Depth matrix:    {depth_out}")
    print(f"  Coverage matrix: {cov_out}")


if __name__ == "__main__":
    main()

