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


USAGE = (
    "Usage: combine_aggregated_stats.py [--shard NAME] [--allow-shrink]\n"
    "         <contig_to_species> <output_dir> <sample.csv> [sample.csv ...]\n"
    "   or: combine_aggregated_stats.py --scan-dir <aggregated_dir> "
    "<contig_to_species> <output_dir>"
)


def main():
    # Flags are stripped before the arity check: the file list is variable
    # length, so the old `len(sys.argv) != 4` equality test cannot apply.

    # Explicit file list is the safe path: the caller (Snakemake) states exactly
    # which samples must be present, and a missing one is an error rather than a
    # silently smaller matrix. --scan-dir reproduces the old globbing behaviour
    # and is deliberately opt-in, because a stale directory silently changes the
    # matrix contents.
    allow_shrink = "--allow-shrink" in sys.argv
    if allow_shrink:
        sys.argv.remove("--allow-shrink")

    # --shard NAME writes depth_matrix_NAME.csv instead of depth_matrix.csv, so a
    # run never overwrites another run's output. Shards are merged separately.
    shard = None
    if "--shard" in sys.argv:
        i = sys.argv.index("--shard")
        shard = sys.argv[i + 1]
        del sys.argv[i:i + 2]

    if len(sys.argv) < 4:
        print(USAGE)
        sys.exit(1)

    if sys.argv[1] == "--scan-dir":
        aggregated_dir = sys.argv[2]
        contig_to_species_path = sys.argv[3]
        output_dir = sys.argv[4]
        explicit_files = None
        print("WARNING: --scan-dir combines whatever is in the directory, which may\n"
              "         include samples from earlier runs and may omit samples that\n"
              "         failed. Prefer passing an explicit file list.", flush=True)
    else:
        contig_to_species_path = sys.argv[1]
        output_dir = sys.argv[2]
        explicit_files = sys.argv[3:]
        aggregated_dir = None
        if not explicit_files:
            print("ERROR: no sample CSVs given")
            sys.exit(1)

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
    if explicit_files is not None:
        input_paths = list(explicit_files)
        missing = [f for f in input_paths if not os.path.exists(f)]
        if missing:
            print(f"ERROR: {len(missing)} declared sample file(s) do not exist, e.g. "
                  f"{missing[:3]}")
            sys.exit(1)
        print(f"Combining {len(input_paths)} sample files (explicit list)")
    else:
        input_paths = [
            os.path.join(aggregated_dir, f) for f in sorted(os.listdir(aggregated_dir))
            if f.endswith("_depth_and_coverage_stats.csv")
            or f.endswith("_aggregated_stats.csv")
        ]
        print(f"Found {len(input_paths)} sample files in {aggregated_dir}")

    if not input_paths:
        print("ERROR: No input files found.")
        sys.exit(1)
    n_expected = len(input_paths)

    # ------------------------------------------------------------------
    # 3. Build column headers
    # ------------------------------------------------------------------
    mean_cols = [f"{sp}_mean" for sp in species_order]
    std_cols = [f"{sp}_std" for sp in species_order]
    header = ["Sample_ID"] + mean_cols + std_cols

    # ------------------------------------------------------------------
    # 4. Process all samples, writing both matrices simultaneously
    # ------------------------------------------------------------------
    suffix = f"_{shard}" if shard else ""
    depth_out = os.path.join(output_dir, f"depth_matrix{suffix}.csv")
    cov_out = os.path.join(output_dir, f"coverage_matrix{suffix}.csv")

    # A shard is expected to be small -- it holds only this run's samples -- so
    # the shrinkage guard applies only when writing an unsharded matrix in place.
    prev_n = None
    if not shard and os.path.exists(depth_out):
        try:
            with open(depth_out) as pf:
                prev_n = max(0, sum(1 for _ in pf) - 1)
            print(f"Existing matrix has {prev_n:,} samples")
        except Exception:
            prev_n = None

    with open(depth_out, "w", newline="") as df, open(cov_out, "w", newline="") as cf:
        dw = csv.writer(df)
        cw = csv.writer(cf)
        dw.writerow(header)
        cw.writerow(header)

        written = []
        for i, filepath in enumerate(input_paths, 1):
            filename = os.path.basename(filepath)
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
                written.append(sample_id)

                if i % 50 == 0 or i == n_expected:
                    print(f"  [{i}/{n_expected}] processed", flush=True)

            except Exception as e:
                # Previously this warned and continued, so an unreadable sample
                # vanished from the matrix without failing the run. A matrix
                # missing samples is worse than no matrix.
                print(f"ERROR reading {filepath}: {e}", flush=True)
                sys.exit(1)

    if len(written) != n_expected:
        print(f"ERROR: wrote {len(written)} rows but {n_expected} were expected")
        sys.exit(1)

    # Shrinkage guard. Running with a reduced samples.tsv -- e.g. only the
    # accessions being backfilled -- rebuilds the matrix from just those, quietly
    # discarding every previously processed sample. The declared-input fix makes
    # that MORE likely, not less, so refuse to shrink an existing matrix unless
    # asked. To backfill, list ALL samples in samples.tsv: Snakemake skips the
    # ones whose aggregated CSV already exists and only runs the missing ones.
    if prev_n is not None and len(written) < prev_n and not allow_shrink:
        print(f"\nERROR: the existing matrix has {prev_n:,} samples but this run "
              f"would write only {len(written):,}.")
        print("       Refusing to overwrite. Either list all samples in samples.tsv")
        print("       (Snakemake will skip the ones already done), or pass")
        print("       --allow-shrink if you really mean to replace the matrix.")
        sys.exit(1)

    manifest = os.path.join(output_dir, f"matrix_manifest{suffix}.tsv")
    with open(manifest, "w") as mf:
        mf.write("Sample_ID\tsource_file\n")
        for sid, path in zip(written, input_paths):
            mf.write(f"{sid}\t{path}\n")

    print(f"\nDone:")
    print(f"  Depth matrix:    {depth_out}")
    print(f"  Coverage matrix: {cov_out}")
    print(f"  Manifest:        {manifest}  ({len(written)} samples)")


if __name__ == "__main__":
    main()

