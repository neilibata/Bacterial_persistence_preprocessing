#!/usr/bin/env python3
"""Merge dated matrix shards into the canonical sample x species matrix.

Each pipeline run writes its own dated shard rather than overwriting a single
matrix. Because column order is fixed by the contig_to_species mapping, shards
are column-identical and can simply be concatenated -- so a backfill of a few
hundred samples never has to touch, or risk, the samples already processed.

    python merge_matrix_shards.py <out_matrix.csv> <shard.csv> [shard.csv ...]

Refuses to merge when it cannot guarantee the result:
  * shard headers differ           -> different species set or order
  * a Sample_ID appears twice      -> ambiguous which run is authoritative
Writes <out_matrix>.provenance.tsv recording which shard each sample came from.
"""
from __future__ import annotations

import csv
import os
import sys


def read_header(path):
    with open(path, newline="") as fh:
        return next(csv.reader(fh), None)


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: merge_matrix_shards.py <out.csv> <shard.csv> [shard.csv ...]")
    out_path = sys.argv[1]
    shards = [s for s in sys.argv[2:] if os.path.exists(s)]
    if not shards:
        sys.exit("ERROR: no shards found")

    # Deterministic order: sort by filename so the merge is reproducible.
    shards = sorted(shards)

    header = read_header(shards[0])
    if header is None:
        sys.exit(f"ERROR: {shards[0]} is empty")

    for s in shards[1:]:
        h = read_header(s)
        if h != header:
            # Column identity is the whole basis for concatenating; if it fails,
            # the shards were built against different contig_to_species files.
            n_diff = sum(1 for a, b in zip(header or [], h or []) if a != b)
            sys.exit(f"ERROR: header mismatch between {shards[0]} and {s}\n"
                     f"       {len(header)} vs {len(h or [])} columns, "
                     f"{n_diff} differing\n"
                     f"       Shards must share a contig_to_species mapping.")

    seen = {}
    rows_written = 0
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(header)
        for s in shards:
            n = 0
            with open(s, newline="") as fh:
                r = csv.reader(fh)
                next(r, None)
                for row in r:
                    if not row:
                        continue
                    sid = row[0].strip()
                    if sid in seen:
                        sys.exit(f"ERROR: {sid} appears in both {seen[sid]} and {s}.\n"
                                 f"       Remove the stale shard, or de-duplicate "
                                 f"before merging.")
                    seen[sid] = s
                    w.writerow(row)
                    n += 1
                    rows_written += 1
            print(f"  {os.path.basename(s):<44} {n:>7,} samples")

    prov = out_path + ".provenance.tsv"
    with open(prov, "w") as pf:
        pf.write("Sample_ID\tshard\n")
        for sid, s in seen.items():
            pf.write(f"{sid}\t{os.path.basename(s)}\n")

    print(f"\n  merged {len(shards)} shard(s) -> {out_path}")
    print(f"  {rows_written:,} samples, {len(header) - 1:,} data columns")
    print(f"  provenance: {prov}")


if __name__ == "__main__":
    main()
