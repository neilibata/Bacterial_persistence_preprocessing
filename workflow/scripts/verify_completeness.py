#!/usr/bin/env python3
"""Verify that every requested accession reached the final matrix.

The pipeline previously had no such check: `combine_aggregated_stats.py` globbed
the aggregated directory, so a sample that failed anywhere upstream simply did
not appear, and the run still reported success. That is how 398 samples went
missing from one reference catalogue without anyone noticing until the analysis
stage, months later.

    python workflow/scripts/verify_completeness.py \
        resources/samples.tsv results/matrices/depth_matrix.csv

Exits non-zero if any requested accession is absent, and writes a report naming
each one so it can be resubmitted.
"""
from __future__ import annotations

import csv
import os
import re
import sys

SRA = re.compile(r"\b([SED]RR\d{6,})\b")


def requested(path):
    """Sample identifiers expected as matrix rows.

    Must mirror `load_samples()` in the Snakefile exactly. For a two-column
    `sample_id <TAB> run1;run2` file the matrix row is the SAMPLE, not its runs
    -- checking for the runs would report every multi-run sample as missing
    while flagging the (correct) sample rows as unexpected.
    """
    ids, seen = [], set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if re.match(r"^sample[_ ]?id\b", line, re.IGNORECASE):
                continue
            parts = re.split(r"[\t,]", line, maxsplit=1)
            if len(parts) == 2 and SRA.search(parts[1]):
                sid = parts[0].strip()
                if sid and sid not in seen:
                    ids.append(sid)
                    seen.add(sid)
                    continue
            for m in SRA.findall(line):
                if m not in seen:
                    ids.append(m)
                    seen.add(m)
    return ids


def in_matrix(path):
    with open(path, newline="") as fh:
        r = csv.reader(fh)
        next(r, None)
        return {row[0].strip() for row in r if row}


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: verify_completeness.py <samples.tsv> <matrix.csv> [report.tsv]")
    want = requested(sys.argv[1])
    got = in_matrix(sys.argv[2])
    report = sys.argv[3] if len(sys.argv) > 3 else "results/matrices/completeness_report.tsv"

    missing = [s for s in want if s not in got]
    extra = sorted(got - set(want))

    os.makedirs(os.path.dirname(report) or ".", exist_ok=True)
    with open(report, "w") as fh:
        fh.write("accession\tstatus\n")
        for s in want:
            fh.write(f"{s}\t{'ok' if s in got else 'MISSING'}\n")
        for s in extra:
            fh.write(f"{s}\tNOT_REQUESTED\n")

    print(f"requested : {len(want):,}")
    print(f"in matrix : {len(got):,}")
    print(f"missing   : {len(missing):,}")
    print(f"unexpected: {len(extra):,}   (present in matrix but not requested)")
    print(f"report    : {report}")

    if extra:
        print(f"\nWARNING: {len(extra)} row(s) were not requested — a stale aggregated/"
              f" directory can carry samples across runs. e.g. {extra[:5]}")
    if missing:
        print(f"\nFAILED: {len(missing)} requested accession(s) never reached the matrix:")
        for s in missing[:20]:
            print(f"  {s}")
        if len(missing) > 20:
            print(f"  ... and {len(missing) - 20} more (see {report})")
        sys.exit(1)
    print("\nPASS - every requested accession is present.")


if __name__ == "__main__":
    main()
