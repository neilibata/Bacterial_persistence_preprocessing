#!/usr/bin/env python3
"""Fetch every run of a sample and concatenate the paired-end ones.

Why this exists rather than a shell loop around fasterq-dump:

`fasterq-dump` reconstructs FASTQs from SRA's "Normalized Format", and whether
it emits `_1`/`_2` depends on what that reconstruction produced. Inferring a
run's layout from which files appeared therefore conflates three different
outcomes -- genuinely single-end, paired but reconstructed oddly, and a failed
download -- and the first version of this step silently skipped all three
alike. A sample could lose most of its data and still exit zero.

ENA publishes `library_layout` per run, so the layout is *looked up*, not
guessed, and ENA serves the FASTQs directly (no reconstruction, md5 published).
SRA remains the fallback when ENA has no record or the transfer fails.

Outcomes are kept distinct:
  paired  -> downloaded and concatenated
  single  -> excluded from both mates, recorded in the skip report
  error   -> the sample FAILS; it is never treated as single-end

    python fetch_runs.py --sample S1 --runs ERR1,ERR2 \\
        --outdir results/fastq --tmpdir results/fastq/.tmp_S1 \\
        --qc-dir results/qc --threads 6
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = "run_accession,library_layout,fastq_ftp,fastq_md5,read_count"
TIMEOUT = 60


def log(msg: str) -> None:
    print(msg, flush=True)


# ─────────────────────────────────────────────────────────────── ENA metadata
def ena_metadata(run: str, retries: int = 3) -> dict | None:
    """Layout and file URLs for a run, or None if ENA has no usable record."""
    url = (f"{ENA_API}?accession={run}&result=read_run"
           f"&fields={ENA_FIELDS}&format=tsv")
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT) as fh:
                text = fh.read().decode()
            break
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            log(f"    ENA query attempt {attempt}/{retries} failed: {exc}")
            if attempt == retries:
                return None
    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    layout = (row.get("library_layout") or "").strip().upper()
    urls = [u for u in (row.get("fastq_ftp") or "").split(";") if u]
    md5s = [m for m in (row.get("fastq_md5") or "").split(";") if m]
    return {"layout": layout, "urls": urls, "md5s": md5s,
            "read_count": row.get("read_count", "")}


# ─────────────────────────────────────────────────────────────── downloading
def md5sum(path: Path, chunk: int = 8 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def count_lines(path: Path) -> int:
    """`wc -l`, falling back to Python only if wc is unavailable."""
    try:
        out = subprocess.run(["wc", "-l", str(path)],
                             capture_output=True, text=True, check=True)
        return int(out.stdout.split()[0])
    except (OSError, subprocess.CalledProcessError, ValueError, IndexError):
        with open(path, "rb") as fh:
            return sum(1 for _ in fh)


def curl(url: str, dest: Path) -> bool:
    """Fetch over HTTPS. ENA's paths are given host-relative."""
    if not url.startswith(("http://", "https://", "ftp://")):
        url = "https://" + url
    cmd = ["curl", "-sSL", "--retry", "5", "--retry-delay", "10",
           "--connect-timeout", "60", "-o", str(dest), url]
    return subprocess.run(cmd).returncode == 0


def fetch_from_ena(run: str, tmp: Path) -> tuple[str, list[Path]]:
    """(status, files) where status is paired | single | unavailable."""
    meta = ena_metadata(run)
    if meta is None:
        log(f"    ENA: no record for {run}")
        return "unavailable", []

    layout, urls, md5s = meta["layout"], meta["urls"], meta["md5s"]
    log(f"    ENA: layout={layout or '?'} files={len(urls)} "
        f"reads={meta['read_count'] or '?'}")

    # Paired requires the two mate files. A PAIRED record with a single file
    # (interleaved submission) is not something this pipeline can consume, so
    # treat it as unavailable and let SRA try rather than silently mis-handling.
    mates = [u for u in urls if u.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
    if layout == "SINGLE" and len(mates) < 2:
        return "single", []
    if len(mates) < 2:
        log(f"    ENA: layout={layout} but {len(mates)} mate file(s); "
            f"falling back to SRA")
        return "unavailable", []

    mates.sort(key=lambda u: u.endswith("_2.fastq.gz"))
    want_md5 = {u: m for u, m in zip(urls, md5s)} if len(md5s) == len(urls) else {}

    out: list[Path] = []
    for url in mates:
        dest = tmp / os.path.basename(url)
        if not curl(url, dest):
            log(f"    ENA: download failed for {url}")
            for p in out:
                p.unlink(missing_ok=True)
            dest.unlink(missing_ok=True)
            return "unavailable", []
        exp = want_md5.get(url)
        if exp:
            got = md5sum(dest)
            if got != exp:
                log(f"    ENA: md5 mismatch for {dest.name} "
                    f"(expected {exp}, got {got})")
                for p in out:
                    p.unlink(missing_ok=True)
                dest.unlink(missing_ok=True)
                return "unavailable", []
        out.append(dest)
    log(f"    ENA: retrieved {len(out)} mate files, md5 verified")
    return "paired", out


def fetch_from_sra(run: str, tmp: Path, threads: int,
                   max_size: str) -> tuple[str, list[Path]]:
    """Fallback path. Returns paired | single | error."""
    if shutil.which("prefetch") is None or shutil.which("fasterq-dump") is None:
        log("    SRA: sra-tools not on PATH")
        return "error", []
    if subprocess.run(["prefetch", "--max-size", max_size, run],
                      cwd=tmp).returncode != 0:
        log(f"    SRA: prefetch failed for {run}")
        return "error", []
    if subprocess.run(["fasterq-dump", run, "--split-3",
                       "--threads", str(threads), "--outdir", str(tmp)],
                      cwd=tmp).returncode != 0:
        log(f"    SRA: fasterq-dump failed for {run}")
        return "error", []

    r1, r2 = tmp / f"{run}_1.fastq", tmp / f"{run}_2.fastq"
    orphan = tmp / f"{run}.fastq"
    shutil.rmtree(tmp / run, ignore_errors=True)

    if r1.exists() and r2.exists():
        # `wc -l` in C, not a Python loop: these files run to tens of millions
        # of lines and iterating them in Python stalls the job for minutes.
        n1, n2 = count_lines(r1), count_lines(r2)
        if n1 != n2:
            log(f"    SRA: mate files differ ({n1 // 4} vs {n2 // 4} reads)")
            return "error", []
        orphan.unlink(missing_ok=True)
        gz = ["pigz", "-p", str(threads)] if shutil.which("pigz") else ["gzip"]
        for p in (r1, r2):
            subprocess.run(gz + [str(p)], check=True)
        log(f"    SRA: retrieved paired output ({n1 // 4} reads/mate)")
        return "paired", [Path(f"{r1}.gz"), Path(f"{r2}.gz")]

    if orphan.exists():
        # Size only -- the exact count is informational and not worth a full
        # pass over a multi-gigabyte file we are about to delete.
        gb = orphan.stat().st_size / 1e9
        orphan.unlink()
        log(f"    SRA: single-end output ({gb:.1f} GB), no pairs")
        return "single", []

    log(f"    SRA: produced no FASTQ output for {run}")
    return "error", []


# ────────────────────────────────────────────────────────────────────── main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--runs", required=True, help="comma- or space-separated")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--tmpdir", required=True, type=Path)
    ap.add_argument("--qc-dir", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--max-size", default="300G")
    ap.add_argument("--prefer", choices=["ena", "sra"], default="ena")
    args = ap.parse_args()

    runs = [r for r in args.runs.replace(",", " ").split() if r]

    # Absolute paths are mandatory: prefetch and fasterq-dump are run with
    # cwd=tmpdir, so a relative --outdir resolves *inside* tmpdir and the output
    # lands somewhere this script never looks. That produced a spurious
    # "produced no FASTQ output" for a run that had downloaded fine.
    args.outdir = args.outdir.resolve()
    args.tmpdir = args.tmpdir.resolve()
    args.qc_dir = args.qc_dir.resolve()

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.qc_dir.mkdir(parents=True, exist_ok=True)
    if args.tmpdir.exists():
        shutil.rmtree(args.tmpdir)
    args.tmpdir.mkdir(parents=True)

    log(f"=== {args.sample}: {len(runs)} run(s) ===")
    log(f"runs: {' '.join(runs)}")

    kept: list[tuple[str, list[Path]]] = []
    skipped: list[tuple[str, str]] = []

    for run in runs:
        log(f"--- {run} ---")
        order = ([fetch_from_ena, fetch_from_sra] if args.prefer == "ena"
                 else [fetch_from_sra, fetch_from_ena])

        # Try every source for PAIRED before concluding single-end. SRA stores
        # some runs with each second mate zero-length, so it reports single-end
        # for a run ENA holds as properly paired (ERR1913111). Stopping at the
        # first "single" would discard real paired data.
        status, files, seen = "error", [], []
        for fn in order:
            if fn is fetch_from_ena:
                status, files = fetch_from_ena(run, args.tmpdir)
            else:
                status, files = fetch_from_sra(run, args.tmpdir,
                                               args.threads, args.max_size)
            seen.append(status)
            if status == "paired":
                break
        if status != "paired":
            status = "single" if "single" in seen else "error"

        if status == "paired":
            kept.append((run, files))
        elif status == "single":
            log(f"    SKIP: {run} is single-end; excluded from both mates")
            skipped.append((run, "single_end"))
        else:
            # Never silently degrade a broken download into "single-end".
            log(f"ERROR: {run} could not be retrieved from ENA or SRA")
            return 1

    skip_report = args.qc_dir / f"{args.sample}.skipped_runs.tsv"
    with open(skip_report, "w") as fh:
        for run, why in skipped:
            fh.write(f"{args.sample}\t{run}\t{why}\n")

    if not kept:
        log(f"ERROR: {args.sample} has no paired-end runs "
            f"({len(skipped)} single-end)")
        return 1

    # Same run set, same order, for both mates -- that is what keeps the
    # concatenated files in register. gzip members concatenate losslessly.
    r1_out = args.outdir / f"{args.sample}_1.fastq.gz"
    r2_out = args.outdir / f"{args.sample}_2.fastq.gz"
    for out, idx in ((r1_out, 0), (r2_out, 1)):
        part = Path(f"{out}.part")
        with open(part, "wb") as dst:
            for _run, files in kept:
                with open(files[idx], "rb") as src:
                    shutil.copyfileobj(src, dst, 8 << 20)
        part.replace(out)

    shutil.rmtree(args.tmpdir, ignore_errors=True)

    if r1_out.stat().st_size == 0 or r2_out.stat().st_size == 0:
        log(f"ERROR: empty output for {args.sample}")
        return 1

    log(f"=== Done {args.sample}: {len(kept)} paired run(s) concatenated, "
        f"{len(skipped)} single-end skipped ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
