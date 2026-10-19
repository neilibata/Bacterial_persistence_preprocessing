#!/usr/bin/env python3
"""Rebuild a reference catalogue from its public sources, byte for byte.

Why this exists
---------------
The four reference catalogues this pipeline aligns against are each tens of
gigabytes, and three of the four are other people's published data. Shipping
copies would mean redistributing their work and would put ~33 GB in a deposit
that nobody can diff. Shipping the *recipe* is smaller, cites them properly,
and is a stronger reproducibility claim: a checksum proves the rebuilt file is
the one the analysis read, where a copy only asserts it.

So the catalogues are not distributed. This script fetches each component from
the archive its authors published it in, renames contigs by the rule recorded
in the manifest, concatenates them in the recorded order, and verifies the
result against per-contig names and lengths. A mismatch stops the build rather
than producing a reference that is subtly not the reported one.

Sources, and what each is
-------------------------
  carter    3,971 genomes, two provenances in one catalogue:
              1,816 GCA_* GenBank assemblies (Carter et al.'s Hadza MAGs),
                    fetched per accession from the NCBI FTP site;
              2,155 MGYG*  MGnify genomes (UHGG), fetched per accession from
                    the EMBL-EBI MGnify API, catalogue human-gut-v2-0-2.
  leviatan  3,594 species representatives, one file:
              Rep_all.fa.gz from figshare 10.6084/m9.figshare.16885261
              (Leviatan et al. 2022, CC BY 4.0). Their 142,912-assembly
              Zenodo deposit is access-restricted by IRB requirement and is
              NOT needed: the representatives are the open release.
  chimp     1,682 dereplicated MAGs:
              ChimpMAGs.tar from Dryad 10.5061/dryad.00000006x
              (Sanders et al. 2023, CC0-1.0).
  combined  5,081 genomes, all of leviatan plus the 1,487 carter genomes that
            dereplication did not fold into a leviatan representative. Built
            from the other two sources; downloads nothing of its own.

Contig naming, which differs per provenance and is not guessable
----------------------------------------------------------------
  MGYG*  used verbatim. MGnify already names contigs MGYG000004406_1,
         _2, ... and the catalogue keeps them, which is why these indices
         start at 1 while the renamed ones start at 0.
  GCA_*  renamed to <accession>_c_<i>, i from 0, in source file order.
  Rep_*  renamed to <Rep_N>_c_<i>, i from 0, in source file order.
  WC*    kept as <SGB>|c_<012d> as the merge wrote them.

Three MGYG genomes (MGYG000001406.1, MGYG000001455.1, MGYG000001456.1) are
missing a few contigs relative to their source, and one more carries a `.1`
version suffix the rest do not. These are not errors to correct: they are what
the reported run aligned against. The manifest lists contigs explicitly for
exactly this reason -- the build selects what the manifest names rather than
re-deriving a rule and silently disagreeing.

Usage
-----
    # prove the method without a 33 GB rebuild: first 25 genomes, then diff
    python build_references.py --ref carter --subset 25 --out proof/

    # the real thing
    python build_references.py --ref combined --out refs/ --index

    # check a reference you already have
    python build_references.py --ref carter --verify-only refs/carter.fna
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
#: Manifests live beside the workflow, not in the deposit: they are a few MB
#: and version with the code that reads them.
MANIFESTS = HERE.parent.parent / "resources" / "reference_manifests"

NCBI_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA"
MGNIFY = "https://www.ebi.ac.uk/metagenomics/api/v1/genomes"
FIGSHARE = "https://api.figshare.com/v2/articles/16885261"
DRYAD_VERSION = "https://datadryad.org/api/v2/versions/228587/files"

#: Sources, for the record and for --cite.
SOURCES = {
    "carter": [
        ("NCBI GenBank assemblies (GCA_*)", NCBI_FTP, "open"),
        ("MGnify human-gut-v2-0-2 (MGYG*)", MGNIFY, "open"),
    ],
    "leviatan": [("figshare 10.6084/m9.figshare.16885261 :: Rep_all.fa.gz",
                  FIGSHARE, "CC BY 4.0")],
    "chimp": [("Dryad 10.5061/dryad.00000006x :: ChimpMAGs.tar",
               DRYAD_VERSION, "CC0-1.0")],
    "combined": [("built from carter and leviatan sources", "-", "-")],
}


# ─────────────────────────────── manifests ────────────────────────────────
def load_manifest(ref: str):
    """(genome rows, [(contig_name, length)]) in reference order."""
    g = MANIFESTS / f"{ref}.genomes.tsv"
    c = MANIFESTS / f"{ref}.contigs.tsv.gz"
    for p in (g, c):
        if not p.exists():
            sys.exit(f"  manifest missing: {p}\n"
                     f"  (expected under {MANIFESTS})")
    genomes = []
    with open(g) as fh:
        next(fh)
        for line in fh:
            order, genome, scheme, n_contigs, total_bp = \
                line.rstrip("\n").split("\t")
            genomes.append({"order": int(order), "genome": genome,
                            "scheme": scheme, "n_contigs": int(n_contigs),
                            "total_bp": int(total_bp)})
    contigs = []
    with gzip.open(c, "rt") as fh:
        next(fh)
        for line in fh:
            name, length = line.rstrip("\n").split("\t")
            contigs.append((name, int(length)))
    return genomes, contigs


# ──────────────────────────────── fetching ────────────────────────────────
def _get(url: str, timeout: int = 180) -> bytes:
    req = urllib.request.Request(url, headers={"Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def fetch_gca(acc: str, cache: Path) -> Path:
    """One GenBank assembly FASTA, by accession.

    The FTP layout splits the numeric part into three-digit directories and
    then appends the submitter's own assembly name, which is not derivable
    from the accession -- so the parent directory is listed to find it. One
    request per genome; the result is cached because 1,816 of these is the
    slow part of a carter rebuild.
    """
    out = cache / f"{acc}.fna.gz"
    if out.exists() and out.stat().st_size:
        return out
    digits = acc.split("_")[1].split(".")[0]
    parent = f"{NCBI_FTP}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/"
    listing = _get(parent, timeout=120).decode("utf8", "replace")
    names = sorted(set(re.findall(
        re.escape(acc) + r"_[^/\"'<>\s]+", listing)))
    if not names:
        raise FileNotFoundError(f"{acc} not found under {parent}")
    name = names[0]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(_get(f"{parent}{name}/{name}_genomic.fna.gz"))
    return out


def fetch_mgyg(acc: str, cache: Path) -> Path:
    """One MGnify genome FASTA, by accession."""
    out = cache / f"{acc}.fna"
    if out.exists() and out.stat().st_size:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    # The accession in the catalogue may carry a version suffix that MGnify's
    # own API does not use; strip it for the request, keep it for the name.
    api_acc = acc.split(".")[0]
    out.write_bytes(_get(f"{MGNIFY}/{api_acc}/downloads/{api_acc}.fna"))
    return out


def fetch_bulk(ref: str, cache: Path) -> Path:
    """The single-archive sources: leviatan from figshare, chimp from Dryad."""
    cache.mkdir(parents=True, exist_ok=True)
    if ref == "leviatan":
        out = cache / "Rep_all.fa.gz"
        if out.exists() and out.stat().st_size > 2_000_000_000:
            return out
        meta = json.loads(_get(FIGSHARE, timeout=120))
        url = next(f["download_url"] for f in meta["files"]
                   if f["name"] == "Rep_all.fa.gz")
    else:
        out = cache / "ChimpMAGs.tar"
        if out.exists() and out.stat().st_size > 8_000_000_000:
            return out
        meta = json.loads(_get(DRYAD_VERSION + "?per_page=100", timeout=120))
        f = next(x for x in meta["_embedded"]["stash:files"]
                 if x["path"] == "ChimpMAGs.tar")
        url = "https://datadryad.org" + \
            f["_links"]["stash:file-download"]["href"]
    print(f"    fetching {out.name} ({url[:70]}...)")
    with urllib.request.urlopen(url, timeout=1800) as r, open(out, "wb") as w:
        while True:
            b = r.read(1 << 22)
            if not b:
                break
            w.write(b)
    return out


# ───────────────────────────── reading sources ────────────────────────────
def read_fasta(path_or_bytes) -> list[tuple[str, str]]:
    """[(name, sequence)] in file order, sequence uppercased and unwrapped."""
    if isinstance(path_or_bytes, (bytes, bytearray)):
        fh = io.StringIO(path_or_bytes.decode("utf8", "replace"))
    elif str(path_or_bytes).endswith(".gz"):
        fh = gzip.open(path_or_bytes, "rt")
    else:
        fh = open(path_or_bytes)
    recs, name, buf = [], None, []
    with fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip().upper())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def source_records(ref: str, g: dict, cache: Path, bulk=None):
    """Renamed contigs for one per-accession genome, in source order.

    Only the per-accession schemes come through here; the single-archive ones
    are streamed by stage_bulk() because holding them in memory does not fit.
    """
    scheme, acc = g["scheme"], g["genome"]
    if scheme == "mgyg":
        # Verbatim: MGnify's own contig names are what the catalogue carries.
        return read_fasta(fetch_mgyg(acc, cache / "mgyg"))
    if scheme == "gca":
        recs = read_fasta(fetch_gca(acc, cache / "gca"))
        return [(f"{acc}_c_{i}", s) for i, (_, s) in enumerate(recs)]
    raise ValueError(f"{scheme!r} is streamed, not fetched per genome")


def _write_records(path: Path, recs) -> None:
    with open(path, "a") as w:
        for name, seq in recs:
            w.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                w.write(seq[i:i + 60] + "\n")


def stage_bulk(ref: str, cache: Path, stage: Path, need: set) -> None:
    """Stream a single-archive source, staging one file per needed genome.

    The archives are 2.7-8.7 GB compressed and expand past 8 GB, so reading
    one into a dict and then keeping the selected contigs as well wanted
    ~18 GB of resident memory -- fine on the machine this was written on and
    not on a 16 GB laptop. Streaming to one small file per genome keeps the
    peak at a single genome, at the cost of one extra copy on disk, which is
    deleted once the catalogue is concatenated.
    """
    stage.mkdir(parents=True, exist_ok=True)
    if ref in ("leviatan", "combined"):
        p = fetch_bulk("leviatan", cache / "bulk")
        print(f"    streaming {p.name} -> {len(need):,} staged genome(s)")
        seen = collections_counter()
        name = None
        buf: list[str] = []
        kept = set()

        def flush():
            if name is None:
                return
            m = re.match(r"^(Rep_\d+)", name)
            if not m or m.group(1) not in need:
                return
            g = m.group(1)
            i = seen[g]
            seen[g] = i + 1
            kept.add(g)
            _write_records(stage / f"{g}.fa", [(f"{g}_c_{i}", "".join(buf))])

        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    flush()
                    name = line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line.strip().upper())
            flush()
        print(f"    staged {len(kept):,} of {len(need):,} needed genome(s)")
        return

    if ref == "chimp":
        # ChimpMAGs.tar holds 3,950 MAGs under two naming conventions -- 225
        # as WC<n>_<n>.fa.contigs.fasta.fa and the rest as
        # WILDCHIMP_<plate>_WC<n>_S<n>__DASTool_bins.<binner>.<n>.fa -- and
        # the catalogue keeps whichever name the MAG arrived with. The 1,682
        # the manifest names are dRep's cluster winners; the other 2,268 are
        # skipped without being read into memory.
        #
        # The contig rule is a prefix, not a renumbering: <MAG>|<contig as it
        # is in the MAG file>. That is why the two conventions leave different
        # contig names in the catalogue -- c_000000000001 for the WC set,
        # k141_<n> for the DASTool set -- and why re-indexing would not
        # reproduce it.
        p = fetch_bulk("chimp", cache / "bulk")
        print(f"    streaming {p.name} -> {len(need):,} staged genome(s)")
        kept = 0
        with tarfile.open(p, "r|") as tf:          # r| : stream, never seek
            for m in tf:
                if not m.isfile():
                    continue
                mag = re.sub(r"\.fa(\.contigs\.fasta\.fa)?$", "",
                             m.name.split("/")[-1])
                if mag not in need:
                    continue
                fo = tf.extractfile(m)
                if fo is None:
                    continue
                _write_records(stage / f"{mag}.fa",
                               [(f"{mag}|{n}", s)
                                for n, s in read_fasta(fo.read())])
                kept += 1
                if kept % 400 == 0:
                    print(f"      {kept:,}/{len(need):,}")
        print(f"    staged {kept:,} of {len(need):,} needed genome(s)")
        return


def collections_counter():
    import collections
    return collections.Counter()



# ──────────────────────────────── building ────────────────────────────────
def build(ref: str, out_dir: Path, cache: Path, subset: int | None,
          do_index: bool) -> int:
    genomes, contigs = load_manifest(ref)
    if subset:
        genomes = genomes[:subset]
        keep = {g["genome"] for g in genomes}
        contigs = [(n, L) for n, L in contigs
                   if _genome_of(n) in keep]
        print(f"  subset: first {len(genomes)} genome(s), "
              f"{len(contigs)} contig(s)")

    want = dict(contigs)
    order = [n for n, _ in contigs]
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / f"{ref}{'.subset' if subset else ''}.fna"

    # One staged file per genome, then concatenated in manifest order. The
    # staging step is what keeps memory flat: nothing larger than a single
    # genome is ever resident, whichever scheme it came from.
    stage = out_dir / f".stage_{ref}{'_subset' if subset else ''}"
    stage.mkdir(parents=True, exist_ok=True)
    bulk_schemes = {g["scheme"] for g in genomes} & {"rep", "wc"}
    if bulk_schemes:
        stage_bulk(ref, cache, stage,
                   {g["genome"] for g in genomes if g["scheme"] in bulk_schemes})

    per_acc = [g for g in genomes if g["scheme"] in ("gca", "mgyg")]
    if per_acc:
        print(f"\n  fetching {len(per_acc):,} per-accession genome(s)")
    for n, g in enumerate(per_acc, 1):
        dst = stage / f"{g['genome']}.fa"
        if dst.exists() and dst.stat().st_size:
            continue
        if n % 250 == 0 or n == len(per_acc):
            print(f"    {n:,}/{len(per_acc):,}")
        try:
            _write_records(dst, source_records(ref, g, cache))
        except Exception as exc:                       # noqa: BLE001
            print(f"    FAILED {g['genome']}: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            return 1

    print(f"\n  concatenating {len(genomes):,} genome(s) in manifest order")
    seen_names: list[str] = []
    bad: list[tuple[str, int, int]] = []
    with open(fasta, "w") as w:
        for g in genomes:
            src = stage / f"{g['genome']}.fa"
            if not src.exists():
                print(f"  staged file missing for {g['genome']}",
                      file=sys.stderr)
                return 1
            name, length, chunks = None, 0, []

            def emit():
                if name is None:
                    return
                if name not in want:
                    return                      # a contig the manifest drops
                seen_names.append(name)
                if length != want[name]:
                    bad.append((name, want[name], length))
                w.write(f">{name}\n")
                for c in chunks:
                    w.write(c + "\n")

            for line in open(src):
                if line.startswith(">"):
                    emit()
                    name = line[1:].split()[0]
                    length, chunks = 0, []
                else:
                    s = line.strip()
                    length += len(s)
                    chunks.append(s)
            emit()

    if bad:
        print(f"\n  {len(bad)} contig(s) differ in length from the manifest:",
              file=sys.stderr)
        for nm, wnt, got_ in bad[:5]:
            print(f"      {nm}: manifest {wnt:,} bp, built {got_:,} bp",
                  file=sys.stderr)
        return 1
    if len(seen_names) != len(order) or seen_names != order:
        missing = [n for n in order if n not in set(seen_names)]
        print(f"\n  built {len(seen_names):,} contigs against "
              f"{len(order):,} in the manifest", file=sys.stderr)
        if missing:
            print(f"  {len(missing)} not produced, e.g. {missing[:4]}",
                  file=sys.stderr)
        else:
            print("  same set, different order -- the manifest order was not "
                  "followed", file=sys.stderr)
        return 1

    import shutil as _sh
    _sh.rmtree(stage, ignore_errors=True)

    h = hashlib.sha256(fasta.read_bytes()).hexdigest()
    print(f"\n  wrote {fasta}  ({fasta.stat().st_size/1e9:.3f} GB)")
    print(f"  {len(order):,} contigs, "
          f"{sum(want.values()):,} bp, all names and lengths match the manifest")
    print(f"  sha256 {h}")
    (out_dir / f"{fasta.name}.sha256").write_text(f"{h}  {fasta.name}\n")

    # The second thing config needs, derived rather than shipped: every
    # catalogue here is one genome per SGB, so a contig's species is the
    # prefix of its own name. Emitting it beside the FASTA means the map
    # cannot drift from the reference it describes -- which shipping a 24 MB
    # copy alongside would eventually allow.
    cmap = out_dir / f"{fasta.with_suffix('').name}.contig_to_species.tsv"
    with open(cmap, "w") as w:
        for name in order:
            w.write(f"{_genome_of(name)}\t{name}\t{want[name]}\n")
    n_sp = len({_genome_of(n) for n in order})
    print(f"  wrote {cmap.name}  ({len(order):,} rows, {n_sp:,} species)")
    print("\n  point config.yaml at these two:")
    print(f"      reference_fasta:    {fasta}")
    print(f"      contig_to_species:  {cmap}")
    print("  the workflow builds the bowtie2 large index itself "
          "(rule build_metagenome_index)")

    if do_index:
        return bowtie2_index(fasta)
    return 0


def _genome_of(contig: str) -> str:
    for rx in (r"^(GCA_\d+\.\d+)_c_\d+$", r"^(MGYG\d+(?:\.\d+)?)_\d+$",
               r"^(Rep_\d+)_c_\d+$", r"^(WC\w+?)\|c_\d+$"):
        m = re.match(rx, contig)
        if m:
            return m.group(1)
    return contig


def bowtie2_index(fasta: Path) -> int:
    """bowtie2-build with --large-index.

    Every one of these catalogues is over the 4 GiB at which bowtie2 needs
    the 64-bit index, so --large-index is not optional -- without it the
    build fails late, after hours of work.
    """
    exe = "bowtie2-build"
    try:
        subprocess.run([exe, "--version"], capture_output=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        print(f"  {exe} not on PATH; skipping the index. Install bowtie2 "
              f"and re-run with --index.", file=sys.stderr)
        return 1
    base = str(fasta.with_suffix(""))
    cmd = [exe, "--large-index", "--threads", str(os.cpu_count() or 4),
           str(fasta), base]
    print(f"\n  {' '.join(cmd)}")
    r = subprocess.run(cmd)
    if r.returncode:
        print("  bowtie2-build failed", file=sys.stderr)
        return r.returncode
    made = sorted(p.name for p in fasta.parent.glob(
        fasta.with_suffix("").name + "*.bt2l"))
    print(f"  built {len(made)} large-index file(s): {made}")
    return 0


def verify_only(ref: str, path: Path) -> int:
    """Check an existing FASTA against the manifest, without rebuilding."""
    _genomes, contigs = load_manifest(ref)
    want = dict(contigs)
    seen, bad = 0, 0
    name, length = None, 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seen += 1
                    if want.get(name) != length:
                        bad += 1
                        if bad <= 5:
                            print(f"    {name}: manifest "
                                  f"{want.get(name, 'absent')}, file {length}")
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        seen += 1
        if want.get(name) != length:
            bad += 1
    partial = seen < len(want)
    print(f"\n  {path.name}: {seen:,} contig(s) checked")
    if bad:
        print(f"  {bad} disagree with the manifest")
    elif partial:
        # A subset build is the expected case here, not a failure: say which
        # it is rather than reporting a count difference as if it were one.
        print(f"  all {seen:,} match; this is a subset of the "
              f"{len(want):,}-contig catalogue")
    elif seen > len(want):
        print(f"  {seen - len(want):,} contig(s) are not in the manifest")
    else:
        print(f"  all {seen:,} match -- complete catalogue")
    return 0 if (not bad and seen <= len(want)) else 1


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", choices=sorted(SOURCES), default=None)
    ap.add_argument("--out", type=Path, default=Path("references"))
    ap.add_argument("--cache", type=Path, default=Path("reference_sources"),
                    help="where fetched sources are kept; re-used across runs")
    ap.add_argument("--subset", type=int, default=None,
                    help="build only the first N genomes, for proving the "
                         "method without a full rebuild")
    ap.add_argument("--index", action="store_true",
                    help="run bowtie2-build --large-index on the result")
    ap.add_argument("--verify-only", type=Path, default=None,
                    help="check this FASTA against the manifest and exit")
    ap.add_argument("--cite", action="store_true",
                    help="print the sources and licences and exit")
    a = ap.parse_args()

    if a.cite:
        for ref, rows in sorted(SOURCES.items()):
            print(f"\n  {ref}")
            for what, where, lic in rows:
                print(f"      {what}")
                print(f"          {where}   [{lic}]")
        return 0
    if not a.ref:
        print("  --ref is required unless you asked for --cite",
              file=sys.stderr)
        return 2
    if a.verify_only:
        return verify_only(a.ref, a.verify_only)
    return build(a.ref, a.out, a.cache, a.subset, a.index)


if __name__ == "__main__":
    raise SystemExit(main())
