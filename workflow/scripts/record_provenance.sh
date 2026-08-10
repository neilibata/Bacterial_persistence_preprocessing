#!/usr/bin/env bash
# Record exactly what produced this run's outputs.
#
# The first version of this pipeline did not capture tool versions, so the
# matrices it produced cannot be tied to a specific bowtie2 or samtools build.
# That is unrecoverable for those outputs; this makes it impossible to repeat.
set -uo pipefail
OUT="${1:?usage: record_provenance.sh <output.txt>}"
mkdir -p "$(dirname "$OUT")"

{
  echo "# Provenance record"
  echo "generated_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host: $(hostname)"
  echo "user: ${USER:-unknown}"
  echo

  echo "## Pipeline revision"
  if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "git_commit: $(git rev-parse HEAD)"
    echo "git_describe: $(git describe --tags --always --dirty 2>/dev/null || echo none)"
    echo "git_dirty: $(test -n "$(git status --porcelain)" && echo yes || echo no)"
  else
    echo "git_commit: (not a git checkout)"
  fi
  echo

  echo "## Tool versions"
  for tool in snakemake bowtie2 samtools trim_galore cutadapt prefetch fasterq-dump pigz python; do
    if command -v "$tool" >/dev/null 2>&1; then
      v=$("$tool" --version 2>&1 | head -1)
      printf '%-14s %s\n' "$tool:" "$v"
    else
      printf '%-14s %s\n' "$tool:" "NOT FOUND"
    fi
  done
  # BBTools reports its version differently
  if command -v repair.sh >/dev/null 2>&1; then
    printf '%-14s %s\n' "bbmap:" "$(repair.sh --version 2>&1 | grep -i version | head -1)"
  else
    printf '%-14s %s\n' "bbmap:" "NOT FOUND"
  fi
  echo

  echo "## Conda environment"
  if command -v conda >/dev/null 2>&1; then
    echo "env: ${CONDA_DEFAULT_ENV:-none}"
    conda list --export 2>/dev/null || echo "(conda list failed)"
  else
    echo "(conda not on PATH)"
  fi
} > "$OUT"

echo "Wrote provenance to $OUT"
