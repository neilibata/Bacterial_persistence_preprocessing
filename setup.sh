#!/bin/bash
# setup.sh
# =========================================================================
# Interactive first-time setup. Detects your cluster environment and writes
# config/config.yaml + profiles/slurm/config.yaml.
#
# Safe to re-run — backs up existing config first.
#
# Usage:
#     bash setup.sh
# =========================================================================

set -euo pipefail
cd "$(dirname "$0")"

echo "============================================================"
echo " Preprocessing for the Bacterial Persistence Pipeline — Setup"
echo "============================================================"
echo ""

# Current value of a config key, or empty if the file or the key is absent.
# Re-running setup.sh must not silently replace a path you already set, so
# every prompt that can be answered from the existing config defaults to it.
cfg_val () {
    local key="$1"
    [[ -f config/config.yaml ]] || return 0
    { grep -E "^${key}:" config/config.yaml || true; } \
        | head -1 | sed -E "s/^${key}:[[:space:]]*//" | tr -d '"'
}

# --- Detect conda location ---
CONDA_SH=""
if command -v conda &>/dev/null; then
    CONDA_BASE=$(conda info --base 2>/dev/null || echo "")
    if [[ -n "$CONDA_BASE" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
        CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
    fi
fi
if [[ -z "$CONDA_SH" ]]; then
    for base in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3"; do
        if [[ -f "$base/etc/profile.d/conda.sh" ]]; then
            CONDA_SH="$base/etc/profile.d/conda.sh"
            break
        fi
    done
fi
if [[ -z "$CONDA_SH" ]]; then
    echo "WARNING: conda not found. You'll need conda to install the pipeline environment."
    echo "Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    echo ""
fi

# --- Check for snakemake ---
if ! command -v snakemake &>/dev/null; then
    echo "Snakemake not found. To install:"
    echo "  conda env create -f environment.yaml"
    echo "  conda activate bact_evo_pipeline"
    echo ""
    echo "Continuing setup anyway (you can install later)..."
    echo ""
fi

# --- Detect SLURM partitions ---
echo "--- Detecting SLURM partitions ---"
if command -v sinfo &>/dev/null; then
    echo "Available partitions on this cluster:"
    sinfo -s --format="%10P %6a %10l %5D %15F" 2>/dev/null || echo "  (could not query sinfo)"
    echo ""
else
    echo "WARNING: sinfo not found — are you on a SLURM cluster?"
    echo ""
fi

# --- Detect SLURM account ---
DETECTED_ACCOUNT=""
if command -v sacctmgr &>/dev/null; then
    DETECTED_ACCOUNT=$(sacctmgr -nP show associations user="$USER" format=account 2>/dev/null | head -1 || echo "")
fi

# --- Interactive prompts ---
echo "--- Configuration ---"
echo "Press Enter to accept defaults shown in [brackets]."
echo ""

read -rp "CPU-heavy partition name (trimming, alignment) [shared]: " PART_CPU
PART_CPU=${PART_CPU:-shared}

read -rp "Shared/light partition name (filtering, coverage) [shared]: " PART_SHARED
PART_SHARED=${PART_SHARED:-shared}

read -rp "Download partition (needs internet) [$PART_CPU]: " PART_DL
PART_DL=${PART_DL:-$PART_CPU}

# Samples above big_sample_mb go to a bigger-memory partition. The shipped
# default names the one this study used, which exists only on that cluster, so
# offer it as a default only if this cluster actually has it.
BIGMEM_DEF=$(cfg_val partition_bigmem)
if [[ -z "$BIGMEM_DEF" ]] || ! sinfo -h -p "$BIGMEM_DEF" -o '%P' 2>/dev/null | grep -q .; then
    BIGMEM_DEF=$PART_CPU
fi
read -rp "Big-memory partition (large alignments, repair) [$BIGMEM_DEF]: " PART_BIGMEM
PART_BIGMEM=${PART_BIGMEM:-$BIGMEM_DEF}

# --- Orchestrator partition and walltime ---------------------------------
# The orchestrator is 1 CPU and 4 GB that mostly sleeps, but it has to stay
# alive for the whole run, and most partitions cap walltime well below that:
# on the cluster this was built for, every general partition stops at 3 days
# while the orchestrator asked for 7. sbatch rejects that outright -- the run
# never starts -- so pick the partition and then clamp the request to whatever
# that partition actually permits.
part_maxtime () {  # minutes, or empty if sinfo cannot say
    command -v sinfo &>/dev/null || return 0
    sinfo -h -p "$1" -o '%l' 2>/dev/null | head -1 | awk '
        /infinite|UNLIMITED/ { print 525600; exit }
        /-/  { split($0,a,"-"); split(a[2],b,":"); print a[1]*1440 + b[1]*60 + b[2]; exit }
        /:/  { n=split($0,b,":"); if (n==3) print b[1]*60+b[2]; else print b[1]; exit }'
}
LONGEST=""
if command -v sinfo &>/dev/null; then
    LONGEST=$(for p in $(sinfo -h -o '%P' 2>/dev/null | tr -d '*' | sort -u); do
                  m=$(part_maxtime "$p"); [[ -n "$m" ]] && echo "$m $p"
              done | sort -rn | head -1 | awk '{print $2}')
fi
read -rp "Orchestrator partition (long walltime, 1 CPU) [${LONGEST:-$PART_SHARED}]: " PART_ORCH
PART_ORCH=${PART_ORCH:-${LONGEST:-$PART_SHARED}}

ORCH_WANT_MIN=10080                        # 7 days, what the pipeline prefers
ORCH_MAX=$(part_maxtime "$PART_ORCH")
if [[ -n "$ORCH_MAX" && "$ORCH_MAX" -gt 0 && "$ORCH_MAX" -lt "$ORCH_WANT_MIN" ]]; then
    ORCH_MIN=$ORCH_MAX
    echo "  note: $PART_ORCH allows at most $((ORCH_MAX/1440))d $(((ORCH_MAX%1440)/60))h;"
    echo "        the orchestrator will request that instead of 7 days."
    echo "        Snakemake resumes where it stopped, so re-submit if it runs out."
else
    ORCH_MIN=$ORCH_WANT_MIN
fi
ORCH_TIME=$(printf '%d-%02d:%02d:00' $((ORCH_MIN/1440)) $(((ORCH_MIN%1440)/60)) $((ORCH_MIN%60)))

DEFAULT_ACCT="${DETECTED_ACCOUNT:-}"
if [[ -n "$DEFAULT_ACCT" ]]; then
    read -rp "SLURM account (detected: $DEFAULT_ACCT) [$DEFAULT_ACCT]: " SLURM_ACCT
    SLURM_ACCT=${SLURM_ACCT:-$DEFAULT_ACCT}
else
    read -rp "SLURM account (leave empty if not required) []: " SLURM_ACCT
    SLURM_ACCT=${SLURM_ACCT:-}
fi

read -rp "Max concurrent SLURM jobs [200]: " MAX_JOBS
MAX_JOBS=${MAX_JOBS:-200}

read -rp "Results directory (use fast scratch) [results]: " RESULTS_DIR
RESULTS_DIR=${RESULTS_DIR:-results}

read -rp "Alignment threads per sample [20]: " ALIGN_T
ALIGN_T=${ALIGN_T:-20}

read -rp "MAPQ threshold [30]: " MAPQ
MAPQ=${MAPQ:-30}

# --- Conda env name ---
read -rp "Conda environment name [bact_evo_pipeline]: " CONDA_ENV
CONDA_ENV=${CONDA_ENV:-bact_evo_pipeline}

echo ""
echo "--- Input files ---"
echo "You provide one file:"
echo "  samples.tsv  — SRA/ENA accessions, one per line"
echo ""
echo "The reference catalogue and its contig map are NOT shipped with this"
echo "repository: three of the four published catalogues are other groups'"
echo "data, and all four are tens of gigabytes. Rebuild the one you want from"
echo "its public source instead -- the result is byte-identical to the"
echo "reference the reported analysis aligned against, and the build verifies"
echo "that itself:"
echo ""
echo "  python workflow/scripts/build_references.py --ref combined --out references/"
echo ""
echo "That writes the FASTA and the contig map, and prints the two paths asked"
echo "for below. --ref takes carter, leviatan, combined or chimp; --cite prints"
echo "every source and its licence. Use --subset 12 first if you want to see"
echo "the whole path work before committing to the full download."
echo ""

REF_FASTA_DEF=$(cfg_val reference_fasta)
REF_FASTA_DEF=${REF_FASTA_DEF:-references/combined.fna}
read -rp "reference_fasta [$REF_FASTA_DEF]: " REF_FASTA
REF_FASTA=${REF_FASTA:-$REF_FASTA_DEF}

CONTIG_MAP_DEF=$(cfg_val contig_to_species)
if [[ -z "$CONTIG_MAP_DEF" ]]; then
    # Derived from the answer above, since build_references.py writes the two
    # beside each other and pairing a map with a different catalogue is silent.
    CONTIG_MAP_DEF="${REF_FASTA%.*}.contig_to_species.tsv"
fi
read -rp "contig_to_species [$CONTIG_MAP_DEF]: " CONTIG_MAP
CONTIG_MAP=${CONTIG_MAP:-$CONTIG_MAP_DEF}

for f in "$REF_FASTA" "$CONTIG_MAP"; do
    [[ -e "$f" ]] || echo "  note: $f does not exist yet — build it before launching."
done
echo ""

# --- Write config ---
CONFIG_FILE="config/config.yaml"
if [[ -f "$CONFIG_FILE" ]]; then
    BACKUP="${CONFIG_FILE}.bak.$(date +%Y%m%d_%H%M%S)"
    cp "$CONFIG_FILE" "$BACKUP"
    echo "Backed up existing config to $BACKUP"
fi

cat > "$CONFIG_FILE" << YAML
# =========================================================================
# config.yaml — Generated by setup.sh on $(date)
# =========================================================================

# --- Input files ---
samples_tsv: "resources/samples.tsv"

# Both of these come from one build_references.py run; neither is shipped.
reference_fasta: "$REF_FASTA"
contig_to_species: "$CONTIG_MAP"

# --- Output ---
results_dir: "$RESULTS_DIR"

# --- SLURM ---
partition_cpu: "$PART_CPU"
partition_shared: "$PART_SHARED"
partition_download: "$PART_DL"
partition_bigmem: "$PART_BIGMEM"
big_sample_mb: 15000
slurm_account: "$SLURM_ACCT"
max_jobs: $MAX_JOBS

# --- Parameters ---
mapq_threshold: $MAPQ
prefetch_max_size: "300G"
align_threads: $ALIGN_T
preprocess_threads: 8

# --- Intermediate file cleanup ---
# The pipeline automatically deletes intermediate files after downstream
# rules succeed, keeping only the important checkpoints:
#
# CLEANUP CHAIN:
#   raw FASTQs         → deleted after trimming
#   trimmed reads      → deleted after PhiX removal
#   PhiX-cleaned       → deleted after human removal
#   human-cleaned      → deleted after repair
#   ★ repaired reads   → KEPT by default (keep_repaired_reads: true)
#   aligned BAMs       → deleted after MAPQ filtering
#   ★ filtered BAMs    → KEPT by default (keep_filtered_bams: true)
#   ★ coverage files   → ALWAYS KEPT
#   ★ aggregated CSVs  → ALWAYS KEPT
#   ★ final matrices   → ALWAYS KEPT
#
# Set to true to keep ALL intermediate files (uses ~5-10x more disk):
keep_all: false

# Repaired reads are useful as a re-alignment checkpoint and for other projects.
# Set to false to delete them after alignment (saves ~4GB per sample):
keep_repaired_reads: true

# Filtered BAMs (high/low MAPQ) are kept by default for downstream analysis.
# Set to false to delete them after coverage extraction (saves disk):
keep_filtered_bams: true
YAML

echo ""
echo "Config written to $CONFIG_FILE"

# --- Write SLURM profile (Snakemake 8/9 format) ---
PROFILE="profiles/slurm/config.yaml"
ACCT_LINE=""
if [[ -n "$SLURM_ACCT" ]]; then
    ACCT_LINE="  - slurm_account=$SLURM_ACCT"
fi

# Build precommand for conda activation in child SLURM jobs
PRECOMMAND_LINE=""
if [[ -n "$CONDA_SH" ]]; then
    PRECOMMAND_LINE="precommand: \"source $CONDA_SH && conda activate $CONDA_ENV\""
    echo ""
    echo "Detected conda at: $CONDA_SH"
    echo "Child SLURM jobs will activate: $CONDA_ENV"
else
    echo ""
    echo "WARNING: Could not auto-detect conda path."
    echo "You'll need to manually add a precommand to $PROFILE"
    echo "Example: precommand: \"source ~/miniconda3/etc/profile.d/conda.sh && conda activate $CONDA_ENV\""
fi

cat > "$PROFILE" << YAML
# Snakemake 8/9 SLURM profile
# Generated by setup.sh on $(date)
# Requires: pip install snakemake-executor-plugin-slurm

executor: slurm
jobs: $MAX_JOBS
latency-wait: 120
keep-going: true
rerun-incomplete: true
printshellcmds: true

# Activate conda in every child SLURM job so tools are in PATH
${PRECOMMAND_LINE:-# precommand: "source ~/miniconda3/etc/profile.d/conda.sh && conda activate $CONDA_ENV"}

default-resources:
  - slurm_partition=$PART_SHARED
  - mem_mb=4000
  - runtime=60
${ACCT_LINE}
YAML

echo "Profile written to $PROFILE"

# --- Update run_pipeline.sh SBATCH headers ---
if [[ -n "$SLURM_ACCT" ]]; then
    sed -i "s/^#SBATCH --account=.*/#SBATCH --account=$SLURM_ACCT/" run_pipeline.sh
fi
sed -i "s/^#SBATCH --partition=.*/#SBATCH --partition=$PART_ORCH/" run_pipeline.sh
sed -i "s/^#SBATCH --time=.*/#SBATCH --time=$ORCH_TIME/" run_pipeline.sh
echo "Orchestrator: $PART_ORCH, --time=$ORCH_TIME"

# --- Create log dirs ---
mkdir -p logs/slurm

# --- Summary ---
echo ""
echo "============================================================"
echo " Setup complete!"
echo "============================================================"
echo ""
echo " Next steps:"
echo "   1. Put your accessions in resources/samples.tsv, and build the"
echo "      catalogue if you have not already:"
echo "        python workflow/scripts/build_references.py --ref combined --out references/"
echo "   2. Install the conda environment (if not done):"
echo "        conda env create -f environment.yaml"
echo "        conda activate $CONDA_ENV"
echo "   3. Preview the pipeline:"
echo "        bash run_pipeline.sh --dry-run"
echo "   4. Launch:"
echo "        sbatch run_pipeline.sh"
echo ""
