#!/bin/bash
#SBATCH --job-name=smk_orchestrator
#SBATCH --output=logs/orchestrator_%j.out
#SBATCH --error=logs/orchestrator_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared
#SBATCH --account=carmody_lab

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    cd "$SLURM_SUBMIT_DIR"
else
    cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

PROFILE="profiles/slurm"
CONDA_ENV="bact_evo_pipeline"
MAX_JOBS=$(grep -oP 'jobs:\s*\K\d+' "$PROFILE/config.yaml" 2>/dev/null || echo 200)

activate_conda() {
    if command -v conda &>/dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate "$CONDA_ENV" 2>/dev/null && return 0
    fi
    for base in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3"; do
        if [[ -f "$base/etc/profile.d/conda.sh" ]]; then
            source "$base/etc/profile.d/conda.sh"
            conda activate "$CONDA_ENV" 2>/dev/null && return 0
        fi
    done
    if command -v module &>/dev/null; then
        module load anaconda 2>/dev/null || module load conda 2>/dev/null || true
        if command -v conda &>/dev/null; then
            eval "$(conda shell.bash hook)"
            conda activate "$CONDA_ENV" 2>/dev/null && return 0
        fi
    fi
    echo "WARNING: Could not activate conda env '$CONDA_ENV'."
}
activate_conda

command -v snakemake &>/dev/null || { echo "ERROR: snakemake not found"; exit 1; }
echo "Snakemake version: $(snakemake --version)"
echo "Working directory: $(pwd)"

mkdir -p logs

ACTION="${1:-run}"

case "$ACTION" in

  run)
    echo "============================================================"
    echo " Pipeline Orchestrator"
    echo " Started:  $(date)"
    echo " Host:     $(hostname)"
    echo " Job ID:   ${SLURM_JOB_ID:-interactive}"
    echo " Max jobs: $MAX_JOBS"
    echo "============================================================"

    LOGFILE="logs/snakemake_$(date +%Y%m%d_%H%M%S).log"

    snakemake \
      --profile "$PROFILE" \
      -j "$MAX_JOBS" \
      --printshellcmds \
      --nocolor \
      2>&1 | tee "$LOGFILE"

    EXIT_CODE=${PIPESTATUS[0]}
    echo "Finished: $(date)  Exit: $EXIT_CODE"
    exit $EXIT_CODE
    ;;

  -n|--dry-run)
    echo "=== DRY RUN ==="
    snakemake --profile "$PROFILE" -j "$MAX_JOBS" -n --printshellcmds
    ;;

  --local)
    CORES="${2:-8}"
    echo "=== LOCAL RUN ($CORES cores) ==="
    snakemake -j "$CORES" --printshellcmds --executor local
    ;;

  --status)
    echo "Pipeline Status — $(date)"
    count() { find "$1" -type f -name "$2" 2>/dev/null | wc -l; }
    printf "  %-35s %s\n" "Downloaded FASTQs:"      "$(count results/fastq '*_1.fastq.gz')"
    printf "  %-35s %s\n" "Repaired reads:"          "$(count results/repaired '*_1_repaired.fq.gz')"
    printf "  %-35s %s\n" "Filtered BAMs (high):"    "$(count results/filtered/high '*_high_mapq.bam')"
    printf "  %-35s %s\n" "Coverage extracted:"       "$(count results/coverage/high '*_coverage.txt')"
    printf "  %-35s %s\n" "Aggregated (per-sample):"  "$(count results/aggregated '*_stats.csv')"
    for mat in results/matrices/depth_matrix.csv results/matrices/coverage_matrix.csv; do
        name=$(basename "$mat")
        [[ -f "$mat" ]] && printf "  %-35s  ✓ done\n" "$name" || printf "  %-35s  ✗ pending\n" "$name"
    done
    for ref in phix human metagenome; do
        s="results/references/${ref}/.index_done"
        [[ -f "$s" ]] && printf "  %-20s  ✓ built\n" "$ref" || printf "  %-20s  ✗ not yet\n" "$ref"
    done
    ;;

  --unlock)
    snakemake --unlock
    echo "Unlocked. Re-run with: sbatch run_pipeline.sh"
    ;;

  -h|--help)
    echo "Usage:"
    echo "  sbatch run_pipeline.sh          Submit as SLURM orchestrator"
    echo "  bash  run_pipeline.sh --dry-run Preview (no submission)"
    echo "  bash  run_pipeline.sh --local 8 Run locally with 8 cores"
    echo "  bash  run_pipeline.sh --status  Check progress"
    echo "  bash  run_pipeline.sh --unlock  Unlock after crash"
    ;;

  *) echo "Unknown: $ACTION"; exit 1 ;;
esac
