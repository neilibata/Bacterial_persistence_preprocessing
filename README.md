# Bacterial Evolution & Health Pipeline

A Snakemake pipeline that processes shotgun metagenomic data from raw SRA accessions to per-species coverage matrices on any SLURM cluster. Built for investigating bacterial evolutionary dynamics across human populations.

## What You Provide

| File | Description |
|---|---|
| `resources/samples.tsv` | SRA accessions (one per line): SRR, ERR, or DRR IDs |
| `resources/reference_metagenome.fna` | Your metagenome reference FASTA (e.g. UHGG subset) |
| `resources/contig_to_species.tsv` | Mapping: `species_id <TAB> contig_name <TAB> length` |

Everything else — downloading PhiX/hg38 genomes, building all Bowtie2 indices, QC, alignment, filtering, coverage extraction, aggregation, and matrix combination — runs automatically.

## Quick Start

```bash
# 1. Clone and enter the repo
git clone https://github.com/YOUR_USER/Bacterial_evolution_health.git
cd Bacterial_evolution_health

# 2. Install the environment (all tools in one conda env)
conda env create -f environment.yaml
conda activate bact_evo_pipeline

# 3. Place your three input files in resources/

# 4. Run interactive setup (detects partitions, conda path, writes config)
bash setup.sh

# 5. Preview what would run
bash run_pipeline.sh --dry-run

# 6. Launch on SLURM
sbatch run_pipeline.sh

# 7. Check progress any time
bash run_pipeline.sh --status
```

## Pipeline Overview

```
                    ┌──────────────────────────────────┐
                    │     ONE-TIME SETUP (automatic)    │
                    │                                    │
                    │  Download PhiX → build bt2 index   │
                    │  Download hg38 → build bt2 index   │
                    │  Your .fna     → build bt2 index   │
                    └──────────────┬───────────────────┘
                                   │
           ┌───────────────────────┼───────────────────────┐
           │                       │                       │
    ┌──────▼──────┐        ┌──────▼──────┐        ┌──────▼──────┐
    │  Sample A   │        │  Sample B   │        │  Sample C   │ ...
    │             │        │             │        │             │
    │  Download   │        │  Download   │        │  Download   │
    │  Trim       │        │  Trim       │        │  Trim       │
    │  PhiX rm    │        │  PhiX rm    │        │  PhiX rm    │
    │  Human rm   │        │  Human rm   │        │  Human rm   │
    │  Repair     │        │  Repair     │        │  Repair     │
    │  Align      │        │  Align      │        │  Align      │
    │  MAPQ split │        │  MAPQ split │        │  MAPQ split │
    │  Coverage   │        │  Coverage   │        │  Coverage   │
    │  Aggregate  │        │  Aggregate  │        │  Aggregate  │
    │     ↓       │        │     ↓       │        │     ↓       │
    │   CSV out   │        │   CSV out   │        │   CSV out   │
    └──────┬──────┘        └──────┬──────┘        └──────┬──────┘
           │                      │                      │
           └──────────────────────┼──────────────────────┘
                                  ▼
                    ┌──────────────────────────┐
                    │   Combine into matrices   │
                    │  depth_matrix.csv         │
                    │  coverage_matrix.csv      │
                    └──────────────────────────┘

         ALL SAMPLES RUN IN PARALLEL on separate SLURM jobs
```

Each box above runs as independent SLURM jobs. With `max_jobs=200`, the pipeline can process ~25-50 samples simultaneously (each sample uses 4-8 job slots across its active rules). Snakemake handles all scheduling — sample A might be aligning while sample B is downloading and sample C is extracting coverage, all concurrently.

## How It Runs on SLURM

```
 ┌───────────────────────────────────────────────────────┐
 │  sbatch run_pipeline.sh                               │
 │  ┌─────────────────────────────────────────────────┐  │
 │  │  Snakemake orchestrator (1 CPU, 4 GB, 7 days)   │  │
 │  │                                                   │  │
 │  │  Submits + monitors child SLURM jobs:            │  │
 │  │    → sbatch trim_galore  (8 CPU, ~16 GB)         │  │
 │  │    → sbatch align        (20 CPU, ~100 GB)       │  │
 │  │    → sbatch filter_mapq  (4 CPU, ~8 GB)          │  │
 │  │    → sbatch coverage     (1 CPU, ~4 GB)          │  │
 │  │    → ...                                          │  │
 │  └─────────────────────────────────────────────────┘  │
 └───────────────────────────────────────────────────────┘
```

The orchestrator itself is a tiny job (1 CPU, 4 GB) that runs for up to 7 days. It uses `sbatch` to submit compute-heavy work as separate SLURM jobs, then polls their status via `sacct`. This means the orchestrator can run on any low-priority partition.

## How Failures Are Handled

Three mechanisms work together:

### 1. Resources scale with input size and double on retry

```
Sample with 500 MB FASTQs:
  Attempt 1: 8 GB RAM,  6 min wall time
  Attempt 2: 16 GB RAM, 12 min
  Attempt 3: 24 GB RAM, 18 min

Sample with 50 GB FASTQs:
  Attempt 1: 100 GB RAM, 10 hr wall time
  Attempt 2: 200 GB RAM, 20 hr
  Attempt 3: 300 GB RAM, 30 hr
```

### 2. Re-run picks up where it left off

```bash
# Some samples failed? Just re-submit:
sbatch run_pipeline.sh
# Completed samples are skipped. Only failures are re-tried.
```

### 3. Safe intermediate cleanup

Intermediate files (raw FASTQs, trimmed reads, aligned BAMs) are automatically deleted **only after** the downstream step succeeds. If a job fails mid-pipeline, the last good intermediate is preserved. See [Disk Space Management](#disk-space-management) below.

## Disk Space Management

The pipeline automatically cleans up intermediate files to save disk space (~5-10x reduction). Two config options control this behavior:

```yaml
# config/config.yaml

# Set to true to keep ALL intermediate files (uses ~5-10x more disk):
keep_all: false

# Keep repaired reads as a re-alignment checkpoint (and for other projects).
# Set to false for large-scale runs to save ~4GB per sample:
keep_repaired_reads: true

# Keep filtered BAMs for downstream analysis.
# Set to false to delete them after coverage extraction:
keep_filtered_bams: true
```

**What gets cleaned up vs. kept:**

| Stage | Files | Default behavior |
|---|---|---|
| Raw FASTQs | `results/fastq/` | Deleted after trimming |
| Trimmed reads | `results/trimmed/` | Deleted after PhiX removal |
| PhiX-cleaned | `results/decontam/*_phiX_*` | Deleted after human removal |
| Human-cleaned | `results/decontam/*_human_*` | Deleted after repair |
| **Repaired reads** | `results/repaired/` | **Kept** (configurable) |
| Aligned BAMs | `results/aligned/` | Deleted after MAPQ filtering |
| **Filtered BAMs** | `results/filtered/` | **Kept** (configurable) |
| **Coverage files** | `results/coverage/` | **Always kept** |
| **Aggregated CSVs** | `results/aggregated/` | **Always kept** |
| **Final matrices** | `results/matrices/` | **Always kept** |

**For large-scale runs (thousands of samples)**, set both `keep_repaired_reads: false` and `keep_filtered_bams: false`, and consider reducing `align_threads` to 8 while increasing `max_jobs` to 990 for better throughput.

Cleanup is safe: Snakemake only deletes an intermediate file **after** the rule consuming it has completed successfully.

## Parallelism Design

The original bash pipeline used `xargs -P 4` to run 4 alignments within a single 80-core SLURM job. This Snakemake version achieves the same throughput differently: it submits 4 separate 20-core alignment jobs. This is better because:

- If one sample fails, only that job fails (not 3 other samples on the same node)
- Resources scale per-sample (a small sample doesn't waste 20 cores for an hour)
- Snakemake can mix stages: while sample A aligns, sample B can trim, sample C can extract coverage

**Concurrency math**: With `max_jobs=200` and `align_threads=20`, the cluster can run:
- ~10 alignments concurrently (the bottleneck stage)
- Plus ~50 lighter jobs (trim, filter, coverage) on shared-partition nodes
- All at the same time, on different samples

## Repository Structure

```
Bacterial_evolution_health/
├── config/
│   └── config.yaml               ← EDIT THIS (or run setup.sh)
├── resources/                     ← YOUR INPUT DATA
│   ├── samples.tsv
│   ├── reference_metagenome.fna
│   └── contig_to_species.tsv
├── workflow/
│   ├── Snakefile                  pipeline rules
│   └── scripts/
│       ├── aggregate_coverage.py  contig → species aggregation
│       └── combine_aggregated_stats.py  per-sample → matrices
├── profiles/
│   └── slurm/
│       ├── config.yaml            SLURM submission settings
│       └── status.py              job status checker (sacct)
├── setup.sh                       interactive first-time setup
├── run_pipeline.sh                orchestrator (sbatch this)
├── environment.yaml               conda env (all tools)
├── README.md
└── .gitignore

results/                           (created by pipeline)
├── references/                    PhiX + hg38 + metagenome indices
├── fastq/                         (temp — auto-deleted)
├── trimmed/                       (temp — auto-deleted)
├── decontam/                      (temp — auto-deleted)
├── repaired/                      ★ clean reads (KEPT by default)
├── aligned/                       (temp — auto-deleted)
├── filtered/{high,low}/           ★ MAPQ-split BAMs (KEPT by default)
├── coverage/{high,low}/           ★ per-contig coverage (KEPT)
├── aggregated/                    ★ per-species CSVs (KEPT)
└── matrices/                      ★ final sample × species matrices ★
```

## Common Operations

| What | Command |
|---|---|
| First-time setup | `bash setup.sh` |
| Submit pipeline | `sbatch run_pipeline.sh` |
| Preview (no submission) | `bash run_pipeline.sh --dry-run` |
| Check progress | `bash run_pipeline.sh --status` |
| After failures | `sbatch run_pipeline.sh` (just re-submit) |
| Unlock after crash | `bash run_pipeline.sh --unlock` |
| Run locally (testing) | `bash run_pipeline.sh --local 8` |
| Force re-run a rule | `snakemake --profile profiles/slurm -j 200 --forcerun align` |
| Add new samples | Add to `resources/samples.tsv`, then `sbatch run_pipeline.sh` |

## Skipping Steps (Using Existing Data)

### Already have FASTQs?

```bash
mkdir -p results/fastq
ln -s /path/to/SRR12345678_1.fastq.gz results/fastq/
ln -s /path/to/SRR12345678_2.fastq.gz results/fastq/
```

### Already have Bowtie2 indices?

```bash
mkdir -p results/references/{phix,human,metagenome}

# Symlink existing index files
ln -s /path/to/phix.*.bt2 results/references/phix/
ln -s /path/to/hg38.*.bt2l results/references/human/
ln -s /path/to/reference.*.bt2l results/references/metagenome/

# Create sentinel files so Snakemake skips the build rules
touch results/references/{phix,human,metagenome}/.index_done
```

## Cluster Portability

The pipeline uses **no hardcoded cluster paths**. Everything is configured through `config/config.yaml`:

| Setting | What it controls |
|---|---|
| `partition_cpu` | Partition for CPU-heavy jobs (trim, align, index build) |
| `partition_shared` | Partition for light jobs (filter, coverage, aggregate) |
| `partition_download` | Partition for jobs needing internet (SRA + ref download) |
| `slurm_account` | SLURM account string (leave `""` if not required) |
| `max_jobs` | Max concurrent SLURM jobs |
| `align_threads` | Cores per alignment job |
| `preprocess_threads` | Cores per preprocessing job |

Run `bash setup.sh` to detect your cluster's partitions and configure automatically. The setup script also auto-detects your conda installation path and writes a `precommand` into the SLURM profile so that child jobs activate the correct environment.

### Clusters without internet on compute nodes

If your compute nodes can't reach the internet, pre-download the data:

```bash
# On a login node (has internet):
# Download PhiX
curl -fSL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" | gunzip > results/references/phix/phix.fna

# Download hg38
curl -fSL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" | gunzip > results/references/human/hg38.fna

# Download SRA data
conda activate bact_evo_pipeline
while read acc; do
    prefetch --max-size 300G "$acc"
    fasterq-dump "$acc" --split-files --outdir results/fastq/
    pigz results/fastq/${acc}*.fastq
done < resources/samples.tsv

# Then submit the pipeline — it will skip all downloads
sbatch run_pipeline.sh
```

## Output

### Per-sample CSVs (`results/aggregated/`)

| Column | Description |
|---|---|
| `Species_ID` | Species/SGB identifier |
| `Weighted_Mean_Depth` | Mean read depth, weighted by contig length |
| `Weighted_Std_Dev_Depth` | Weighted standard deviation of depth |
| `Weighted_Mean_Coverage` | Breadth of coverage (0–1), weighted by contig length |
| `Weighted_Std_Dev_Coverage` | Weighted standard deviation of coverage |

### Combined matrices (`results/matrices/`)

Two sample × species matrices are produced as the final deliverables:

- **`depth_matrix.csv`** — rows = samples, columns = `Sample_ID` + `{species}_mean` + `{species}_std` for depth
- **`coverage_matrix.csv`** — same layout for breadth-of-coverage

Species column order is fixed by `contig_to_species.tsv`, so all matrices are directly comparable.

## Software

All installed via `environment.yaml`:

| Tool | Stage |
|---|---|
| Snakemake ≥ 8.0 | Orchestration |
| snakemake-executor-plugin-slurm | SLURM integration |
| SRA Toolkit + pigz | Download |
| Trim Galore | Preprocessing |
| Bowtie2 | Preprocessing + Alignment |
| BBTools (repair.sh) | Preprocessing |
| samtools | Alignment + Filtering + Coverage |
| Python + pandas | Aggregation + Matrix combination |
