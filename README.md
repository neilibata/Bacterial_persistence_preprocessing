# Bacterial Evolution & Health Pipeline

> ### Versions
>
> | Tag | What it is |
> | --- | --- |
> | **`v1.0-as-published`** | The pipeline exactly as it produced the matrices used in the manuscript. Preserved unchanged for provenance. It has known defects (below) — use it to reproduce the published matrices, not for new work. |
> | **`v2.0`** and later | Same analysis, with the defects fixed. Use this for new processing and backfills. |
>
> **Defects in v1.0, all silent:**
> 1. `combine_aggregated_stats.py` globbed the aggregated directory instead of
>    using the file list Snakemake declared, so the matrix contained whatever was
>    on disk — omitting samples that failed upstream, and including samples left
>    over from earlier runs.
> 2. An unreadable per-sample CSV was warned about and skipped; the matrix was
>    written anyway.
> 3. No check that every requested accession reached the matrix. 398 samples were
>    absent from one reference catalogue and this was not discovered until
>    analysis, months later.
> 4. `repair.sh` (a JVM tool) ran without `-Xmx`, so it sized its heap from node
>    memory and was killed by the SLURM cgroup at the requested limit.
> 5. Alignment memory was modelled as `input_size x 5`, but bowtie2's footprint is
>    dominated by the resident index and is essentially flat in sample size — so
>    large samples were under-provisioned and OOM-killed.
> 6. `bowtie2 -p N` plus `samtools --threads 4` twice requested more threads than
>    the allocation had.
> 7. No record of tool versions, so v1.0 outputs cannot be tied to specific builds.
> 8. Every accession was treated as its own sample, so a sample sequenced across
>    several runs produced one matrix row *per run* instead of one per sample.
>    16% of this cohort is multi-run. The damage is not repairable after the
>    fact — breadth is a percentage of genome covered, so per-run breadths
>    cannot be pooled without the original alignments; 89 samples had to be
>    reprocessed from raw reads. *(Fixed in v3.0.)*
>
> v2.0 fixes 1–7; v3.0 fixes 8. Outputs are written as **dated shards**
> (`depth_matrix_YYYYMMDD.csv`) and merged into the canonical matrix, so a
> backfill never rewrites samples that are already done.


A Snakemake pipeline that processes shotgun metagenomic data from raw SRA accessions to per-species coverage matrices on any SLURM cluster. Built for investigating bacterial evolutionary dynamics across human populations.

## What You Provide

| File | Description |
|---|---|
| `resources/samples.tsv` | Sample list — see **Sample list formats** below |

That is the only file you supply. The reference catalogue and its contig map
are **rebuilt from their public sources** by
`workflow/scripts/build_references.py` — see
**[Reference catalogues](#reference-catalogues)** — and are deliberately not
committed: each is 3.6–11.8 GB, and three of the four are other groups'
published data.

Everything else — downloading PhiX/hg38 genomes, building all Bowtie2 indices, QC, alignment, filtering, coverage extraction, aggregation, and matrix combination — runs automatically.

### Sample list formats

**One accession per line** — each accession is treated as its own sample:

```
ERR7745335
SRR25477699
```

**`sample_id <TAB> run1;run2;…`** — for samples sequenced across several runs:

```
sample_id	runs
FAT_006-22-0-0	ERR1297836;ERR1297797;ERR1297781;ERR1297773
FAT_008-22-0-0	ERR1297645;ERR1297705
```

Both formats can be mixed in one file, and a `sample_id` header is optional.

Use the second form whenever one biological sample has multiple runs. The runs
are downloaded, checked for paired-end output, and **concatenated in listed
order before trimming**, so the sample produces exactly one matrix row.

This matters more than it looks. If multi-run samples are listed one-accession-
per-line instead, the pipeline emits one matrix row *per run*, and those rows
cannot be merged afterwards — breadth is a percentage of genome covered, so
per-run breadths cannot be pooled without the original alignments. In this
project 16% of samples are multi-run (up to 35 runs each), and getting this
wrong for one reference catalogue cost 89 samples that had to be reprocessed
from scratch.

Any run that does not yield both `_1` and `_2` FASTQs aborts the sample rather
than being silently dropped, since a missing mate would desynchronise the
read pairing of every run concatenated after it.

## Quick Start

```bash
# 1. Clone and enter the repo
git clone https://github.com/YOUR_USER/Bacterial_evolution_health.git
cd Bacterial_evolution_health

# 2. Install the environment (all tools in one conda env)
conda env create -f environment.yaml
conda activate bact_evo_pipeline

# 3. Put your sample list at resources/samples.tsv

# 4. Rebuild the reference catalogue you want (see Reference catalogues).
#    --subset 12 first if you want to watch the whole path work cheaply.
python workflow/scripts/build_references.py --ref combined --out references/

# 5. Run interactive setup (detects partitions, conda path, writes config)
bash setup.sh

# 6. Preview what would run
bash run_pipeline.sh --dry-run

# 7. Launch on SLURM
sbatch run_pipeline.sh

# 8. Check progress any time
bash run_pipeline.sh --status
```

## Reference catalogues

The pipeline aligns against one of four catalogues. None is committed here:
each is 3.6–11.8 GB, and three of the four are other groups' published data,
so this repository ships the **recipe** rather than copies. Rebuilding gives a
file that is byte-identical to the one the reported analysis used, and the
build proves that rather than assuming it.

```bash
python workflow/scripts/build_references.py --ref combined --out references/
```

| `--ref` | Genomes | Size | Source | Licence |
|---|---|---|---|---|
| `carter` | 3,971 | 9.1 GB | NCBI GenBank (1,816 `GCA_*`) + MGnify `human-gut-v2-0-2` (2,155 `MGYG*`) | open |
| `leviatan` | 3,594 | 8.9 GB | figshare [`10.6084/m9.figshare.16885261`](https://doi.org/10.6084/m9.figshare.16885261) → `Rep_all.fa.gz` | CC BY 4.0 |
| `combined` | 5,081 | 11.8 GB | all of `leviatan`, plus the 1,487 `carter` genomes dereplication kept | — |
| `chimp` | 1,682 | 3.7 GB | Dryad [`10.5061/dryad.00000006x`](https://doi.org/10.5061/dryad.00000006x) → `ChimpMAGs.tar` | CC0 1.0 |

`combined` is the reference the manuscript reports; `chimp` is the negative
control. `--cite` prints every source and licence without downloading
anything.

**Cite the catalogues, not just this repository:** Carter et al. (*Cell* 2023;
*Nature* 2026), Leviatan et al. (*Nat Commun* 13, 3863, 2022 — whose own
construction code is at [erans99/GutReferenceSet](https://github.com/erans99/GutReferenceSet)),
and Sanders et al. (*Nat Microbiol* 8, 1039–1050, 2023).

### What the build produces

Two files, which are the two the pipeline needs:

- `<ref>.fna` — the catalogue, plus a `.sha256` beside it
- `<ref>.contig_to_species.tsv` — `species <TAB> contig <TAB> length`, derived
  from the reference rather than shipped, so the two cannot drift apart

The Bowtie2 index is **not** built here: `rule build_metagenome_index` does it
inside the workflow, with `--large-index`, which every one of these catalogues
needs. Run `build_references.py --index` only if you want the index outside a
workflow run.

Put the two paths into `config.yaml` — the build prints them:

```yaml
reference_fasta:   references/combined.fna
contig_to_species: references/combined.contig_to_species.tsv
```

### Proving it works before committing to the download

```bash
python workflow/scripts/build_references.py --ref carter --subset 12 --out proof/
```

That fetches 12 genomes (~16 MB), renames and concatenates them, and checks
every contig name and length against the manifest. It exercises both of
`carter`'s fetch routes and the full rename-and-verify path in under a minute.
`--verify-only <fasta>` checks a catalogue you already have.

### Cost of a full rebuild

| `--ref` | Downloads | Notes |
|---|---|---|
| `carter` | 3,971 requests | the slow one: one listing plus one fetch per GenBank accession. Hours, and resumable — the source cache is re-used. |
| `leviatan` | 1 file, 2.7 GB | fast |
| `chimp` | 1 file, 8.7 GB | fast |
| `combined` | re-uses `carter` + `leviatan` sources | build `carter` first and the fetches are already cached |

Budget roughly 15 GB for the source cache on top of the catalogue itself, and
a further ~22 GB if you let the workflow build the `combined` index.

### Why the manifests exist

`resources/reference_manifests/` holds, per catalogue, the genome order and
every contig's name and length (12 MB in total). The build follows them rather
than re-deriving a naming rule, because the rules are not uniform and the
differences are not guessable:

- `MGYG*` contigs are used **verbatim** — MGnify already names them
  `MGYG000004406_1`, which is why these indices start at 1 and the others at 0
- `GCA_*` and `Rep_*` are renamed `<genome>_c_<i>`, `i` from 0, in source order
- chimp contigs are prefixed, not renumbered: `<MAG>|<original contig>`, which
  is why `c_000000000001` and `k141_<n>` names sit side by side
- three MGYG genomes are short a few contigs against their source
  (`MGYG000001406.1`, `MGYG000001455.1`, `MGYG000001456.1`), and a fourth
  carries a `.1` suffix the rest do not

Those are properties of the catalogue the analysis actually read, not defects
to correct, so the manifest records them and the build reproduces them. A
mismatch in any contig name or length stops the build rather than producing a
reference that is subtly not the reported one.

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
 │  │  Snakemake orchestrator (1 CPU, 4 GB, long)     │  │
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

The orchestrator itself is a tiny job (1 CPU, 4 GB) that submits compute-heavy
work as separate SLURM jobs via `sbatch` and polls them with `sacct`. Because
it only sleeps and polls, it can sit on any low-priority partition — but it
has to **stay alive for the whole run**, and that is where walltime limits
bite.

Most partitions cap well below a long run. On Cannon every general partition
stops at 3 days:

| partition | max walltime |
|---|---|
| `shared`, `serial_requeue`, `sapphire`, `bigmem` | 3 days |
| `intermediate`, `bigmem_intermediate` | 14 days |
| `unrestricted` | 365 days |

`sbatch` **rejects** an over-limit request rather than trimming it, so asking
7 days on `shared` means the run never starts. `setup.sh` therefore asks which
partition to put the orchestrator on, defaults to the longest-walltime
partition it can find in `sinfo`, and clamps `--time` to what that partition
actually allows. The committed default is a safe 3 days on `shared`.

If the orchestrator does run out of time, re-submit: Snakemake skips completed
samples and continues.

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
