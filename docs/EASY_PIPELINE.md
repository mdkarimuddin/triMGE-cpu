# Easy pipeline (recommended)

One command to:

1. **Combine** FASTA files (optional)
2. **Classify** with k-mer k5 + SVM (3-class + P–P subtypes)
3. **Write automated reports** — tables and figures for both stages

No HyenaDNA. No GPU.

## What you need

- A folder of nucleotide FASTAs (`.fna` / `.fa` / `.fasta`), **or** one combined FASTA
- Exported `models/` directory (see [INSTALLATION.md](INSTALLATION.md))
- Package installed with report extras: `pip install -e ".[report]"`

## Quick start

### Folder of FASTAs

```bash
module load python-data   # on CSC
cd /scratch/project_xxxxxxx/triMGE-cpu

trimge-cpu-pipeline \
  --input-dir /path/to/fastas \
  --workdir /path/to/my_run \
  --models-dir models
```

### Single combined FASTA

```bash
trimge-cpu-pipeline \
  --fasta /path/to/all_sequences.fna \
  --workdir /path/to/my_run \
  --models-dir models
```

### Explicit list of files

```bash
trimge-cpu-pipeline \
  -i sample1.fna sample2.fna sample3.fa \
  --workdir /path/to/my_run \
  --models-dir models
```

## Output layout

```text
my_run/
  01_combined/
    all_sequences.fna
    all_sequences.fna.manifest.tsv   # maps combined IDs → source files
  02_predictions/
    predictions.tsv                  # 3-class + PP subtype columns
  03_report/
    01_three_class/
      tables/                        # counts, confidence, length summaries
      figures/                       # PNG + PDF plots
    02_pp_subtype/
      tables/                        # subtype counts (phage-plasmids only)
      figures/                       # subtype composition plots
    pipeline_summary.json
    README.md
  pipeline_summary.json
```

### Predictions columns (`02_predictions/predictions.tsv`)

| Column | Meaning |
|--------|---------|
| `sequence_id` | Sequence header (may include `filename\|id` prefix) |
| `length` | Sequence length (bp) |
| `mge_class` | `phage` / `phage_plasmid` / `plasmid` |
| `mge_class_display` | Display label |
| `mge_confidence` | Max class probability |
| `pp_subtype` | Subtype name (only if `phage_plasmid`) |
| `pp_subtype_confidence` | Max subtype probability |

## Useful options

| Flag | Effect |
|------|--------|
| `--recursive` | Search FASTAs under subfolders of `--input-dir` |
| `--no-prefix` | Do not prefix IDs with source filename when combining |
| `--min-length N` | Drop sequences shorter than N bp after cleaning (default: 100) |
| `--skip-classify` | Reuse existing `02_predictions/predictions.tsv` |
| `--skip-report` | Stop after classification (no figures) |
| `--report-steps 1` | Only 3-class report |
| `--report-steps 2` | Only P–P subtype report |
| `--report-steps 1,2` | Both (default) |
| `--verbose` | Debug logging |

### Examples

Classify only (no report):

```bash
trimge-cpu-pipeline \
  --input-dir /path/to/fastas \
  --workdir /path/to/my_run \
  --models-dir models \
  --skip-report
```

Rebuild reports from an existing run:

```bash
trimge-cpu-pipeline \
  --fasta /path/to/my_run/01_combined/all_sequences.fna \
  --workdir /path/to/my_run \
  --models-dir models \
  --skip-classify
```

Or call the report CLI directly:

```bash
trimge-cpu-report \
  -i /path/to/my_run/02_predictions/predictions.tsv \
  -o /path/to/my_run/03_report
```

## What the reports cover

- **Step 1 (3-class):** class counts, confidence / length summaries, bar and distribution plots  
- **Step 2 (P–P subtypes):** runs only on sequences predicted as `phage_plasmid`; subtype counts and composition figures  

Details: [REPORTS.md](REPORTS.md).

## Notes

- Sequences shorter than `--min-length` (default 100 bp) are skipped.
- Contigs ≥ ~10 kb work best (see manuscript fragment analysis).
- This package is the CPU baseline (~91.4% 3-class). For max accuracy (+0.9 pp) use the hybrid `triMGE` package with HyenaDNA + GPU.
