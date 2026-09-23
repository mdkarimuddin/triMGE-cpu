# triMGE-cpu

**triMGE-cpu** — CPU-only sequence classification of mobile genetic elements.

No HyenaDNA. No GPU. No PyTorch.

Uses the manuscript's best **k-mer-only** models:

1. **3-class MGE** — phage / phage-plasmid / plasmid  
   Model: **k-mer k5 + SVM** (~91.4% held-out accuracy)

2. **P-P subtyping** — only for predicted phage-plasmids  
   Model: **k-mer k5 + SVM** (15 subtypes, ~92.2% accuracy)

```text
FASTA input
   │
   └─ k-mer k5 features (1024-d)
         │
         ▼
   3-class SVM  ──► phage / phage-plasmid / plasmid
         │
         └─ if phage-plasmid ──► k5 + SVM ──► 15 P-P subtypes
```

Compared with full `triMGE` (HyenaDNA + k5 RF = 92.3%), this package trades ~0.9 pp of three-class accuracy for a light, portable, CPU-only install.

## Install

```bash
cd /scratch/project_xxxxxxx/triMGE-cpu
pip install -e ".[report]"
```

## Export models (one-time)

```bash
python scripts/export_models.py \
  --data-dir /scratch/project_xxxxxxx/unified_mge_analysis/data \
  --s1-table /scratch/project_xxxxxxx/phage-plasmid-analysis/Phage-plasmids/Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx \
  --output-dir models \
  --fit-on train_val
```

## Classify

```bash
trimge-cpu -i sequences.fasta -o predictions.tsv --models-dir models
```

JSON:

```bash
trimge-cpu -i sequences.fasta -o predictions.json --format json --models-dir models
```

## Easy pipeline (recommended)

Combine FASTAs → classify → **automated 3-class + P–P subtype reports**:

```bash
trimge-cpu-pipeline \
  --input-dir /path/to/fastas \
  --workdir /path/to/run_output \
  --models-dir models
```

Full guide: **[docs/EASY_PIPELINE.md](docs/EASY_PIPELINE.md)**

## Commands

| Command | Role |
|---------|------|
| `trimge-cpu` | Classify FASTA |
| `trimge-cpu-combine` | Merge FASTAs |
| `trimge-cpu-report` | Tables/figures from predictions TSV |
| `trimge-cpu-pipeline` | Combine → classify → report |

## Documentation

| Guide | Description |
|-------|-------------|
| [docs/README.md](docs/README.md) | Documentation index |
| [Installation](docs/INSTALLATION.md) | Install + export models |
| [Easy pipeline](docs/EASY_PIPELINE.md) | End-to-end usage (start here) |
| [Usage](docs/USAGE.md) | Individual CLI options |
| [Reports](docs/REPORTS.md) | Automated 3-class and P–P report outputs |
| [Models](docs/MODELS.md) | Model details |

## Related

- GPU / hybrid package: `/scratch/project_xxxxxxx/triMGE`
- Manuscript: `unified_mge_analysis/Updated_analyses_manuscript/manuscript_v8.md`

## License

MIT
