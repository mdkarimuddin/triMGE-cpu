# triMGE-cpu

triMGE-cpu classifies mobile genetic element sequences on CPU with k-mer models. It does not use HyenaDNA, a GPU, or PyTorch. A run first labels each sequence as phage, phage-plasmid, or plasmid, then assigns a phage-plasmid subtype to the sequences called phage-plasmid.

The default pipeline writes both stages:

1. **3-class classification** — CSV tables and PDF figures in `03_report/01_three_class/`
2. **P-P classification** — CSV tables and PDF figures in `03_report/02_pp_subtype/`

Step 2 figures are written when at least one sequence is predicted as phage-plasmid. Trained models are in `models/`. Install the report extra so the PDFs are produced.

```bash
pip install -e ".[report]"

trimge-cpu-pipeline \
  --input-dir /path/to/fastas \
  --workdir /path/to/run_output \
  --models-dir models
```

Outputs: `run_output/01_combined`, `02_predictions`, and `03_report`. Details: [docs/EASY_PIPELINE.md](docs/EASY_PIPELINE.md).

Classify only, without reports:

```bash
trimge-cpu -i sequences.fasta -o predictions.tsv --models-dir models
```

Rebuild reports from an existing predictions table:

```bash
trimge-cpu-report -i predictions.tsv -o report_out
```

| Command | Role |
|---------|------|
| `trimge-cpu` | Classify a FASTA |
| `trimge-cpu-combine` | Merge FASTAs |
| `trimge-cpu-report` | CSV and PDF reports for both steps |
| `trimge-cpu-pipeline` | Combine, classify, then report |

| Guide | Description |
|-------|-------------|
| [docs/README.md](docs/README.md) | Documentation index |
| [Installation](docs/INSTALLATION.md) | Install |
| [Easy pipeline](docs/EASY_PIPELINE.md) | End-to-end run |
| [Usage](docs/USAGE.md) | CLI options |
| [Reports](docs/REPORTS.md) | CSV and PDF outputs for both steps |
| [Models](docs/MODELS.md) | Classifiers included in `models/` |

GPU package: [triMGE-gpu](https://version.helsinki.fi/mduddin/triMGE-gpu).

## License

MIT
