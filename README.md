# triMGE-cpu

triMGE-cpu classifies mobile genetic element sequences with k-mer models. A run first labels each sequence as phage, phage-plasmid, or plasmid, then assigns a phage-plasmid subtype to the sequences called phage-plasmid.

The default pipeline writes both stages:

1. **3-class classification** — CSV tables and PDF figures in `03_report/01_three_class/`
2. **P-P classification** — CSV tables and PDF figures in `03_report/02_pp_subtype/`

Step 2 figures are written when at least one sequence is predicted as phage-plasmid. Trained models are already in `models/`. You need Python 3.10 or newer. Input files are nucleotide FASTA (`.fa`, `.fna`, or `.fasta`).

## Beginner steps

**1. Clone the repository and enter it**

```bash
git clone https://version.helsinki.fi/mduddin/triMGE-cpu.git
cd triMGE-cpu
```

**2. Install once**

On CSC, load Python 3.12 first. The default `python3` on a login node can be 3.9, which cannot install this package.

```bash
module load python-data
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[report]"
```

The `[report]` extra installs Matplotlib, which writes the PDF figures. Do not skip it.

**3. Put your sequences in one folder**

```text
my_fastas/
  sample_a.fasta
  sample_b.fasta
```

One file is enough. Sequences shorter than 100 bp are skipped.

**4. Run**

```bash
trimge-cpu-pipeline \
  --input-dir /path/to/my_fastas \
  --workdir /path/to/my_run \
  --models-dir models
```

If you already have a single FASTA, use `--fasta /path/to/sequences.fasta` instead of `--input-dir`.

Use a new `--workdir` for each new set of sequences. Keep the virtual environment active (`source .venv/bin/activate`) in later sessions.

**5. Open the results**

```text
my_run/02_predictions/predictions.tsv
my_run/03_report/01_three_class/     phage / phage-plasmid / plasmid
my_run/03_report/02_pp_subtype/      phage-plasmid subtypes
```

In `predictions.tsv`, `mge_class` is the 3-class label and `mge_confidence` is its confidence. `pp_subtype` and `pp_subtype_confidence` are filled only for phage-plasmid sequences.

## Other commands

| Command | Role |
|---------|------|
| `trimge-cpu` | Classify a FASTA |
| `trimge-cpu-combine` | Merge FASTAs |
| `trimge-cpu-report` | CSV and PDF reports for both steps |
| `trimge-cpu-pipeline` | Combine, classify, then report |

Classify only, without reports:

```bash
trimge-cpu -i sequences.fasta -o predictions.tsv --models-dir models
```

Rebuild reports from an existing predictions table:

```bash
trimge-cpu-report -i predictions.tsv -o report_out
```

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

triMGE-cpu is released under the [MIT License](LICENSE) (copyright 2026 Md Karim Uddin). You may use, copy, modify, and distribute this package, including the trained classifiers in `models/`, for academic or commercial work. Keep the copyright notice and the license text with any copy.

The package is provided as is, without warranty.
