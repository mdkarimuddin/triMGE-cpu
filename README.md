# triMGE-cpu

[![tests](https://github.com/mdkarimuddin/triMGE-cpu/actions/workflows/tests.yml/badge.svg)](https://github.com/mdkarimuddin/triMGE-cpu/actions/workflows/tests.yml)

triMGE-cpu classifies mobile genetic element sequences with k-mer models. A run first labels each sequence as phage, phage-plasmid, or plasmid, then assigns a phage-plasmid subtype to the sequences called phage-plasmid.

The default pipeline writes both stages:

1. **3-class classification** — CSV tables and PDF figures in `03_report/01_three_class/`
2. **P-P classification** — CSV tables and PDF figures in `03_report/02_pp_subtype/`

Step 2 figures are written when at least one sequence is predicted as phage-plasmid. Trained models are already in `models/`. You need Python 3.11 or newer. Input files are nucleotide FASTA (`.fa`, `.fna`, or `.fasta`).

## Results

The manuscript numbers are the models fit on the training split only. The weights in `models/` were fit on the training and validation sequences and then scored on the same 638 test sequences.

| Stage | Classifier | Manuscript | Shipped weights |
|-------|------------|------------|-----------------|
| Phage / phage-plasmid / plasmid | 5-mer SVM | 91.4% accuracy and macro F1 | 91.9% accuracy and macro F1 |
| Phage-plasmid subtype | 5-mer SVM | 92.2% accuracy, 92.8% macro F1 | See [docs/MODELS.md](docs/MODELS.md) |

The subtype percentage is the manuscript result on true phage-plasmid sequences. Accuracy stays high on contigs of about 10 kb or longer and is 69.9% at 1 kb. This package still classifies sequences down to 100 bp.

## Beginner steps

**1. Clone the repository and enter it**

```bash
git clone https://github.com/mdkarimuddin/triMGE-cpu.git
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

The GPU package, triMGE-gpu, is distributed separately.

## Public demo

https://trimge-cpu-wxkbzdebrh3nejpcjdhjum.streamlit.app/

`app.py` is that page. It accepts a pasted FASTA, an uploaded file, or the included example. One public run is limited to 100 sequences and 2,000,000 bp. Keep scikit-learn in the 1.8 series so the shipped `.joblib` files load. The page needs Streamlit 1.51 or newer.

```bash
pip install -r requirements.txt
streamlit run app.py
```

## License

triMGE-cpu is released under the [MIT License](LICENSE) (copyright 2026 Md Karim Uddin). You may use, copy, modify, and distribute this package, including the trained classifiers in `models/`, for academic or commercial work. Keep the copyright notice and the license text with any copy.

The package is provided as is, without warranty.
