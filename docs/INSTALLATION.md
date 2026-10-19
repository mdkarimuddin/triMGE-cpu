# Installation

## Requirements

- Python ≥ 3.11
- Dependencies: `numpy`, `scikit-learn`, `joblib`, `pandas`, `openpyxl`
- Optional for figures: `matplotlib` (install with `[report]`)

On CSC, load the Python 3.12 module before creating a virtual environment. The default `python3` on a login node can be 3.9, which cannot install scikit-learn 1.8.

```bash
module load python-data
```

## Install the package

```bash
cd /scratch/project_xxxxxxx/triMGE-cpu
pip install -e ".[report]"
```

This registers:

| Command | Purpose |
|---------|---------|
| `trimge-cpu` | Classify one FASTA |
| `trimge-cpu-combine` | Merge FASTAs |
| `trimge-cpu-report` | Build tables/figures from predictions |
| `trimge-cpu-pipeline` | Full easy pipeline |

If a console script has a broken shebang (old env), run via module:

```bash
python -m trimgecpu -i sequences.fasta -o predictions.tsv --models-dir models
```

## Export models (one-time)

Models are already present under `models/` if you used the project export. To rebuild:

```bash
python scripts/export_models.py \
  --data-dir /scratch/project_xxxxxxx/unified_mge_analysis/data \
  --s1-table /scratch/project_xxxxxxx/phage-plasmid-analysis/Phage-plasmids/Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx \
  --output-dir models \
  --fit-on train_val
```

Expected files in `models/`:

- `mge_3class_svm.joblib`
- `mge_3class_scaler.joblib`
- `pp_subtype_svm.joblib`
- `pp_subtype_scaler.joblib`
- `pp_label_encoder.joblib`
- `class_metadata.json`

See [MODELS.md](MODELS.md).

## Verify

```bash
python -c "from trimgecpu.models import load_model_bundle; print(load_model_bundle('models').metadata['tool'])"
# → triMGE-cpu

pytest -q
```
