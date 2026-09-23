# Models (triMGE-cpu)

## Deployed models

| Task | Features | Classifier | Manuscript metric |
|------|----------|------------|-------------------|
| 3-class MGE | k-mer k5 (1024-d) | SVM (RBF, C=1.0) | 91.4% held-out (train-only fit in manuscript; train+val export ≈91.9%) |
| P-P subtyping | k-mer k5 (1024-d) | SVM (RBF, C=10) | 92.2% accuracy |

Both stages use `StandardScaler` before SVM. No HyenaDNA. No GPU.

## Artifacts

| File | Contents |
|------|----------|
| `mge_3class_svm.joblib` | Fitted 3-class SVM |
| `mge_3class_scaler.joblib` | Scaler for 3-class features |
| `pp_subtype_svm.joblib` | Fitted P-P SVM |
| `pp_subtype_scaler.joblib` | Scaler for P-P features |
| `pp_label_encoder.joblib` | Subtype name encoder |
| `class_metadata.json` | Metadata |

## Export

```bash
python scripts/export_models.py \
  --data-dir /scratch/project_xxxxxxx/unified_mge_analysis/data \
  --s1-table /scratch/project_xxxxxxx/phage-plasmid-analysis/Phage-plasmids/Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx \
  --output-dir models \
  --fit-on train_val
```

## Why not HyenaDNA?

The hybrid HyenaDNA+k5 RF reaches 92.3% (+0.9 pp). For portable CPU deployment, k-mer SVM is the recommended default.
