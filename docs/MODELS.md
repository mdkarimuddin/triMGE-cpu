# Models (triMGE-cpu)

## Deployed models

| Task | Features | Classifier | Manuscript test | Shipped weights on that test |
|------|----------|------------|-----------------|------------------------------|
| 3-class MGE | k-mer k5 (1024-d) | SVM (RBF, C=1.0) | 91.4% accuracy and macro F1 | 91.9% accuracy and macro F1 |
| P-P subtyping | k-mer k5 (1024-d) | SVM (RBF, C=10) | 92.2% accuracy, 92.8% macro F1 | See below |

Both stages use `StandardScaler` before SVM.

The manuscript percentages are the models fit on the training split only. The weights in `models/` were fit on train plus validation. Scored on the same 638 test sequences, those shipped weights reach 91.9% accuracy and 91.9% macro F1.

The 92.2% subtype accuracy is the manuscript result on true phage-plasmid sequences. Stage 2 runs only after stage 1 calls a sequence phage-plasmid. The manuscript length check for this k5 SVM stays high on contigs of at least about 10 kb and falls to 69.9% at 1 kb. This package still classifies sequences down to 100 bp.

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

The benchmark tables are not in this repository. The command below is how `models/` was written. Another machine cannot rerun it without those tables.

```bash
python scripts/export_models.py \
  --data-dir /scratch/project_xxxxxxx/unified_mge_analysis/data \
  --s1-table /scratch/project_xxxxxxx/phage-plasmid-analysis/Phage-plasmids/Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx \
  --output-dir models \
  --fit-on train_val
```

