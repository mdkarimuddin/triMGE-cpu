# Models (triMGE-cpu)

After `python scripts/export_models.py`, this directory contains:

| File | Contents |
|------|----------|
| `mge_3class_svm.joblib` | k5 + SVM for 3-class MGE |
| `mge_3class_scaler.joblib` | StandardScaler for 3-class features |
| `pp_subtype_svm.joblib` | k5 + SVM for 15 P-P subtypes |
| `pp_subtype_scaler.joblib` | StandardScaler for P-P features |
| `pp_label_encoder.joblib` | LabelEncoder → subtype names |
| `class_metadata.json` | Model names, classes, hyperparameters |

All models are included in this directory. Load them with scikit-learn 1.8.x, the range pinned in `pyproject.toml`. A newer scikit-learn can refuse to read these `.joblib` files.
