"""Load serialized sklearn models and metadata (CPU / k-mer only)."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import joblib
import numpy as np
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC

from .constants import (
    MGE_CLASS_NAMES,
    MGE_SVM_PARAMS,
    PP_SUBTYPE_CLASSES,
    PP_SVM_PARAMS,
)


@dataclass
class ModelBundle:
    classifier_3class: SVC
    mge_scaler: StandardScaler
    pp_classifier: SVC
    pp_scaler: StandardScaler
    pp_label_encoder: LabelEncoder
    metadata: dict


def default_models_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "models"


def load_model_bundle(models_dir: str | Path | None = None) -> ModelBundle:
    """Load the exported CPU 3-class and P-P subtype models."""
    root = Path(models_dir or default_models_dir())
    metadata_path = root / "class_metadata.json"
    mge_path = root / "mge_3class_svm.joblib"
    mge_scaler_path = root / "mge_3class_scaler.joblib"
    pp_model_path = root / "pp_subtype_svm.joblib"
    pp_scaler_path = root / "pp_subtype_scaler.joblib"
    pp_encoder_path = root / "pp_label_encoder.joblib"

    required = [
        metadata_path,
        mge_path,
        mge_scaler_path,
        pp_model_path,
        pp_scaler_path,
        pp_encoder_path,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing model artifacts:\n  - "
            + "\n  - ".join(missing)
            + "\nRun `python scripts/export_models.py` first."
        )

    with open(metadata_path) as handle:
        metadata = json.load(handle)

    return ModelBundle(
        classifier_3class=joblib.load(mge_path),
        mge_scaler=joblib.load(mge_scaler_path),
        pp_classifier=joblib.load(pp_model_path),
        pp_scaler=joblib.load(pp_scaler_path),
        pp_label_encoder=joblib.load(pp_encoder_path),
        metadata=metadata,
    )


def save_model_bundle(bundle: ModelBundle, models_dir: str | Path) -> None:
    root = Path(models_dir)
    root.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle.classifier_3class, root / "mge_3class_svm.joblib")
    joblib.dump(bundle.mge_scaler, root / "mge_3class_scaler.joblib")
    joblib.dump(bundle.pp_classifier, root / "pp_subtype_svm.joblib")
    joblib.dump(bundle.pp_scaler, root / "pp_subtype_scaler.joblib")
    joblib.dump(bundle.pp_label_encoder, root / "pp_label_encoder.joblib")
    with open(root / "class_metadata.json", "w") as handle:
        json.dump(bundle.metadata, handle, indent=2)


def build_default_metadata() -> dict:
    return {
        "tool": "triMGE-cpu",
        "tool_version": "0.1.0",
        "mge_model": {
            "name": "k5 + SVM",
            "features": ["kmer_k5"],
            "class_names": MGE_CLASS_NAMES,
            "hyperparameters": MGE_SVM_PARAMS,
            "requires_gpu": False,
        },
        "pp_subtype_model": {
            "name": "k5 + SVM",
            "features": ["kmer_k5"],
            "class_names": PP_SUBTYPE_CLASSES,
            "hyperparameters": PP_SVM_PARAMS,
            "requires_gpu": False,
        },
    }


def predict_mge(
    bundle: ModelBundle,
    kmer_features: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Predict 3-class labels and probabilities from k-mer features."""
    scaled = bundle.mge_scaler.transform(kmer_features)
    labels = bundle.classifier_3class.predict(scaled)
    probabilities = bundle.classifier_3class.predict_proba(scaled)
    class_names = [MGE_CLASS_NAMES[int(label)] for label in labels]
    return labels, probabilities, class_names


def predict_pp_subtypes(
    bundle: ModelBundle,
    kmer_features: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Predict P-P subtype labels and probabilities."""
    scaled = bundle.pp_scaler.transform(kmer_features)
    encoded = bundle.pp_classifier.predict(scaled)
    probabilities = bundle.pp_classifier.predict_proba(scaled)
    subtype_names = list(bundle.pp_label_encoder.inverse_transform(encoded))
    return encoded, probabilities, subtype_names
