"""End-to-end CPU MGE classification pipeline (k-mer only)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .constants import KMER_K, MGE_DISPLAY_NAMES
from .fasta import FastaRecord
from .kmer import batch_kmer_features
from .models import ModelBundle, load_model_bundle, predict_mge, predict_pp_subtypes


@dataclass
class SequencePrediction:
    sequence_id: str
    length: int
    mge_class: str
    mge_class_display: str
    mge_confidence: float
    mge_probabilities: dict[str, float]
    pp_subtype: str | None = None
    pp_subtype_confidence: float | None = None
    pp_subtype_probabilities: dict[str, float] | None = None


class TriMGECPUClassifier:
    """
    triMGE-cpu two-stage classifier (no HyenaDNA / GPU):
      1) k5 + SVM -> phage / phage-plasmid / plasmid
      2) k5 + SVM on predicted phage-plasmids -> 15 subtypes
    """

    def __init__(self, models_dir: str | None = None) -> None:
        self.bundle = load_model_bundle(models_dir)

    def predict_records(self, records: Sequence[FastaRecord]) -> list[SequencePrediction]:
        if not records:
            return []

        sequences = [record.sequence for record in records]
        kmer_features = batch_kmer_features(sequences, k=KMER_K)

        labels, probabilities, class_names = predict_mge(self.bundle, kmer_features)
        predictions: list[SequencePrediction] = []

        pp_indices = [idx for idx, cls in enumerate(class_names) if cls == "phage_plasmid"]
        pp_subtype_names: dict[int, str] = {}
        pp_subtype_probs: dict[int, dict[str, float]] = {}
        pp_subtype_conf: dict[int, float] = {}

        if pp_indices:
            pp_kmer = kmer_features[pp_indices]
            _, pp_proba, pp_names = predict_pp_subtypes(self.bundle, pp_kmer)
            encoder_classes = list(self.bundle.pp_label_encoder.classes_)
            for local_idx, global_idx in enumerate(pp_indices):
                pp_subtype_names[global_idx] = pp_names[local_idx]
                pp_subtype_conf[global_idx] = float(np.max(pp_proba[local_idx]))
                pp_subtype_probs[global_idx] = {
                    encoder_classes[i]: float(pp_proba[local_idx][i])
                    for i in range(len(encoder_classes))
                }

        class_name_list = self.bundle.metadata["mge_model"]["class_names"]
        for idx, record in enumerate(records):
            class_name = class_names[idx]
            class_probs = {
                class_name_list[i]: float(probabilities[idx][i])
                for i in range(len(class_name_list))
            }
            prediction = SequencePrediction(
                sequence_id=record.sequence_id,
                length=len(record.sequence),
                mge_class=class_name,
                mge_class_display=MGE_DISPLAY_NAMES[class_name],
                mge_confidence=float(np.max(probabilities[idx])),
                mge_probabilities=class_probs,
            )
            if idx in pp_subtype_names:
                prediction.pp_subtype = pp_subtype_names[idx]
                prediction.pp_subtype_confidence = pp_subtype_conf[idx]
                prediction.pp_subtype_probabilities = pp_subtype_probs[idx]
            predictions.append(prediction)

        return predictions
