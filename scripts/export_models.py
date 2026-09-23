#!/usr/bin/env python3
"""
Export k-mer k5 SVM models from the training artifacts.

Example:
  python scripts/export_models.py \\
    --data-dir /scratch/project_xxxxxxx/unified_mge_analysis/data \\
    --s1-table /scratch/project_xxxxxxx/phage-plasmid-analysis/Phage-plasmids/Publication_related_data/Supplementary_data/S1_All_MGE_information_table.xlsx \\
    --output-dir models \\
    --fit-on train_val
"""

from __future__ import annotations

import argparse
import logging
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from trimgecpu.constants import (  # noqa: E402
    MGE_CLASS_NAMES,
    MGE_SVM_PARAMS,
    PP_SVM_PARAMS,
)
from trimgecpu.models import ModelBundle, build_default_metadata, save_model_bundle  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def get_broad_class(mge_type: str | None) -> str:
    if not mge_type:
        return "Unknown"
    text = str(mge_type)
    mapping = [
        (("P1_1", "P1_2", "P1 community"), "P1"),
        (("cp32",), "cp32"),
        (("N15",), "N15"),
        (("pKpn",), "pKpn"),
        (("pSLy3",), "pSLy3"),
        (("SSU5",), "SSU5_pHCM2"),
        (("pMT1",), "pMT1"),
        (("AB_1", "AB_2", "AB community"), "AB"),
        (("PiSa",), "PiSa"),
        (("IEBH",), "IEBH-Waukesha"),
        (("VP882",), "VP882"),
        (("p48",), "p48"),
        (("Actino",), "Actinophage_A"),
        (("BigBertha",), "BigBertha"),
    ]
    for needles, label in mapping:
        if any(needle in text for needle in needles):
            return label
    if "Plasmid" in text:
        return "Misclassified"
    return "Singleton"


def load_pp_labels(s1_path: Path, metadata_path: Path, min_samples: int = 20):
    import openpyxl

    workbook = openpyxl.load_workbook(s1_path, read_only=True)
    worksheet = workbook["Sheet1"]
    accession_to_type = {}
    for row_idx, row in enumerate(worksheet.iter_rows(values_only=True)):
        if row_idx == 0:
            continue
        accession_to_type[row[0]] = row[7]
    workbook.close()

    metadata = pd.read_csv(metadata_path)
    pp_df = metadata[metadata["label"] == 1].copy()
    pp_df["accession"] = pp_df["sequence_id"].apply(
        lambda value: value.rsplit(".", 1)[0] if "." in value else value
    )
    pp_df["mge_detail"] = pp_df["accession"].map(accession_to_type)
    missing = pp_df["mge_detail"].isna()
    if missing.any():
        pp_df.loc[missing, "mge_detail"] = pp_df.loc[missing, "sequence_id"].map(accession_to_type)
    pp_df["broad_class"] = pp_df["mge_detail"].apply(get_broad_class)

    counts = pp_df["broad_class"].value_counts()
    valid_classes = [
        cls
        for cls, count in counts.items()
        if count >= min_samples and cls not in {"Misclassified", "Unknown"}
    ]
    pp_filtered = pp_df[pp_df["broad_class"].isin(valid_classes)].copy()
    return pp_filtered.index.values, pp_filtered["broad_class"].values, valid_classes


def select_fit_indices(split_data: dict, fit_on: str, n_total: int) -> np.ndarray:
    if fit_on == "train":
        return split_data["train_idx"]
    if fit_on == "train_val":
        return np.concatenate([split_data["train_idx"], split_data["val_idx"]])
    if fit_on == "all":
        return np.arange(n_total)
    raise ValueError(f"Unknown fit_on value: {fit_on}")


def export_models(
    data_dir: Path,
    s1_table: Path,
    output_dir: Path,
    fit_on: str = "train_val",
) -> None:
    with open(data_dir / "split_indices.pkl", "rb") as handle:
        split_data = pickle.load(handle)
    with open(data_dir / "sequences_all.pkl", "rb") as handle:
        seq_data = pickle.load(handle)

    kmer_all = np.load(data_dir / "kmer_k5_all.npy")
    labels_all = seq_data["labels"]

    fit_idx = select_fit_indices(split_data, fit_on, len(labels_all))
    mge_scaler = StandardScaler()
    x_fit = mge_scaler.fit_transform(kmer_all[fit_idx])
    mge_model = SVC(**MGE_SVM_PARAMS)
    mge_model.fit(x_fit, labels_all[fit_idx])
    logger.info("Trained 3-class SVM (k5) on %s samples (%s)", len(fit_idx), fit_on)

    # Quick held-out check when possible
    if fit_on != "all" and "test_idx" in split_data:
        test_idx = split_data["test_idx"]
        x_te = mge_scaler.transform(kmer_all[test_idx])
        acc = float((mge_model.predict(x_te) == labels_all[test_idx]).mean())
        logger.info("Held-out test accuracy (3-class k5 SVM): %.4f", acc)

    pp_indices, pp_labels, pp_classes = load_pp_labels(s1_table, data_dir / "metadata_all.csv")
    pp_kmer = kmer_all[pp_indices]
    pp_encoder = LabelEncoder()
    pp_y = pp_encoder.fit_transform(pp_labels)

    pp_scaler = StandardScaler()
    pp_x = pp_scaler.fit_transform(pp_kmer)
    pp_model = SVC(**PP_SVM_PARAMS)
    pp_model.fit(pp_x, pp_y)
    logger.info(
        "Trained P-P SVM on %s phage-plasmid sequences (%s classes)",
        len(pp_y),
        len(pp_classes),
    )

    metadata = build_default_metadata()
    metadata["export"] = {
        "fit_on": fit_on,
        "n_mge_samples": int(len(fit_idx)),
        "n_pp_samples": int(len(pp_y)),
        "pp_classes": list(pp_encoder.classes_),
        "mge_class_names": MGE_CLASS_NAMES,
        "note": "k-mer k5 SVM models shipped with the package",
    }

    bundle = ModelBundle(
        classifier_3class=mge_model,
        mge_scaler=mge_scaler,
        pp_classifier=pp_model,
        pp_scaler=pp_scaler,
        pp_label_encoder=pp_encoder,
        metadata=metadata,
    )
    save_model_bundle(bundle, output_dir)
    logger.info("Saved model artifacts to %s", output_dir)


def main() -> int:
    parser = argparse.ArgumentParser(description="Export CPU inference models for triMGE-cpu")
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--s1-table", required=True, type=Path)
    parser.add_argument("--output-dir", default=REPO_ROOT / "models", type=Path)
    parser.add_argument(
        "--fit-on",
        choices=["train", "train_val", "all"],
        default="train_val",
        help="Which unified split partition(s) to fit on (default: train_val)",
    )
    args = parser.parse_args()

    export_models(
        data_dir=args.data_dir,
        s1_table=args.s1_table,
        output_dir=args.output_dir,
        fit_on=args.fit_on,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
