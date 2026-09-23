"""Load and validate triMGE prediction tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REQUIRED_COLS = [
    "sequence_id",
    "length",
    "mge_class",
    "mge_confidence",
]


def load_predictions(path: Path | str) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Predictions missing columns {missing}; got {list(df.columns)}")
    df = df.copy()
    df["mge_class"] = df["mge_class"].astype(str)
    df["length"] = pd.to_numeric(df["length"], errors="coerce")
    df["mge_confidence"] = pd.to_numeric(df["mge_confidence"], errors="coerce")
    if "pp_subtype" not in df.columns:
        df["pp_subtype"] = ""
    if "pp_subtype_confidence" not in df.columns:
        df["pp_subtype_confidence"] = float("nan")
    else:
        df["pp_subtype_confidence"] = pd.to_numeric(df["pp_subtype_confidence"], errors="coerce")
    df["pp_subtype"] = df["pp_subtype"].fillna("").astype(str)
    return df


def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)
