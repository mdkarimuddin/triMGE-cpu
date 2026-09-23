"""Step 1: 3-class summary tables + visualizations."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .io_utils import ensure_dirs
from .style import CLASS_COLORS, CLASS_DISPLAY, CLASS_ORDER, apply_mpl_style


def _save(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(f"{stem}.png", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def generate_tables(df: pd.DataFrame, tables_dir: Path) -> dict[str, Path]:
    ensure_dirs(tables_dir)
    out: dict[str, Path] = {}

    counts = (
        df["mge_class"]
        .value_counts()
        .reindex(CLASS_ORDER)
        .fillna(0)
        .astype(int)
        .rename_axis("mge_class")
        .reset_index(name="n")
    )
    counts["pct"] = (100 * counts["n"] / counts["n"].sum()).round(2)
    counts["mge_class_display"] = counts["mge_class"].map(CLASS_DISPLAY)
    p = tables_dir / "three_class_counts.csv"
    counts.to_csv(p, index=False)
    out["counts"] = p

    conf = (
        df.groupby("mge_class", as_index=False)["mge_confidence"]
        .agg(n="count", mean="mean", median="median", std="std", min="min", max="max")
    )
    conf = conf.set_index("mge_class").reindex(CLASS_ORDER).reset_index()
    conf["mge_class_display"] = conf["mge_class"].map(CLASS_DISPLAY)
    p = tables_dir / "three_class_confidence_summary.csv"
    conf.to_csv(p, index=False)
    out["confidence"] = p

    length = (
        df.groupby("mge_class", as_index=False)["length"]
        .agg(n="count", mean="mean", median="median", std="std", min="min", max="max")
    )
    length = length.set_index("mge_class").reindex(CLASS_ORDER).reset_index()
    length["mge_class_display"] = length["mge_class"].map(CLASS_DISPLAY)
    p = tables_dir / "three_class_length_summary.csv"
    length.to_csv(p, index=False)
    out["length"] = p

    # per-sequence table (normalized copy for downstream)
    keep = [
        c
        for c in [
            "sequence_id",
            "length",
            "mge_class",
            "mge_class_display",
            "mge_confidence",
            "pp_subtype",
            "pp_subtype_confidence",
        ]
        if c in df.columns
    ]
    seq_path = tables_dir / "three_class_predictions.tsv"
    df[keep].to_csv(seq_path, sep="\t", index=False)
    out["predictions"] = seq_path
    return out


def generate_figures(df: pd.DataFrame, figures_dir: Path) -> list[Path]:
    ensure_dirs(figures_dir)
    apply_mpl_style()
    written: list[Path] = []

    counts = df["mge_class"].value_counts().reindex(CLASS_ORDER).fillna(0)
    fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
    x = np.arange(len(CLASS_ORDER))
    ax.bar(
        x,
        counts.values,
        color=[CLASS_COLORS[c] for c in CLASS_ORDER],
        edgecolor="white",
        width=0.7,
    )
    for i, v in enumerate(counts.values):
        ax.text(i, v + max(counts.values) * 0.02, f"{int(v)}", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels([CLASS_DISPLAY[c] for c in CLASS_ORDER])
    ax.set_ylabel("Sequences")
    ax.set_title(f"triMGE 3-class calls (n={len(df)})", loc="left", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    stem = figures_dir / "01_three_class_counts"
    _save(fig, stem)
    written.append(Path(f"{stem}.png"))

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)
    data, labels, colors = [], [], []
    for c in CLASS_ORDER:
        vals = df.loc[df["mge_class"] == c, "mge_confidence"].dropna().values
        if len(vals) == 0:
            continue
        data.append(vals)
        labels.append(f"{CLASS_DISPLAY[c]}\n(n={len(vals)})")
        colors.append(CLASS_COLORS[c])
    bp = axes[0].boxplot(data, tick_labels=labels, patch_artist=True, showfliers=False)
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.8)
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("Confidence")
    axes[0].set_title("A  Confidence by class", loc="left", fontweight="bold")
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    data, labels, colors = [], [], []
    for c in CLASS_ORDER:
        vals = df.loc[df["mge_class"] == c, "length"].dropna().values
        if len(vals) == 0:
            continue
        data.append(vals)
        labels.append(f"{CLASS_DISPLAY[c]}\n(n={len(vals)})")
        colors.append(CLASS_COLORS[c])
    bp = axes[1].boxplot(data, tick_labels=labels, patch_artist=True, showfliers=False)
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.8)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Length (bp, log)")
    axes[1].set_title("B  Length by class", loc="left", fontweight="bold")
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)
    stem = figures_dir / "02_three_class_confidence_length"
    _save(fig, stem)
    written.append(Path(f"{stem}.png"))

    # pie / donut
    fig, ax = plt.subplots(figsize=(5.2, 5.0), constrained_layout=True)
    sizes = [int(counts[c]) for c in CLASS_ORDER]
    ax.pie(
        sizes,
        labels=[CLASS_DISPLAY[c] for c in CLASS_ORDER],
        colors=[CLASS_COLORS[c] for c in CLASS_ORDER],
        autopct=lambda p: f"{p:.1f}%" if p > 0 else "",
        startangle=90,
        wedgeprops=dict(width=0.42, edgecolor="white"),
        pctdistance=0.75,
    )
    ax.set_title(f"3-class composition (n={len(df)})", fontweight="bold")
    stem = figures_dir / "03_three_class_composition"
    _save(fig, stem)
    written.append(Path(f"{stem}.png"))
    return written


def run_step1(df: pd.DataFrame, step_dir: Path) -> dict:
    tables_dir = step_dir / "tables"
    figures_dir = step_dir / "figures"
    ensure_dirs(step_dir, tables_dir, figures_dir)
    tables = generate_tables(df, tables_dir)
    figures = generate_figures(df, figures_dir)
    summary = {
        "n_sequences": int(len(df)),
        "class_counts": df["mge_class"].value_counts().to_dict(),
        "tables": {k: str(v) for k, v in tables.items()},
        "figures": [str(f) for f in figures],
    }
    (step_dir / "step1_summary.json").write_text(
        __import__("json").dumps(summary, indent=2) + "\n"
    )
    return summary
