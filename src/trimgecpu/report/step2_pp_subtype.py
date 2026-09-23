"""Step 2: PP-subtype summary tables + visualizations (phage_plasmid only)."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .io_utils import ensure_dirs
from .style import SUBTYPE_COLORS, apply_mpl_style, subtype_color


def _save(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(f"{stem}.png", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def filter_pp(df: pd.DataFrame) -> pd.DataFrame:
    pp = df[df["mge_class"] == "phage_plasmid"].copy()
    pp["pp_subtype"] = pp["pp_subtype"].replace("", "Unknown").fillna("Unknown")
    # empty string after strip
    pp.loc[pp["pp_subtype"].astype(str).str.strip().eq(""), "pp_subtype"] = "Unknown"
    return pp


def _subtype_order(pp: pd.DataFrame) -> list[str]:
    preferred = [s for s in SUBTYPE_COLORS if s != "Unknown"]
    present = list(pp["pp_subtype"].value_counts().index)
    ordered = [s for s in preferred if s in present]
    ordered += [s for s in present if s not in ordered]
    return ordered


def generate_tables(pp: pd.DataFrame, tables_dir: Path) -> dict[str, Path]:
    ensure_dirs(tables_dir)
    out: dict[str, Path] = {}

    pred_path = tables_dir / "pp_subtype_predictions.tsv"
    pp.to_csv(pred_path, sep="\t", index=False)
    out["predictions"] = pred_path

    counts = (
        pp["pp_subtype"]
        .value_counts()
        .rename_axis("pp_subtype")
        .reset_index(name="n")
    )
    counts["pct"] = (100 * counts["n"] / counts["n"].sum()).round(2)
    p = tables_dir / "pp_subtype_counts.csv"
    counts.to_csv(p, index=False)
    out["counts"] = p

    conf = (
        pp.groupby("pp_subtype", as_index=False)
        .agg(
            n=("sequence_id", "count"),
            mean_mge_confidence=("mge_confidence", "mean"),
            median_mge_confidence=("mge_confidence", "median"),
            mean_pp_confidence=("pp_subtype_confidence", "mean"),
            median_pp_confidence=("pp_subtype_confidence", "median"),
            median_length=("length", "median"),
            mean_length=("length", "mean"),
        )
        .sort_values("n", ascending=False)
    )
    p = tables_dir / "pp_subtype_summary.csv"
    conf.to_csv(p, index=False)
    out["summary"] = p
    return out


def generate_figures(pp: pd.DataFrame, figures_dir: Path) -> list[Path]:
    ensure_dirs(figures_dir)
    apply_mpl_style()
    written: list[Path] = []
    order = _subtype_order(pp)
    counts = pp["pp_subtype"].value_counts().reindex(order).fillna(0)

    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    x = np.arange(len(order))
    ax.bar(x, counts.values, color=[subtype_color(s) for s in order], edgecolor="white", width=0.72)
    for i, v in enumerate(counts.values):
        ax.text(i, v + max(counts.max(), 1) * 0.02, f"{int(v)}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Sequences")
    ax.set_title(f"PP subtypes among phage-plasmid calls (n={len(pp)})", loc="left", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    stem = figures_dir / "01_pp_subtype_counts"
    _save(fig, stem)
    written.append(Path(f"{stem}.png"))

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.9), constrained_layout=True)
    panels = [
        ("mge_confidence", "A  3-class confidence", False),
        ("pp_subtype_confidence", "B  Subtype confidence", False),
        ("length", "C  Length", True),
    ]
    for ax, (col, title, logy) in zip(axes, panels):
        data, labels, colors = [], [], []
        for s in order:
            vals = pp.loc[pp["pp_subtype"] == s, col].dropna().values
            if len(vals) == 0:
                continue
            data.append(vals)
            labels.append(f"{s}\n(n={len(vals)})")
            colors.append(subtype_color(s))
        if not data:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title, loc="left", fontweight="bold")
            continue
        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, showfliers=False)
        for patch, colr in zip(bp["boxes"], colors):
            patch.set_facecolor(colr)
            patch.set_alpha(0.8)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.tick_params(axis="x", labelsize=7)
        if logy:
            ax.set_yscale("log")
            ax.set_ylabel("bp (log)")
        else:
            ax.set_ylim(0, 1.05)
            ax.set_ylabel("confidence")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    stem = figures_dir / "02_pp_subtype_confidence_length"
    _save(fig, stem)
    written.append(Path(f"{stem}.png"))

    fig, ax = plt.subplots(figsize=(5.5, 5.0), constrained_layout=True)
    sizes = [int(counts[s]) for s in order]
    ax.pie(
        sizes,
        labels=order,
        colors=[subtype_color(s) for s in order],
        autopct=lambda p: f"{p:.1f}%" if p > 0 else "",
        startangle=90,
        wedgeprops=dict(width=0.42, edgecolor="white"),
        pctdistance=0.75,
        textprops={"fontsize": 8},
    )
    ax.set_title(f"PP subtype composition (n={len(pp)})", fontweight="bold")
    stem = figures_dir / "03_pp_subtype_composition"
    _save(fig, stem)
    written.append(Path(f"{stem}.png"))
    return written


def run_step2(df: pd.DataFrame, step_dir: Path) -> dict:
    tables_dir = step_dir / "tables"
    figures_dir = step_dir / "figures"
    ensure_dirs(step_dir, tables_dir, figures_dir)

    pp = filter_pp(df)
    if pp.empty:
        summary = {
            "n_phage_plasmid": 0,
            "warning": "No phage_plasmid predictions; step 2 produced empty tables only.",
            "tables": {},
            "figures": [],
        }
        empty = tables_dir / "pp_subtype_predictions.tsv"
        pp.to_csv(empty, sep="\t", index=False)
        (step_dir / "step2_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        return summary

    tables = generate_tables(pp, tables_dir)
    figures = generate_figures(pp, figures_dir)
    summary = {
        "n_phage_plasmid": int(len(pp)),
        "subtype_counts": pp["pp_subtype"].value_counts().to_dict(),
        "tables": {k: str(v) for k, v in tables.items()},
        "figures": [str(f) for f in figures],
    }
    (step_dir / "step2_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
