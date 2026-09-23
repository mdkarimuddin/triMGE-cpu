"""Shared plotting style for report figures."""

from __future__ import annotations

import matplotlib as mpl

CLASS_ORDER = ["phage", "phage_plasmid", "plasmid"]
CLASS_DISPLAY = {
    "phage": "phage",
    "phage_plasmid": "phage-plasmid",
    "plasmid": "plasmid",
}
CLASS_COLORS = {
    "phage": "#2C7BB6",
    "phage_plasmid": "#FDAE61",
    "plasmid": "#1A9850",
}

# Common subtypes seen in practice; extras appended at runtime
SUBTYPE_COLORS = {
    "Singleton": "#7B7B7B",
    "N15": "#2C7BB6",
    "VP882": "#FDAE61",
    "IEBH-Waukesha": "#D7191C",
    "Unknown": "#BDBDBD",
}


def apply_mpl_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.linewidth": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.dpi": 200,
            "pdf.fonttype": 42,
        }
    )


def subtype_color(name: str) -> str:
    return SUBTYPE_COLORS.get(str(name), "#9E9E9E")
