"""Public Streamlit demo for the triMGE-cpu classifiers."""

from __future__ import annotations

import csv
import html
import io
import sys
from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from trimgecpu.constants import (  # noqa: E402
    DEMO_MAX_BASES,
    DEMO_MAX_SEQUENCES,
    MIN_SEQUENCE_LENGTH,
    PP_SUBTYPE_CLASSES,
    demo_rejection,
)
from trimgecpu.fasta import parse_fasta_text  # noqa: E402
from trimgecpu.pipeline import TriMGECPUClassifier  # noqa: E402

EXAMPLE_FASTA = ROOT / "tests" / "fixtures" / "example_mge.fasta"
REPO_URL = "https://github.com/mdkarimuddin/triMGE-cpu"

st.set_page_config(
    page_title="triMGE-cpu · MGE classifier",
    page_icon=":material/genetics:",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# --------------------------------------------------------------------------
# Class palette. One place, used by chips, cards, table and charts.
# --------------------------------------------------------------------------
CLASSES = [
    # key, display, colour, tint
    ("phage", "Phage", "#1F4E79", "#E4EDF6"),
    ("phage_plasmid", "Phage-Plasmid", "#C47B2B", "#F8ECDD"),
    ("plasmid", "Plasmid", "#2F6F4E", "#E2F0E8"),
]
CLASS_ORDER = [display for _, display, _, _ in CLASSES]
CLASS_COLORS = [color for _, _, color, _ in CLASSES]
COLOR_BY_DISPLAY = {display: color for _, display, color, _ in CLASSES}
TINT_BY_DISPLAY = {display: tint for _, display, _, tint in CLASSES}


@st.cache_resource
def load_classifier() -> TriMGECPUClassifier:
    return TriMGECPUClassifier(models_dir=ROOT / "models")


# --------------------------------------------------------------------------
# Data helpers (unchanged behaviour)
# --------------------------------------------------------------------------
def predictions_to_tsv(predictions) -> str:
    buffer = io.StringIO()
    fieldnames = [
        "sequence_id",
        "length",
        "mge_class",
        "mge_class_display",
        "mge_confidence",
        "pp_subtype",
        "pp_subtype_confidence",
    ]
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for pred in predictions:
        writer.writerow(
            {
                "sequence_id": pred.sequence_id,
                "length": pred.length,
                "mge_class": pred.mge_class,
                "mge_class_display": pred.mge_class_display,
                "mge_confidence": f"{pred.mge_confidence:.4f}",
                "pp_subtype": pred.pp_subtype or "",
                "pp_subtype_confidence": (
                    f"{pred.pp_subtype_confidence:.4f}"
                    if pred.pp_subtype_confidence is not None
                    else ""
                ),
            }
        )
    return buffer.getvalue()


def call_probability(pred) -> float:
    """Probability of the class that was called (not the largest class probability)."""
    return pred.mge_probabilities[pred.mge_class]


def top_probability_class(pred) -> str:
    """Display name of the class with the largest probability."""
    top_key = max(pred.mge_probabilities, key=pred.mge_probabilities.get)
    return next(display for key, display, _, _ in CLASSES if key == top_key)


def label_differs_from_top(pred) -> bool:
    return top_probability_class(pred) != pred.mge_class_display


def result_frame(predictions) -> pd.DataFrame:
    rows = []
    for pred in predictions:
        rows.append(
            {
                "Sequence": pred.sequence_id,
                "Length (bp)": pred.length,
                "Class": pred.mge_class_display,
                "Call probability": call_probability(pred),
                "Highest-probability class": top_probability_class(pred),
                "Phage-plasmid subtype": pred.pp_subtype or "—",
                "Subtype confidence": pred.pp_subtype_confidence,
            }
        )
    return pd.DataFrame(rows)


def class_probability_frame(predictions) -> pd.DataFrame:
    rows = []
    for pred in predictions:
        for key, label, _, _ in CLASSES:
            rows.append(
                {
                    "Sequence": pred.sequence_id,
                    "Class": label,
                    "Probability": pred.mge_probabilities[key],
                }
            )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Charts
# --------------------------------------------------------------------------
def _axis_config(chart: alt.Chart) -> alt.Chart:
    return (
        chart.configure_axis(
            labelFont="IBM Plex Sans",
            titleFont="IBM Plex Sans",
            labelColor="#4A5560",
            titleColor="#4A5560",
            gridColor="#ECE8E0",
            domainColor="#D9D4C9",
            tickColor="#D9D4C9",
        )
        .configure_legend(labelFont="IBM Plex Sans", orient="top", direction="horizontal")
        .configure_view(stroke=None)
    )


def class_probability_chart(predictions) -> alt.Chart:
    frame = class_probability_frame(predictions)
    height = max(140, 34 * len(predictions))
    chart = (
        alt.Chart(frame)
        .mark_bar(height={"band": 0.72})
        .encode(
            y=alt.Y(
                "Sequence:N",
                sort=[pred.sequence_id for pred in predictions],
                title=None,
                axis=alt.Axis(labelFont="IBM Plex Mono", labelLimit=180),
            ),
            x=alt.X(
                "Probability:Q",
                stack="normalize",
                title="Class probability",
                axis=alt.Axis(format="%", tickCount=5),
            ),
            color=alt.Color(
                "Class:N",
                sort=CLASS_ORDER,
                scale=alt.Scale(domain=CLASS_ORDER, range=CLASS_COLORS),
                title=None,
            ),
            order=alt.Order("class_rank:Q"),
            tooltip=["Sequence", "Class", alt.Tooltip("Probability:Q", format=".1%")],
        )
        .transform_calculate(
            class_rank="indexof(['Phage','Phage-Plasmid','Plasmid'], datum.Class)"
        )
        .properties(height=height)
    )
    return _axis_config(chart)


def subtype_probability_chart(pred) -> alt.Chart:
    frame = pd.DataFrame(
        [
            {"Subtype": subtype, "Probability": probability}
            for subtype, probability in pred.pp_subtype_probabilities.items()
        ]
    ).sort_values("Probability", ascending=False)
    frame["Top"] = frame["Subtype"] == pred.pp_subtype
    chart = (
        alt.Chart(frame)
        .mark_bar(cornerRadiusEnd=3, height={"band": 0.7})
        .encode(
            y=alt.Y("Subtype:N", sort=frame["Subtype"].tolist(), title=None),
            x=alt.X(
                "Probability:Q",
                scale=alt.Scale(domain=[0, 1]),
                title="Subtype probability",
                axis=alt.Axis(format="%", tickCount=5),
            ),
            color=alt.condition(
                alt.datum.Top, alt.value("#C47B2B"), alt.value("#E6CFB2")
            ),
            tooltip=["Subtype", alt.Tooltip("Probability:Q", format=".1%")],
        )
        .properties(height=24 * len(frame))
    )
    return _axis_config(chart)


# --------------------------------------------------------------------------
# HTML building blocks
# --------------------------------------------------------------------------
def esc(value) -> str:
    return html.escape(str(value), quote=True)


def format_length(n: int) -> str:
    if n >= 1_000_000:
        return f"{n / 1_000_000:.2f} Mb"
    if n >= 10_000:
        return f"{n / 1_000:.1f} kb"
    return f"{n:,} bp"


def probability_bar(pred) -> str:
    segments = []
    for key, label, color, _ in CLASSES:
        p = pred.mge_probabilities[key]
        if p <= 0:
            continue
        text = f"{p:.0%}" if p >= 0.12 else ""
        segments.append(
            f'<span class="tm-seg" style="width:{p * 100:.2f}%;background:{color}" '
            f'title="{label}: {p:.1%}">{text}</span>'
        )
    return f'<div class="tm-bar">{"".join(segments)}</div>'


def result_card(pred) -> str:
    label = pred.mge_class_display
    color = COLOR_BY_DISPLAY[label]
    tint = TINT_BY_DISPLAY[label]
    subtype = ""
    if pred.pp_subtype:
        subtype = (
            '<div class="tm-subtype">'
            '<span class="tm-subtype-label">Subtype</span>'
            f'<span class="tm-subtype-value">{esc(pred.pp_subtype)}</span>'
            f'<span class="tm-subtype-conf">{pred.pp_subtype_confidence:.0%}</span>'
            "</div>"
        )
    short_warning = ""
    if pred.length < 10_000:
        short_warning = '<span class="tm-flag" title="Accuracy is highest on contigs of about 10 kb or longer">short contig</span>'
    differs_flag = ""
    if label_differs_from_top(pred):
        top = top_probability_class(pred)
        top_p = max(pred.mge_probabilities.values())
        differs_flag = (
            '<span class="tm-flag" title="The SVM label and its probability estimates disagree. '
            'Treat this call as uncertain.">'
            f"highest probability: {esc(top)} {top_p:.0%}</span>"
        )
    return f"""
    <div class="tm-card" style="--accent:{color};--tint:{tint}">
      <div class="tm-card-head">
        <div class="tm-id" title="{esc(pred.sequence_id)}">{esc(pred.sequence_id)}</div>
        <span class="tm-badge">{esc(label)}</span>
      </div>
      <div class="tm-meta">
        <span>{format_length(pred.length)}</span>
        <span class="tm-sep">·</span>
        <span>call probability <b>{call_probability(pred):.0%}</b></span>
        {short_warning}
        {differs_flag}
      </div>
      {probability_bar(pred)}
      {subtype}
    </div>
    """


def composition_strip(predictions) -> str:
    counts = {label: 0 for label in CLASS_ORDER}
    for pred in predictions:
        counts[pred.mge_class_display] += 1
    total = max(1, len(predictions))
    segs, legend = [], []
    for label in CLASS_ORDER:
        n = counts[label]
        color = COLOR_BY_DISPLAY[label]
        if n:
            segs.append(f'<span style="width:{n / total * 100:.2f}%;background:{color}"></span>')
        legend.append(
            f'<span class="tm-legend-item"><i style="background:{color}"></i>'
            f"{label} <b>{n}</b></span>"
        )
    return (
        f'<div class="tm-comp">{"".join(segs)}</div>'
        f'<div class="tm-legend">{"".join(legend)}</div>'
    )


def stat(label: str, value: str) -> str:
    return f'<div class="tm-stat"><div class="tm-stat-v">{value}</div><div class="tm-stat-l">{label}</div></div>'


# --------------------------------------------------------------------------
# Styles
# --------------------------------------------------------------------------
def inject_styles() -> None:
    st.html(
        """
        <style>
        :root {
          --tm-ink: #14202B;
          --tm-muted: #5B6670;
          --tm-line: #E3DFD6;
          --tm-paper: #FFFFFF;
          --tm-navy: #0F3B5F;
          --tm-mono: "IBM Plex Mono", ui-monospace, monospace;
        }
        .block-container { max-width: 1120px; padding-top: 2.2rem; padding-bottom: 3rem; }
        [data-testid="stHeader"] { background: transparent; }
        footer { visibility: hidden; }

        /* ---------- hero ---------- */
        .tm-hero {
          position: relative; overflow: hidden;
          border-radius: 18px; padding: 2.3rem 2.4rem 2rem;
          background: radial-gradient(120% 140% at 100% 0%, #1D5A86 0%, #0F3B5F 45%, #0A2740 100%);
          color: #F3F6F9; margin-bottom: 1.4rem;
        }
        .tm-hero::after {
          content: "ACGTTGCAATGCGTACCGTA GGCATTACGTAGCTTAGCCA TTAGCCGATACGGATCCGTA";
          position: absolute; right: -2rem; bottom: -0.8rem; width: 44%;
          font-family: var(--tm-mono); font-size: 2.1rem; line-height: 1.05;
          letter-spacing: 0.32em; color: rgba(255,255,255,0.045);
          word-break: break-all; pointer-events: none;
        }
        .tm-kicker {
          display: inline-flex; align-items: center; gap: .5rem;
          font-size: .74rem; font-weight: 600; letter-spacing: .14em; text-transform: uppercase;
          color: #9FC3E3; margin-bottom: .6rem;
        }
        .tm-kicker i { width: 1.6rem; height: 2px; background: #E0A15A; display: inline-block; }
        .tm-title {
          font-size: 3rem; font-weight: 700; letter-spacing: -0.035em; line-height: 1;
          margin: 0 0 .7rem 0; color: #FFFFFF;
        }
        .tm-title span { color: #E0A15A; font-weight: 500; }
        .tm-lede { max-width: 40rem; font-size: 1.08rem; line-height: 1.55; color: #D3DEE8; margin: 0; }
        .tm-stats { display: flex; flex-wrap: wrap; gap: 2.2rem; margin-top: 1.8rem; position: relative; z-index: 1; }
        .tm-stat-v { font-size: 1.55rem; font-weight: 600; color: #FFFFFF; letter-spacing: -0.02em; }
        .tm-stat-l { font-size: .78rem; color: #9FB4C8; margin-top: .1rem; }
        .tm-links { margin-top: 1.4rem; display: flex; gap: .6rem; flex-wrap: wrap; position: relative; z-index: 1; }
        .tm-links a {
          font-size: .84rem; font-weight: 500; color: #F3F6F9 !important; text-decoration: none;
          border: 1px solid rgba(255,255,255,.28); border-radius: 999px; padding: .32rem .85rem;
        }
        .tm-links a:hover { background: rgba(255,255,255,.1); }

        /* ---------- section headings ---------- */
        .tm-h { display: flex; align-items: baseline; gap: .7rem; margin: .2rem 0 .7rem; }
        .tm-h b { font-size: 1.05rem; font-weight: 600; color: var(--tm-ink); }
        .tm-h span { font-size: .85rem; color: var(--tm-muted); }
        .tm-step {
          font-family: var(--tm-mono); font-size: .72rem; color: var(--tm-navy);
          background: #E6EEF5; border-radius: 6px; padding: .1rem .42rem;
        }

        /* ---------- how-it-works panel ---------- */
        .tm-flow { display: flex; flex-direction: column; gap: .55rem; }
        .tm-flow-step { display: grid; grid-template-columns: 1.7rem 1fr; gap: .65rem; align-items: start; }
        .tm-flow-n {
          width: 1.7rem; height: 1.7rem; border-radius: 50%; background: #0F3B5F; color: #fff;
          font-size: .78rem; font-weight: 600; display: grid; place-items: center;
        }
        .tm-flow-t { font-weight: 600; font-size: .92rem; color: var(--tm-ink); }
        .tm-flow-d { font-size: .84rem; color: var(--tm-muted); line-height: 1.45; }
        .tm-chips { display: flex; flex-wrap: wrap; gap: .4rem; margin-top: .35rem; }
        .tm-chip {
          display: inline-flex; align-items: center; gap: .35rem; font-size: .78rem; font-weight: 500;
          padding: .16rem .55rem; border-radius: 999px; background: var(--tint); color: var(--accent);
        }
        .tm-chip i { width: .5rem; height: .5rem; border-radius: 50%; background: var(--accent); }
        .tm-note { font-size: .78rem; color: var(--tm-muted); border-top: 1px dashed var(--tm-line); padding-top: .65rem; margin-top: .3rem; }

        /* ---------- results ---------- */
        .tm-summary { display: grid; grid-template-columns: repeat(4, 1fr); gap: .8rem; margin-bottom: .9rem; }
        .tm-kpi { background: var(--tm-paper); border: 1px solid var(--tm-line); border-radius: 12px; padding: .85rem 1rem; }
        .tm-kpi-v { font-size: 1.5rem; font-weight: 600; letter-spacing: -0.02em; color: var(--tm-ink); }
        .tm-kpi-l { font-size: .78rem; color: var(--tm-muted); }
        .tm-comp { display: flex; height: 10px; border-radius: 999px; overflow: hidden; background: #ECE8E0; }
        .tm-comp span { display: block; height: 100%; }
        .tm-legend { display: flex; gap: 1.1rem; flex-wrap: wrap; margin: .5rem 0 1rem; font-size: .84rem; color: var(--tm-muted); }
        .tm-legend-item { display: inline-flex; align-items: center; gap: .4rem; }
        .tm-legend-item i { width: .6rem; height: .6rem; border-radius: 3px; display: inline-block; }
        .tm-legend-item b { color: var(--tm-ink); }

        .tm-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: .8rem; }
        .tm-card {
          background: var(--tm-paper); border: 1px solid var(--tm-line); border-left: 4px solid var(--accent);
          border-radius: 12px; padding: .9rem 1rem .95rem; display: flex; flex-direction: column; gap: .55rem;
        }
        .tm-card-head { display: flex; justify-content: space-between; align-items: center; gap: .6rem; }
        .tm-id { font-family: var(--tm-mono); font-size: .9rem; font-weight: 500; color: var(--tm-ink);
                 white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
        .tm-badge { flex: none; font-size: .76rem; font-weight: 600; padding: .18rem .6rem; border-radius: 999px;
                    background: var(--tint); color: var(--accent); }
        .tm-meta { font-size: .8rem; color: var(--tm-muted); display: flex; align-items: center; gap: .4rem; flex-wrap: wrap; }
        .tm-meta b { color: var(--tm-ink); font-weight: 600; }
        .tm-sep { color: #C9C3B6; }
        .tm-flag { font-size: .7rem; color: #8A5A1A; background: #FBF1E2; border-radius: 4px; padding: .05rem .4rem; }
        .tm-bar { display: flex; height: 18px; border-radius: 6px; overflow: hidden; background: #F0ECE4; }
        .tm-seg { display: flex; align-items: center; justify-content: center; height: 100%;
                  font-size: .66rem; font-weight: 600; color: rgba(255,255,255,.92); white-space: nowrap; }
        .tm-subtype { display: flex; align-items: center; gap: .5rem; font-size: .82rem;
                      background: #FBF5EC; border-radius: 8px; padding: .35rem .6rem; }
        .tm-subtype-label { color: #8A5A1A; font-size: .72rem; text-transform: uppercase; letter-spacing: .08em; font-weight: 600; }
        .tm-subtype-value { font-family: var(--tm-mono); font-weight: 500; color: var(--tm-ink); }
        .tm-subtype-conf { margin-left: auto; font-weight: 600; color: #C47B2B; }

        .tm-footer { display: flex; justify-content: space-between; flex-wrap: wrap; gap: .6rem;
                     font-size: .8rem; color: var(--tm-muted); border-top: 1px solid var(--tm-line);
                     padding-top: 1rem; margin-top: 2.2rem; }
        .tm-footer a { color: var(--tm-navy) !important; }

        /* Streamlit widget polish */
        [data-testid="stVerticalBlockBorderWrapper"] { background: #FFFFFF; }
        textarea { font-family: var(--tm-mono) !important; font-size: .82rem !important; }
        .stButton > button[kind="primary"], .stDownloadButton > button { font-weight: 600; }

        @media (max-width: 760px) {
          .tm-hero { padding: 1.6rem 1.3rem; }
          .tm-title { font-size: 2.2rem; }
          .tm-stats { gap: 1.3rem; }
          .tm-summary { grid-template-columns: repeat(2, 1fr); }
        }
        </style>
        """
    )


def section(step: str, title: str, hint: str = "") -> None:
    st.html(
        f'<div class="tm-h"><span class="tm-step">{step}</span><b>{title}</b>'
        f"<span>{hint}</span></div>"
    )


# --------------------------------------------------------------------------
# Page
# --------------------------------------------------------------------------
def render_hero() -> None:
    stats = "".join(
        [
            stat("MGE classes", "3"),
            stat("phage-plasmid subtypes", str(len(PP_SUBTYPE_CLASSES))),
            stat("3-class accuracy, 638 test sequences", "91.9%"),
            stat("features", "5-mer"),
        ]
    )
    st.html(
        f"""
        <div class="tm-hero">
          <div class="tm-kicker"><i></i>Mobile genetic element classifier</div>
          <h1 class="tm-title">triMGE<span>-cpu</span></h1>
          <p class="tm-lede">Tell phages, phage-plasmids and plasmids apart from sequence alone.
          Paste a contig. Get a class and its probabilities. Phage-plasmids also get a subtype.</p>
          <div class="tm-stats">{stats}</div>
          <div class="tm-links">
            <a href="{REPO_URL}" target="_blank">GitHub repository ↗</a>
            <a href="{REPO_URL}#readme" target="_blank">Run it locally ↗</a>
          </div>
        </div>
        """
    )


def render_how_it_works() -> None:
    chips = "".join(
        f'<span class="tm-chip" style="--accent:{c};--tint:{t}"><i></i>{d}</span>'
        for _, d, c, t in CLASSES
    )
    st.html(
        f"""
        <div class="tm-flow">
          <div class="tm-flow-step"><div class="tm-flow-n">1</div><div>
            <div class="tm-flow-t">Count 5-mers</div>
            <div class="tm-flow-d">Each sequence becomes a 1,024-dimension k-mer profile.</div>
          </div></div>
          <div class="tm-flow-step"><div class="tm-flow-n">2</div><div>
            <div class="tm-flow-t">Three-class SVM</div>
            <div class="tm-flow-d">Assigns one of three classes.</div>
            <div class="tm-chips">{chips}</div>
          </div></div>
          <div class="tm-flow-step"><div class="tm-flow-n">3</div><div>
            <div class="tm-flow-t">Subtype SVM</div>
            <div class="tm-flow-d">Phage-plasmid calls get one of {len(PP_SUBTYPE_CLASSES)} subtypes, such as P1, N15 or SSU5.</div>
          </div></div>
          <div class="tm-note">Accuracy: 91.9% for these shipped weights on 638 test sequences;
          the manuscript model, trained on the training split only, reaches 91.4%.
          Best on contigs ≥ 10 kb. Sequences under {MIN_SEQUENCE_LENGTH} bp are skipped.
          Demo limit: {DEMO_MAX_SEQUENCES} sequences or {DEMO_MAX_BASES / 1e6:.0f} Mb per run.</div>
        </div>
        """
    )


def render_input() -> tuple[str, bool]:
    source = st.segmented_control(
        "Input",
        ["Paste FASTA", "Upload file", "Try example"],
        default="Paste FASTA",
        label_visibility="collapsed",
    ) or "Paste FASTA"
    fasta_text = ""
    if source == "Paste FASTA":
        fasta_text = st.text_area(
            "FASTA",
            height=210,
            placeholder=">contig_1\nATGCGTACCGTAGGCATTACGTAGCTTAGCC...\n>contig_2\n...",
            label_visibility="collapsed",
        )
    elif source == "Upload file":
        uploaded = st.file_uploader(
            "FASTA file",
            type=["fa", "fna", "fasta", "txt"],
            label_visibility="collapsed",
            help="Up to 5 MB. .fa, .fna, .fasta or .txt",
        )
        if uploaded is not None:
            fasta_text = uploaded.getvalue().decode("utf-8", errors="replace")
    else:
        fasta_text = EXAMPLE_FASTA.read_text()
        n = fasta_text.count(">")
        st.code(
            "\n".join(
                line if line.startswith(">") else (line[:60] + "…" if len(line) > 60 else line)
                for line in fasta_text.splitlines()[:6]
            ),
            language=None,
        )
        st.caption(f"{n} sequences from the bundled example: three short synthetic repeats and two real contigs.")
    left, right = st.columns([1, 2], vertical_alignment="center")
    classify = left.button(
        "Classify sequences",
        type="primary",
        icon=":material/play_arrow:",
        disabled=not fasta_text.strip(),
        width="stretch",
    )
    if fasta_text.strip():
        right.caption(f"{fasta_text.count('>')} header(s) · {len(fasta_text):,} characters")
    return fasta_text, classify


def render_results(predictions, skipped: int, total_bases: int) -> None:
    st.html('<div style="height:.6rem"></div>')
    section("02", "Results", f"{len(predictions)} sequence(s) classified")
    if skipped:
        st.info(f"Skipped {skipped} sequence(s) shorter than {MIN_SEQUENCE_LENGTH} bp.", icon=":material/info:")

    phage_plasmids = sum(1 for pred in predictions if pred.pp_subtype)
    mean_conf = sum(call_probability(p) for p in predictions) / len(predictions)
    n_differ = sum(1 for p in predictions if label_differs_from_top(p))
    st.html(
        f"""
        <div class="tm-summary">
          <div class="tm-kpi"><div class="tm-kpi-v">{len(predictions)}</div><div class="tm-kpi-l">sequences</div></div>
          <div class="tm-kpi"><div class="tm-kpi-v">{format_length(total_bases)}</div><div class="tm-kpi-l">total length</div></div>
          <div class="tm-kpi"><div class="tm-kpi-v">{phage_plasmids}</div><div class="tm-kpi-l">phage-plasmid calls</div></div>
          <div class="tm-kpi"><div class="tm-kpi-v">{mean_conf:.0%}</div><div class="tm-kpi-l">mean call probability</div></div>
        </div>
        {composition_strip(predictions)}
        """
    )

    cards_tab, table_tab, chart_tab, subtype_tab = st.tabs(
        [
            ":material/view_agenda: Cards",
            ":material/table: Table",
            ":material/bar_chart: Probabilities",
            ":material/category: Subtypes",
        ]
    )
    with cards_tab:
        st.html(f'<div class="tm-grid">{"".join(result_card(p) for p in predictions)}</div>')
        caption = (
            "The label is the SVM decision. The bar shows its class probabilities, "
            "and the call probability is the share given to the called class. "
            "The largest segment can differ from the label."
        )
        if n_differ:
            caption += f" {n_differ} of {len(predictions)} call(s) here differ from the highest probability."
        st.caption(caption)
    with table_tab:
        table = result_frame(predictions)
        st.dataframe(
            table,
            width="stretch",
            hide_index=True,
            column_config={
                "Sequence": st.column_config.TextColumn(width="medium"),
                "Length (bp)": st.column_config.NumberColumn(format="localized"),
                "Class": st.column_config.TextColumn(),
                "Call probability": st.column_config.ProgressColumn(
                    min_value=0, max_value=1, format="percent"
                ),
                "Highest-probability class": st.column_config.TextColumn(),
                "Subtype confidence": st.column_config.ProgressColumn(
                    min_value=0, max_value=1, format="percent"
                ),
            },
        )
    with chart_tab:
        st.altair_chart(class_probability_chart(predictions), width="stretch")
    with subtype_tab:
        typed = [p for p in predictions if p.pp_subtype_probabilities]
        if not typed:
            st.info("No sequence was called phage-plasmid, so no subtype was assigned.", icon=":material/info:")
        else:
            chosen = typed[0]
            if len(typed) > 1:
                ids = [p.sequence_id for p in typed]
                pick = st.selectbox("Phage-plasmid sequence", ids)
                chosen = typed[ids.index(pick)]
            st.markdown(
                f"`{chosen.sequence_id}` → **{chosen.pp_subtype}** "
                f"({chosen.pp_subtype_confidence:.0%})"
            )
            st.altair_chart(subtype_probability_chart(chosen), width="stretch")

    st.download_button(
        "Download results (TSV)",
        data=predictions_to_tsv(predictions),
        file_name="trimge_cpu_predictions.tsv",
        mime="text/tab-separated-values",
        icon=":material/download:",
    )


def run_classification(fasta_text: str) -> None:
    if demo_rejection(len(fasta_text), 0, 0) == "chars":
        st.error("This paste is too large for the public demo. Use a FASTA under 5 million characters.")
        return
    records, skipped = parse_fasta_text(fasta_text)
    total_bases = sum(len(record.sequence) for record in records)
    rejection = demo_rejection(len(fasta_text), len(records), total_bases)
    if not records:
        st.warning(f"No sequence of at least {MIN_SEQUENCE_LENGTH} bp was found.")
    elif rejection == "batch":
        st.error(
            f"The public demo accepts up to {DEMO_MAX_SEQUENCES} sequences "
            f"and {DEMO_MAX_BASES:,} bp in one run. "
            f"This input has {len(records)} sequences and {total_bases:,} bp. "
            "For a larger set, run trimge-cpu locally."
        )
    else:
        with st.spinner("Counting k-mers and classifying…"):
            predictions = load_classifier().predict_records(records)
        st.session_state["results"] = {
            "predictions": predictions,
            "skipped": skipped,
            "total_bases": total_bases,
        }


inject_styles()
render_hero()

input_col, info_col = st.columns([1.7, 1], gap="medium")
with input_col:
    with st.container(border=True):
        section("01", "Input", "FASTA, one or more sequences")
        fasta_text, classify = render_input()
with info_col:
    with st.container(border=True):
        section("··", "How it works")
        render_how_it_works()

if classify:
    run_classification(fasta_text)

if "results" in st.session_state:
    saved = st.session_state["results"]
    render_results(saved["predictions"], saved["skipped"], saved["total_bases"])

st.html(
    f"""
    <div class="tm-footer">
      <span>triMGE-cpu · research classifier · MIT License · Md Karim Uddin</span>
      <span>A call is a sequence label, not a biological confirmation. ·
      <a href="{REPO_URL}" target="_blank">Source</a></span>
    </div>
    """
)
