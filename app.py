"""Public Streamlit demo for the triMGE-cpu classifiers."""

from __future__ import annotations

import csv
import io
import sys
from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from trimgecpu.constants import MIN_SEQUENCE_LENGTH
from trimgecpu.fasta import FastaRecord, clean_sequence
from trimgecpu.pipeline import TriMGECPUClassifier

MAX_SEQUENCES = 100
MAX_BASES = 2_000_000
MAX_CHARS = 5_000_000
EXAMPLE_FASTA = ROOT / "tests" / "fixtures" / "example_mge.fasta"

st.set_page_config(
    page_title="triMGE-cpu",
    page_icon=":material/genetics:",
    layout="wide",
)


@st.cache_resource
def load_classifier() -> TriMGECPUClassifier:
    return TriMGECPUClassifier(models_dir=ROOT / "models")


def records_from_text(text: str) -> tuple[list[FastaRecord], int]:
    """Parse pasted or uploaded FASTA. A bare sequence is accepted as one record."""
    records: list[FastaRecord] = []
    skipped = 0
    current_id: str | None = None
    chunks: list[str] = []
    bare_index = 0

    def flush() -> None:
        nonlocal skipped
        if current_id is None:
            return
        sequence = clean_sequence("".join(chunks))
        if len(sequence) >= MIN_SEQUENCE_LENGTH:
            records.append(FastaRecord(current_id, sequence))
        elif sequence:
            skipped += 1

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            current_id = line[1:].split()[0] or f"sequence_{len(records) + 1}"
            chunks = []
            continue
        if current_id is None:
            bare_index += 1
            current_id = f"sequence_{bare_index}"
        chunks.append(line)
    flush()
    return records, skipped


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


def result_rows(predictions) -> list[dict[str, str | int | float]]:
    rows = []
    for pred in predictions:
        rows.append(
            {
                "Sequence": pred.sequence_id,
                "Length (bp)": pred.length,
                "Class": pred.mge_class_display,
                "Confidence": round(pred.mge_confidence, 4),
                "Phage-plasmid subtype": pred.pp_subtype or "",
                "Subtype confidence": (
                    round(pred.pp_subtype_confidence, 4)
                    if pred.pp_subtype_confidence is not None
                    else ""
                ),
            }
        )
    return rows


CLASS_ORDER = ["Phage", "Phage-Plasmid", "Plasmid"]
CLASS_COLORS = ["#1F4E79", "#C47B2B", "#2F6F4E"]


def class_probability_frame(predictions) -> pd.DataFrame:
    rows = []
    for pred in predictions:
        for key, label in (
            ("phage", "Phage"),
            ("phage_plasmid", "Phage-Plasmid"),
            ("plasmid", "Plasmid"),
        ):
            rows.append(
                {
                    "Sequence": pred.sequence_id,
                    "Class": label,
                    "Probability": pred.mge_probabilities[key],
                }
            )
    return pd.DataFrame(rows)


def class_probability_chart(predictions) -> alt.Chart:
    frame = class_probability_frame(predictions)
    height = max(180, 52 * len(predictions))
    return (
        alt.Chart(frame)
        .mark_bar()
        .encode(
            y=alt.Y("Sequence:N", sort=[pred.sequence_id for pred in predictions], title=None),
            x=alt.X("Probability:Q", scale=alt.Scale(domain=[0, 1]), title="Probability"),
            color=alt.Color(
                "Class:N",
                sort=CLASS_ORDER,
                scale=alt.Scale(domain=CLASS_ORDER, range=CLASS_COLORS),
                title=None,
            ),
            yOffset=alt.YOffset("Class:N", sort=CLASS_ORDER),
            tooltip=[
                "Sequence",
                "Class",
                alt.Tooltip("Probability:Q", format=".3f"),
            ],
        )
        .properties(height=height)
    )


def subtype_probability_chart(predictions) -> alt.Chart | None:
    typed = [pred for pred in predictions if pred.pp_subtype_probabilities]
    if not typed:
        return None
    rows = []
    for pred in typed:
        for subtype, probability in pred.pp_subtype_probabilities.items():
            rows.append(
                {
                    "Sequence": pred.sequence_id,
                    "Subtype": subtype,
                    "Probability": probability,
                }
            )
    frame = pd.DataFrame(rows)
    subtype_order = (
        frame.groupby("Subtype")["Probability"].max().sort_values(ascending=False).index.tolist()
    )
    if len(typed) == 1:
        return (
            alt.Chart(frame)
            .mark_bar(color="#C47B2B")
            .encode(
                y=alt.Y("Subtype:N", sort=subtype_order, title=None),
                x=alt.X("Probability:Q", scale=alt.Scale(domain=[0, 1]), title="Probability"),
                tooltip=["Subtype", alt.Tooltip("Probability:Q", format=".3f")],
            )
            .properties(height=max(320, 22 * len(subtype_order)))
        )
    return (
        alt.Chart(frame)
        .mark_bar()
        .encode(
            y=alt.Y("Subtype:N", sort=subtype_order, title=None),
            x=alt.X("Probability:Q", scale=alt.Scale(domain=[0, 1]), title="Probability"),
            color=alt.Color("Sequence:N", title=None),
            yOffset="Sequence:N",
            tooltip=[
                "Sequence",
                "Subtype",
                alt.Tooltip("Probability:Q", format=".3f"),
            ],
        )
        .properties(height=max(360, 28 * len(subtype_order)))
    )


def paint_class(value: str) -> str:
    color = {
        "Phage": "#1F4E79",
        "Phage-Plasmid": "#C47B2B",
        "Plasmid": "#2F6F4E",
    }.get(value)
    if not color:
        return ""
    return f"color: {color}; font-weight: 600"


def render_results(predictions, skipped: int, total_bases: int) -> None:
    st.markdown("### Results")
    if skipped:
        st.info(f"Skipped {skipped} sequence(s) shorter than 100 bp.")
    phage_plasmids = sum(1 for pred in predictions if pred.pp_subtype)
    with st.container(border=True):
        left, middle, right = st.columns(3)
        left.metric("Sequences", len(predictions))
        middle.metric("Phage-plasmid calls", phage_plasmids)
        right.metric("Bases", f"{total_bases:,}")
        table = pd.DataFrame(result_rows(predictions))
        st.dataframe(
            table.style.map(paint_class, subset=["Class"]),
            use_container_width=True,
            hide_index=True,
        )
        st.download_button(
            "Download TSV",
            data=predictions_to_tsv(predictions),
            file_name="trimge_cpu_predictions.tsv",
            mime="text/tab-separated-values",
        )
    st.markdown("### Class probabilities")
    st.altair_chart(class_probability_chart(predictions), use_container_width=True)
    st.caption(
        "The label is the SVM decision. The bars are calibrated probabilities, "
        "so the tallest bar can differ from the call."
    )
    subtype_chart = subtype_probability_chart(predictions)
    if subtype_chart is not None:
        st.markdown("### Phage-plasmid subtype probabilities")
        st.altair_chart(subtype_chart, use_container_width=True)


def inject_styles() -> None:
    st.markdown(
        """
        <style>
        .block-container { max-width: 980px; padding-top: 1.4rem; }
        .trimge-kicker {
            margin: 0;
            color: #1F4E79;
            font-size: 0.78rem;
            font-weight: 600;
            letter-spacing: 0.08em;
            text-transform: uppercase;
        }
        .trimge-title {
            margin: 0.15rem 0 0.35rem 0;
            font-size: 2.4rem;
            line-height: 1.1;
            letter-spacing: -0.03em;
        }
        .trimge-lede { margin: 0 0 0.9rem 0; font-size: 1.05rem; color: #243040; }
        .trimge-legend { display: flex; flex-wrap: wrap; gap: 0.55rem; margin: 0 0 1.1rem 0; }
        .trimge-chip {
            display: inline-flex;
            align-items: center;
            gap: 0.4rem;
            padding: 0.28rem 0.65rem;
            border-radius: 999px;
            background: #FFFFFF;
            border: 1px solid #E2E6EC;
            font-size: 0.86rem;
        }
        .trimge-dot { width: 0.55rem; height: 0.55rem; border-radius: 50%; display: inline-block; }
        footer { visibility: hidden; }
        </style>
        """,
        unsafe_allow_html=True,
    )


inject_styles()
st.markdown(
    """
    <p class="trimge-kicker">Mobile genetic element classifier</p>
    <h1 class="trimge-title">triMGE-cpu</h1>
    <p class="trimge-lede">Classify a nucleotide sequence as phage, phage-plasmid, or plasmid.
    A phage-plasmid call also receives one of 15 subtypes.</p>
    <div class="trimge-legend">
      <span class="trimge-chip"><span class="trimge-dot" style="background:#1F4E79"></span>Phage</span>
      <span class="trimge-chip"><span class="trimge-dot" style="background:#C47B2B"></span>Phage-plasmid</span>
      <span class="trimge-chip"><span class="trimge-dot" style="background:#2F6F4E"></span>Plasmid</span>
    </div>
    """,
    unsafe_allow_html=True,
)

with st.container(border=True):
    source = st.radio(
        "Input",
        ["Paste FASTA", "Upload a file", "Example sequences"],
        horizontal=True,
    )
    fasta_text = ""
    if source == "Paste FASTA":
        fasta_text = st.text_area(
            "FASTA",
            height=200,
            placeholder=">sequence_1\nACGTACGT...",
            label_visibility="collapsed",
        )
    elif source == "Upload a file":
        uploaded = st.file_uploader("FASTA file", type=["fa", "fna", "fasta", "txt"])
        if uploaded is not None:
            fasta_text = uploaded.getvalue().decode("utf-8", errors="replace")
    else:
        fasta_text = EXAMPLE_FASTA.read_text()
        st.caption("Five sequences from the included example.")
    classify = st.button("Classify", type="primary", disabled=not fasta_text.strip())

with st.expander("About this classifier"):
    st.markdown(
        "The app runs the 5-mer SVM weights shipped with triMGE-cpu. "
        "On the manuscript test set those weights reach 91.9% accuracy. "
        "The manuscript reports 91.4% for the model fit on the training split only. "
        "Sequences shorter than 100 bp are skipped. The same classifier is strongest "
        "on contigs of about 10 kb or longer."
    )

if classify:
    if len(fasta_text) > MAX_CHARS:
        st.error("This paste is too large for the public demo. Use a FASTA under 5 million characters.")
    else:
        records, skipped = records_from_text(fasta_text)
        total_bases = sum(len(record.sequence) for record in records)
        if not records:
            st.warning("No sequence of at least 100 bp was found.")
        elif len(records) > MAX_SEQUENCES or total_bases > MAX_BASES:
            st.error(
                f"The public demo accepts up to {MAX_SEQUENCES} sequences "
                f"and {MAX_BASES:,} bp in one run. "
                f"This input has {len(records)} sequences and {total_bases:,} bp. "
                "For a larger set, run trimge-cpu locally."
            )
        else:
            with st.spinner("Classifying"):
                predictions = load_classifier().predict_records(records)
            st.session_state["results"] = {
                "predictions": predictions,
                "skipped": skipped,
                "total_bases": total_bases,
            }

if "results" in st.session_state:
    saved = st.session_state["results"]
    render_results(saved["predictions"], saved["skipped"], saved["total_bases"])

st.divider()
st.caption(
    "triMGE-cpu is a research classifier (MIT License, Md Karim Uddin). "
    "A call is a sequence label, not a biological confirmation."
)
