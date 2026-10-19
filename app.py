"""Public Streamlit demo for the triMGE-cpu classifiers."""

from __future__ import annotations

import csv
import io
import sys
from pathlib import Path

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


st.title("triMGE-cpu")
st.markdown(
    "Classify a nucleotide sequence as **phage**, **phage-plasmid**, or **plasmid**. "
    "A sequence called phage-plasmid also receives one of 15 phage-plasmid subtypes."
)
st.caption(
    "The app runs the 5-mer SVM weights shipped with triMGE-cpu. "
    "On the manuscript test set those weights reach 91.9% accuracy. "
    "Sequences shorter than 100 bp are skipped. The same classifier is strongest "
    "on contigs of about 10 kb or longer."
)

source = st.radio(
    "Input",
    ["Paste FASTA", "Upload a file", "Example sequences"],
    horizontal=True,
)

fasta_text = ""
if source == "Paste FASTA":
    fasta_text = st.text_area(
        "FASTA",
        height=220,
        placeholder=">sequence_1\nACGTACGT...",
        label_visibility="collapsed",
    )
elif source == "Upload a file":
    uploaded = st.file_uploader("FASTA file", type=["fa", "fna", "fasta", "txt"])
    if uploaded is not None:
        fasta_text = uploaded.getvalue().decode("utf-8", errors="replace")
else:
    fasta_text = EXAMPLE_FASTA.read_text()
    st.caption(f"Using {EXAMPLE_FASTA.name}.")

classify = st.button("Classify", type="primary", disabled=not fasta_text.strip())

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
            if skipped:
                st.info(f"Skipped {skipped} sequence(s) shorter than 100 bp.")
            phage_plasmids = sum(1 for pred in predictions if pred.pp_subtype)
            left, middle, right = st.columns(3)
            left.metric("Sequences", len(predictions))
            middle.metric("Phage-plasmid calls", phage_plasmids)
            right.metric("Bases", f"{total_bases:,}")
            st.dataframe(result_rows(predictions), use_container_width=True, hide_index=True)
            st.download_button(
                "Download TSV",
                data=predictions_to_tsv(predictions),
                file_name="trimge_cpu_predictions.tsv",
                mime="text/tab-separated-values",
            )

st.divider()
st.caption(
    "triMGE-cpu is a research classifier (MIT License, Md Karim Uddin). "
    "A call is a sequence label, not a biological confirmation."
)
