"""Pin the shipped models on the included example and the public-demo parser."""

from pathlib import Path

from trimgecpu.constants import DEMO_MAX_BASES, DEMO_MAX_CHARS, DEMO_MAX_SEQUENCES, demo_rejection
from trimgecpu.fasta import parse_fasta, parse_fasta_text
from trimgecpu.pipeline import TriMGECPUClassifier

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "tests" / "fixtures" / "example_mge.fasta"

# Labels and confidences from the shipped 5-mer SVMs. Confidence is the maximum
# class probability, rounded to 4 decimals, matching the TSV writer.
EXPECTED = [
    ("example_seq_a", 112, "plasmid", 0.7531, None),
    ("example_seq_b", 112, "plasmid", 0.7531, None),
    ("example_seq_c", 112, "plasmid", 0.7531, None),
    ("k141_3005", 31539, "phage_plasmid", 0.6749, "Singleton"),
    ("k141_603", 90738, "plasmid", 0.6051, None),
]


def test_example_predictions_stay_stable():
    classifier = TriMGECPUClassifier(REPO / "models")
    assert classifier.bundle.classifier_3class.probability
    predictions = classifier.predict_records(parse_fasta(EXAMPLE))
    observed = [
        (
            pred.sequence_id,
            pred.length,
            pred.mge_class,
            round(pred.mge_confidence, 4),
            pred.pp_subtype,
        )
        for pred in predictions
    ]
    assert observed == EXPECTED


def test_parse_fasta_text_accepts_a_bare_sequence_and_skips_short_ones():
    long_sequence = "ACGT" * 30
    text = f">kept extra text\n{long_sequence}\n>short\nACGT\n"
    records, skipped = parse_fasta_text(text)
    assert skipped == 1
    assert [record.sequence_id for record in records] == ["kept"]
    bare, bare_skipped = parse_fasta_text(long_sequence)
    assert bare_skipped == 0
    assert bare[0].sequence_id == "sequence_1"
    assert len(bare[0].sequence) == len(long_sequence)


def test_public_demo_limits():
    assert demo_rejection(DEMO_MAX_CHARS, DEMO_MAX_SEQUENCES, DEMO_MAX_BASES) is None
    assert demo_rejection(DEMO_MAX_CHARS + 1, 1, 100) == "chars"
    assert demo_rejection(100, DEMO_MAX_SEQUENCES + 1, 100) == "batch"
    assert demo_rejection(100, 1, DEMO_MAX_BASES + 1) == "batch"
