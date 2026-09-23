"""Command-line interface for triMGE-cpu."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from pathlib import Path

from .fasta import parse_fasta
from .pipeline import TriMGECPUClassifier

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="trimge-cpu",
        description=(
            "triMGE-cpu: CPU-only classification of mobile genetic elements. "
            "Stage 1: phage / phage-plasmid / plasmid (k-mer k5 + SVM). "
            "Stage 2: P-P subtype for predicted phage-plasmids (k-mer k5 + SVM)."
        ),
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV or JSON file")
    parser.add_argument(
        "--format",
        choices=["tsv", "json"],
        default="tsv",
        help="Output format (default: tsv)",
    )
    parser.add_argument(
        "--models-dir",
        default=None,
        help="Directory containing exported model artifacts (default: ./models)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    return parser


def write_tsv(predictions, output_path: Path) -> None:
    fieldnames = [
        "sequence_id",
        "length",
        "mge_class",
        "mge_class_display",
        "mge_confidence",
        "pp_subtype",
        "pp_subtype_confidence",
    ]
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
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


def write_json(predictions, output_path: Path) -> None:
    payload = []
    for pred in predictions:
        item = {
            "sequence_id": pred.sequence_id,
            "length": pred.length,
            "mge_class": pred.mge_class,
            "mge_class_display": pred.mge_class_display,
            "mge_confidence": pred.mge_confidence,
            "mge_probabilities": pred.mge_probabilities,
            "pp_subtype": pred.pp_subtype,
            "pp_subtype_confidence": pred.pp_subtype_confidence,
            "pp_subtype_probabilities": pred.pp_subtype_probabilities,
        }
        payload.append(item)
    output_path.write_text(json.dumps(payload, indent=2))


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.exists():
        logger.error("Input FASTA not found: %s", input_path)
        return 1

    records = parse_fasta(input_path)
    if not records:
        logger.error("No sequences retained after length filtering")
        return 1

    logger.info("Loaded %s sequences from %s", len(records), input_path)
    classifier = TriMGECPUClassifier(models_dir=args.models_dir)
    predictions = classifier.predict_records(records)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if args.format == "json" or output_path.suffix.lower() == ".json":
        write_json(predictions, output_path)
    else:
        write_tsv(predictions, output_path)

    n_pp = sum(1 for p in predictions if p.pp_subtype)
    logger.info(
        "Wrote %s predictions (%s phage-plasmids typed) → %s",
        len(predictions),
        n_pp,
        output_path,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
