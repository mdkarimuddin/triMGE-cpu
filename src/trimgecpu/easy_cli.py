"""CLI for the easy end-to-end CPU pipeline."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from .constants import MIN_SEQUENCE_LENGTH
from .easy_pipeline import run_easy_pipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="trimge-cpu-pipeline",
        description=(
            "Easy triMGE-cpu pipeline: "
            "combine FASTAs → classify (k-mer SVM, CPU-only) → tables/figures."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--input-dir",
        help="Directory containing FASTA files to combine (*.fna/*.fa/*.fasta)",
    )
    src.add_argument(
        "-i",
        "--input",
        nargs="+",
        help="One or more FASTA files to combine",
    )
    src.add_argument(
        "--fasta",
        help="Already-combined single FASTA",
    )
    p.add_argument(
        "-w",
        "--workdir",
        required=True,
        help="Output working directory (combined + predictions + report)",
    )
    p.add_argument(
        "--models-dir",
        default=None,
        help="Directory with exported triMGE-cpu models (default: ./models)",
    )
    p.add_argument("--recursive", action="store_true", help="Recurse under --input-dir")
    p.add_argument(
        "--no-prefix",
        action="store_true",
        help="Do not prefix sequence IDs with source filename when combining",
    )
    p.add_argument("--min-length", type=int, default=MIN_SEQUENCE_LENGTH)
    p.add_argument(
        "--skip-classify",
        action="store_true",
        help="Reuse workdir/02_predictions/predictions.tsv",
    )
    p.add_argument(
        "--skip-report",
        action="store_true",
        help="Stop after classification (no figures)",
    )
    p.add_argument(
        "--report-steps",
        default="1,2",
        help="Report steps to run (1=3-class, 2=PP subtype)",
    )
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    steps = tuple(int(x.strip()) for x in args.report_steps.split(",") if x.strip())
    models_dir = Path(args.models_dir) if args.models_dir else None
    try:
        summary = run_easy_pipeline(
            workdir=Path(args.workdir),
            input_files=[Path(p) for p in args.input] if args.input else None,
            input_dir=Path(args.input_dir) if args.input_dir else None,
            recursive=args.recursive,
            single_fasta=Path(args.fasta) if args.fasta else None,
            models_dir=models_dir,
            prefix_with_filename=not args.no_prefix,
            min_length=args.min_length,
            skip_classify=args.skip_classify,
            skip_report=args.skip_report,
            steps_report=steps,
            verbose=args.verbose,
        )
    except Exception as exc:
        logging.getLogger(__name__).error("%s", exc)
        return 1
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
