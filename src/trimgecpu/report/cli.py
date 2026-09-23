"""CLI for the triMGE report pipeline."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .pipeline import run_report_pipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="trimge-cpu-report",
        description=(
            "triMGE-cpu reporting pipeline: "
            "Step 1 = 3-class tables/figures; "
            "Step 2 = phage-plasmid subtype tables/figures."
        ),
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="triMGE-cpu predictions TSV (output of `trimge-cpu`)",
    )
    p.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory for report tables and figures",
    )
    p.add_argument(
        "--steps",
        default="1,2",
        help="Comma-separated steps to run (default: 1,2). Example: --steps 1",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose logging",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    log = logging.getLogger("trimgecpu.report")

    steps = tuple(int(x.strip()) for x in args.steps.split(",") if x.strip())
    for s in steps:
        if s not in (1, 2):
            log.error("Invalid step %s (allowed: 1, 2)", s)
            return 2

    log.info("Running report pipeline steps=%s", steps)
    result = run_report_pipeline(
        predictions=Path(args.input),
        output_dir=Path(args.output_dir),
        steps=steps,
    )
    log.info("Done. n=%s → %s", result["n_sequences"], result["output_dir"])
    if "step1" in result:
        log.info("Step 1 classes: %s", result["step1"].get("class_counts"))
    if "step2" in result:
        log.info(
            "Step 2 PP n=%s subtypes=%s",
            result["step2"].get("n_phage_plasmid"),
            result["step2"].get("subtype_counts"),
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
