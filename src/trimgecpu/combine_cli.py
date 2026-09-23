"""CLI: combine multiple FASTA files."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from .combine import combine_fastas, discover_fastas
from .constants import MIN_SEQUENCE_LENGTH


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="trimge-cpu-combine",
        description="Combine multiple FASTA files into one (optional filename prefixes).",
    )
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("-i", "--input", nargs="+", help="One or more FASTA files")
    g.add_argument("--input-dir", help="Directory of FASTA files (*.fna, *.fa, *.fasta)")
    p.add_argument("-o", "--output", required=True, help="Combined FASTA output path")
    p.add_argument("--recursive", action="store_true", help="With --input-dir, search recursively")
    p.add_argument(
        "--no-prefix",
        action="store_true",
        help="Do not prefix IDs with source filename (may collide)",
    )
    p.add_argument(
        "--separator",
        default="|",
        help="Separator between filename prefix and original ID (default: |)",
    )
    p.add_argument(
        "--min-length",
        type=int,
        default=MIN_SEQUENCE_LENGTH,
        help=f"Skip sequences shorter than this after cleaning (default: {MIN_SEQUENCE_LENGTH})",
    )
    p.add_argument(
        "--manifest",
        default=None,
        help="Manifest TSV path (default: <output>.manifest.tsv)",
    )
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    log = logging.getLogger("trimgecpu.combine")

    if args.input_dir:
        inputs = discover_fastas(Path(args.input_dir), recursive=args.recursive)
        if not inputs:
            log.error("No FASTA files found in %s", args.input_dir)
            return 1
    else:
        inputs = [Path(x) for x in args.input]
        missing = [str(p) for p in inputs if not p.exists()]
        if missing:
            log.error("Missing files: %s", ", ".join(missing))
            return 1

    log.info("Combining %d FASTA file(s) → %s", len(inputs), args.output)
    result = combine_fastas(
        inputs,
        Path(args.output),
        prefix_with_filename=not args.no_prefix,
        separator=args.separator,
        min_length=args.min_length,
        manifest_path=Path(args.manifest) if args.manifest else None,
    )
    log.info("%s", json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
