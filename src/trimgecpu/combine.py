"""Combine multiple FASTA files into one, with collision-safe IDs."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from .constants import MIN_SEQUENCE_LENGTH
from .fasta import clean_sequence

FASTA_GLOBS = ("*.fna", "*.fa", "*.fasta", "*.fas")


def _safe_prefix(path: Path) -> str:
    """Filesystem-safe prefix from filename (no extension)."""
    stem = path.stem
    stem = re.sub(r"[^\w.\-]+", "_", stem)
    return stem or "file"


def discover_fastas(input_dir: Path, recursive: bool = False) -> list[Path]:
    files: list[Path] = []
    patterns = FASTA_GLOBS
    if recursive:
        for pat in patterns:
            files.extend(input_dir.rglob(pat))
    else:
        for pat in patterns:
            files.extend(input_dir.glob(pat))
    # unique, sorted
    uniq = sorted({p.resolve() for p in files if p.is_file()})
    return [Path(p) for p in uniq]


def combine_fastas(
    inputs: list[Path],
    output_fasta: Path,
    *,
    prefix_with_filename: bool = True,
    separator: str = "|",
    min_length: int = MIN_SEQUENCE_LENGTH,
    manifest_path: Path | None = None,
) -> dict:
    """
    Concatenate FASTA files into ``output_fasta``.

    When ``prefix_with_filename`` is True, headers become
    ``{filename_stem}{separator}{original_id}`` to avoid collisions.
    """
    if not inputs:
        raise ValueError("No FASTA input files provided")

    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path is None:
        manifest_path = output_fasta.with_suffix(output_fasta.suffix + ".manifest.tsv")

    n_in = 0
    n_written = 0
    n_skipped_short = 0
    rows: list[dict] = []

    with output_fasta.open("w") as out, manifest_path.open("w", newline="") as mf:
        writer = csv.DictWriter(
            mf,
            fieldnames=[
                "source_file",
                "source_prefix",
                "original_id",
                "combined_id",
                "length",
                "kept",
            ],
            delimiter="\t",
        )
        writer.writeheader()

        for fp in inputs:
            prefix = _safe_prefix(fp)
            current_id: str | None = None
            chunks: list[str] = []

            def flush() -> None:
                nonlocal n_in, n_written, n_skipped_short, current_id, chunks
                if current_id is None:
                    return
                n_in += 1
                seq = clean_sequence("".join(chunks))
                combined = (
                    f"{prefix}{separator}{current_id}" if prefix_with_filename else current_id
                )
                kept = len(seq) >= min_length
                writer.writerow(
                    {
                        "source_file": str(fp),
                        "source_prefix": prefix if prefix_with_filename else "",
                        "original_id": current_id,
                        "combined_id": combined,
                        "length": len(seq),
                        "kept": kept,
                    }
                )
                if kept:
                    out.write(f">{combined}\n")
                    for i in range(0, len(seq), 80):
                        out.write(seq[i : i + 80] + "\n")
                    n_written += 1
                else:
                    n_skipped_short += 1
                current_id = None
                chunks = []

            with fp.open() as handle:
                for raw in handle:
                    line = raw.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        flush()
                        current_id = line[1:].split()[0]
                        chunks = []
                    else:
                        chunks.append(line)
                flush()

    return {
        "n_input_files": len(inputs),
        "n_sequences_seen": n_in,
        "n_sequences_written": n_written,
        "n_skipped_short": n_skipped_short,
        "output_fasta": str(output_fasta.resolve()),
        "manifest": str(manifest_path.resolve()),
        "prefix_with_filename": prefix_with_filename,
        "min_length": min_length,
    }
