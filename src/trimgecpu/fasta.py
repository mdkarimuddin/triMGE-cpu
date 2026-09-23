"""FASTA parsing and sequence cleaning."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from .constants import MIN_SEQUENCE_LENGTH


@dataclass(frozen=True)
class FastaRecord:
    sequence_id: str
    sequence: str


def clean_sequence(sequence: str) -> str:
    """Keep only ACGT bases and uppercase."""
    return "".join(base for base in sequence.upper() if base in "ACGT")


def parse_fasta(path: str | Path) -> list[FastaRecord]:
    """Parse a FASTA file into cleaned records."""
    records: list[FastaRecord] = []
    current_id: str | None = None
    chunks: list[str] = []

    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seq = clean_sequence("".join(chunks))
                    if len(seq) >= MIN_SEQUENCE_LENGTH:
                        records.append(FastaRecord(current_id, seq))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if current_id is not None:
        seq = clean_sequence("".join(chunks))
        if len(seq) >= MIN_SEQUENCE_LENGTH:
            records.append(FastaRecord(current_id, seq))

    return records


def iter_fasta_records(path: str | Path) -> Iterator[FastaRecord]:
    """Yield cleaned FASTA records one at a time."""
    yield from parse_fasta(path)
