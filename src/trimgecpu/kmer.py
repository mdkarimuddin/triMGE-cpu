"""K-mer frequency feature extraction."""

from __future__ import annotations

from itertools import product

import numpy as np


def kmer_vocabulary(k: int) -> tuple[list[str], dict[str, int]]:
    """Return ordered k-mer list and lookup table."""
    kmers = ["".join(chars) for chars in product("ACGT", repeat=k)]
    return kmers, {km: idx for idx, km in enumerate(kmers)}


def compute_kmer_frequencies(sequence: str, k: int = 5) -> np.ndarray:
    """
    Compute normalized k-mer frequency vector.

  Matches `compute_kmer_frequencies()` in the unified pipeline.
    """
    _, kmer_to_idx = kmer_vocabulary(k)
    counts = np.zeros(len(kmer_to_idx), dtype=np.float64)

    if len(sequence) < k:
        return counts

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        idx = kmer_to_idx.get(kmer)
        if idx is not None:
            counts[idx] += 1

    total = counts.sum()
    if total > 0:
        counts /= total
    return counts


def batch_kmer_features(sequences: list[str], k: int = 5) -> np.ndarray:
    """Build a (n_sequences, 4^k) feature matrix."""
    return np.vstack([compute_kmer_frequencies(seq, k=k) for seq in sequences])
