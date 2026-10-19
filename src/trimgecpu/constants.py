"""Shared constants for CPU-only MGE classification."""

MGE_CLASS_NAMES = ["phage", "phage_plasmid", "plasmid"]
MGE_CLASS_LABELS = {
    0: "phage",
    1: "phage_plasmid",
    2: "plasmid",
}
MGE_DISPLAY_NAMES = {
    "phage": "Phage",
    "phage_plasmid": "Phage-Plasmid",
    "plasmid": "Plasmid",
}

# 15 P-P subtypes from the manuscript benchmark (alphabetical LabelEncoder order).
PP_SUBTYPE_CLASSES = [
    "AB",
    "Actinophage_A",
    "BigBertha",
    "IEBH-Waukesha",
    "N15",
    "P1",
    "PiSa",
    "SSU5_pHCM2",
    "Singleton",
    "VP882",
    "cp32",
    "p48",
    "pKpn",
    "pMT1",
    "pSLy3",
]

KMER_K = 5
MIN_SEQUENCE_LENGTH = 100

# Limits for one run of the public demo.
DEMO_MAX_SEQUENCES = 100
DEMO_MAX_BASES = 2_000_000
DEMO_MAX_CHARS = 5_000_000


def demo_rejection(n_chars: int, n_sequences: int, n_bases: int) -> str | None:
    """Return why a public-demo submission is too large, or None when it is allowed."""
    if n_chars > DEMO_MAX_CHARS:
        return "chars"
    if n_sequences > DEMO_MAX_SEQUENCES or n_bases > DEMO_MAX_BASES:
        return "batch"
    return None

# Best CPU 3-class model: k-mer k5 + SVM (91.4% held-out; manuscript baseline).
# Matches phage-plasmid-analysis/puhti_scripts/scripts/02_train_kmer_models.py
MGE_SVM_PARAMS = {
    "kernel": "rbf",
    "C": 1.0,
    "probability": True,
    "random_state": 42,
}

# Best P-P subtype model: k-mer k5 + SVM (92.2%).
PP_SVM_PARAMS = {
    "kernel": "rbf",
    "C": 10,
    "gamma": "scale",
    "probability": True,
    "random_state": 42,
}
