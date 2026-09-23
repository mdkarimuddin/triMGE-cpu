"""triMGE-cpu — CPU-only three-way MGE classification (k-mer SVM)."""

__version__ = "0.1.0"

__all__ = ["TriMGECPUClassifier", "SequencePrediction", "__version__"]


def __getattr__(name: str):
    if name == "TriMGECPUClassifier":
        from .pipeline import TriMGECPUClassifier

        return TriMGECPUClassifier
    if name == "SequencePrediction":
        from .pipeline import SequencePrediction

        return SequencePrediction
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
