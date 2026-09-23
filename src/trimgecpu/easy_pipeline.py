"""Easy end-to-end CPU pipeline: combine FASTAs → classify → report."""

from __future__ import annotations

import json
import logging
from pathlib import Path

from .combine import combine_fastas, discover_fastas
from .constants import MIN_SEQUENCE_LENGTH
from .report.pipeline import run_report_pipeline

logger = logging.getLogger(__name__)


def _run_classify(
    input_fasta: Path,
    output_tsv: Path,
    *,
    models_dir: Path | None,
    verbose: bool,
) -> None:
    from .cli import main as classify_main

    argv = ["-i", str(input_fasta), "-o", str(output_tsv)]
    if models_dir is not None:
        argv.extend(["--models-dir", str(models_dir)])
    if verbose:
        argv.append("--verbose")
    rc = classify_main(argv)
    if rc not in (0, None):
        raise RuntimeError(f"trimge-cpu classify failed with exit code {rc}")


def run_easy_pipeline(
    *,
    workdir: Path,
    input_files: list[Path] | None = None,
    input_dir: Path | None = None,
    recursive: bool = False,
    single_fasta: Path | None = None,
    models_dir: Path | None = None,
    prefix_with_filename: bool = True,
    min_length: int = MIN_SEQUENCE_LENGTH,
    skip_classify: bool = False,
    skip_report: bool = False,
    steps_report: tuple[int, ...] = (1, 2),
    verbose: bool = False,
) -> dict:
    """
    End-to-end CPU user pipeline.

    Layout under ``workdir``::

        01_combined/all_sequences.fna (+ manifest)
        02_predictions/predictions.tsv
        03_report/
        pipeline_summary.json
    """
    workdir = Path(workdir)
    combined_dir = workdir / "01_combined"
    pred_dir = workdir / "02_predictions"
    report_dir = workdir / "03_report"
    for d in (workdir, combined_dir, pred_dir, report_dir):
        d.mkdir(parents=True, exist_ok=True)

    summary: dict = {"workdir": str(workdir.resolve()), "tool": "triMGE-cpu"}

    combined_fasta = combined_dir / "all_sequences.fna"
    if single_fasta is not None:
        single_fasta = Path(single_fasta)
        if not single_fasta.exists():
            raise FileNotFoundError(single_fasta)
        if single_fasta.resolve() != combined_fasta.resolve():
            combined_fasta.write_bytes(single_fasta.read_bytes())
        summary["combine"] = {
            "mode": "single_fasta",
            "source": str(single_fasta.resolve()),
            "output_fasta": str(combined_fasta.resolve()),
        }
        logger.info("Using single FASTA → %s", combined_fasta)
    else:
        if input_dir is not None:
            files = discover_fastas(Path(input_dir), recursive=recursive)
        else:
            files = [Path(p) for p in (input_files or [])]
        if not files:
            raise FileNotFoundError("No FASTA files found to combine")
        manifest = combine_fastas(
            files,
            combined_fasta,
            prefix_with_filename=prefix_with_filename,
            min_length=min_length,
        )
        summary["combine"] = {
            "mode": "combine",
            "n_files": len(files),
            "output_fasta": str(combined_fasta.resolve()),
            "manifest_rows": len(manifest),
        }
        logger.info("Combined %s FASTAs → %s", len(files), combined_fasta)

    pred_tsv = pred_dir / "predictions.tsv"
    if skip_classify and pred_tsv.exists():
        logger.info("Reusing predictions: %s", pred_tsv)
        summary["classify"] = {"skipped": True, "predictions": str(pred_tsv.resolve())}
    else:
        logger.info("Classifying with triMGE-cpu (k-mer SVM)…")
        _run_classify(
            combined_fasta,
            pred_tsv,
            models_dir=models_dir,
            verbose=verbose,
        )
        summary["classify"] = {
            "skipped": False,
            "predictions": str(pred_tsv.resolve()),
            "models_dir": str(models_dir) if models_dir else "default",
        }

    if skip_report:
        summary["report"] = {"skipped": True}
    else:
        logger.info("Building report…")
        run_report_pipeline(
            predictions=pred_tsv,
            output_dir=report_dir,
            steps=steps_report,
        )
        summary["report"] = {
            "skipped": False,
            "output_dir": str(report_dir.resolve()),
            "steps": list(steps_report),
        }

    summary_path = workdir / "pipeline_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Pipeline summary → %s", summary_path)
    return summary
