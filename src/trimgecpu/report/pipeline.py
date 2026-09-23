"""Orchestrate Step 1 (3-class) then Step 2 (PP subtype) reporting."""

from __future__ import annotations

import json
from pathlib import Path

from .io_utils import ensure_dirs, load_predictions
from .step1_three_class import run_step1
from .step2_pp_subtype import run_step2


def run_report_pipeline(
    predictions: Path | str,
    output_dir: Path | str,
    steps: tuple[int, ...] = (1, 2),
) -> dict:
    """
    Run post-classification reporting.

    Parameters
    ----------
    predictions:
        triMGE predictions TSV (from `trimge -i ... -o ...`).
    output_dir:
        Root directory for report outputs.
    steps:
        1 = three-class tables/figures; 2 = PP-subtype tables/figures.
    """
    predictions = Path(predictions)
    output_dir = Path(output_dir)
    ensure_dirs(output_dir)

    df = load_predictions(predictions)
    result: dict = {
        "input_predictions": str(predictions.resolve()),
        "output_dir": str(output_dir.resolve()),
        "n_sequences": int(len(df)),
        "steps_requested": list(steps),
    }

    if 1 in steps:
        step1_dir = output_dir / "01_three_class"
        result["step1"] = run_step1(df, step1_dir)

    if 2 in steps:
        step2_dir = output_dir / "02_pp_subtype"
        result["step2"] = run_step2(df, step2_dir)

    (output_dir / "pipeline_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    readme = output_dir / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# triMGE-cpu report pipeline output",
                "",
                f"- Input: `{predictions}`",
                f"- Sequences: **{len(df)}**",
                "",
                "## Layout",
                "",
                "- `01_three_class/tables` — class counts, confidence/length summaries, prediction copy",
                "- `01_three_class/figures` — 3-class bar/box/donut plots",
                "- `02_pp_subtype/tables` — PP-only predictions + subtype summaries",
                "- `02_pp_subtype/figures` — subtype composition and distributions",
                "- `pipeline_summary.json` — machine-readable run summary",
                "",
                "Regenerate:",
                "```bash",
                f"trimge-cpu-report -i {predictions} -o {output_dir}",
                "```",
                "",
            ]
        )
    )
    return result
