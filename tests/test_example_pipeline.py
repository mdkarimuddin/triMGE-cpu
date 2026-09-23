"""One-command check that the shipped models write both report steps."""

from pathlib import Path

import pytest

from trimgecpu.easy_pipeline import run_easy_pipeline

REPO = Path(__file__).resolve().parents[1]
FASTA = REPO / "tests" / "fixtures" / "example_mge.fasta"


def test_example_fasta_writes_both_report_steps(tmp_path):
    pytest.importorskip("matplotlib")
    pytest.importorskip("sklearn")
    work = tmp_path / "run"
    run_easy_pipeline(
        workdir=work,
        single_fasta=FASTA,
        models_dir=REPO / "models",
        min_length=100,
    )
    assert (work / "02_predictions" / "predictions.tsv").is_file()
    assert (work / "03_report" / "01_three_class" / "tables" / "three_class_counts.csv").is_file()
    assert list((work / "03_report" / "01_three_class" / "figures").glob("*.pdf"))
    assert (work / "03_report" / "02_pp_subtype").is_dir()
    step2_tables = work / "03_report" / "02_pp_subtype" / "tables"
    assert step2_tables.is_dir()
    assert list(step2_tables.glob("*"))
