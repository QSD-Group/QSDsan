from pathlib import Path

import yaml

EVAL = Path(__file__).resolve().parents[1] / "eval" / "questions.yaml"


def test_questions_file_is_well_formed():
    items = yaml.safe_load(EVAL.read_text())
    assert len(items) >= 4
    for item in items:
        assert item["question"]
        assert item["expect"] in {"cite_substring", "refuse"}
        if item["expect"] == "cite_substring":
            assert item["value"]


def test_seeded_cases_present():
    items = yaml.safe_load(EVAL.read_text())
    questions = " ".join(i["question"].lower() for i in items)
    assert "wastestream" in questions
    assert "pitfalls" in questions
    assert "bsm1" in questions
    assert "scada" in questions or "modbus" in questions


def test_has_code_adapter_regression_case():
    # A question that only the code adapter can answer: its expected citation is
    # a GitHub blob URL into QSDsan source (doctests/API), not a readthedocs page.
    items = yaml.safe_load(EVAL.read_text())
    code_values = " ".join(
        i.get("value", "") for i in items if i["expect"] == "cite_substring"
    )
    assert "github.com/QSD-Group/QSDsan/blob" in code_values
