import runpy
import textwrap
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / ".github" / "scripts" / "sync_dependency_lower_bounds.py"


def _load_script():
    return runpy.run_path(str(SCRIPT))


def test_sync_bounds_updates_minor_and_major_changes(tmp_path):
    script = _load_script()
    pyproject = tmp_path / "pyproject.toml"
    requirements = tmp_path / "requirements-release.txt"
    pyproject.write_text(textwrap.dedent("""
    [project]
    dependencies = [
        "thermosteam>=0.53.4",
        "biosteam>=2.53.10",
        "SALib>=1.4.5",
    ]
    """).strip())
    requirements.write_text("thermosteam==0.54.0\nbiosteam==3.0.0\n")

    changed = script["sync_lower_bounds"](pyproject, requirements)

    assert changed is True
    text = pyproject.read_text()
    assert '"thermosteam>=0.54.0"' in text
    assert '"biosteam>=3.0.0"' in text
    assert '"SALib>=1.4.5"' in text


def test_sync_bounds_ignores_patch_changes(tmp_path):
    script = _load_script()
    pyproject = tmp_path / "pyproject.toml"
    requirements = tmp_path / "requirements-release.txt"
    original = textwrap.dedent("""
    [project]
    dependencies = [
        "thermosteam>=0.53.4",
        "biosteam>=2.53.10",
    ]
    """).strip()
    pyproject.write_text(original)
    requirements.write_text("thermosteam==0.53.5\nbiosteam==2.53.11\n")

    changed = script["sync_lower_bounds"](pyproject, requirements)

    assert changed is False
    assert pyproject.read_text() == original
