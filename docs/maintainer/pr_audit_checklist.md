# QSDsan PR Audit Checklist

Use this checklist when reviewing QSDsan pull requests that touch BioSTEAM/Thermosteam imports, stream/unit APIs, unit registries, or stream/unit documentation.

## First Pass

1. Inspect the PR diff before editing:
   - `git diff --name-only <base>...HEAD`
   - `git diff <base>...HEAD -- qsdsan tests docs`
2. Identify whether the PR touches:
   - BioSTEAM/Thermosteam imports
   - private upstream APIs
   - unit registry/unit-of-measure logic
   - `Stream`, `SanStream`, `WasteStream`, `SanUnit`, or BioSTEAM unit wrappers
   - stream/unit tutorials or API docs
3. Keep changes narrow. Do not turn an audit into a broad style sweep.

## Import Boundary Rules

Private upstream APIs must be centralized:

- Allowed location: `qsdsan/_compat.py`
- Flag elsewhere: `biosteam._*`, `thermosteam._*`, `bst._*`, `tmo._*`, `from biosteam._...`, `from thermosteam._...`

Search:

```powershell
rg -n "biosteam\._|thermosteam\._|bst\._|tmo\._|from biosteam\._|from thermosteam\._" qsdsan -g "*.py"
```

Public BioSTEAM/Thermosteam APIs are fine when they express a true upstream dependency:

- Keep direct imports for specific BioSTEAM units, facilities, decorators, design tools, exceptions, reactions, and settings helpers when those APIs are the thing being wrapped or subclassed.
- Prefer QSDsan facade imports for generic shared objects already exposed by `qsdsan`, especially `Stream`, `Unit`, `System`, `HeatUtility`, `PowerUtility`, `Scope`, `Model`, `Metric`, and `Parameter`, when import cycles allow it.

Current rule for sanunit modules:

- Do not use `from biosteam import Stream` in `qsdsan/sanunits`; use `from .. import Stream`.
- Do not replace specific imports such as `bst.units.Flash`, `biosteam.units.Pump`, design tools, or decorators just for cosmetic consistency.

## Unit Registry Rules

Avoid private Pint internals such as `ureg._units`. Prefer public parsing or conversion APIs with exception handling.

For unit-definition idempotency, verify repeated definitions still skip already-defined aliases and define missing aliases only when needed.

Relevant test:

```powershell
.\.venv\Scripts\python.exe -m pytest tests\test_units_of_measure.py
```

## Stream Taxonomy

Preserve this public meaning:

- `Stream`: thermosteam/BioSTEAM material stream behavior only.
- `SanStream`: `Stream` plus QSDsan stream-level LCA/impact functionality.
- `WasteStream`: `SanStream` plus wastewater-specific aggregate properties and influent characterization models.

Do not rename these casually and do not make `qsdsan.Stream` mean `WasteStream`.

When documenting:

- Use `:class:`~.SanStream`` and `:class:`~.WasteStream`` for QSDsan classes in Sphinx prose.
- Explain that `WasteStream` does not require literal disposal as waste; it means wastewater-modeling capabilities.

## Tests To Choose From

Run focused tests based on the touched area:

```powershell
.\.venv\Scripts\python.exe -m pytest tests\test_compat.py tests\test_import_consistency.py
.\.venv\Scripts\python.exe -m pytest tests\test_units_of_measure.py
.\.venv\Scripts\python.exe -m pytest tests\test_waste_stream.py tests\test_sanunit.py tests\test_bst_units.py
.\.venv\Scripts\python.exe -m pytest tests\test_dyn_sys.py tests\test_process.py tests\test_component.py
```

For docs-only stream taxonomy changes, at minimum run stream-related tests and, if practical, a Sphinx build. If the full docs build fails because of pre-existing warnings/errors, report that clearly and identify whether failures are unrelated.

## Review Output

When reviewing a PR, lead with concrete findings:

- Severity
- File and line
- Why it matters for QSDsan's BioSTEAM/Thermosteam boundary or public API
- Minimal suggested fix

If no issues are found, say that directly and mention the focused tests or scans that were run.
