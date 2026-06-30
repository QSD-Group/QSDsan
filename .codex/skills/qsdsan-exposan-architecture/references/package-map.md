# QSDsan/EXPOsan Package Map

Use this reference when an agent needs the package boundary, dependency policy, or unit/process namespace layout.

## Package Roles

QSDsan is the main package. It provides the reusable engine for building sanitation and resource recovery systems:

- `Component`, `Components`, `SanStream`, `WasteStream`
- `SanUnit` and the `unit_operations` package (`bst`/`static`/`dynamic`)
- process models and dynamic simulation interfaces
- `System`, TEA, LCA, construction, equipment, utilities
- public APIs intended for downstream projects

EXPOsan is the system repository. It contains example, research, benchmark, and user-built systems that compose QSDsan:

- modules such as `adm`, `asm`, `bsm1`, `bsm2`, `htl`, `metab`, `werf`
- system assembly functions, scenario analyses, validation data, and datasets
- project-specific units only when they are not generally reusable
- tests that protect system-level behavior and published benchmarks

When in doubt, ask: "Would a different sanitation modeler reuse this independent of one EXPOsan case study?" If yes, it probably belongs in QSDsan. If no, keep it in EXPOsan.

## Dependency Boundary

Local development should install the sibling clones in editable mode:

```powershell
uv pip install -e ".\QSDsan[dev]"
uv pip install -e ".\EXPOsan[dev]"
```

CI integration may intentionally use GitHub sibling dependencies:

```bash
pip install --no-cache-dir -e ".[ci,integration]"
```

Do not add EXPOsan as a normal QSDsan runtime dependency. QSDsan may have an integration extra for testing against EXPOsan. EXPOsan depends on QSDsan because EXPOsan systems are built from QSDsan primitives.

## Unit & Process Namespaces

Built-in units are organized by behavior under `qsdsan.unit_operations`:

```text
qsdsan.unit_operations
|-- bst
|-- static
`-- dynamic
```

Process models live in `qsdsan.process_models`. `qsdsan.sanunits` and `qsdsan.processes` remain as back-compat aliases (resolving to `unit_operations` and `process_models`); use the canonical names in new code. When adding a unit:

- Classify it by behavior into `bst`/`static`/`dynamic` (see the heuristics below), not by file location.
- Keep `qsdsan.unit_operations.<ClassName>` importable and add a test that the namespace exports it.

## Classification Heuristics

Use `bst` for BioSTEAM-inherited unit operations with QSDsan-added behavior:

- mixers, splitters, pumps, flash/distillation wrappers, heat exchangers, storage tanks, process-water helpers
- units that mainly adapt stream initialization or QSDsan naming
- no sanitation-specific design assumptions beyond integration glue

Use `static` for steady-state QSDsan units:

- sanitation fixtures, treatment beds, septic systems, sludge handling, screening, sedimentation, crop application
- hydrothermal, hydroprocessing, electrochemical, membrane, gas-extraction, resource-recovery units
- units with design/cost/construction/LCA logic but no primary dynamic state model

Use `dynamic` for explicit dynamic-model units:

- dynamic influent, hydraulic delay, dynamic mixers/splitters where applicable
- CSTR/PFR/suspended-growth reactors, ASM/ADM/PM2-coupled units, dynamic junctions
- classes with `isdynamic`, `state`, `dstate`, `_init_state`, `_update_state`, `_update_dstate`, `_compile_AE`, or `_compile_ODE`

If a class supports both steady-state and dynamic operation, classify by its defining contract. If dynamic state behavior is central to the class, prefer `dynamic` and document static use as a mode.

## EXPOsan Import Conventions

EXPOsan should consume QSDsan through public APIs:

- Prefer public imports such as `from qsdsan import unit_operations as su` (and `process_models as pc`).
- Use `from qsdsan.unit_operations import bst, static, dynamic` when code benefits from explicit grouping.
- Avoid imports from private QSDsan implementation modules such as `_pumping.py`.
- Keep EXPOsan systems focused on system assembly and project-specific analyses.
- If a project-specific unit becomes broadly reusable, promote it to QSDsan (with QSDsan tests) before consuming it from EXPOsan.

## EXPOsan Module Template

Each EXPOsan system or system family should live in one top-level module. Use `system.py` when the module exposes one system, and `systems.py` when it exposes multiple systems or configurations. Keep existing modules stable unless there is a focused reason to refactor.

Preferred shape for new modules:

```text
exposan/<module>/
  __init__.py
  system.py or systems.py
  model.py or models.py
  _components.py           optional
  _units.py                optional, for project-specific units
  _process_settings.py     optional
  _tea.py or _lca.py       optional
  data/                    input data only
  results/                 created only when writing results
  figures/                 created only when writing figures
```

Keep `__init__.py` thin: paths, light public imports, `load()`, and lazy access errors. Build systems in `system.py` or `systems.py`; build uncertainty/sensitivity models in `model.py` or `models.py`; keep large constants in data files or focused settings modules.

Use `None` rather than mutable default dictionaries in public functions. For multi-system dispatch, prefer an explicit registry such as `SYSTEM_CREATORS = {'A': create_systemA}` over `globals()` lookup.

## Testing Guidance

For QSDsan changes:

- test public namespace exports
- test legacy import compatibility
- test dynamic state contracts with focused unit tests
- keep EXPOsan-dependent tests behind integration coverage

For EXPOsan changes:

- test affected systems or modules directly
- preserve benchmark expectations when available
- use `state_reset_hook='reset_cache'` for dynamic simulations that already rely on it
- avoid widening tolerances unless the numerical reason is understood and documented

Useful local checks:

```powershell
uv run python -c "import qsdsan, exposan; print(qsdsan.__version__, exposan.__version__)"
uv run pytest .\QSDsan\tests
uv run pytest .\EXPOsan\tests
```
