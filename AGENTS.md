# Agent Context

QSDsan is the reusable modeling engine for sanitation and resource-recovery systems. Put generic components, streams, units, processes, system behavior, TEA/LCA utilities, and public APIs here.

For cross-repo structure, dependency, import, unit/process namespace, or tutorial/docs authoring work, use the architecture skill for the agent runtime in use:

```text
.codex/skills/qsdsan-exposan-architecture/SKILL.md
.claude/skills/qsdsan-exposan-architecture/SKILL.md
```

Keep EXPOsan out of QSDsan runtime dependencies. EXPOsan-dependent checks belong in explicit integration coverage.

## Changelog & releases

- The changelog is `CHANGELOG.rst` at the repo root (surfaced in the docs via `docs/source/CHANGELOG.rst`, which just `.. include::`s it). EXPOsan has no changelog of its own, so record **major** EXPOsan changes here; minor ones such as typo fixes need no entry.
- Versions are static in each repo's `pyproject.toml` (`version = "X.Y.Z"`). EXPOsan also pins an exact `qsdsan==X.Y.Z`, which must be bumped to match for coordinated releases.
- For a release, add a `` `X.Y.Z`_ `` section at the top of `CHANGELOG.rst` with a matching footer link target (`.. _X.Y.Z: https://github.com/QSD-Group/QSDsan/releases/tag/vX.Y.Z`), then bump the `pyproject` version(s).

## Packaging data files

- Any non-`.py` file QSDsan reads at runtime (e.g., `units_of_measure.txt`, files under `data/`) **must** be declared in `[tool.setuptools.package-data]` in `pyproject.toml`. Editable installs read the source tree, so a missing entry still works locally — but the built wheel/sdist omits the file, which breaks `import qsdsan` for non-editable installs (and fails every EXPOsan test at collection, since importing qsdsan is the first thing they do).
- After adding or relocating such a file, verify it ships: build a wheel (`uv build --wheel`) and confirm the file is inside it. Local (editable) tests passing is **not** sufficient to catch this class of bug.

## Testing conventions

- For `qsdsan.unit_operations.bst.*` (BioSTEAM-inherited wrappers), tests layer three things but deliberately skip frozen numeric baselines: (1) **smoke** — instantiate + `simulate()`, (2) **parity** — compare against the raw `bst.units.X` equivalent (flows, `installed_cost`, `utility_cost`, `power_utility.rate`), (3) **add-on** — assert the SanUnit mixin surface (`add_OPEX`, `lifetime`, `construction`, etc.) persists through `__init__`/`simulate()`. Numeric-drift detection is EXPOsan's job (system-level integration tests), not QSDsan's — duplicating it here just adds churn on every BioSTEAM bump. QSDsan-native units (static/dynamic) may warrant tighter numeric coverage.
- Any new bst wrapper that inherits BioSTEAM-first in its MRO (or overrides `__init__`) must call `self._init_sanunit_addons(**addon_kwargs)` after the BioSTEAM parent `__init__`, or the LCA/add-on attribute surface silently won't exist on instances. Simpler alternative: inherit `(SanUnit, BSTX)` and call `SanUnit.__init__` explicitly (see `_abstract.py`/`_tank.py`/`_pumping.py`).
- `--doctest-modules` (set in `pytest.ini`) only collects doctests from modules **within the path pytest is told to scan** — a workflow invocation like `pytest tests` silently skips every `qsdsan/**/*.py` doctest. Run bare `pytest` (no path arg) so both `qsdsan/` doctests and `tests/` unit tests are collected.
- Tutorial notebooks are execution-tested by `.github/workflows/build-notebooks.yml` (`pytest --nbval-lax --current-env`), not by the docs build (which only re-executes notebooks that ship zero stored outputs) and not by the main pytest job (`docs/` is in `norecursedirs`). `--nbval-lax` only fails on unexpected exceptions; deliberate teaching-error cells must be tagged `raises-exception` in cell metadata (don't use try/except-print for those).
- Do not proactively add new test files for a bug fix — verify inline (reproduce, confirm the fix, run existing tests) and let the user decide whether a regression test is wanted.
- Always update docstrings/parameter docs/inline examples in the same change as any renamed/removed parameter, default, or behavior — stale docs are treated as bugs, not follow-ups.

## Known upstream/engineering gotchas

- On **released** biosteam (2.53.11 and earlier, what QSDsan currently pins), `System.set_tolerance(rT=...)` writes to `temperature_tolerance` (the absolute tolerance), not `relative_temperature_tolerance` — both `T=` and `rT=` set the same slot in `_system.py`. Set the relative-temperature tolerance via the attribute directly instead. **Fixed on BioSTEAM master** (`T=`/`rT=` now assign separate slots, confirmed in `BioSTEAM-platform/biosteam/biosteam/_system.py`) — re-check whether a released PyPI version already carries the fix before applying this workaround, and drop the workaround once QSDsan's biosteam floor is bumped past the fixed release.
- `System` solver defaults (biosteam 2.53.x / qsdsan 1.5.x): `molar_tolerance`=1.0 kmol/hr, `relative_molar_tolerance`=0.01, `temperature_tolerance`=0.1 K, `relative_temperature_tolerance`=0.001, `maxiter`=200, `method`='aitken'. Non-convergence raises `RuntimeError`.
- thermosteam 0.53.5 (vs 0.53.4) started estimating and storing `Hf`/`HHV`/`LHV`/`combustion` for any chemical with a formula but no measured Hf (a Dulong/Boie fuel correlation applied out-of-domain to inorganics). QSDsan guards this in `qsdsan/_component.py` (`_has_measured_Hf` + `Component._clear_estimated_energetics`): these fields are forced back to `None` unless the chemical is organic, has a measured Hf, or the user supplied one explicitly. Keep this guard regardless of the thermosteam version installed — never pin thermosteam's own missing-data estimation or display output in a QSDsan doctest/test, it differs by version.
- When a `_design` method produces absurd geometry (extreme aspect ratio, huge area), do not "fix" it by swapping which flow attribute the sizing call uses (e.g. `_mixed.F_vol` → `ins[0].F_vol`) — the constraint is usually correctly flagging an unphysical operating point (e.g. METAB's UASB+M ~400× sidestream recirculation for degassing). The right fix is a geometric sanity bound plus a warning, not hiding the physics violation by re-deriving Q.
- Testing against unreleased BioSTEAM/thermosteam master: master and the last PyPI release can share the same version string, so distinguish by `biosteam.__file__` not `__version__`. A CI failure involving numba `cache=True` functions (FUG/flexsolve/`_tea.py`) that won't reproduce locally at the exact same commits is a stale JIT cache first, a real regression second.

## Public API surface

QSDsan 1.5.3+ re-exports select BioSTEAM/thermosteam names so most users don't need `import biosteam`/`import thermosteam` directly (`qs.settings`, `qs.Thermo`, `qs.PowerUtility`, etc. — see `docs/source/api/public_api.rst` and `tests/test_public_api.py`, which asserts each re-export `is` its source so an upstream rename fails CI). The only biosteam global that needed a settable-accessor wrapper is `bst.CE` (exposed as `qs.CEPCI`) — anything else is either a class attribute or a singleton mutated in place, so a plain re-export already writes through to shared state. `qsdsan.utils.cost` is a thin wrapper (not a plain re-export): it also accepts `CEPCI=` as an alias for biosteam's `CE=`.

## Changelog discipline

Keep every `CHANGELOG.rst` entry to 1-3 sentences: what changed, the user-facing effect, one notable gotcha if any. Do not enumerate every affected class/file — name 1-2 representative examples instead of a manifest. Tutorial/docs changes collapse into a single umbrella bullet describing topics, with no section numbers (`§4`, `§5.6`) since tutorials get renumbered and those references go stale.

## Tutorial & docs prose conventions

- Address the audience as "users", not "students" (tutorials reach practitioners/researchers too, not just a classroom).
- Minimize em-dashes; avoid the word "knobs" (use "settings"/"parameters"/"controls" instead).
- Notes/call-outs are `<div class="alert alert-info">` blocks, not italicized asides.
- Never put inline `` `code` `` inside `*italic*` spans, inside markdown link text, or a markdown link inside an italic span — nbsphinx/pandoc (RTD-only, not GitHub/Jupyter) renders these broken (literal backticks, raw RST, or a swallowed space before the link).
- Keep the `<a class="anchor" id="sN">` / `<a href="#sN">` navigation markers intact in every tutorial except `1_Helpful_Basics.ipynb` — this is for visual consistency across the series, not real HTML navigation.
- Figures (diagrams, system schematics) are hand-authored SVG **pairs**, `<name>_light.svg` / `<name>_dark.svg` under `docs/source/images/<topic>/`, wired via Furo's `only-light`/`only-dark` classes; captions go in the surrounding markdown, never inside the SVG.
- Deliberate teaching errors: let the cell raise and show the traceback, tag it `raises-exception` in cell metadata. Do not wrap in try/except-print.
- Re-execute a touched notebook with the venv kernel before considering it done (see the dev-environment notes in the platform-root `AGENTS.md`).
