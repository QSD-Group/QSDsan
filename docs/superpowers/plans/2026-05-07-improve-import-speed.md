# Improve QSDsan Import Speed Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `import qsdsan` reliably complete quickly, then reduce QSDsan's avoidable import-time work without breaking the public API.

**Architecture:** Treat this as two layers: first stabilize upstream Numba/BioSTEAM import behavior in the local environment, then slim QSDsan's eager import graph while preserving legacy top-level names. The safest first code changes are lazy imports for optional/heavy utilities (`stats`, docs examples, plotting/SALib/SymPy helpers), followed by a more cautious public API compatibility pass for `sanunits`.

**Tech Stack:** Python 3.12, QSDsan 1.4.4, BioSTEAM 2.53.10, Thermosteam, Numba, pytest, `python -X importtime`.

---

## Baseline Findings

Measurements were taken on Windows from `C:\Users\Yalin\Documents\Coding\QSDsan-unit-refactor\QSDsan` using `..\.venv\Scripts\python.exe` with the repo installed editable.

Important evidence:

- `import qsdsan` without a configured `NUMBA_CACHE_DIR` did not finish after 6 minutes.
- `faulthandler.dump_traceback_later(60, exit=True)` showed the import blocked at:
  - `qsdsan/__init__.py:26`, importing `biosteam`
  - `biosteam/__init__.py:37`, applying `@numba.njit(cache=True)` to `f_dummy`
  - `numba/core/caching.py`, inside `ensure_cache_path`
  - `tempfile.NamedTemporaryFile` / `_mkstemp_inner`
- Setting `NUMBA_CACHE_DIR` to a repo-local cache made upstream imports complete:
  - `import thermosteam`: about 7.6 s
  - `import biosteam`: about 29.0 s
  - `import qsdsan`: about 20.6 s in a subsequent run
- With `NUMBA_CACHE_DIR=$PWD\.numba_cache`, `python -X importtime -c "import qsdsan"` showed:
  - `qsdsan`: about 18.8 s cumulative
  - `biosteam`: about 18.0 s cumulative
  - `thermosteam`: about 5.7 s cumulative inside `biosteam`
  - `qsdsan.utils`: about 0.44 s, mostly `qsdsan.utils.parsing` pulling full `sympy`
  - `qsdsan.sanunits`: about 0.037 s after upstream imports are already loaded
  - `qsdsan.stats`: about 0.22 s, pulling `seaborn`, `SALib`, and `scipy.signal`

Interpretation:

- The immediate "import never finishes" problem is probably environmental/upstream Numba cache behavior, not QSDsan application code.
- The dominant successful-import cost is BioSTEAM/Thermosteam, which QSDsan imports eagerly at top level.
- QSDsan still has cleanup opportunities: `qsdsan.__init__` eagerly imports `stats`, `sanunits`, `processes`, `equipments`, docs examples, plotting utilities, SALib utilities, and full SymPy helpers.

## File Map

- Modify: `qsdsan/__init__.py`
  - Owns top-level API exports and eager imports.
  - Likely place for PEP 562 `__getattr__` lazy exports, if adopted.
- Modify: `qsdsan/utils/__init__.py`
  - Currently imports every utility module eagerly, including plotting and dynamics helpers.
- Modify: `qsdsan/utils/parsing.py`
  - Imports full `sympy`; likely can narrow imports and avoid `solve` unless needed.
- Modify: `qsdsan/stats.py`
  - Imports plotting/SALib modules eagerly; can move many imports into functions.
- Modify: `qsdsan/utils/doc_examples.py`
  - Imports `chaospy.distributions` eagerly; can move inside `create_example_model`.
- Modify cautiously: `qsdsan/sanunits/__init__.py`
  - Imports all unit modules and creates one cached Numba function at import.
  - This is public API-heavy and should be changed only after `stats`/utils wins are verified.
- Test: `tests/test_import_consistency.py`
  - Extend to assert lazy public names still resolve.
- Test: new `tests/test_import_speed.py`
  - Add behavioral guardrails without making CI flaky by default.
- Docs: `CONTRIBUTING.rst` or `docs/source` developer page
  - Document `NUMBA_CACHE_DIR` workaround if reproducible.

## Chunk 1: Stabilize Baseline and Add Import Benchmark Harness

### Task 1: Add a local benchmark script

**Files:**
- Create: `benchmarks/import_time.py`
- Test: manual command only

- [ ] **Step 1: Create benchmark directory**

Run:

```powershell
New-Item -ItemType Directory -Force benchmarks
```

Expected: `benchmarks` exists.

- [ ] **Step 2: Add a small subprocess benchmark**

Create `benchmarks/import_time.py`:

```python
from __future__ import annotations

import argparse
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path


def run_import(module: str, repeat: int, timeout: float, numba_cache_dir: str | None) -> list[float]:
    env = os.environ.copy()
    if numba_cache_dir:
        env["NUMBA_CACHE_DIR"] = str(Path(numba_cache_dir).resolve())

    times = []
    for _ in range(repeat):
        start = time.perf_counter()
        subprocess.run(
            [sys.executable, "-c", f"import {module}"],
            check=True,
            timeout=timeout,
            env=env,
        )
        times.append(time.perf_counter() - start)
    return times


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("module", nargs="?", default="qsdsan")
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--numba-cache-dir", default=None)
    args = parser.parse_args()

    times = run_import(args.module, args.repeat, args.timeout, args.numba_cache_dir)
    print(f"module={args.module}")
    print(f"times={', '.join(f'{i:.3f}s' for i in times)}")
    print(f"min={min(times):.3f}s mean={statistics.mean(times):.3f}s")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run with local Numba cache**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe benchmarks\import_time.py qsdsan --repeat 3 --numba-cache-dir .\.numba_cache
```

Expected: three timings are printed and all complete under 120 s.

- [ ] **Step 4: Capture import graph**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -X importtime -c "import qsdsan" 2> importtime_qsdsan_before.txt
```

Expected: `importtime_qsdsan_before.txt` exists for local review only. Do not commit it.

- [ ] **Step 5: Commit**

```powershell
git add benchmarks/import_time.py
git commit -m "test: add import time benchmark helper"
```

## Chunk 2: Fix the Local Import Hang Before Optimizing Code

### Task 2: Reproduce and document the Numba cache issue

**Files:**
- Modify: `CONTRIBUTING.rst`
- Optional Modify: `.gitignore`

- [ ] **Step 1: Verify hang stack**

Run:

```powershell
..\.venv\Scripts\python.exe -c "import faulthandler; faulthandler.dump_traceback_later(60, exit=True); import qsdsan"
```

Expected on the affected machine: stack shows `biosteam/__init__.py` applying `@numba.njit(cache=True)` and Numba cache setup inside `tempfile.NamedTemporaryFile`.

- [ ] **Step 2: Verify workaround**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -c "import qsdsan; print(qsdsan.__version__)"
```

Expected: import completes and prints `1.4.4`.

- [ ] **Step 3: Ignore local cache artifacts**

If `.numba_cache/` is created in the repo root, add this to `.gitignore`:

```gitignore
.numba_cache/
```

- [ ] **Step 4: Document developer workaround**

Add a short note to `CONTRIBUTING.rst`:

```rst
Import-time troubleshooting
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If ``import qsdsan`` appears to hang on Windows while importing BioSTEAM/Numba,
try setting a local Numba cache directory before running tests or profiling:

.. code-block:: powershell

   $env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'

This avoids slow or blocked cache checks in the default temporary directory.
```

- [ ] **Step 5: Verify docs-only change**

Run:

```powershell
git diff -- CONTRIBUTING.rst .gitignore
```

Expected: only the troubleshooting note and optional `.numba_cache/` ignore entry changed.

- [ ] **Step 6: Commit**

```powershell
git add CONTRIBUTING.rst .gitignore
git commit -m "docs: document numba cache workaround for imports"
```

## Chunk 3: Add Import API Guardrails

### Task 3: Strengthen import compatibility tests

**Files:**
- Modify: `tests/test_import_consistency.py`
- Create: `tests/test_import_speed.py`

- [ ] **Step 1: Inspect current consistency test**

Run:

```powershell
Get-Content tests\test_import_consistency.py
```

Expected: understand what names are already covered before adding more.

- [ ] **Step 2: Add representative public API checks**

Add tests that assert these names remain available after `import qsdsan as qs`:

```python
def test_top_level_public_api_representatives():
    import qsdsan as qs

    for name in (
        "Component",
        "Components",
        "SanStream",
        "WasteStream",
        "Process",
        "Processes",
        "ImpactIndicator",
        "ImpactItem",
        "Construction",
        "Equipment",
        "SanUnit",
        "TEA",
        "LCA",
        "stats",
        "sanunits",
        "processes",
        "equipments",
    ):
        assert hasattr(qs, name), name
```

- [ ] **Step 3: Add import smoke subprocess test**

Create `tests/test_import_speed.py`:

```python
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def test_import_qsdsan_completes_with_local_numba_cache(tmp_path):
    env = os.environ.copy()
    env.setdefault("NUMBA_CACHE_DIR", str(tmp_path / "numba-cache"))

    subprocess.run(
        [sys.executable, "-c", "import qsdsan"],
        check=True,
        timeout=120,
        env=env,
        cwd=Path(__file__).resolve().parents[1],
    )
```

Keep this as a completion guard, not a strict speed threshold; hard performance thresholds are usually flaky on CI.

- [ ] **Step 4: Run new tests**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -m pytest tests\test_import_consistency.py tests\test_import_speed.py -q
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```powershell
git add tests/test_import_consistency.py tests/test_import_speed.py
git commit -m "test: guard qsdsan import compatibility"
```

## Chunk 4: Lazily Import Heavy Stats Dependencies

### Task 4: Defer SALib, seaborn, matplotlib, and scipy-heavy imports

**Files:**
- Modify: `qsdsan/stats.py`
- Test: `tests/test_import_consistency.py`, targeted stats tests if present, import benchmark

- [ ] **Step 1: Identify which functions need each heavy import**

Search:

```powershell
rg "sns|plt|morris_sampler|fast_sampler|rbd_sampler|sobol_sampler|morris|fast|rbd_fast|sobol|sa_plt_morris|plot_spearman" qsdsan\stats.py
```

Expected: map each imported object to the functions that use it.

- [ ] **Step 2: Replace top-level heavy imports with lightweight placeholders**

Change the top of `qsdsan/stats.py` from eager imports:

```python
import seaborn as sns
from matplotlib import pyplot as plt
from SALib.sample import (...)
from SALib.analyze import morris, fast, rbd_fast, sobol
from SALib.plotting import morris as sa_plt_morris
from biosteam.plots import plot_spearman
```

to local imports inside the functions that use them.

Keep `numpy`, `pandas`, `biosteam`, `Iterable`, and `warn` top-level unless profiling shows they matter enough to justify more churn.

- [ ] **Step 3: Use local imports near first use**

Example pattern:

```python
def generate_samples(...):
    from SALib.sample import (
        morris as morris_sampler,
        fast_sampler,
        latin as rbd_sampler,
        saltelli as sobol_sampler,
    )
    ...
```

For plotting functions:

```python
def plot_uncertainties(...):
    import seaborn as sns
    from matplotlib import pyplot as plt
    ...
```

- [ ] **Step 4: Run import benchmark**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe benchmarks\import_time.py qsdsan --repeat 3 --numba-cache-dir .\.numba_cache
```

Expected: `qsdsan.stats` no longer pulls `seaborn`, `SALib`, or `scipy.signal` during plain import.

- [ ] **Step 5: Run tests**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -m pytest tests\test_import_consistency.py -q
```

Expected: pass.

- [ ] **Step 6: Commit**

```powershell
git add qsdsan/stats.py
git commit -m "perf: defer stats plotting and sensitivity imports"
```

## Chunk 5: Lazily Import Documentation Example Dependencies

### Task 5: Move `chaospy` import into example model creation

**Files:**
- Modify: `qsdsan/utils/doc_examples.py`
- Test: `tests/test_import_consistency.py`

- [ ] **Step 1: Move `chaospy` import**

Change:

```python
from chaospy import distributions as shape
```

to a local import inside `create_example_model` before the first use of `shape`:

```python
def create_example_model(...):
    from chaospy import distributions as shape
    ...
```

- [ ] **Step 2: Run targeted import check**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -X importtime -c "import qsdsan" 2> importtime_qsdsan_after_doc_examples.txt
Select-String -Path importtime_qsdsan_after_doc_examples.txt -Pattern "chaospy"
```

Expected: no `chaospy` import during plain `import qsdsan`.

- [ ] **Step 3: Run doc example smoke check**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -c "from qsdsan.utils import create_example_model; m=create_example_model(); print(type(m).__name__)"
```

Expected: prints `Model`.

- [ ] **Step 4: Commit**

```powershell
git add qsdsan/utils/doc_examples.py
git commit -m "perf: defer chaospy import for doc examples"
```

## Chunk 6: Narrow SymPy Imports From Utility Import Path

### Task 6: Reduce `qsdsan.utils.parsing` import cost

**Files:**
- Modify: `qsdsan/utils/parsing.py`
- Test: `tests/test_process.py`, `tests/test_import_consistency.py`

- [ ] **Step 1: Check actual SymPy symbol usage**

Run:

```powershell
rg "symbols|sympify|simplify|Matrix|solve|parse_expr" qsdsan\utils\parsing.py
```

Expected: identify exact functions needing each import.

- [ ] **Step 2: Move rarely used SymPy imports into functions**

Keep only imports required by functions used on the common path. If `solve`, `Matrix`, or `parse_expr` are only used in one function, import them inside that function.

Example:

```python
def get_stoichiometric_coeff(...):
    from sympy import Matrix, simplify, solve, symbols
    from sympy.parsing.sympy_parser import parse_expr
    ...
```

- [ ] **Step 3: Re-run importtime**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -X importtime -c "import qsdsan" 2> importtime_qsdsan_after_sympy.txt
Select-String -Path importtime_qsdsan_after_sympy.txt -Pattern "sympy|qsdsan.utils.parsing"
```

Expected: plain import no longer imports the full `sympy` package through `qsdsan.utils.parsing`, or the cumulative time is materially lower.

- [ ] **Step 4: Run process tests**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -m pytest tests\test_process.py tests\test_import_consistency.py -q
```

Expected: pass.

- [ ] **Step 5: Commit**

```powershell
git add qsdsan/utils/parsing.py
git commit -m "perf: defer sympy imports in parsing utilities"
```

## Chunk 7: Consider Lazy Utility Module Exports

### Task 7: Decide whether to add PEP 562 lazy exports to `qsdsan.utils`

**Files:**
- Modify: `qsdsan/utils/__init__.py`
- Test: `tests/test_import_consistency.py`, any utility-specific tests

- [ ] **Step 1: Re-measure after Chunks 4-6**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -X importtime -c "import qsdsan" 2> importtime_qsdsan_after_utils.txt
Select-String -Path importtime_qsdsan_after_utils.txt -Pattern "qsdsan.utils"
```

Expected: decide whether `qsdsan.utils` remains worth changing.

- [ ] **Step 2: If worthwhile, introduce lazy module map**

Use PEP 562 `__getattr__` only for utility names from modules that are not required during core import. Keep `ureg`, `auom`, `ruom`, `copy_attr`, unit parsing, and other core dependencies eager if core modules need them.

Example skeleton:

```python
_LAZY_MODULES = {
    "palettes": ".colors",
    "colormaps": ".colors",
    "ExogenousDynamicVariable": ".dynamics",
}


def __getattr__(name):
    try:
        module_name = _LAZY_MODULES[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    from importlib import import_module
    value = getattr(import_module(module_name, __name__), name)
    globals()[name] = value
    return value
```

- [ ] **Step 3: Preserve `from qsdsan.utils import *`**

If `__all__` includes lazy names, verify star import resolves them:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -c "from qsdsan.utils import *; print('utils star ok')"
```

Expected: no missing names.

- [ ] **Step 4: Run tests**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -m pytest tests\test_import_consistency.py tests\test_units_of_measure.py tests\test_process.py -q
```

Expected: pass.

- [ ] **Step 5: Commit only if benchmark improves**

```powershell
git add qsdsan/utils/__init__.py
git commit -m "perf: lazily expose selected utility helpers"
```

## Chunk 8: Evaluate Top-Level `qsdsan` Lazy Exports

### Task 8: Defer optional subpackages from top-level import

**Files:**
- Modify: `qsdsan/__init__.py`
- Test: `tests/test_import_consistency.py`, full focused test set

- [ ] **Step 1: Identify top-level subpackages not needed for core classes**

Candidates already imported at the end of `qsdsan/__init__.py`:

```python
equipments,
processes,
sanunits,
stats,
```

Do not start by lazily importing `_component`, `_components`, `_sanstream`, `_waste_stream`, `_process`, `_impact_*`, `_construction`, `_equipment`, `_transportation`, `_sanunit`, `_tea`, or `_lca`; those define common top-level public classes.

- [ ] **Step 2: Replace optional package imports with lazy attributes**

Add:

```python
_LAZY_SUBMODULES = {
    "equipments": ".equipments",
    "processes": ".processes",
    "sanunits": ".sanunits",
    "stats": ".stats",
}


def __getattr__(name):
    if name in _LAZY_SUBMODULES:
        from importlib import import_module
        module = import_module(_LAZY_SUBMODULES[name], __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
```

Only do this after confirming no core module relies on `qsdsan.sanunits` being imported during `qsdsan.__init__`.

- [ ] **Step 3: Preserve `utils.wwtpump` behavior**

Current code does:

```python
from .sanunits import wwtpump
utils.__all__ = (*utils.__all__, 'wwtpump')
setattr(utils, 'wwtpump', wwtpump)
```

This forces `sanunits` import. Replace with a tiny lazy wrapper or lazy utility attribute.

Possible wrapper in `qsdsan/__init__.py`:

```python
def _get_wwtpump():
    from .sanunits import wwtpump
    return wwtpump
```

Then implement the actual lazy exposure in `qsdsan.utils` so `qsdsan.utils.wwtpump` imports `qsdsan.sanunits` only when accessed.

- [ ] **Step 4: Run compatibility checks**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -c "import qsdsan as qs; print(qs.sanunits.WWTpump); print(qs.stats.define_inputs)"
..\.venv\Scripts\python.exe -c "from qsdsan import sanunits, stats; print(sanunits.WWTpump, stats.define_inputs)"
..\.venv\Scripts\python.exe -c "from qsdsan.utils import wwtpump; print(wwtpump)"
```

Expected: all commands succeed.

- [ ] **Step 5: Run tests**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -m pytest tests\test_import_consistency.py tests\test_sanunit.py tests\test_bst_units.py -q
```

Expected: pass.

- [ ] **Step 6: Benchmark**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe benchmarks\import_time.py qsdsan --repeat 5 --numba-cache-dir .\.numba_cache
```

Expected: measurable reduction versus Chunk 1 baseline. If reduction is under ~5%, reconsider whether this complexity is worth keeping.

- [ ] **Step 7: Commit**

```powershell
git add qsdsan/__init__.py qsdsan/utils/__init__.py
git commit -m "perf: lazily import optional qsdsan subpackages"
```

## Chunk 9: Optional, High-Risk `sanunits` Lazy Export Refactor

### Task 9: Split `qsdsan.sanunits` public API into lazy exports only if needed

**Files:**
- Modify: `qsdsan/sanunits/__init__.py`
- Test: `tests/test_sanunit.py`, `tests/test_bst_units.py`, EXPOsan smoke tests if available

- [ ] **Step 1: Benchmark whether this is worth it**

Current evidence with upstream imports already loaded showed `qsdsan.sanunits` around 0.037 s cumulative. This is probably not the first place to spend time.

Proceed only if:

- `sanunits` becomes a large fraction after upstream and stats/utils changes, or
- the project specifically wants `import qsdsan` not to import all unit classes.

- [ ] **Step 2: Build a name-to-module map**

For each unit module, map exported names from each module's `__all__` to the module path.

Do this manually or with a helper script, then encode:

```python
_LAZY_EXPORTS = {
    "WWTpump": "._pumping",
    "MixTank": "._abstract",
    ...
}
```

- [ ] **Step 3: Keep `dydt_cstr` behavior explicit**

`dydt_cstr` is currently defined in `qsdsan/sanunits/__init__.py` and decorated with `@njit(cache=True)`.

Options:

- Leave it eager if many unit modules depend on it.
- Move it to a small `qsdsan/sanunits/_dynamics.py` module and import only where needed.

- [ ] **Step 4: Add `__getattr__` and preserve star import**

Use the same PEP 562 pattern as earlier. Verify:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -c "from qsdsan.sanunits import *; print(WWTpump)"
```

- [ ] **Step 5: Run broad tests**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -m pytest tests\test_sanunit.py tests\test_bst_units.py tests\test_junctions.py tests\test_dyn_sys.py -q
```

Expected: pass.

- [ ] **Step 6: Commit only if payoff justifies risk**

```powershell
git add qsdsan/sanunits/__init__.py
git commit -m "perf: lazily expose sanunit classes"
```

## Chunk 10: Final Verification and Reporting

### Task 10: Compare before/after and run focused tests

**Files:**
- No required modifications

- [ ] **Step 1: Run final benchmark**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe benchmarks\import_time.py qsdsan --repeat 7 --numba-cache-dir .\.numba_cache
```

Expected: record min and mean.

- [ ] **Step 2: Capture final import graph**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -X importtime -c "import qsdsan" 2> importtime_qsdsan_after.txt
```

Expected: `stats`, `chaospy`, full `sympy`, `seaborn`, and `SALib` are absent from the plain import path unless intentionally retained.

- [ ] **Step 3: Run focused tests**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -m pytest tests\test_import_consistency.py tests\test_import_speed.py tests\test_units_of_measure.py tests\test_component.py tests\test_process.py tests\test_sanunit.py -q
```

Expected: pass.

- [ ] **Step 4: Run full QSDsan tests if time allows**

Run:

```powershell
$env:NUMBA_CACHE_DIR=(Resolve-Path .).Path + '\.numba_cache'
..\.venv\Scripts\python.exe -m pytest tests -q
```

Expected: pass, or document unrelated failures separately.

- [ ] **Step 5: Clean profiling artifacts**

Remove local-only files:

```powershell
Remove-Item importtime_qsdsan_before.txt, importtime_qsdsan_after*.txt -ErrorAction SilentlyContinue
```

Do not remove `.numba_cache` if it is useful locally; keep it ignored.

- [ ] **Step 6: Final summary**

Report:

- Original behavior without `NUMBA_CACHE_DIR`: hung beyond 6 minutes with Numba cache stack.
- Baseline with local `NUMBA_CACHE_DIR`: mean/min import time.
- Final with local `NUMBA_CACHE_DIR`: mean/min import time.
- Import graph removals, especially `seaborn`, `SALib`, `chaospy`, and `sympy`.
- Tests run and results.
