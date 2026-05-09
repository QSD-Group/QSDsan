# QSDsan Unit Namespaces Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add clear `qsdsan.sanunits.bst`, `qsdsan.sanunits.static`, and `qsdsan.sanunits.dynamic` namespaces while preserving existing `qsdsan.sanunits.<ClassName>` imports.

**Architecture:** Start with re-export modules so the public structure changes without moving implementation files. Tests lock down import behavior and representative class placement before deeper refactors. Documentation explains the three unit groups and leaves source-file movement as a later, lower-risk migration.

**Tech Stack:** Python package modules, `pytest`, existing QSDsan/Sphinx documentation.

---

## File Structure

- Create `qsdsan/sanunits/bst.py`: public re-export namespace for BioSTEAM-compatible QSDsan wrappers.
- Create `qsdsan/sanunits/static.py`: public re-export namespace for steady-state QSDsan sanitation/resource-recovery units.
- Create `qsdsan/sanunits/dynamic.py`: public re-export namespace for dynamic QSDsan units.
- Modify `qsdsan/sanunits/__init__.py`: import the three namespaces and add them to `__all__`.
- Create `tests/test_sanunit_namespaces.py`: verify imports, representative classes, and backward compatibility.
- Modify `docs/source/api/sanunits/_index.rst`: mention the three namespaces and link their API pages.
- Create `docs/source/api/sanunits/bst.rst`, `docs/source/api/sanunits/static.rst`, and `docs/source/api/sanunits/dynamic.rst`: namespace API pages.

## Chunk 1: Public Namespace Modules

### Task 1: Add Failing Namespace Import Tests

**Files:**
- Create: `tests/test_sanunit_namespaces.py`

- [ ] **Step 1: Write tests for namespace availability**

```python
import qsdsan as qs


def test_sanunit_namespaces_are_public():
    from qsdsan.sanunits import bst, static, dynamic

    assert qs.sanunits.bst is bst
    assert qs.sanunits.static is static
    assert qs.sanunits.dynamic is dynamic
```

- [ ] **Step 2: Run the new test and verify it fails**

Run: `pytest tests/test_sanunit_namespaces.py::test_sanunit_namespaces_are_public -v`

Expected: FAIL because `qsdsan.sanunits.bst`, `static`, and `dynamic` do not exist yet.

### Task 2: Add Empty Namespace Modules

**Files:**
- Create: `qsdsan/sanunits/bst.py`
- Create: `qsdsan/sanunits/static.py`
- Create: `qsdsan/sanunits/dynamic.py`
- Modify: `qsdsan/sanunits/__init__.py`

- [ ] **Step 1: Create empty modules**

Each new module should start with a short docstring and an empty `__all__`:

```python
"""BioSTEAM-compatible unit operations for QSDsan."""

__all__ = ()
```

- [ ] **Step 2: Import modules from `qsdsan.sanunits.__init__`**

Add this near the bottom of `qsdsan/sanunits/__init__.py`, after implementation modules are imported:

```python
from . import bst, static, dynamic
```

Add namespace names to `__all__`:

```python
__all__ = (
    ...
    'bst',
    'static',
    'dynamic',
)
```

- [ ] **Step 3: Run namespace availability test**

Run: `pytest tests/test_sanunit_namespaces.py::test_sanunit_namespaces_are_public -v`

Expected: PASS.

- [ ] **Step 4: Commit Chunk 1**

```bash
git add qsdsan/sanunits/__init__.py qsdsan/sanunits/bst.py qsdsan/sanunits/static.py qsdsan/sanunits/dynamic.py tests/test_sanunit_namespaces.py
git commit -m "Add sanunit namespace modules"
```

## Chunk 2: Class Classification By Re-Export

### Task 3: Lock Down Representative Classes

**Files:**
- Modify: `tests/test_sanunit_namespaces.py`

- [ ] **Step 1: Add representative class tests**

```python
def test_bst_namespace_reexports_compatible_wrappers():
    from qsdsan.sanunits import bst

    assert bst.Mixer is qs.sanunits.Mixer
    assert bst.Splitter is qs.sanunits.Splitter
    assert bst.Pump is qs.sanunits.Pump
    assert bst.Flash is qs.sanunits.Flash
    assert bst.BinaryDistillation is qs.sanunits.BinaryDistillation
    assert bst.HXutility is qs.sanunits.HXutility
    assert bst.StorageTank is qs.sanunits.StorageTank
    assert bst.IsothermalCompressor is qs.sanunits.IsothermalCompressor


def test_static_namespace_reexports_qsdsan_static_units():
    from qsdsan.sanunits import static

    assert static.Excretion is qs.sanunits.Excretion
    assert static.SepticTank is qs.sanunits.SepticTank
    assert static.Sedimentation is qs.sanunits.Sedimentation
    assert static.Screening is qs.sanunits.Screening
    assert static.SludgePasteurization is qs.sanunits.SludgePasteurization


def test_dynamic_namespace_reexports_dynamic_units():
    from qsdsan.sanunits import dynamic

    assert dynamic.DynamicInfluent is qs.sanunits.DynamicInfluent
    assert dynamic.HydraulicDelay is qs.sanunits.HydraulicDelay
    assert dynamic.CSTR is qs.sanunits.CSTR
    assert dynamic.PFR is qs.sanunits.PFR
```

- [ ] **Step 2: Run tests and verify they fail**

Run: `pytest tests/test_sanunit_namespaces.py -v`

Expected: FAIL because the namespace modules do not re-export classes yet.

### Task 4: Populate `bst` Re-Exports

**Files:**
- Modify: `qsdsan/sanunits/bst.py`

- [ ] **Step 1: Add BioSTEAM-compatible re-exports**

Start with classes already tested in `tests/test_bst_units.py` and obvious thin wrappers:

```python
"""BioSTEAM-compatible unit operations for QSDsan."""

from ._abstract import Mixer, Splitter, FakeSplitter, ReversedSplitter
from ._compressor import IsothermalCompressor
from ._distillation import (
    BinaryDistillation,
    ShortcutColumn,
    MESHDistillation,
    AdiabaticMultiStageVLEColumn,
)
from ._facilities import ProcessWaterCenter
from ._flash import Flash
from ._heat_exchanging import HeatExchangerNetwork, HXprocess, HXutility
from ._pumping import Pump
from ._tank import Tank, MixTank, StorageTank

__all__ = (
    'Mixer',
    'Splitter',
    'FakeSplitter',
    'ReversedSplitter',
    'IsothermalCompressor',
    'BinaryDistillation',
    'ShortcutColumn',
    'MESHDistillation',
    'AdiabaticMultiStageVLEColumn',
    'ProcessWaterCenter',
    'Flash',
    'HeatExchangerNetwork',
    'HXprocess',
    'HXutility',
    'Pump',
    'Tank',
    'MixTank',
    'StorageTank',
)
```

- [ ] **Step 2: Run `bst` namespace tests**

Run: `pytest tests/test_sanunit_namespaces.py::test_bst_namespace_reexports_compatible_wrappers tests/test_bst_units.py -v`

Expected: PASS.

### Task 5: Populate `static` Re-Exports

**Files:**
- Modify: `qsdsan/sanunits/static.py`

- [ ] **Step 1: Re-export steady-state QSDsan unit groups**

Use explicit imports from existing implementation modules. Include classes whose primary contract is static simulation, sanitation/resource-recovery design, cost, LCA, or system-specific behavior.

Start with the classes needed by the test, then expand with the full classification table after reviewing each module-level `__all__`.

- [ ] **Step 2: Run static namespace tests**

Run: `pytest tests/test_sanunit_namespaces.py::test_static_namespace_reexports_qsdsan_static_units -v`

Expected: PASS.

### Task 6: Populate `dynamic` Re-Exports

**Files:**
- Modify: `qsdsan/sanunits/dynamic.py`

- [ ] **Step 1: Re-export dynamic units**

Use explicit imports for units with `_compile_AE`, `_compile_ODE`, `isdynamic`, process-model coupling, or dedicated dynamic behavior.

Start with:

```python
from ._dynamic_influent import DynamicInfluent
from ._pumping import HydraulicDelay
from ._suspended_growth_bioreactor import CSTR, PFR
```

Then add other dynamic process units after checking their implementation and tests.

- [ ] **Step 2: Run dynamic namespace tests**

Run: `pytest tests/test_sanunit_namespaces.py::test_dynamic_namespace_reexports_dynamic_units -v`

Expected: PASS.

- [ ] **Step 3: Run all namespace tests**

Run: `pytest tests/test_sanunit_namespaces.py tests/test_bst_units.py -v`

Expected: PASS.

- [ ] **Step 4: Commit Chunk 2**

```bash
git add qsdsan/sanunits/bst.py qsdsan/sanunits/static.py qsdsan/sanunits/dynamic.py tests/test_sanunit_namespaces.py
git commit -m "Classify sanunits into public namespaces"
```

## Chunk 3: Documentation

### Task 7: Add Namespace API Pages

**Files:**
- Create: `docs/source/api/sanunits/bst.rst`
- Create: `docs/source/api/sanunits/static.rst`
- Create: `docs/source/api/sanunits/dynamic.rst`
- Modify: `docs/source/api/sanunits/_index.rst`

- [ ] **Step 1: Add `bst` API page**

```rst
BioSTEAM-Compatible Units
=========================

.. automodule:: qsdsan.sanunits.bst
   :members:
   :undoc-members:
   :show-inheritance:
```

- [ ] **Step 2: Add `static` API page**

```rst
Static QSDsan Units
===================

.. automodule:: qsdsan.sanunits.static
   :members:
   :undoc-members:
   :show-inheritance:
```

- [ ] **Step 3: Add `dynamic` API page**

```rst
Dynamic QSDsan Units
====================

.. automodule:: qsdsan.sanunits.dynamic
   :members:
   :undoc-members:
   :show-inheritance:
```

- [ ] **Step 4: Link pages from the sanunits index**

Add the three namespace pages near the top of the `toctree` in `docs/source/api/sanunits/_index.rst`.

- [ ] **Step 5: Run documentation-adjacent import check**

Run: `python -c "import qsdsan as qs; from qsdsan.sanunits import bst, static, dynamic; print(bst.__all__[:3], static.__all__[:3], dynamic.__all__[:3])"`

Expected: command exits successfully and prints representative class names.

- [ ] **Step 6: Commit Chunk 3**

```bash
git add docs/source/api/sanunits/_index.rst docs/source/api/sanunits/bst.rst docs/source/api/sanunits/static.rst docs/source/api/sanunits/dynamic.rst
git commit -m "Document sanunit namespaces"
```

## Chunk 4: Classification Audit

### Task 8: Audit Full `sanunits.__all__`

**Files:**
- Modify: `qsdsan/sanunits/bst.py`
- Modify: `qsdsan/sanunits/static.py`
- Modify: `qsdsan/sanunits/dynamic.py`
- Modify: `tests/test_sanunit_namespaces.py`

- [ ] **Step 1: Generate the current class inventory**

Run: `python -c "import qsdsan as qs; print('\\n'.join(qs.sanunits.__all__))"`

Expected: prints all public sanunit names.

- [ ] **Step 2: Add coverage test for classified public classes**

Add a test that checks every public class-like name in `qs.sanunits.__all__` is either intentionally classified or listed in an `UNCLASSIFIED` tuple in the test.

- [ ] **Step 3: Classify or explicitly defer each public class**

Update namespace `__all__` tuples until the audit test passes. If a class is ambiguous, defer it in the test with a comment and add it to the design note's open decisions later.

- [ ] **Step 4: Run full namespace test suite**

Run: `pytest tests/test_sanunit_namespaces.py tests/test_bst_units.py tests/test_sanunit.py -v`

Expected: PASS.

- [ ] **Step 5: Commit Chunk 4**

```bash
git add qsdsan/sanunits/bst.py qsdsan/sanunits/static.py qsdsan/sanunits/dynamic.py tests/test_sanunit_namespaces.py
git commit -m "Audit sanunit namespace classification"
```

## Chunk 5: Final Verification

### Task 9: Run Focused Regression Tests

**Files:**
- No source edits expected.

- [ ] **Step 1: Run focused tests**

Run: `pytest tests/test_sanunit_namespaces.py tests/test_bst_units.py tests/test_sanunit.py tests/test_dyn_sys.py -v`

Expected: PASS.

- [ ] **Step 2: Check git status**

Run: `git status --short`

Expected: only unrelated pre-existing files should remain modified.

- [ ] **Step 3: Commit any missed documentation or test cleanup**

Only commit files related to this namespace work.

```bash
git add <relevant files>
git commit -m "Finalize sanunit namespace migration"
```

## Notes And Open Questions

- Do not move existing implementation files in this first pass.
- Do not remove existing `qsdsan.sanunits.<ClassName>` imports.
- Be conservative with the `dynamic` namespace: include only units with explicit dynamic behavior.
- Some classes may reasonably appear in both `static` and `dynamic` if they support both modes. Prefer one canonical namespace first, then revisit after user feedback.
- The docs should describe `dynamic` as QSDsan's distinctive modeling layer, not merely a utility bucket.
