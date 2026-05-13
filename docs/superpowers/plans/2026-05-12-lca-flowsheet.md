# LCA Flowsheet Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the global LCA registries (`ImpactItem`, `Construction`, `Transportation`) with per-flowsheet registries by subclassing BioSTEAM's `Flowsheet`, so switching between systems requires no manual `clear_lca_registries()` calls.

**Architecture:** `SanFlowsheet` subclasses `bst.Flowsheet` and adds three per-flowsheet `Registry` instances (`item`, `construction`, `transportation`). `SanMainFlowsheet` subclasses both `SanFlowsheet` and `bst.MainFlowsheet`, overriding `set_flowsheet()` to also swap the three LCA class-level registries. `ImpactIndicator` remains a global singleton (shared across all flowsheets — indicators must be consistent across systems). `SanUnit` gains a `_construction_specs` class attribute for declaring default material usage as pure data; `LCA._update_system()` resolves specs to `Construction` objects at LCA creation time, never at unit creation time.

**Tech Stack:** Python 3, BioSTEAM (`bst.Flowsheet`, `bst.MainFlowsheet`), thermosteam (`Registry`, `@registered`), QSDsan (`ImpactItem`, `Construction`, `Transportation`, `SanUnit`, `LCA`).

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `qsdsan/_flowsheet.py` | **Create** | `SanFlowsheet` + `SanMainFlowsheet` class definitions only — no singletons, no qsdsan imports at module level |
| `qsdsan/__init__.py` | **Modify** | Import `_flowsheet`; create `qs_main_flowsheet` singleton after all LCA classes are loaded; replace `Flowsheet` and `main_flowsheet` exports |
| `qsdsan/_sanunit.py` | **Modify** | Add `_construction_specs = ()` class attribute |
| `qsdsan/_lca.py` | **Modify** | Add `_resolve_construction_specs()` called from `_update_system()`; remove `clear_lca_registries = clear_lca_registries` alias |
| `qsdsan/utils/misc.py` | **Modify** | Add `DeprecationWarning` to `clear_lca_registries` |
| `tests/test_flowsheet.py` | **Create** | All tests for this feature |

**Not changed:** `_construction.py`, `_transportation.py`, `_impact_indicator.py`, `_impact_item.py` — their `@registered` mechanism works unchanged; `set_flowsheet()` swaps their class-level `registry` attribute from outside.

---

## Key Design Constraints

**`object.__setattr__` for new registry attributes:** `Flowsheet.__setattr__` raises `AttributeError` once a flowsheet is registered in `FlowsheetRegistry`. Since `super().__new__(cls, ID)` registers the flowsheet before returning, any attribute assignment in `SanFlowsheet.__new__` after the `super()` call must use `object.__setattr__(self, key, value)` to bypass the protection.

**No qsdsan imports at module level in `_flowsheet.py`:** `_flowsheet.py` is imported early in `__init__.py`, before `ImpactItem`, `Construction`, `Transportation` are defined. All qsdsan imports inside `_flowsheet.py` must be lazy (inside method bodies).

**Singleton created late in `__init__.py`:** The `qs_main_flowsheet` singleton calls `set_flowsheet()`, which lazily imports `ImpactItem`, `Construction`, `Transportation`. This import must succeed, so the singleton is created at the bottom of `__init__.py` after all LCA classes are loaded.

**`ImpactIndicator` registry is never swapped:** It uses the thermosteam global default forever. `SanFlowsheet` has no `indicator` registry attribute. `set_flowsheet()` does not touch `ImpactIndicator.registry`.

**`_construction_specs` resolution never mutates unit state across LCA creations:** When `LCA._update_system()` resolves specs into `Construction` objects, it checks if a construction for that item ID already exists on the unit (deduplication). Spec-derived `Construction` objects are registered in the current flowsheet's `Construction.registry`, so `flowsheet.clear()` cleans them up automatically via the detach loop.

---

## Task 1: Write failing tests — flowsheet isolation

**Files:**
- Create: `QSDsan/tests/test_flowsheet.py`

- [ ] **Step 1: Create the test file with all flowsheet tests**

```python
# tests/test_flowsheet.py
__all__ = ('test_flowsheet',)

def test_flowsheet():
    import qsdsan as qs
    from qsdsan._flowsheet import SanFlowsheet, SanMainFlowsheet

    # ── 1. SanFlowsheet has the three LCA registries ──────────────────────
    fs = qs.Flowsheet('_test_lca_1')
    assert hasattr(fs, 'item'),          "flowsheet missing 'item' registry"
    assert hasattr(fs, 'construction'),  "flowsheet missing 'construction' registry"
    assert hasattr(fs, 'transportation'),"flowsheet missing 'transportation' registry"
    # ImpactIndicator is global — no 'indicator' registry on flowsheet
    assert not hasattr(fs, 'indicator'), "flowsheet should NOT have 'indicator' registry"

    # ── 2. qs.Flowsheet IS SanFlowsheet ───────────────────────────────────
    assert qs.Flowsheet is SanFlowsheet, "qs.Flowsheet should be SanFlowsheet"

    # ── 3. qs.main_flowsheet IS SanMainFlowsheet ──────────────────────────
    assert isinstance(qs.main_flowsheet, SanMainFlowsheet), \
        "qs.main_flowsheet should be a SanMainFlowsheet instance"

    # ── 4. Item registry is isolated per flowsheet ───────────────────────────
    GWP = qs.ImpactIndicator('GlobalWarming_t', unit='kg CO2-eq')  # global

    with qs.Flowsheet('_test_sys_a') as fs_a:
        steel_a = qs.ImpactItem('Steel_t', 'kg', GlobalWarming_t=2.55)
        assert qs.ImpactItem.get_item('Steel_t') is steel_a

    # After exiting sys_a, we are back to the previous flowsheet
    assert qs.ImpactItem.get_item('Steel_t') is None, \
        "Steel_t should not be visible outside its flowsheet"

    with qs.Flowsheet('_test_sys_b') as fs_b:
        assert qs.ImpactItem.get_item('Steel_t') is None, \
            "Steel_t from sys_a should not leak into sys_b"
        steel_b = qs.ImpactItem('Steel_t', 'kg', GlobalWarming_t=3.0)
        assert qs.ImpactItem.get_item('Steel_t') is steel_b

    # ── 5. ImpactIndicators are shared across all flowsheets ─────────────────
    with qs.Flowsheet('_test_sys_c') as fs_c:
        assert qs.ImpactIndicator.get_indicator('GlobalWarming_t') is GWP, \
            "ImpactIndicator should be visible in any flowsheet"

    assert qs.ImpactIndicator.get_indicator('GlobalWarming_t') is GWP, \
        "ImpactIndicator should persist after leaving flowsheet"

    # ── 6. clear() detaches StreamImpactItem from its stream ─────────────────
    components = qs.Components.load_default()
    qs.set_thermo(components)

    with qs.Flowsheet('_test_sys_d') as fs_d:
        s = qs.SanStream('_test_stream_d', H2O=100, units='kg/hr')
        item_d = qs.StreamImpactItem(linked_stream=s, GlobalWarming_t=1.5)
        assert s.stream_impact_item is item_d
        fs_d.clear()
        assert s.stream_impact_item is None, \
            "clear() should unlink StreamImpactItem from stream"
        assert qs.ImpactItem.get_item(item_d.ID) is None, \
            "clear() should remove item from registry"

    # ── 7. clear() detaches Construction from its unit ───────────────────────
    with qs.Flowsheet('_test_sys_e') as fs_e:
        concrete = qs.ImpactItem('Concrete_t', 'kg', GlobalWarming_t=0.1)
        M = qs.unit_operations.MixTank('_test_M_e', ins=qs.WasteStream(H2O=100, units='kg/hr'))
        c = qs.Construction(item=concrete, quantity=50, quantity_unit='kg',
                            linked_unit=M)
        M.construction = [c]
        assert len(M.construction) == 1
        fs_e.clear()
        assert M.construction == [], \
            "clear() should empty unit.construction"
        assert qs.ImpactItem.get_item('Concrete_t') is None

    # Cleanup global indicator used in this test
    GWP.deregister()


if __name__ == '__main__':
    test_flowsheet()
```

- [ ] **Step 2: Run the tests to confirm they fail for the right reason**

```
cd QSDsan
python -m pytest tests/test_flowsheet.py -v 2>&1 | head -40
```

Expected: `ImportError: cannot import name 'SanFlowsheet' from 'qsdsan._flowsheet'` (module does not exist yet).

---

## Task 2: Implement `_flowsheet.py`

**Files:**
- Create: `QSDsan/qsdsan/_flowsheet.py`

- [ ] **Step 1: Create the file**

```python
# qsdsan/_flowsheet.py
import biosteam as bst
from thermosteam import AbstractStream
from thermosteam.utils import Registry
from biosteam._unit import AbstractUnit
from biosteam._system import System

__all__ = ('SanFlowsheet', 'SanMainFlowsheet')


class SanFlowsheet(bst.Flowsheet):
    """
    Flowsheet that additionally owns per-flowsheet registries for
    ImpactItem, Construction, and Transportation.

    ImpactIndicator is intentionally NOT included — it is global and
    shared across all flowsheets so that indicator definitions stay
    consistent across systems.
    """

    def __new__(cls, ID=None):
        # super().__new__ registers self in FlowsheetRegistry, triggering
        # __setattr__ protection. Use object.__setattr__ for new attrs.
        self = super().__new__(cls, ID)
        object.__setattr__(self, 'item',           Registry())
        object.__setattr__(self, 'construction',   Registry())
        object.__setattr__(self, 'transportation', Registry())
        return self

    @classmethod
    def from_registries(cls, ID, stream, unit, system,
                        item=None, construction=None, transportation=None):
        # Call parent which uses object.__new__ + manual attr setting,
        # then adds self to FlowsheetRegistry at the end.
        self = bst.Flowsheet.from_registries(ID, stream, unit, system)
        # self is now registered, so use object.__setattr__.
        object.__setattr__(self, 'item',           item          or Registry())
        object.__setattr__(self, 'construction',   construction  or Registry())
        object.__setattr__(self, 'transportation', transportation or Registry())
        # Fix the class since parent used bst.Flowsheet, not cls.
        object.__setattr__(self, '__class__', cls)
        return self

    @property
    def registries(self):
        return (self.stream, self.unit, self.system,
                self.item, self.construction, self.transportation)

    def clear(self, reset_ticket_numbers=True):
        # Detach StreamImpactItems from their streams before clearing.
        for lca_item in list(self.item.data.values()):
            linked = getattr(lca_item, '_linked_stream', None)
            if linked is not None:
                linked._stream_impact_item = None
                object.__setattr__(lca_item, '_linked_stream', None)
        # Detach Construction objects from their units.
        for c in list(self.construction.data.values()):
            if c._linked_unit is not None:
                c._linked_unit._construction = []
        # Detach Transportation objects from their units.
        for t in list(self.transportation.data.values()):
            if t._linked_unit is not None:
                t._linked_unit._transportation = []
        super().clear(reset_ticket_numbers)
        self.item.clear()
        self.construction.clear()
        self.transportation.clear()


class SanMainFlowsheet(SanFlowsheet, bst.MainFlowsheet):
    """
    The active QSDsan flowsheet. Swaps ImpactItem, Construction, and
    Transportation class-level registries in addition to the BioSTEAM
    stream/unit/system registries whenever the active flowsheet changes.
    """

    __slots__ = ()

    def set_flowsheet(self, flowsheet, new=False):
        # Resolve string to a SanFlowsheet (not plain bst.Flowsheet).
        if isinstance(flowsheet, str):
            if not new and flowsheet in self.flowsheet:
                flowsheet = self.flowsheet[flowsheet]
            else:
                flowsheet = SanFlowsheet(flowsheet)
        # Delegate stream/unit/system registry swap to BioSTEAM.
        bst.MainFlowsheet.set_flowsheet(self, flowsheet)
        # After the call above, self.__dict__ IS flowsheet.__dict__.
        dct = self.__dict__
        if 'item' not in dct:
            return  # plain bst.Flowsheet was passed — skip LCA swap
        # Lazy imports: _flowsheet.py is loaded before these classes exist.
        from . import ImpactItem, Construction, Transportation
        ImpactItem.registry     = dct['item']
        Construction.registry   = dct['construction']
        Transportation.registry = dct['transportation']
```

- [ ] **Step 2: Run the tests — expect a different failure now**

```
python -m pytest tests/test_flowsheet.py -v 2>&1 | head -40
```

Expected: tests import `SanFlowsheet` successfully but fail on `qs.Flowsheet is SanFlowsheet` (not yet wired up in `__init__.py`).

---

## Task 3: Wire up the new flowsheet in `__init__.py`

**Files:**
- Modify: `QSDsan/qsdsan/__init__.py`

- [ ] **Step 1: Read the current `__init__.py` to understand the exact lines to change**

Current relevant lines (approximately):
```python
Flowsheet = _bst.Flowsheet          # line ~45
main_flowsheet = _bst.main_flowsheet # line ~46
...
def default():
    _bst.default()
    utils.clear_lca_registries()    # line ~117
```

- [ ] **Step 2: Replace `Flowsheet` and `main_flowsheet` assignments and update `default()`**

Remove:
```python
Flowsheet = _bst.Flowsheet
main_flowsheet = _bst.main_flowsheet
```

Add near the top (after `import biosteam as _bst`, before `from ._component import *`):
```python
from ._flowsheet import SanFlowsheet, SanMainFlowsheet
```

At the **bottom** of `__init__.py`, after all LCA class imports (`from ._lca import *` etc.), add the singleton creation block:

```python
# ── QSDsan main flowsheet (must be after all LCA classes are imported) ──────
from biosteam._unit import AbstractUnit as _AbstractUnit
from thermosteam import AbstractStream as _AbstractStream

_qs_default_fs = SanFlowsheet.from_registries(
    'default',
    _AbstractStream.registry,
    _AbstractUnit.registry,
    _bst.System.registry,
    # item/construction/transportation start empty — LCA registries are clean
)
main_flowsheet = object.__new__(SanMainFlowsheet)
_bst.MainFlowsheet.set_flowsheet(main_flowsheet, _qs_default_fs)
# Now swap LCA registries for the default flowsheet.
from . import ImpactItem as _ImpactItem, Construction as _Construction
from . import Transportation as _Transportation
_ImpactItem.registry     = _qs_default_fs.item
_Construction.registry   = _qs_default_fs.construction
_Transportation.registry = _qs_default_fs.transportation

Flowsheet = SanFlowsheet
F = main_flowsheet
del _qs_default_fs, _AbstractUnit, _AbstractStream
del _ImpactItem, _Construction, _Transportation
```

Update `default()`:
```python
def default():
    _bst.default()          # resets thermo; switches bst.main_flowsheet to 'default'
    main_flowsheet.set_flowsheet('default', new=True)  # fresh default flowsheet
```

Update `__all__` to include `'Flowsheet'`, `'main_flowsheet'`, `'F'` if not already there (check first).

- [ ] **Step 3: Run the flowsheet tests — all should pass now**

```
python -m pytest tests/test_flowsheet.py -v
```

Expected: all tests PASS.

- [ ] **Step 4: Run the existing test suite to check for regressions**

```
python -m pytest tests/ -v --ignore=tests/test_exposan.py -x 2>&1 | tail -30
```

Expected: all existing tests PASS (no regressions from flowsheet change).

- [ ] **Step 5: Commit**

```bash
git add qsdsan/_flowsheet.py qsdsan/__init__.py tests/test_flowsheet.py
git commit -m "feat: add SanFlowsheet with per-flowsheet LCA registries

ImpactItem, Construction, and Transportation now register to the active
flowsheet instead of a global singleton. ImpactIndicator remains global.
Switching systems via 'with qs.Flowsheet(\"sysB\"):' automatically isolates
all LCA state without manual clear_lca_registries() calls.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Task 4: Write failing tests for `_construction_specs`

**Files:**
- Modify: `QSDsan/tests/test_flowsheet.py`

- [ ] **Step 1: Append `test_construction_specs` to `test_flowsheet.py`**

```python
def test_construction_specs():
    import qsdsan as qs
    from qsdsan.unit_operations import MixTank

    components = qs.Components.load_default()
    qs.set_thermo(components)

    # ── Define a unit subclass with default construction specs ────────────────
    class SpecTank(MixTank):
        # Each spec: dict with keys item, quantity, quantity_unit,
        # and optionally lifetime, lifetime_unit.
        _construction_specs = (
            dict(item='SpecSteel', quantity=10., quantity_unit='kg',
                 lifetime=20., lifetime_unit='yr'),
            dict(item='SpecConcrete', quantity=5., quantity_unit='m3'),
        )

    # ── 1. LCA raises clearly when required items are not loaded ─────────────
    with qs.Flowsheet('_test_spec_a') as fs_a:
        feed = qs.WasteStream('_spec_feed_a', H2O=1000, units='kg/hr')
        tank = SpecTank('_spec_tank_a', ins=feed)
        sys = qs.System('_spec_sys_a', path=(tank,))
        sys.simulate()

        GWP = qs.ImpactIndicator('GWP_spec', unit='kg CO2-eq')
        # Item NOT loaded — LCA should raise
        import pytest
        with pytest.raises(RuntimeError, match="SpecSteel"):
            lca = qs.LCA(sys, lifetime=20, indicators=(GWP,),
                         simulate_system=False)

    # ── 2. LCA resolves specs when items are loaded ───────────────────────────
    with qs.Flowsheet('_test_spec_b') as fs_b:
        feed = qs.WasteStream('_spec_feed_b', H2O=1000, units='kg/hr')
        tank = SpecTank('_spec_tank_b', ins=feed)
        sys = qs.System('_spec_sys_b', path=(tank,))
        sys.simulate()

        GWP = qs.ImpactIndicator('GWP_spec', unit='kg CO2-eq')
        steel   = qs.ImpactItem('SpecSteel',   'kg', GWP_spec=2.55)
        concrete= qs.ImpactItem('SpecConcrete','m3', GWP_spec=100.)

        lca = qs.LCA(sys, lifetime=20, indicators=(GWP,),
                     simulate_system=False)

        # Spec-derived constructions appear in LCA construction units
        assert tank in lca.construction_units, \
            "unit with _construction_specs should appear in lca.construction_units"

        impacts = lca.get_construction_impacts()
        # 10 kg steel × 2.55 kg CO2-eq/kg = 25.5; 5 m3 concrete × 100 = 500
        from numpy.testing import assert_allclose
        assert_allclose(impacts['GWP_spec'], 25.5 + 500., rtol=1e-6)

    # ── 3. Explicit unit.construction overrides spec for same item ────────────
    with qs.Flowsheet('_test_spec_c') as fs_c:
        feed = qs.WasteStream('_spec_feed_c', H2O=1000, units='kg/hr')
        tank = SpecTank('_spec_tank_c', ins=feed)

        GWP    = qs.ImpactIndicator('GWP_spec', unit='kg CO2-eq')
        steel  = qs.ImpactItem('SpecSteel',   'kg', GWP_spec=2.55)
        conc   = qs.ImpactItem('SpecConcrete','m3', GWP_spec=100.)

        # User provides explicit construction for SpecSteel → spec is skipped
        explicit = qs.Construction(item=steel, quantity=99., quantity_unit='kg',
                                   linked_unit=tank)
        tank.construction = [explicit]

        sys = qs.System('_spec_sys_c', path=(tank,))
        sys.simulate()
        lca = qs.LCA(sys, lifetime=20, indicators=(GWP,), simulate_system=False)

        impacts = lca.get_construction_impacts()
        # 99 kg steel × 2.55 = 252.45; 5 m3 concrete × 100 = 500
        assert_allclose(impacts['GWP_spec'], 252.45 + 500., rtol=1e-6)

    # ── 4. _construction_specs defaults to empty tuple on base SanUnit ────────
    assert qs.SanUnit._construction_specs == (), \
        "SanUnit._construction_specs should default to empty tuple"

    # Cleanup
    qs.ImpactIndicator.get_indicator('GWP_spec').deregister()


if __name__ == '__main__':
    test_flowsheet()
    test_construction_specs()
```

- [ ] **Step 2: Run to confirm it fails for the right reason**

```
python -m pytest tests/test_flowsheet.py::test_construction_specs -v 2>&1 | head -20
```

Expected: `AttributeError: type object 'SanUnit' has no attribute '_construction_specs'`

---

## Task 5: Implement `_construction_specs` on `SanUnit` and resolution in `LCA`

**Files:**
- Modify: `QSDsan/qsdsan/_sanunit.py:173`
- Modify: `QSDsan/qsdsan/_lca.py:254`

- [ ] **Step 1: Add `_construction_specs` class attribute to `SanUnit`**

In [_sanunit.py](QSDsan/qsdsan/_sanunit.py), find the line:
```python
    _init_lca = AbstractMethod
```
(currently around line 173). Add the class attribute immediately above it:

```python
    # Subclasses declare default material usage as a tuple of dicts.
    # Keys: item (str, ImpactItem ID), quantity (float), quantity_unit (str),
    #       lifetime (float, optional), lifetime_unit (str, optional, default 'yr').
    # LCA resolves these at LCA creation time — items must exist in the active
    # flowsheet by then. Override at instance level to change quantities,
    # or set unit.construction directly to bypass specs entirely.
    _construction_specs = ()

    _init_lca = AbstractMethod
```

- [ ] **Step 2: Add `_resolve_construction_specs()` to `LCA` and call it from `_update_system()`**

In [_lca.py](QSDsan/qsdsan/_lca.py), find `_update_system()` (around line 254). Add a new method `_resolve_construction_specs` right after `_update_system`, then call it from `_update_system`.

The full replacement for `_update_system` and the new method:

```python
    def _update_system(self, system):
        for u in system.units:
            if not isinstance(u, SanUnit):
                continue
            if getattr(u, 'construction', []):
                self._construction_units.add(u)
            if getattr(u, 'transportation', []):
                self._transportation_units.add(u)
        self._resolve_construction_specs(system)  # new call
        self._construction_units = sorted(self._construction_units,
                                          key=lambda u: u.ID)
        self._transportation_units = sorted(self._transportation_units,
                                            key=lambda u: u.ID)
        for s in (i for i in system.feeds+system.products):
            if not hasattr(s, 'stream_impact_item'):
                continue
            if s.stream_impact_item:
                self._lca_streams.add(s)
        self._lca_streams = sorted(self._lca_streams, key=lambda s: s.ID)
        self._system = system
        try:
            system._LCA = self
        except AttributeError:
            pass

    def _resolve_construction_specs(self, system):
        """Create Construction objects from SanUnit._construction_specs.

        Called during LCA initialisation. Items must already exist in the
        active flowsheet's ImpactItem registry; raises RuntimeError if not.
        Skips any spec whose item ID already appears in unit.construction
        so that explicit user-defined constructions take precedence.
        """
        from . import Construction as _Construction
        for u in system.units:
            if not isinstance(u, SanUnit):
                continue
            specs = getattr(u, '_construction_specs', ())
            if not specs:
                continue
            existing_ids = {c.item.ID for c in u.construction if c.item}
            added = []
            for spec in specs:
                item_id = spec['item']
                if item_id in existing_ids:
                    continue  # user-defined construction takes precedence
                item = ImpactItem.get_item(item_id)
                if item is None:
                    raise RuntimeError(
                        f"ImpactItem '{item_id}' is listed in "
                        f"{type(u).__name__}._construction_specs but was not "
                        f"found in the current flowsheet. Load it before "
                        f"creating the LCA, or set {u.ID}.construction "
                        f"explicitly to override the spec."
                    )
                c = _Construction(
                    ID=f'{u.ID}_{item_id}',
                    linked_unit=u,
                    item=item,
                    quantity=spec['quantity'],
                    quantity_unit=spec.get('quantity_unit', ''),
                    lifetime=spec.get('lifetime'),
                    lifetime_unit=spec.get('lifetime_unit', 'yr'),
                )
                added.append(c)
            if added:
                u.construction = list(u.construction) + added
                self._construction_units.add(u)
```

- [ ] **Step 3: Remove the stale `clear_lca_registries` alias from `LCA`**

In [_lca.py](QSDsan/qsdsan/_lca.py), find and remove:
```python
    clear_lca_registries = clear_lca_registries
```
(currently around line 252, just before `_update_system`). Remove this line entirely — it assigns a module-level function as a class attribute, which is no longer meaningful given the new flowsheet-based approach.

- [ ] **Step 4: Run the construction specs tests**

```
python -m pytest tests/test_flowsheet.py::test_construction_specs -v
```

Expected: all assertions PASS.

- [ ] **Step 5: Run the full test suite**

```
python -m pytest tests/ -v --ignore=tests/test_exposan.py -x 2>&1 | tail -30
```

Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add qsdsan/_sanunit.py qsdsan/_lca.py tests/test_flowsheet.py
git commit -m "feat: add _construction_specs to SanUnit; resolve in LCA._update_system

Units declare default material usage as class-level _construction_specs
(tuple of dicts). LCA resolves specs at creation time, not at unit init
time, so ImpactItem objects don't need to exist when a unit is defined.
Explicit unit.construction entries always take precedence over specs.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Task 6: Deprecate `clear_lca_registries`

**Files:**
- Modify: `QSDsan/qsdsan/utils/misc.py:28`

- [ ] **Step 1: Add `DeprecationWarning` to `clear_lca_registries`**

Replace the body of `clear_lca_registries` in [utils/misc.py](QSDsan/qsdsan/utils/misc.py):

```python
def clear_lca_registries(print_msg=False):
    '''
    Clear registries related to LCA.

    .. deprecated::
        Use ``with qs.Flowsheet('new_system'):`` to isolate LCA state per
        system. Switching flowsheets via the context manager automatically
        clears ImpactItem, Construction, and Transportation registries without
        needing this function. Call ``flowsheet.clear()`` to explicitly reset
        a flowsheet's LCA state.
    '''
    import warnings
    warnings.warn(
        "clear_lca_registries() is deprecated. Use 'with qs.Flowsheet(\"name\"):' "
        "to isolate LCA state per system, or call flowsheet.clear() explicitly.",
        DeprecationWarning,
        stacklevel=2,
    )
    from qsdsan import ImpactIndicator, ImpactItem, Construction, Transportation
    for lca_cls in (ImpactIndicator, ImpactItem, Construction, Transportation):
        lca_cls.clear_registry(print_msg)
```

- [ ] **Step 2: Run the doctests that currently call `clear_lca_registries()` to confirm they still pass** (they call the function so they now emit warnings, but should still pass)

```
python -m pytest --doctest-modules qsdsan/_impact_indicator.py qsdsan/_impact_item.py qsdsan/_construction.py qsdsan/_lca.py -v 2>&1 | tail -20
```

Expected: PASS (warnings are emitted but doctests don't fail on warnings by default).

- [ ] **Step 3: Update doctests that call `clear_lca_registries()` to suppress the deprecation warning**

In each doctest that ends with:
```python
>>> from qsdsan.utils import clear_lca_registries
>>> clear_lca_registries()
```

Replace with (using `# doctest: +SKIP` since these are just cleanup calls, not part of the tested API):
```python
>>> # doctest cleanup — normally handled by flowsheet context manager
>>> import warnings; warnings.filterwarnings('ignore', category=DeprecationWarning)
>>> from qsdsan.utils import clear_lca_registries
>>> clear_lca_registries()
```

Files to update: [_impact_indicator.py](QSDsan/qsdsan/_impact_indicator.py), [_impact_item.py](QSDsan/qsdsan/_impact_item.py), [_construction.py](QSDsan/qsdsan/_construction.py), [_lca.py](QSDsan/qsdsan/_lca.py), [utils/misc.py](QSDsan/qsdsan/utils/misc.py).

- [ ] **Step 4: Run doctests again**

```
python -m pytest --doctest-modules qsdsan/_impact_indicator.py qsdsan/_impact_item.py qsdsan/_construction.py qsdsan/_lca.py qsdsan/utils/misc.py -v 2>&1 | tail -20
```

Expected: all PASS with no unexpected output.

- [ ] **Step 5: Commit**

```bash
git add qsdsan/utils/misc.py qsdsan/_impact_indicator.py qsdsan/_impact_item.py qsdsan/_construction.py qsdsan/_lca.py
git commit -m "deprecate: clear_lca_registries() — use Flowsheet context manager instead

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Task 7: Integration test — switching between two full systems

**Files:**
- Modify: `QSDsan/tests/test_flowsheet.py`

- [ ] **Step 1: Append `test_two_system_switch` to `test_flowsheet.py`**

```python
def test_two_system_switch():
    """Verify that switching between two LCA-enabled systems produces no
    cross-contamination in registries or unit state."""
    import qsdsan as qs
    from numpy.testing import assert_allclose

    components = qs.Components.load_default()
    qs.set_thermo(components)

    GWP = qs.ImpactIndicator('GWP_switch', unit='kg CO2-eq')

    # ── System A ──────────────────────────────────────────────────────────────
    with qs.Flowsheet('_switch_sys_a') as fs_a:
        feed_a = qs.WasteStream('_sw_feed_a', H2O=1000, units='kg/hr')
        M_a    = qs.unit_operations.MixTank('_sw_M_a', ins=feed_a)
        sys_a  = qs.System('_sw_sys_a', path=(M_a,))
        sys_a.simulate()

        steel_a   = qs.ImpactItem('SwitchSteel',   'kg', GWP_switch=2.0)
        constr_a  = qs.Construction(item=steel_a, quantity=10., quantity_unit='kg',
                                    linked_unit=M_a)
        M_a.construction = [constr_a]

        lca_a = qs.LCA(sys_a, lifetime=10, indicators=(GWP,),
                       simulate_system=False)
        impacts_a = lca_a.get_construction_impacts()
        # 10 kg × 2.0 = 20
        assert_allclose(impacts_a['GWP_switch'], 20., rtol=1e-6)

    # ── After exiting sys_a: SwitchSteel should not be visible ───────────────
    assert qs.ImpactItem.get_item('SwitchSteel') is None, \
        "SwitchSteel from sys_a should not be visible after exiting its flowsheet"

    # ── System B (same item name, different CF) ───────────────────────────────
    with qs.Flowsheet('_switch_sys_b') as fs_b:
        assert qs.ImpactItem.get_item('SwitchSteel') is None, \
            "SwitchSteel should not leak from sys_a into sys_b"

        feed_b = qs.WasteStream('_sw_feed_b', H2O=500, units='kg/hr')
        M_b    = qs.unit_operations.MixTank('_sw_M_b', ins=feed_b)
        sys_b  = qs.System('_sw_sys_b', path=(M_b,))
        sys_b.simulate()

        # Same name, different CF — no conflict warning expected
        steel_b  = qs.ImpactItem('SwitchSteel', 'kg', GWP_switch=5.0)
        constr_b = qs.Construction(item=steel_b, quantity=3., quantity_unit='kg',
                                   linked_unit=M_b)
        M_b.construction = [constr_b]

        lca_b = qs.LCA(sys_b, lifetime=10, indicators=(GWP,),
                       simulate_system=False)
        impacts_b = lca_b.get_construction_impacts()
        # 3 kg × 5.0 = 15
        assert_allclose(impacts_b['GWP_switch'], 15., rtol=1e-6)

    # ── sys_a impacts are unchanged ───────────────────────────────────────────
    # Re-enter sys_a flowsheet and verify its state is intact
    with qs.Flowsheet('_switch_sys_a'):
        assert_allclose(
            lca_a.get_construction_impacts()['GWP_switch'], 20., rtol=1e-6,
            err_msg="sys_a LCA impacts should be unchanged after sys_b was created"
        )

    # Cleanup
    GWP.deregister()


if __name__ == '__main__':
    test_flowsheet()
    test_construction_specs()
    test_two_system_switch()
```

- [ ] **Step 2: Run the integration test**

```
python -m pytest tests/test_flowsheet.py::test_two_system_switch -v
```

Expected: PASS.

- [ ] **Step 3: Run the complete test suite one final time**

```
python -m pytest tests/ -v --ignore=tests/test_exposan.py 2>&1 | tail -40
```

Expected: all PASS with no regressions.

- [ ] **Step 4: Final commit**

```bash
git add tests/test_flowsheet.py
git commit -m "test: add integration test for two-system LCA switch

Verifies that switching between two systems using flowsheet context managers
produces no cross-contamination in ImpactItem registries or unit construction
state, and that the same item name can be used in different systems with
different CFs without conflict.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Self-Review

### Spec coverage check

| Requirement | Covered by |
|---|---|
| Per-flowsheet ImpactItem/Construction/Transportation registries | Tasks 1–3 |
| ImpactIndicator remains global across flowsheets | Task 1 test 5, Task 2 (no indicator registry) |
| Context manager isolates LCA state | Task 1 test 4, Task 7 |
| `clear()` detaches objects before clearing | Task 1 tests 6–7, Task 2 `SanFlowsheet.clear()` |
| `_construction_specs` on `SanUnit` | Tasks 4–5 |
| Spec resolution raises clearly when items missing | Task 4 test 1 |
| Explicit `unit.construction` overrides specs | Task 4 test 3 |
| `clear_lca_registries` deprecated | Task 6 |
| No regressions in existing tests | Tasks 3, 5, 7 |
| `qs.Flowsheet` and `qs.main_flowsheet` updated | Task 3 |

### Placeholder scan

None found — all steps contain complete code.

### Type consistency

- `SanFlowsheet` defined in Task 2, referenced as `qs.Flowsheet` from Task 3 onward ✓
- `SanMainFlowsheet` defined in Task 2, checked via `isinstance` in Task 1 ✓
- `_construction_specs` added to `SanUnit` in Task 5, used via `getattr(u, '_construction_specs', ())` in the same task ✓
- `lca.construction_units` used in Task 4 test 2 — confirm property name in `_lca.py` (it is `construction_units` as a property returning `self._construction_units`) ✓
- `lca.get_construction_impacts()` returns `{indicator_ID: float}` dict — used consistently ✓
