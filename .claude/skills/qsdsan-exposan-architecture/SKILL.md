---
name: qsdsan-exposan-architecture
description: Use when changing QSDsan or EXPOsan structure, dependencies, imports, sanunit namespaces, dynamic unit behavior, examples, tests, or CI across the paired QSDsan/EXPOsan repositories.
---

# QSDsan/EXPOsan Architecture

QSDsan is the reusable modeling engine: components, streams, units, processes, systems, TEA, and LCA. EXPOsan is the applied systems catalog: published or user-built sanitation/resource-recovery systems that compose QSDsan primitives.

Read `references/package-map.md` before making cross-repository changes, dependency changes, `sanunits` namespace changes, or EXPOsan import migrations.

## Working Rules

| Situation | Rule |
| --- | --- |
| Generic modeling capability | Put it in QSDsan. |
| A specific case study, scenario, flowsheet, analysis, or dataset | Keep it in EXPOsan. |
| EXPOsan needs a reusable primitive | Add or expose the primitive in QSDsan, then update EXPOsan to consume the public API. |
| QSDsan tests need EXPOsan systems | Keep them in explicit integration coverage; do not make core QSDsan tests depend on EXPOsan. |
| Local editable development | Use sibling checkouts, not GitHub package dependencies. |
| CI integration | GitHub sibling dependencies are allowed through integration extras. |

## Unit Namespace Refactor

Classify `qsdsan.sanunits` by behavior, not by file location:

- `bst`: BioSTEAM-inherited unit operations with QSDsan-added behavior.
- `static`: steady-state sanitation/resource-recovery units with QSDsan design, costing, construction, TEA, LCA, or WasteStream assumptions.
- `dynamic`: units with explicit dynamic state contracts such as `isdynamic`, `state`, `dstate`, `_init_state`, `_compile_AE`, or `_compile_ODE`.

Preserve legacy imports like `qsdsan.sanunits.Pump` while adding clearer namespaces. Prefer re-export modules first; move implementation files only after exports and tests are stable.

## QSDsan Tutorial Docs

Do not rely on `docs/source/templates`; the old tutorial template was removed because it drifted from the maintained notebooks. When adding or revising tutorials, use the current notebooks in `docs/source/tutorials` as the source of truth.

Keep the common tutorial conventions unless there is a focused reason to differ: a top anchor, a short contents list with valid anchors when useful, and the maintained notebook style for contributor attribution. When adding a version check, use `qsdsan.__version__`. After tutorial edits, build the docs and check for `toc.not_included` and undefined-label warnings.

## LCA Registry Architecture

### Per-flowsheet vs. global registries

| Class | Registry scope | Swapped on `set_flowsheet`? |
|---|---|---|
| `ImpactIndicator` | Per-flowsheet (`flowsheet.indicator`) | Yes |
| `ImpactItem` | Per-flowsheet (`flowsheet.item`) | Yes |
| `Construction` | Per-flowsheet (`flowsheet.construction`) | Yes |
| `Transportation` | Per-flowsheet (`flowsheet.transportation`) | Yes |

`qs.Flowsheet` is `SanFlowsheet` (a subclass of BioSTEAM's `Flowsheet`). It adds four extra `Registry` attributes (`indicator`, `item`, `construction`, `transportation`) alongside the existing `stream`, `unit`, `system` ones.

`qs.main_flowsheet` is a `SanMainFlowsheet` instance. Its `set_flowsheet()` override swaps all seven registries (BioSTEAM's three plus the four LCA ones) atomically.

### Flowsheet context manager — the correct isolation pattern

```python
with qs.Flowsheet('sysA') as fs_a:
    steel = qs.ImpactItem('Steel', 'kg', GWP=2.55)   # lives in fs_a.item
    # build units, system, LCA here

# After exiting: qs.ImpactItem.get_item('Steel') is None
# Steel is only accessible while fs_a is the active flowsheet.

with qs.Flowsheet('sysB') as fs_b:
    # fs_b.item starts empty — no cross-contamination from sysA
    steel_b = qs.ImpactItem('Steel', 'kg', GWP=5.0)
    ...
```

### `_construction_specs` — declarative default materials

Units that always require certain construction materials can declare them as a class attribute instead of creating `Construction` objects in `__init__`. This avoids requiring `ImpactItem` objects to exist at unit-creation time.

```python
class ConcreteReactor(SanUnit):
    _construction_specs = (
        dict(item='Concrete', quantity=5., quantity_unit='m3'),
        dict(item='Steel',    quantity=10., quantity_unit='kg',
             lifetime=20.,    lifetime_unit='yr'),
    )
```

`LCA.__init__` resolves specs lazily: for each spec whose `item` ID is not already in `unit._construction`, it looks up the item in the current flowsheet and creates a `Construction` object. If the item is missing it raises `RuntimeError` with the missing item name. Explicit `unit.construction` entries take precedence over specs for the same item ID.

### `clear_lca_registries()` — deprecated

`qs.utils.clear_lca_registries()` now emits `DeprecationWarning`. The replacement patterns are:

| Old pattern | New pattern |
|---|---|
| `clear_lca_registries()` + global items | `with qs.Flowsheet('name'):` |
| `flowsheet.clear(); clear_lca_registries()` | `flowsheet.clear()` (handles LCA items too) |

### EXPOsan migration pattern

EXPOsan modules that create multiple sub-systems (A, B, C…) with separate flowsheets must reload LCA items in each flowsheet's context. The standard pattern:

```python
def create_system(system_ID='A', flowsheet=None):
    if flowsheet is None:
        flowsheet_ID = f'br{system_ID}'
        if hasattr(main_flowsheet.flowsheet, flowsheet_ID):
            getattr(main_flowsheet.flowsheet, flowsheet_ID).clear()
        flowsheet = Flowsheet(flowsheet_ID)
        main_flowsheet.set_flowsheet(flowsheet)
        reload_lca = True   # ← always True: items are per-flowsheet
    else:
        reload_lca = False
    _load_lca_data(reload_lca)
    ...
```

Remove any `clear_lca_registries()` calls that preceded the `flowsheet.clear()` call — they are now redundant because `SanFlowsheet.clear()` already detaches `StreamImpactItem` objects from their streams and clears `item`, `construction`, and `transportation` registries.

## Release Conventions

- Always bump QSDsan and EXPOsan versions together — they are released as a paired set.
- Update `CHANGELOG.rst` in both repos before tagging. The tag triggers the release workflow; there is no post-tag opportunity to amend the changelog.
- Tag format: `v*.*.*` (e.g. `v1.4.5`). The workflow verifies the tag matches `pyproject.toml` version and will fail if they differ.

## Change Checklist

1. Check both repos for consumers before changing public names.
2. Add or update QSDsan tests for engine behavior and namespace exports.
3. Add or update EXPOsan tests only for system-level behavior affected by the change.
4. Keep import migrations mechanical and public: prefer `from qsdsan import sanunits as su` or `from qsdsan.sanunits import static, dynamic, bst`.
5. Do not copy BioSTEAM/ThermoSTEAM implementations into QSDsan unless QSDsan adds meaningful behavior.
6. Record unresolved classification or migration decisions in the relevant design/spec doc.
7. When adding LCA objects to a new unit class, use `_construction_specs` instead of creating `ImpactItem`/`Construction` objects in `__init__`. This avoids requiring items to be pre-loaded at unit creation time.
8. **Always update documentation alongside code changes.** For every changed parameter, behavior, or default: update the class/function docstring `Parameters` section, any inline examples that reference old values, and property docstrings. Stale docs (e.g., old default values, removed IDs) are treated as bugs.
