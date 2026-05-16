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

## Change Checklist

1. Check both repos for consumers before changing public names.
2. Add or update QSDsan tests for engine behavior and namespace exports.
3. Add or update EXPOsan tests only for system-level behavior affected by the change.
4. Keep import migrations mechanical and public: prefer `from qsdsan import sanunits as su` or `from qsdsan.sanunits import static, dynamic, bst`.
5. Do not copy BioSTEAM/ThermoSTEAM implementations into QSDsan unless QSDsan adds meaningful behavior.
6. Record unresolved classification or migration decisions in the relevant design/spec doc.
7. When adding LCA objects to a new unit class, use `_construction_specs` instead of creating `ImpactItem`/`Construction` objects in `__init__`. This avoids requiring items to be pre-loaded at unit creation time.
8. **Always update documentation alongside code changes.** For every changed parameter, behavior, or default: update the class/function docstring `Parameters` section, any inline examples that reference old values, and property docstrings. Stale docs (e.g., old default values, removed IDs) are treated as bugs.
