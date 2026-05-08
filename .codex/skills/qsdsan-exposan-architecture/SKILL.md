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

- `bst`: BioSTEAM-compatible wrappers with minimal QSDsan integration.
- `static`: steady-state sanitation/resource-recovery units with QSDsan design, costing, construction, TEA, LCA, or WasteStream assumptions.
- `dynamic`: units with explicit dynamic state contracts such as `isdynamic`, `state`, `dstate`, `_init_state`, `_compile_AE`, or `_compile_ODE`.

Preserve legacy imports like `qsdsan.sanunits.Pump` while adding clearer namespaces. Prefer re-export modules first; move implementation files only after exports and tests are stable.

## Change Checklist

1. Check both repos for consumers before changing public names.
2. Add or update QSDsan tests for engine behavior and namespace exports.
3. Add or update EXPOsan tests only for system-level behavior affected by the change.
4. Keep import migrations mechanical and public: prefer `from qsdsan import sanunits as su` or `from qsdsan.sanunits import static, dynamic, bst`.
5. Do not copy BioSTEAM/ThermoSTEAM implementations into QSDsan unless QSDsan adds meaningful behavior.
6. Record unresolved classification or migration decisions in the relevant design/spec doc.
