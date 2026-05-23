---
name: qsdsan-exposan-architecture
description: Use when changing QSDsan or EXPOsan structure, dependencies, imports, unit/process namespaces, dynamic unit behavior, tutorial/docs authoring, examples, tests, or CI across the paired QSDsan/EXPOsan repositories.
---

# QSDsan/EXPOsan Architecture

QSDsan is the reusable modeling engine: components, streams, units, processes, systems, TEA, and LCA. EXPOsan is the applied systems catalog: published or user-built sanitation/resource-recovery systems that compose QSDsan primitives.

Read `references/package-map.md` before making cross-repository changes, dependency changes, unit/process namespace work, or EXPOsan import changes.

## Working Rules

| Situation | Rule |
| --- | --- |
| Generic modeling capability | Put it in QSDsan. |
| A specific case study, scenario, flowsheet, analysis, or dataset | Keep it in EXPOsan. |
| EXPOsan needs a reusable primitive | Add or expose the primitive in QSDsan, then update EXPOsan to consume the public API. |
| QSDsan tests need EXPOsan systems | Keep them in explicit integration coverage; do not make core QSDsan tests depend on EXPOsan. |
| Local editable development | Use sibling checkouts, not GitHub package dependencies. |
| CI integration | GitHub sibling dependencies are allowed through integration extras. |

## Unit & Process Namespaces

Built-in units live in `qsdsan.unit_operations`, split into three behavior-based sub-namespaces:

- `bst`: BioSTEAM-inherited unit operations with QSDsan-added behavior.
- `static`: steady-state sanitation/resource-recovery units with QSDsan design, costing, construction, TEA, LCA, or WasteStream assumptions.
- `dynamic`: units with explicit dynamic state contracts such as `isdynamic`, `state`, `dstate`, `_init_state`, `_compile_AE`, or `_compile_ODE`.

Process models live in `qsdsan.process_models`. `qsdsan.sanunits` and `qsdsan.processes` remain only as back-compat aliases (resolving to `unit_operations` and `process_models`); use the canonical names in new code. Classify a new unit by behavior, not file location, into `bst`/`static`/`dynamic`, and add a test that the namespace exports it.

## QSDsan Tutorial Docs

The notebooks in `docs/source/tutorials` are the source of truth — match their conventions when adding or revising one. (There is no separate tutorial template; do not add one.)

### Notebook skeleton (cells in order)

1. **Title cell:** `# <Title> <a class="anchor" id="top"></a>`, then the Binder badge, the Colab note, **Prepared by**, **Learning objectives** (~3 bullets), **Prerequisites** (links to the tutorials it builds on), **Covered topics** (a short TOC of `<a href="#sN">N. ...</a>` links), and a **Companion video** blockquote stating which `QSDsan` version it was recorded against.
2. **Setup:** a markdown cell beginning with `<!-- tutorial-setup-section -->` and a `## Setup` heading, then a code cell: `import qsdsan as qs` + `print(f'This tutorial was made with qsdsan v{qs.__version__}.')`. Tutorials that also use EXPOsan print both versions (`import qsdsan as qs, exposan` -> `... v{qs.__version__} and exposan v{exposan.__version__}`).
3. **Body sections.**
4. **Nav footer (last cell):** `<!-- tutorial-nav-footer -->`, a `---` rule, then `<a href="#top">↑ Back to top</a>`.

### Headings & navigation

- Sections are `## N. <Title> <a class="anchor" id="sN"></a>`; subsections `### N.M.`; sub-subsections `#### N.M.K.`. Numbering must be contiguous and consistent — never leave an unnumbered `####` among numbered siblings.
- **Sentence-case** all headings.
- The only "back to top" link is the nav-footer one; do not sprinkle per-section `Back to top` links (furo already provides a sidebar TOC and a floating back-to-top button). No horizontal rules between headings.

### Prose & style

- Address the reader as **"users"**, not "students".
- **Minimize em-dashes**; prefer parentheses or colons.
- Use `<div class="alert alert-info">` admonitions for **Note**/**Tip**/**Heads up** callouts, and collapsible `<details><summary><i>Python Aside: ... (click to expand)</i></summary>` blocks for Python-language tangents (the tutorials index advertises these).
- Don't paste a class's full signature/parameter docstring into markdown — it silently drifts. Describe the few relevant parameters and point readers to `?ClassName` or the API docs.
- Give every stream/unit an explicit ID matching its variable name (`ww1 = qs.WasteStream('ww1', ...)`); reusing an ID across cells triggers "replaced in registry" warnings, and an auto-assigned ID (e.g. `ws1`) won't match the variable a reader sees.
- When first showing `.diagram()`, use the default format and mention that `format='html'` gives an interactive diagram with hover-able stream/unit info.

### Current API names

Use the canonical names, not the legacy aliases: `unit_operations` (not `sanunits`), `process_models` (not `processes`), `TEA` (not `SimpleTEA`), `qs.CEPCI` / `qs.CEPCI_by_year`. Where it aids understanding, note the alias once (e.g., "`qsdsan.sanunits` is a legacy alias for `qsdsan.unit_operations`"). Prefer the `qs.unit_operations.<Name>` / `qs.process_models.<Name>` namespaces.

### Execution, outputs & anchors

- `nbsphinx_execute` defaults to `'auto'`: a notebook with **no** stored outputs is re-executed at build time (so it must run end-to-end), while one **with** stored outputs is shown as-is (its "made with v..." line and numbers reflect whichever version last ran it — re-run to refresh). `nbsphinx_allow_errors = True`, so deliberate teaching errors are fine, but every *unintentional* error must be cleared.
- Cross-tutorial anchors: nbsphinx slugs preserve case and dots, turn spaces into hyphens, and drop backticks (e.g. `4_SanUnit_basic.html#1.1.-SanUnit-and-unit_operations`); RST pages (e.g. `index.rst`) use lowercase Sphinx slugs. Verify a target against the built HTML rather than guessing.

### Editing & verification

- Edit notebooks with `NotebookEdit` or a content-matched script that re-dumps JSON with `indent=1, ensure_ascii=False` and a trailing newline; clear `outputs`/`execution_count` on code cells whose source changed.
- After edits, run the affected notebook(s) end-to-end (only intentional teaching errors should remain) and build the docs, checking for `toc.not_included` and undefined-label warnings.

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

### EXPOsan per-flowsheet LCA pattern

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
4. Prefer public imports: `from qsdsan import unit_operations as su` (and `process_models as pc`) or `from qsdsan.unit_operations import bst, static, dynamic`.
5. Do not copy BioSTEAM/ThermoSTEAM implementations into QSDsan unless QSDsan adds meaningful behavior.
6. Record unresolved classification decisions in the relevant design/spec doc.
7. When adding LCA objects to a new unit class, use `_construction_specs` instead of creating `ImpactItem`/`Construction` objects in `__init__`. This avoids requiring items to be pre-loaded at unit creation time.
8. **Always update documentation alongside code changes.** For every changed parameter, behavior, or default: update the class/function docstring `Parameters` section, any inline examples that reference old values, and property docstrings. Stale docs (e.g., old default values, removed IDs) are treated as bugs.
