# Modeling Notes & Pitfalls + FAQ restructure — design spec

**Date:** 2026-05-29
**Scope:** QSDsan docs PR #1 of a sequenced series. Subsequent specs will cover Process Specifications, Operational Flexibility (AgileSystem), Comparative TEA/LCA, and Convergence tutorials.

## Goal

Give QSDsan users a runnable reference for modeling surprises (the new notebook) and a cleaner non-modeling reference area (the restructured FAQ). Both delivered together in one PR because the FAQ split is required to make a clean home for the AI-Assisted Coding page, which currently lives stranded as tutorial 14.

## Deliverables

1. **`QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb`** — new notebook, 12 entries across 4 groups.
2. **`QSDsan/docs/source/tutorials/5_SanUnit_advanced.ipynb`** §1.1 — new subsection on `simulate()` / `_summary` mechanics; supports notebook entry 4.3.
3. **`QSDsan/docs/source/faq/`** — new subdirectory replacing the single `FAQ.rst`:
   - `index.rst` (umbrella with toctree)
   - `errors.rst`
   - `tips.rst`
   - `styling.rst`
   - `ai_assisted_coding.rst` (moved from `tutorials/14_AI_Assisted_Development.rst`)
4. **`QSDsan/docs/source/FAQ.rst`** — deleted (content moves into `faq/`).
5. **`QSDsan/docs/source/tutorials/14_AI_Assisted_Development.rst`** — deleted (content moves to `faq/ai_assisted_coding.rst`).
6. **`QSDsan/docs/source/tutorials/index.rst`** — drop `14_AI_Assisted_Development`, add `14_Modeling_Notes_and_Pitfalls`.
7. **`QSDsan/docs/source/index.rst`** — update top-level toctree references (single `FAQ` entry becomes `faq/index`).

Out of scope for this PR (tracked as follow-ups, not blocking): the other four tutorials (Process Specifications, AgileSystem, Comparative TEA/LCA, Convergence), the SanUnit-docstring-examples sweep, the deepening of tutorial 5 cost-decorator/auxiliary-unit content.

## Conventions

All notebook work follows the QSDsan tutorial skeleton from the `qsdsan-exposan-architecture` skill:

- Title cell with Binder badge, Colab note, Prepared by, Learning objectives (~3 bullets), Prerequisites, and Covered topics TOC. No Companion video blockquote for this notebook (per user direction).
- Setup cell with `<!-- tutorial-setup-section -->` marker + `qs.__version__` print line.
- Nav footer with `<!-- tutorial-nav-footer -->` and a single back-to-top link. No per-section back-to-top links.
- Sentence-case headings. Section anchors `sN`, subsection `sN.M`. Contiguous numbering.
- Address audience as "users".
- Minimize em-dashes.
- Alert-info divs for notes; never italic-wrapped inline code.
- No `[`code` text](url)` link form; use `the [code documentation](url)`.
- Canonical names in code: `qs.unit_operations`, `qs.process_models`, `qs.TEA`, `qs.CEPCI`, `qs.CEPCI_by_year`. Aliases (`sanunits`, `processes`) mentioned once if at all.
- Explicit stream/unit IDs matching variable names.
- `qsdsan` only in code cells (no `exposan` import unless an entry demands it).

Notebook entries use the agreed **Symptom → Why → Fix** template, with an optional fourth section `<div class="alert alert-info"><b>Coming from BioSTEAM?</b>...` only when the diff itself is the source of surprise.

## Notebook `14_Modeling_Notes_and_Pitfalls.ipynb` — structure

**Title:** "Modeling Notes & Pitfalls"
**Learning objectives:**
- Recognize common QSDsan modeling surprises by symptom.
- Diagnose whether a surprising result is a misconfiguration, a unit-design signal, or expected platform behavior.
- Apply the fix patterns for stream/component, unit-design, TEA/LCA, and BioSTEAM-inherited gotchas.

**Prerequisites:** notebooks 2 (Component), 3 (WasteStream), 4 (SanUnit basic), 6 (System). Group 2.3 cross-references 11 (Dynamic Simulation); Group 3 cross-references 7 (TEA) and 8 (LCA).

**Covered topics TOC:** four group links, `s1` through `s4`.

### Entry list — 4 groups × 3 entries

**§1. Streams and components** (anchor `s1`)

| # | Title | Symptom shown | Fix |
|---|---|---|---|
| 1.1 | `WasteStream` vs `SanStream` vs `Stream` | A unit silently drops composition info after passing a plain `Stream` into a WW path. | Rule of thumb: anything downstream of WW characterization needs `WasteStream`; pure mass/utility streams can stay `Stream`. |
| 1.2 | `Component` is not `Chemical` | A `Component` setter that "should" work (e.g., a Thermosteam property) raises. | Use `Component`-native attributes; `particulate=True` changes phase behavior; some Thermosteam tricks don't apply. |
| 1.3 | Stream IDs and the registry | Reusing an ID across cells produces `replaced in registry`; a later `System` resolves to the wrong stream. | Always pass an explicit ID matching the variable name. Auto-IDs (`ws1`) drift from the variable a reader sees. |

**§2. Unit design and simulation** (anchor `s2`)

| # | Title | Symptom shown | Fix |
|---|---|---|---|
| 2.1 | Absurd geometry from `_design` is a signal | A reactor sized as a "pancake" (very wide, very shallow) or a pipe with 1 cm diameter. | The design algorithm is flagging an unphysical operating point. Rebalance flow/loading, don't swap the reference flow. (METAB UASB+M is the canonical example but the entry doesn't require it.) |
| 2.2 | Recycle convergence with biokinetic models | Steady-state ASM/ADM loop stalls at tight tolerance. | Loosen tolerance before doubling iterations; check initial-guess plausibility. Links to `6_System.ipynb`. |
| 2.3 | Dynamic simulation timing pitfalls | Sim runs but every state stays at its initial value. | Trace to a missing `_init_state` override or an incorrect `t_span`/`atol`/`rtol`. Links to `11_Dynamic_Simulation.ipynb`. |

**§3. TEA and LCA** (anchor `s3`)

| # | Title | Symptom shown | Fix |
|---|---|---|---|
| 3.1 | CEPCI year drifts costs silently | A benchmark unit's installed cost is ~20% off from a paper. | Set `qs.CEPCI = qs.CEPCI_by_year[YEAR]` at the top of every analysis; assert it. Default index is fixed at import. |
| 3.2 | Purchase vs installed vs total capital cost | "My CAPEX doesn't add up." | Walk through the three fields, Lang factor effect, what `_cost` populates vs what TEA sums. Most issues are reading the wrong field. |
| 3.3 | LCA functional unit and dynamic vs static inventory | Same system, ~10× different impact under different FU. | Pick the FU deliberately (kg COD removed / m³ treated / person-year served). Note when dynamic-LCA diverges from static. Links to `8_LCA.ipynb`. |

**§4. Behavior inherited from BioSTEAM** (anchor `s4`)

| # | Title | Symptom shown | Fix |
|---|---|---|---|
| 4.1 | Pressure and energy balance idealizations | A pumping unit reports zero energy cost. | Pressure defaults to ambient, T changes ignored. Usually fine for WW; show one case where it isn't and how to override. Includes "Coming from BioSTEAM?" sidenote. |
| 4.2 | Degrees of freedom and over-constraint | A System with N+1 specs on N DOF: solver thrashes or silently picks one. | Brief preview only; full treatment in the future Process Specifications tutorial. |
| 4.3 | `simulate()` runs more than `_run` | A parameter sweep where capital cost never changes across iterations. | Calling `_run` alone leaves `_design` and `_cost` stale. Use `unit.simulate()` or `system.simulate()`. Links to the new tutorial 5 §1.1 subsection. |

## Tutorial 5 §1.1 — new `_summary` / `simulate()` subsection

Added inside the existing `### 1.1. Fundamental methods` section of `5_SanUnit_advanced.ipynb`, after the existing coverage of `_run`, `_design`, `_cost`.

Format: short markdown paragraph + a Python Aside `<details><summary>` block for the call-graph detail.

Content:
- `Unit.simulate(run=True)` calls `self.run()` (which dispatches to `_run`), then `self._summary(...)` which in turn calls `_design` and `_cost`.
- Calling `_run` directly leaves design and cost stale. Calling `simulate()` is the safe entrypoint.
- Tie to the existing comment `M2.simulate() # don't forget this!` already in tutorial 4 — that comment now points readers to this subsection (one-line update to tutorial 4 in the same PR).

Tutorial 5 's anchor IDs renumber if the section is inserted as a new sub-subsection; alternative is to keep it inline within 1.1 without new numbering. Decision: inline within 1.1 to avoid renumbering existing anchors.

## FAQ restructure

### Before

```
docs/source/
  FAQ.rst                  (single file, ~230 lines, 3 sections: errors / tips / styling)
  tutorials/
    14_AI_Assisted_Development.rst   (514 lines, stranded as "tutorial 14")
```

### After

```
docs/source/
  faq/
    index.rst              (TOC + one-paragraph intro pointing at the 4 sub-pages)
    errors.rst             (existing "Common Errors" + new: Version compatibility)
    tips.rst               (existing "Tips" + new: Dev environment setup, EXPOsan orientation)
    styling.rst            (existing "Styling" + new: Citing QSDsan)
    ai_assisted_coding.rst (moved from tutorials/14_AI_Assisted_Development.rst, unmodified)
  tutorials/
    14_Modeling_Notes_and_Pitfalls.ipynb  (new)
```

### Sub-page contents

| Page | Existing content | New entries |
|---|---|---|
| `errors.rst` | Graphviz install, ModuleNotFoundError, numba "underlying object vanished", UnicodeDecodeError | Version compatibility (Python 3.12+, current BioSTEAM/Thermosteam pins, pyproject.toml as source of truth) |
| `tips.rst` | Archive Branch, Pickle Protocol, Private Fork, Upgrade Python | Dev environment setup (`pip install -e .[dev]`, kernelspec for `.venv`, Windows path note); EXPOsan orientation (what it is, where systems live, how to run `exposan.bsm1`) |
| `styling.rst` | QSDsan vs qsdsan capitalization | Citing QSDsan (paragraph + BibTeX) |
| `ai_assisted_coding.rst` | (full move from current tutorials/14) | none — content moves verbatim |

### Cross-links

- **`faq/errors.rst`** top: one-line callout pointing to the new notebook for modeling-time surprises (number wrong / unit sized weirdly / sim won't converge).
- **`14_Modeling_Notes_and_Pitfalls.ipynb`** intro: one-line callout pointing to `faq/errors.rst` for install/env errors.
- **Tutorial 4** comment "don't forget this!": link to tutorial 5 §1.1 new subsection.

## Sphinx toctree updates

- `docs/source/index.rst`: replace the single `FAQ` entry with `faq/index`. Verify any other top-level toctree (e.g., in a `getting_started.rst` if present) still resolves.
- `docs/source/tutorials/index.rst`: remove `14_AI_Assisted_Development`, add `14_Modeling_Notes_and_Pitfalls`. Update the tutorial table (anchors, titles) per the existing layout.
- Verify nbsphinx anchor slugs for the new notebook follow the case-preserving / dot-preserving / hyphenated-space pattern (skill calls this out).

## Verification

Per `superpowers:verification-before-completion`, before claiming done:

1. **Notebook executes end-to-end.** Run `14_Modeling_Notes_and_Pitfalls.ipynb` top-to-bottom with cleared outputs in the project `.venv`. Only intentional teaching errors remain. (`nbsphinx_allow_errors = True` permits these; we still want a clean run for cells that aren't deliberate errors.)
2. **Tutorial 5 still executes.** Verify the new `_summary` subsection didn't introduce execution issues.
3. **Sphinx build is clean.** Run the docs build, confirm zero `toc.not_included` warnings and zero undefined-label warnings for the new and renamed pages.
4. **Cross-link resolution.** Open the built HTML for the notebook intro callout, FAQ errors callout, and tutorial 4 comment link — each must resolve to the right page.
5. **FAQ search affordance.** Confirm sphinx search still surfaces "FAQ" as a hit for queries like "graphviz" and "private fork" — the URL change from `FAQ.html#section` to `faq/errors.html#section` shouldn't break the search index, but verify.

## Risks and tradeoffs

- **FAQ link rot.** Anyone with a bookmark to `FAQ.html#graphviz-installation` lands on 404 after this PR. Default mitigation: ship a single-line `FAQ.rst` stub with `:orphan:` and a meta-refresh to `faq/index.html`, so external bookmarks land on the new index. Drop the stub later if it's never hit.
- **Notebook execution time.** Twelve entries with minimal toy snippets should run in under 30 s; the `_summary` parameter-sweep symptom for 4.3 may push longer. If it does, the symptom switches from "watch cost not change" to a static screenshot reference.
- **AI-assisted coding visibility.** Moving from `tutorials/14` to `faq/ai_assisted_coding` reduces tutorial-list prominence. Acceptable per user direction; the FAQ TOC and an explicit link from `tutorials/index.rst` ("Looking for AI-assisted coding guidance? See `faq/ai_assisted_coding`") keep it discoverable.

## Implementation order (for the plan)

1. Tutorial 5 §1.1 `_summary` subsection + tutorial 4 one-line link update.
2. FAQ directory scaffolding (`faq/index.rst`, empty sub-pages).
3. Split existing `FAQ.rst` into `errors.rst`, `tips.rst`, `styling.rst`.
4. Move `tutorials/14_AI_Assisted_Development.rst` to `faq/ai_assisted_coding.rst`.
5. Add the 4 new FAQ entries.
6. Update `docs/source/index.rst` toctree.
7. Build a skeleton of `14_Modeling_Notes_and_Pitfalls.ipynb` — title cell, setup cell, 4 group headers, nav footer.
8. Write entries §1.1 → §4.3 in order, executing each as it's added.
9. Add cross-links (FAQ ↔ notebook, tutorial 4 → tutorial 5).
10. Final docs build + search-index check + verification checklist.
