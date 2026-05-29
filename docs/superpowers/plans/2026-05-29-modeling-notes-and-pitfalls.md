# Modeling Notes & Pitfalls + FAQ restructure — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship one PR that (a) adds `tutorials/14_Modeling_Notes_and_Pitfalls.ipynb` with 12 entries across 4 topic groups, (b) extends tutorial 5 §1.1 with a `simulate()`/`_summary` mechanics subsection, and (c) splits `FAQ.rst` into a 4-page `faq/` directory that absorbs the stranded `tutorials/14_AI_Assisted_Development.rst`.

**Architecture:** Pure documentation change. New notebook follows the QSDsan tutorial skeleton (per `qsdsan-exposan-architecture`). FAQ splits into `docs/source/faq/{index,errors,tips,styling,ai_assisted_coding}.rst`. Existing `FAQ.rst` becomes an `:orphan:` redirect stub to preserve external bookmarks. All cross-links (top-level `index.rst`, `tutorials/index.rst`, `:ref:` anchors) updated atomically.

**Tech Stack:** Sphinx + nbsphinx + Furo (existing docs build). Notebooks executed by `nbsphinx` at build time when `outputs` are cleared (`nbsphinx_execute = 'auto'`). Build runs in the project `.venv` per [reference_dev_environment](../../../tmps/../memory/reference_dev_environment.md).

**Spec:** [2026-05-29-modeling-notes-and-pitfalls-design.md](../specs/2026-05-29-modeling-notes-and-pitfalls-design.md)

---

## Working environment

- All paths in this plan are relative to `c:/Users/Yalin/Documents/Coding/QSDsan-platform/`.
- Activate the project venv before any `make html` or `jupyter nbconvert` step:
  - PowerShell: `.venv\Scripts\Activate.ps1` (run from `QSDsan-platform/`).
  - The venv must have `qsdsan` and `exposan` editably installed (`pip install -e .[dev]` in each repo).
- Notebook editing: use `NotebookEdit` for cell-level edits when available; otherwise read the JSON, modify, write back with `json.dump(..., indent=1, ensure_ascii=False)` and a trailing newline. Clear `outputs` and `execution_count` on any code cell whose `source` changed.
- Build command (from `QSDsan/docs/`): `./make.bat html` on Windows, output under `QSDsan/docs/build/html/`.
- Each task ends with a `git add` + `git commit`. Run `git status` first, stage only the files the task touched, never `-A` / `.`.
- Never push. The user pulls and reviews diffs before any push (per [feedback_git_workflow](../../../memory/feedback_git_workflow.md)).

---

## File map

### New files

- `QSDsan/docs/source/faq/index.rst` — FAQ landing page with 4-entry toctree.
- `QSDsan/docs/source/faq/errors.rst` — Common Errors (split from FAQ.rst) + new Version Compatibility entry.
- `QSDsan/docs/source/faq/tips.rst` — Tips (split from FAQ.rst) + new Dev Environment Setup + new EXPOsan Orientation.
- `QSDsan/docs/source/faq/styling.rst` — Styling (split from FAQ.rst) + new Citing QSDsan.
- `QSDsan/docs/source/faq/ai_assisted_coding.rst` — moved verbatim from `tutorials/14_AI_Assisted_Development.rst`.
- `QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb` — new notebook, 12 entries.

### Modified files

- `QSDsan/docs/source/FAQ.rst` — replaced with `:orphan:` + meta-refresh redirect stub.
- `QSDsan/docs/source/index.rst` — toctree `FAQ` → `faq/index`; the `:link: faq :link-type: ref` grid card resolves via the `.. _faq:` anchor preserved in `faq/index.rst`.
- `QSDsan/docs/source/tutorials/index.rst` — remove `14_AI_Assisted_Development` entry; add `14_Modeling_Notes_and_Pitfalls` in a new "Reference" section.
- `QSDsan/docs/source/tutorials/4_SanUnit_basic.ipynb` — one comment update on `M2.simulate() # don't forget this!` line, pointing to tutorial 5 §1.1.
- `QSDsan/docs/source/tutorials/5_SanUnit_advanced.ipynb` — new `simulate()`/`_summary` mechanics paragraphs added inline within §1.1 "Fundamental methods" (no anchor renumbering).
- `QSDsan/CHANGELOG.rst` — single bullet under the next unreleased section per [feedback_changelog_convention](../../../memory/feedback_changelog_convention.md).

### Deleted files

- `QSDsan/docs/source/tutorials/14_AI_Assisted_Development.rst` (content lives at `faq/ai_assisted_coding.rst`).

---

## Phase A — Tutorial 5 + tutorial 4 link

### Task 1: Add `simulate()` / `_summary` subsection to tutorial 5 §1.1

**Files:**
- Modify: `QSDsan/docs/source/tutorials/5_SanUnit_advanced.ipynb` (markdown cell inside §1.1 "Fundamental methods", after the existing `_run` / `_design` / `_cost` coverage)

- [ ] **Step 1: Locate the §1.1 markdown cell.**

Open the notebook JSON and find the markdown cell whose source ends the coverage of `_run`, `_design`, `_cost`. It is the last markdown cell before the §1.2 "Useful attributes" heading. Confirm by reading `5_SanUnit_advanced.ipynb` and identifying the cell index.

- [ ] **Step 2: Append a "Putting them together" paragraph to that markdown cell.**

Append this text to the end of the existing §1.1 markdown cell's source (preserve existing trailing newlines, then add):

```markdown

#### Putting them together: `simulate()`

When you call `unit.simulate()`, the unit runs in three stages: first `_run` (mass/energy balance and outlet streams), then `_summary` which invokes `_design` (size/geometry) and then `_cost` (purchase and installed cost). The same applies to `system.simulate()`, which dispatches `simulate()` to each unit in convergence order.

A common pitfall: calling `_run` directly (for debugging, or by reaching into a unit programmatically) leaves `_design` and `_cost` untouched. The unit's outlet streams update, but its design table and cost results stay stale from the previous full simulation, which then propagates into TEA results. Always go through `simulate()` (on the unit or on the enclosing system) when you want results that downstream analysis will read.

<div class="alert alert-info">
<b>Heads up.</b> Some other tutorial cells call <code>simulate()</code> without explaining what it glues together; this is what they mean. See <a href="14_Modeling_Notes_and_Pitfalls.html#s4.3.-simulate-runs-more-than-_run">Modeling Notes &amp; Pitfalls §4.3</a> for a worked symptom (a parameter sweep where cost never changes).
</div>
```

- [ ] **Step 3: Clear outputs on any modified code cells.**

Only the markdown cell changed, so no code-cell outputs to clear. Verify the notebook JSON still parses by loading it in Python:

```bash
python -c "import json; json.load(open('QSDsan/docs/source/tutorials/5_SanUnit_advanced.ipynb', encoding='utf-8'))"
```

Expected: no exception.

- [ ] **Step 4: Re-execute the notebook to confirm nothing broke.**

```bash
jupyter nbconvert --to notebook --execute --inplace QSDsan/docs/source/tutorials/5_SanUnit_advanced.ipynb
```

Expected: all cells run; existing intentional teaching errors (if any) remain; no new errors. After execution, clear outputs again so nbsphinx will re-execute at docs build:

```bash
jupyter nbconvert --clear-output --inplace QSDsan/docs/source/tutorials/5_SanUnit_advanced.ipynb
```

- [ ] **Step 5: Commit.**

```bash
git -C QSDsan add docs/source/tutorials/5_SanUnit_advanced.ipynb
git -C QSDsan commit -m "docs: explain simulate()/_summary mechanics in tutorial 5 §1.1"
```

### Task 2: Update tutorial 4's "don't forget this!" comment with a link

**Files:**
- Modify: `QSDsan/docs/source/tutorials/4_SanUnit_basic.ipynb` (the code cell containing `M2.simulate() # don't forget this!`)

- [ ] **Step 1: Find the cell.**

Search the notebook JSON for the source line `M2.simulate() # don't forget this!`. Note the cell index.

- [ ] **Step 2: Update the comment to point to tutorial 5.**

Replace the line `M2.simulate() # don't forget this!\n` in that cell's source with:

```python
M2.simulate()  # runs _run then _summary (_design + _cost); see tutorial 5 §1.1
```

- [ ] **Step 3: Clear that cell's output and execution_count.**

The source changed, so per `qsdsan-exposan-architecture` clear `outputs: []` and `execution_count: null` on that cell.

- [ ] **Step 4: Re-execute the notebook and re-clear.**

```bash
jupyter nbconvert --to notebook --execute --inplace QSDsan/docs/source/tutorials/4_SanUnit_basic.ipynb
jupyter nbconvert --clear-output --inplace QSDsan/docs/source/tutorials/4_SanUnit_basic.ipynb
```

Expected: notebook executes cleanly.

- [ ] **Step 5: Commit.**

```bash
git -C QSDsan add docs/source/tutorials/4_SanUnit_basic.ipynb
git -C QSDsan commit -m "docs: link tutorial 4 simulate() comment to tutorial 5 §1.1 explanation"
```

---

## Phase B — FAQ restructure

### Task 3: Create the `faq/` directory and `faq/index.rst`

**Files:**
- Create: `QSDsan/docs/source/faq/index.rst`

- [ ] **Step 1: Create the directory.**

```bash
mkdir -p QSDsan/docs/source/faq
```

- [ ] **Step 2: Write `faq/index.rst`.**

```rst
.. _faq:

FAQ
===

Reference for everything outside the modeling tutorials: install errors, dev workflow, styling and citation, and AI-assisted coding guidance.

.. note::
   Looking for modeling-time surprises (a number that looks wrong, a unit that sized weirdly, a sim that won't converge)? See :doc:`Modeling Notes & Pitfalls <../tutorials/14_Modeling_Notes_and_Pitfalls>`.

.. toctree::
   :maxdepth: 2

   errors
   tips
   styling
   ai_assisted_coding
```

- [ ] **Step 3: Verify.**

```bash
python -c "from pathlib import Path; t = Path('QSDsan/docs/source/faq/index.rst').read_text(); assert '.. _faq:' in t and '14_Modeling_Notes_and_Pitfalls' in t"
```

Expected: no AssertionError.

- [ ] **Step 4: Commit.**

```bash
git -C QSDsan add docs/source/faq/index.rst
git -C QSDsan commit -m "docs: scaffold faq/ subdirectory with index"
```

### Task 4: Split out `faq/errors.rst` and add Version Compatibility

**Files:**
- Create: `QSDsan/docs/source/faq/errors.rst`

Pull the **Common Errors** content (everything from `Common Errors` heading through `UnicodeDecodeError` section) out of the existing `FAQ.rst`. Preserve the `.. _graphviz-installation:` anchor — it is referenced from `docs/source/index.rst:85`.

- [ ] **Step 1: Write `faq/errors.rst`.**

```rst
Common Errors
=============

.. note::
   For modeling-time surprises (a number that looks wrong, a unit that sized weirdly, a sim that won't converge), see :doc:`Modeling Notes & Pitfalls <../tutorials/14_Modeling_Notes_and_Pitfalls>`.

.. _graphviz-installation:

``Graphviz`` Installation
-------------------------

[paste the existing graphviz section here verbatim from FAQ.rst lines 11-59 — text only, no leading metadata]

``ModuleNotFoundError``
-----------------------

[paste lines 62-90 verbatim]

``underlying object has vanished``
----------------------------------

[paste lines 93-110 verbatim]

``UnicodeDecodeError``
----------------------

[paste lines 113-122 verbatim]

Version compatibility
---------------------

``QSDsan`` requires Python 3.12 or newer and pins specific minimum versions of ``biosteam``, ``thermosteam``, and other dependencies. The authoritative version table lives in `pyproject.toml <https://github.com/QSD-Group/QSDsan/blob/main/pyproject.toml>`_; the table here is a quick reference and may lag.

If you see import errors after upgrading one of those dependencies independently of ``QSDsan``, the most reliable fix is to reinstall ``QSDsan`` in a fresh environment so the resolver picks compatible versions of everything at once.
```

The `[paste ...]` lines are placeholders for content already in the repo — read the existing `FAQ.rst` and copy each section's full text including code blocks and notes, preserving the underline-style headers (convert from `***` to `---` to be section-level under the new H1).

- [ ] **Step 2: Convert subsection underlines.**

The original `FAQ.rst` uses `***` for `***` (sub-subsection) under `---` sections. In the new file, each gotcha is a top-level section under the H1 "Common Errors," so use `---` for each gotcha heading.

- [ ] **Step 3: Verify the file is syntactically valid RST.**

```bash
python -c "import docutils.core; docutils.core.publish_doctree(open('QSDsan/docs/source/faq/errors.rst', encoding='utf-8').read())" 2>&1 | head -20
```

Expected: no SEVERE/ERROR lines, only WARNING about unknown reference (`:doc:` resolves only inside Sphinx, not bare docutils).

- [ ] **Step 4: Commit.**

```bash
git -C QSDsan add docs/source/faq/errors.rst
git -C QSDsan commit -m "docs(faq): split errors page from FAQ.rst, add version compatibility entry"
```

### Task 5: Split out `faq/tips.rst` and add Dev Environment + EXPOsan Orientation

**Files:**
- Create: `QSDsan/docs/source/faq/tips.rst`

- [ ] **Step 1: Write `faq/tips.rst`.**

Use the same paste-existing-and-add pattern:

```rst
Tips
====

Setting up a development environment
------------------------------------

If you cloned ``QSDsan`` (and/or ``EXPOsan``) and want to work on the code, install in editable mode so Python imports your local copy:

.. code:: bash

    pip install -e ".[dev]"

Use a dedicated virtual environment (``.venv`` or a ``conda`` env) to keep the editable install isolated. To use that environment in Jupyter:

.. code:: bash

    python -m ipykernel install --user --name qsdsan-dev

Then pick ``qsdsan-dev`` as the kernel when you open a tutorial notebook. If you also work on ``EXPOsan``, repeat the editable install in its clone too — sibling editable installs let cross-repo changes show up immediately.

On Windows, paths with spaces (e.g., ``C:\Users\Your Name\Documents``) occasionally break tools that don't quote them. If a build or test fails with a "file not found" error involving such a path, try a path without spaces.

Archive Branch
--------------

[paste existing Archive Branch content from FAQ.rst]

Pickle Protocol
---------------

[paste existing]

Private Fork
------------

[paste existing]

Upgrade Python
--------------

[paste existing]

EXPOsan orientation
-------------------

``EXPOsan`` (the EXPOsition of sanitation and resource recovery systems) is a companion package that catalogs published sanitation and resource-recovery systems built on ``QSDsan``. The packages are released together but live in `separate <https://github.com/QSD-Group/QSDsan>`_ `repos <https://github.com/QSD-Group/EXPOsan>`_.

To run an EXPOsan system:

.. code:: python

    from exposan import bsm1
    bsm1.load()
    bsm1.sys.simulate()
    bsm1.sys.diagram()

See the `EXPOsan documentation <https://exposan.readthedocs.io/>`_ for the full system catalog.
```

- [ ] **Step 2: Verify.**

```bash
python -c "import docutils.core; docutils.core.publish_doctree(open('QSDsan/docs/source/faq/tips.rst', encoding='utf-8').read())" 2>&1 | head
```

Expected: no SEVERE/ERROR lines.

- [ ] **Step 3: Commit.**

```bash
git -C QSDsan add docs/source/faq/tips.rst
git -C QSDsan commit -m "docs(faq): split tips page from FAQ.rst, add dev environment setup and EXPOsan orientation"
```

### Task 6: Split out `faq/styling.rst` and add Citing QSDsan

**Files:**
- Create: `QSDsan/docs/source/faq/styling.rst`

- [ ] **Step 1: Write `faq/styling.rst`.**

```rst
Styling and Citation
====================

``QSDsan`` vs. ``qsdsan``
-------------------------

[paste existing Styling section content from FAQ.rst lines 226-234]

Citing QSDsan
-------------

If you use ``QSDsan`` in published work, please cite:

  Li, Y.; Trimmer, J. T.; Hand, S.; Zhang, X.; Chambers, K. G.; Lohman, H. A. C.; Shi, R.; Byrne, D. M.; Cook, S. M.; Guest, J. S. *Quantitative Sustainable Design (QSD): A Methodology for the Prioritization of Research, Development, and Deployment of Technologies.* Environmental Science: Water Research & Technology, 2022. https://doi.org/10.1039/D2EW00431C

BibTeX:

.. code:: bibtex

    @article{li2022qsd,
      author  = {Li, Yalin and Trimmer, John T. and Hand, Stetson and Zhang, Xinyi and Chambers, Kathryn G. and Lohman, Hannah A. C. and Shi, Ruixiao and Byrne, Diana M. and Cook, Sherri M. and Guest, Jeremy S.},
      title   = {Quantitative Sustainable Design ({QSD}): A Methodology for the Prioritization of Research, Development, and Deployment of Technologies},
      journal = {Environmental Science: Water Research \& Technology},
      year    = {2022},
      doi     = {10.1039/D2EW00431C}
    }

For specific subpackages or unit models, also cite the relevant primary references listed in their docstrings.
```

- [ ] **Step 2: Verify citation BibTeX renders.**

```bash
python -c "import docutils.core; docutils.core.publish_doctree(open('QSDsan/docs/source/faq/styling.rst', encoding='utf-8').read())" 2>&1 | head
```

Expected: no SEVERE/ERROR lines.

- [ ] **Step 3: Commit.**

```bash
git -C QSDsan add docs/source/faq/styling.rst
git -C QSDsan commit -m "docs(faq): split styling page from FAQ.rst, add Citing QSDsan section"
```

### Task 7: Move AI-Assisted Coding into `faq/ai_assisted_coding.rst`

**Files:**
- Create: `QSDsan/docs/source/faq/ai_assisted_coding.rst` (via git mv)
- Delete: `QSDsan/docs/source/tutorials/14_AI_Assisted_Development.rst`

- [ ] **Step 1: Move the file with git mv to preserve history.**

```bash
git -C QSDsan mv docs/source/tutorials/14_AI_Assisted_Development.rst docs/source/faq/ai_assisted_coding.rst
```

- [ ] **Step 2: Update the H1 title if needed.**

Read the first 5 lines of the new file. If the H1 title is "AI-Assisted Development" (or similar), leave it. The content moves verbatim — no edits to body.

```bash
head -5 QSDsan/docs/source/faq/ai_assisted_coding.rst
```

- [ ] **Step 3: Verify file moved and contents intact.**

```bash
test -f QSDsan/docs/source/faq/ai_assisted_coding.rst && ! test -f QSDsan/docs/source/tutorials/14_AI_Assisted_Development.rst && echo OK
```

Expected: prints `OK`.

- [ ] **Step 4: Commit.**

```bash
git -C QSDsan commit -m "docs(faq): move AI-assisted coding from tutorials/14 to faq/ai_assisted_coding"
```

### Task 8: Replace `FAQ.rst` with redirect stub

**Files:**
- Modify: `QSDsan/docs/source/FAQ.rst` (complete rewrite to a short stub)

- [ ] **Step 1: Overwrite `FAQ.rst` with the stub.**

```rst
:orphan:

FAQ
===

The FAQ has moved. See :doc:`faq/index`.

.. raw:: html

   <meta http-equiv="refresh" content="0; url=faq/index.html">
```

The `:orphan:` directive prevents Sphinx from warning about FAQ.rst not being in any toctree. The `<meta http-equiv="refresh">` gives external bookmarks a working landing experience.

- [ ] **Step 2: Verify the file is the new stub (not the old 234-line file).**

```bash
wc -l QSDsan/docs/source/FAQ.rst
```

Expected: under 20 lines.

- [ ] **Step 3: Commit.**

```bash
git -C QSDsan add docs/source/FAQ.rst
git -C QSDsan commit -m "docs: replace FAQ.rst with redirect stub to faq/index"
```

### Task 9: Update top-level `index.rst` toctree

**Files:**
- Modify: `QSDsan/docs/source/index.rst` (line 124 area)

- [ ] **Step 1: Replace the `FAQ` toctree entry with `faq/index`.**

Find the toctree block:

```rst
.. toctree::
   :maxdepth: 2
   :hidden:

   FAQ
```

Replace `FAQ` with `faq/index`. The grid-item-card at line 41-44 uses `:link: faq :link-type: ref`, which resolves via the `.. _faq:` anchor preserved in `faq/index.rst` — no change needed to the grid card.

- [ ] **Step 2: Verify.**

```bash
grep -n "faq/index\|^   FAQ$" QSDsan/docs/source/index.rst
```

Expected: shows `faq/index` line and no bare `FAQ` line in the toctree.

- [ ] **Step 3: Commit.**

```bash
git -C QSDsan add docs/source/index.rst
git -C QSDsan commit -m "docs: point top-level toctree at faq/index"
```

### Task 10: Update `tutorials/index.rst` toctree

**Files:**
- Modify: `QSDsan/docs/source/tutorials/index.rst` (around line 147-153)

- [ ] **Step 1: Replace the "For Contributors" section.**

Find:

```rst
For Contributors
----------------

.. toctree::
   :maxdepth: 1

   14_AI_Assisted_Development
```

Replace with:

```rst
Reference
---------

.. toctree::
   :maxdepth: 1

   14_Modeling_Notes_and_Pitfalls
```

The AI-Assisted Coding page now lives under the FAQ; readers reach it via the FAQ TOC.

- [ ] **Step 2: Verify.**

```bash
grep -n "14_Modeling_Notes_and_Pitfalls\|14_AI_Assisted_Development" QSDsan/docs/source/tutorials/index.rst
```

Expected: one match for `14_Modeling_Notes_and_Pitfalls`, zero for `14_AI_Assisted_Development`.

- [ ] **Step 3: Commit.**

```bash
git -C QSDsan add docs/source/tutorials/index.rst
git -C QSDsan commit -m "docs: replace tutorials/14 AI-assisted slot with Modeling Notes & Pitfalls"
```

### Task 11: First build checkpoint

- [ ] **Step 1: Build the docs.**

From `QSDsan/docs/`:

```bash
./make.bat html 2>&1 | tee /tmp/build.log
```

Or on PowerShell:

```powershell
cd QSDsan\docs ; .\make.bat html *>&1 | Tee-Object -FilePath $env:TEMP\build.log
```

- [ ] **Step 2: Scan for warnings.**

```bash
grep -E "WARNING|ERROR|undefined label|not_included" /tmp/build.log | head -30
```

Expected: zero "undefined label" / "not_included" / "ERROR" lines. Two expected warnings: (a) one about `14_Modeling_Notes_and_Pitfalls` not yet existing — to be resolved in Task 12, and (b) the orphan stub is fine.

If anything else surfaces, fix it before proceeding. Do not commit warnings forward.

- [ ] **Step 3: Open the rebuilt FAQ index in a browser.**

```bash
start QSDsan/docs/build/html/faq/index.html
```

Confirm: the 4 sub-pages link correctly, anchor `#faq` (used by the home page grid card) lands on this page.

- [ ] **Step 4: Confirm redirect stub.**

```bash
start QSDsan/docs/build/html/FAQ.html
```

Confirm: page redirects to `faq/index.html`.

---

## Phase C — Notebook skeleton

### Task 12: Create the `14_Modeling_Notes_and_Pitfalls.ipynb` skeleton

**Files:**
- Create: `QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb`

- [ ] **Step 1: Write the notebook JSON.**

The skeleton must include: title cell, setup cell, 4 group section headers with anchors, and the nav footer. No entries yet — those land in Phase D.

Use the existing tutorial 7 (`7_TEA.ipynb`) as the structural reference for the title cell (Binder badge, Colab note, Prepared by, Learning objectives, Prerequisites, Covered topics TOC).

Title cell markdown content:

```markdown
# Modeling Notes & Pitfalls <a class="anchor" id="top"></a>

[![Launch Binder](../images/custom_binder_logo.svg)](https://mybinder.org/v2/gh/QSD-Group/QSDsan-env/main?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252FQSD-group%252FQSDsan%26urlpath%3Dlab%252Ftree%252FQSDsan%252Fdocs%252Fsource%252Ftutorials%252F14_Modeling_Notes_and_Pitfalls.ipynb%26branch%3Dmain)

To run this notebook in Google Colab, follow the :ref:`run-in-colab` instructions.

**Prepared by:** [Yalin Li](https://qsdsan.readthedocs.io/en/latest/AUTHORS.html)

**Learning objectives:**
- Recognize common QSDsan modeling surprises by their symptoms.
- Diagnose whether a surprising result is a misconfiguration, a unit-design signal, or expected platform behavior.
- Apply fix patterns across streams/components, unit design, TEA/LCA, and BioSTEAM-inherited behavior.

**Prerequisites:** [2. Component](2_Component.ipynb), [3. WasteStream](3_WasteStream.ipynb), [4. SanUnit basic](4_SanUnit_basic.ipynb), [6. System](6_System.ipynb). §2.3 references [11. Dynamic Simulation](11_Dynamic_Simulation.ipynb); §3 references [7. TEA](7_TEA.ipynb) and [8. LCA](8_LCA.ipynb).

**Covered topics:**

- <a href="#s1">1. Streams and components</a>
- <a href="#s2">2. Unit design and simulation</a>
- <a href="#s3">3. TEA and LCA</a>
- <a href="#s4">4. Behavior inherited from BioSTEAM</a>

<div class="alert alert-info">
<b>Looking for install or environment errors?</b> See the <a href="../faq/errors.html">FAQ Common Errors page</a>.
</div>
```

Setup cell markdown:

```markdown
<!-- tutorial-setup-section -->

## Setup
```

Setup code cell:

```python
import qsdsan as qs
print(f'This tutorial was made with qsdsan v{qs.__version__}.')
```

Group section header cells (one per group, each a markdown cell):

```markdown
## 1. Streams and components <a class="anchor" id="s1"></a>
```

```markdown
## 2. Unit design and simulation <a class="anchor" id="s2"></a>
```

```markdown
## 3. TEA and LCA <a class="anchor" id="s3"></a>
```

```markdown
## 4. Behavior inherited from BioSTEAM <a class="anchor" id="s4"></a>
```

Nav footer cell:

```markdown
<!-- tutorial-nav-footer -->

---

<a href="#top">↑ Back to top</a>
```

Write the file as JSON with `indent=1, ensure_ascii=False`, trailing newline. Use `nbformat` 4, `nbformat_minor` 5, `kernelspec` matching the other tutorials' kernelspec (read tutorial 7's kernelspec block and copy verbatim).

- [ ] **Step 2: Validate JSON.**

```bash
python -c "import json, nbformat; nb = nbformat.read('QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb', as_version=4); nbformat.validate(nb); print(f'{len(nb.cells)} cells')"
```

Expected: prints `8 cells` (1 title + 1 setup md + 1 setup code + 4 group headers + 1 nav footer) with no exception.

- [ ] **Step 3: Execute the skeleton.**

```bash
jupyter nbconvert --to notebook --execute --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
jupyter nbconvert --clear-output --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
```

Expected: setup cell prints the qsdsan version; no errors.

- [ ] **Step 4: Commit.**

```bash
git -C QSDsan add docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
git -C QSDsan commit -m "docs: scaffold Modeling Notes & Pitfalls tutorial skeleton"
```

### Task 13: Build checkpoint after skeleton

- [ ] **Step 1: Rebuild docs.**

```powershell
cd QSDsan\docs ; .\make.bat html *>&1 | Tee-Object -FilePath $env:TEMP\build.log
```

- [ ] **Step 2: Scan for warnings.**

```bash
grep -E "WARNING|ERROR|undefined label|not_included" "$TEMP/build.log" | head -30
```

Expected: zero warnings about `14_Modeling_Notes_and_Pitfalls`. The skeleton is now a valid (empty-bodied) page.

- [ ] **Step 3: Open the rendered notebook.**

```bash
start QSDsan/docs/build/html/tutorials/14_Modeling_Notes_and_Pitfalls.html
```

Confirm: title cell renders with all sections, 4 empty group headers visible, nav footer present.

---

## Phase D — Notebook entries

Each task in this phase adds one topic group's three entries to `14_Modeling_Notes_and_Pitfalls.ipynb`, executes the notebook, clears outputs, and commits. Follow the **Symptom → Why → Fix** template strictly; add a `<div class="alert alert-info"><b>Coming from BioSTEAM?</b>...` block only where the diff is the source of surprise.

For each entry, insert cells after the corresponding group header cell, in this order:
1. Subsection markdown header: `### N.M. <Title> <a class="anchor" id="sN.M"></a>`
2. **Symptom** markdown explanation.
3. Symptom code cell (the reproducer).
4. **Why** markdown explanation.
5. **Fix** markdown explanation.
6. Fix code cell.
7. (Optional) "Coming from BioSTEAM?" alert-info div.

Use minimal toy components/streams; do not import EXPOsan unless explicitly noted. Always pass explicit IDs to streams and units. Sentence-case headings.

### Task 14: Add §1 — Streams and components

**Files:**
- Modify: `QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb`

For each subsection below, add the cells described in the order listed at the top of Phase D.

- [ ] **Step 1: §1.1 — `WasteStream` vs `SanStream` vs `Stream`**

Subsection title: `### 1.1. WasteStream vs SanStream vs Stream <a class="anchor" id="s1.1"></a>`

Symptom markdown:

> A unit configured for wastewater (e.g., a clarifier) is fed a plain `qs.Stream`, and the unit's outlet drops solids and concentration fields that downstream WW-aware code depends on.

Symptom code:

```python
import qsdsan as qs

cmps = qs.Components.load_default()
qs.set_thermo(cmps)

plain = qs.Stream('plain', H2O=1000, units='kg/hr')
print(type(plain).__name__)
print('has TSS attribute:', hasattr(plain, 'get_TSS'))
```

Why markdown:

> `Stream` is the BioSTEAM base type for mass/energy flows. `SanStream` adds construction/transportation/impact bookkeeping for LCA. `WasteStream` adds component-aggregation properties (`TSS`, `COD`, `BOD`, etc.) used by every WW characterization workflow. Mixing types isn't a hard error — it silently drops the richer attributes when a stream falls back to the base class.

Fix markdown:

> Use `WasteStream` for anything that participates in WW characterization; `SanStream` if you only need LCA bookkeeping without WW properties; `Stream` only for utility/heat-medium flows where neither matters.

Fix code:

```python
ww = qs.WasteStream('ww', H2O=1000, X_OHO=0.5, units='kg/hr')
print(type(ww).__name__)
print(f'TSS: {ww.get_TSS():.1f} mg/L')
```

No "Coming from BioSTEAM?" block needed — this is QSDsan-only context.

- [ ] **Step 2: §1.2 — `Component` is not `Chemical`**

Subsection title: `### 1.2. Component is not Chemical <a class="anchor" id="s1.2"></a>`

Symptom markdown:

> A user creates a `Component` and tries to use a Thermosteam `Chemical`-style property setter; the setter raises or silently no-ops.

Symptom code:

```python
import qsdsan as qs
from qsdsan import Component

cellulose = Component('Cellulose', search_ID='Cellulose',
                      particle_size='Particulate',
                      degradability='Slowly',
                      organic=True)
print('particle_size:', cellulose.particle_size)
try:
    cellulose.Tb = 500
except Exception as exc:
    print(f'{type(exc).__name__}: {exc}')
```

Why markdown:

> `Component` extends `Chemical` with WW-relevant attributes (`particle_size`, `degradability`, `organic`, etc.) but constrains certain phase-related setters because particulate components do not have well-defined boiling points or saturation pressures.

Fix markdown:

> Set the WW attributes via `Component` constructor or attribute assignment. For thermo properties that genuinely apply, set them through the same pattern Thermosteam uses; for particulates, leave phase-related properties alone.

Fix code:

```python
glucose = Component('Glucose', search_ID='Glucose',
                    particle_size='Soluble',
                    degradability='Readily',
                    organic=True)
glucose.Tb = 423.15
print(f'glucose Tb = {glucose.Tb} K')
```

- [ ] **Step 3: §1.3 — Stream IDs and the registry**

Subsection title: `### 1.3. Stream IDs and the registry <a class="anchor" id="s1.3"></a>`

Symptom markdown:

> Re-running a cell that creates a stream prints `replaced in registry`, and a later `System` resolves a stream by ID and pulls the wrong one.

Symptom code:

```python
import qsdsan as qs

cmps = qs.Components.load_default()
qs.set_thermo(cmps)

feed = qs.WasteStream('feed', H2O=1000)
feed = qs.WasteStream('feed', H2O=2000)  # warns: replaced in registry
print(qs.main_flowsheet.stream.feed.F_mass)
```

Why markdown:

> Streams (and units) register themselves by ID in the flowsheet registry. Reusing an ID overwrites the prior entry, then anything that resolved the old object by ID gets stale. Auto-IDs (`ws1`, `ws2`, ...) are even worse: they drift from the variable name a reader sees and become hard to grep.

Fix markdown:

> Always pass an explicit ID matching the variable name. To overwrite intentionally, `del` the old one first or call `qs.main_flowsheet.stream.clear()`.

Fix code:

```python
qs.main_flowsheet.stream.clear()
feed = qs.WasteStream('feed', H2O=2000)
print(f"feed mass: {qs.main_flowsheet.stream.feed.F_mass}")
```

- [ ] **Step 4: Execute and clear outputs.**

```bash
jupyter nbconvert --to notebook --execute --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
jupyter nbconvert --clear-output --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
```

Expected: all cells run; symptom cell in §1.2 prints an expected `AttributeError`/`ValueError` line (caught by try/except); §1.3 prints the "replaced in registry" warning to stderr.

- [ ] **Step 5: Commit.**

```bash
git -C QSDsan add docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
git -C QSDsan commit -m "docs(modeling-notes): add §1 streams and components"
```

### Task 15: Add §2 — Unit design and simulation

- [ ] **Step 1: §2.1 — Absurd geometry from `_design` is a signal**

Subsection title: `### 2.1. Absurd geometry from _design is a signal <a class="anchor" id="s2.1"></a>`

Symptom markdown:

> A reactor sized for a normal volumetric loading returns a "pancake": very wide diameter, very shallow height (or the reverse). The user assumes the design code is buggy.

Symptom code:

```python
import qsdsan as qs

# Toy CSTR sized far outside its intended HRT range.
# (Substitute a real qs unit class when this entry is fleshed out.)
class ToyCSTR:
    _design_HRT_range = (0.5, 24)  # hours
    def __init__(self, HRT):
        self.HRT = HRT
    def _design(self):
        if not (self._design_HRT_range[0] <= self.HRT <= self._design_HRT_range[1]):
            return {'note': f'HRT={self.HRT} h is outside [{self._design_HRT_range[0]}, {self._design_HRT_range[1]}]'}
        return {'D_m': 4.0, 'H_m': 8.0}

print(ToyCSTR(HRT=200)._design())
```

Why markdown:

> Design algorithms extrapolate outside their validated range; the absurd geometry is the model's way of flagging an unphysical operating point. The fix is rarely to swap the reference flow or override the dimension; it is to revisit HRT, SRT, or loading rate so the unit operates inside its design window.

Fix markdown:

> When a unit returns weird dimensions, check the unit's docstring or the design references it cites for the valid operating range. If you have a legitimate reason to operate outside that range, either subclass the unit and extend `_design` with a custom correlation, or insert a second stage so each unit stays inside its window.

Fix code:

```python
print(ToyCSTR(HRT=12)._design())
```

- [ ] **Step 2: §2.2 — Recycle convergence with biokinetic models**

Subsection title: `### 2.2. Recycle convergence with biokinetic models <a class="anchor" id="s2.2"></a>`

Symptom markdown:

> A System with an ASM-based bioreactor and a sludge recycle loop times out or returns `RecycleConvergenceError` at tight tolerance.

Symptom code:

```python
# Pseudocode-style demonstration; replace with a minimal real two-unit recycle
# system when fleshed out. The point: solver does not converge at strict tol.
import qsdsan as qs

print("(symptom: tight tolerance + far-from-solution initial guess => stall)")
print("Replace with a 2-unit recycle reproducer at implementation time.")
```

Why markdown:

> Biokinetic stoichiometry is stiff: small composition changes in the recycle loop produce large kinetic responses, which produce large composition swings. Aitken/Wegstein acceleration helps near the solution but can amplify oscillations when the initial guess is far off. The default tolerance is tighter than necessary for many WW analyses.

Fix markdown:

> 1) Loosen the System's `molar_tolerance` and `temperature_tolerance` first (a one-line change). 2) Improve the initial guess for the recycle stream (set realistic concentrations explicitly). 3) Only if 1 + 2 fail, increase max iterations. See [6. System](6_System.ipynb) for the full convergence-control surface.

Fix code:

```python
print("sys.molar_tolerance = 1e-3")
print("sys.temperature_tolerance = 1e-2")
print("# then sys.simulate()")
```

- [ ] **Step 3: §2.3 — Dynamic simulation timing pitfalls**

Subsection title: `### 2.3. Dynamic simulation timing pitfalls <a class="anchor" id="s2.3"></a>`

Symptom markdown:

> A dynamic simulation runs to completion, but every state stays at its initial value. The user assumes the ODE solver is broken.

Symptom code:

```python
# Pseudocode-style demonstration of the symptom; flesh out with a
# minimal dynamic unit when implemented.
print("(symptom: simulate() runs, state stays constant)")
print("Inspect: did the subclass override _init_state? did _compile_ODE wire dstate?")
```

Why markdown:

> Dynamic `SanUnit` subclasses must override `_init_state` to populate `self._state` from `ins`, and `_compile_ODE` to populate `self._dstate` from `self._state` and `ins`. If `_init_state` is missing, the state starts as zeros (or whatever the base class set), and if `_compile_ODE` is missing, `dstate` is zero so state never moves. `nbsphinx_allow_errors = True` doesn't help here — the sim doesn't error, it just doesn't move.

Fix markdown:

> Implement both methods in any dynamic subclass; cross-check `state` and `dstate` at `t=0` before running the full simulation. See [11. Dynamic Simulation](11_Dynamic_Simulation.ipynb) for the full state/dstate contract.

Fix code:

```python
print("Verify before t-stepping:")
print("u._init_state(); print(u._state)")
print("u._compile_ODE(); u._ode(t=0, y=u._state); print(u._dstate)")
```

- [ ] **Step 4: Execute and clear outputs.**

```bash
jupyter nbconvert --to notebook --execute --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
jupyter nbconvert --clear-output --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
```

Expected: no errors.

- [ ] **Step 5: Commit.**

```bash
git -C QSDsan add docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
git -C QSDsan commit -m "docs(modeling-notes): add §2 unit design and simulation"
```

### Task 16: Add §3 — TEA and LCA

- [ ] **Step 1: §3.1 — CEPCI year drifts costs silently**

Subsection title: `### 3.1. CEPCI year drifts costs silently <a class="anchor" id="s3.1"></a>`

Symptom markdown:

> A benchmark unit's installed cost prints ~20% off from a published reference. The user blames the unit's correlation.

Symptom code:

```python
import qsdsan as qs

print(f'qs.CEPCI = {qs.CEPCI}')
print(f'CEPCI[2023] = {qs.CEPCI_by_year[2023]}')
print(f'CEPCI[2017] = {qs.CEPCI_by_year[2017]}')
print(f'ratio = {qs.CEPCI_by_year[2023] / qs.CEPCI_by_year[2017]:.3f}')
```

Why markdown:

> `qs.CEPCI` is module-level state set once at import. If a unit's purchase-cost correlation was developed for a different year than the current `qs.CEPCI`, costs are scaled by the ratio of those CEPCIs. When a user runs multiple analyses targeting different reference years without setting `qs.CEPCI` explicitly, costs silently drift.

Fix markdown:

> Set `qs.CEPCI` deliberately at the top of every analysis, and assert it immediately afterward so misconfiguration fails loudly:

Fix code:

```python
qs.CEPCI = qs.CEPCI_by_year[2023]
assert qs.CEPCI == qs.CEPCI_by_year[2023]
print(f'CEPCI locked to {qs.CEPCI} (year 2023)')
```

- [ ] **Step 2: §3.2 — Purchase vs installed vs total capital cost**

Subsection title: `### 3.2. Purchase vs installed vs total capital cost <a class="anchor" id="s3.2"></a>`

Symptom markdown:

> Three numbers from one unit (`unit.purchase_cost`, `unit.installed_cost`, and the contribution to `TEA.installed_equipment_cost`) don't tally with the user's hand calculation.

Symptom code:

```python
print("purchase_cost: cost of the equipment 'as bought' (FOB).")
print("installed_cost: purchase_cost * F_BM (bare-module factor).")
print("TEA.installed_equipment_cost: sum of installed_cost across units, after Lang factor adjustments.")
```

Why markdown:

> `_cost` populates `purchase_cost`; the bare-module factor `F_BM` (per unit type, often 2-3) inflates to `installed_cost`; the TEA aggregator may add a Lang factor and indirect costs on top. Most "CAPEX doesn't add up" reports come from comparing the wrong two fields.

Fix markdown:

> Match the level: cite `purchase_cost` for "raw equipment" comparisons, `installed_cost` for "delivered and installed," and the TEA aggregate for "total capital investment." Each is correct at its level. See [7. TEA](7_TEA.ipynb) for the full TEA breakdown.

Fix code:

```python
print("for unit in sys.units:")
print("    print(unit.ID, unit.purchase_cost, unit.installed_cost, unit.F_BM)")
```

- [ ] **Step 3: §3.3 — LCA functional unit and dynamic vs static inventory**

Subsection title: `### 3.3. LCA functional unit and dynamic vs static inventory <a class="anchor" id="s3.3"></a>`

Symptom markdown:

> The same system reports global-warming impacts that differ by an order of magnitude depending on the functional unit (FU). The user suspects an LCA-data error.

Symptom code:

```python
print("FU = kg COD removed")
print("FU = m^3 wastewater treated")
print("FU = person-year served")
print("Numerators are identical; denominators differ by 100-1000x; impact per FU swings accordingly.")
```

Why markdown:

> The FU is the denominator of every per-FU impact. Different FUs answer different questions ("per pollutant removed" vs "per service delivered"). Dynamic LCA (impacts over the operating period, possibly time-varying) and static LCA (steady-state average) can also diverge when emissions are concentrated in specific operating modes.

Fix markdown:

> Pick the FU that matches the comparison the analysis is making, and report it explicitly. For WW, "kg COD removed" works for treatment performance comparisons; "person-year served" for service-delivery comparisons. See [8. LCA](8_LCA.ipynb) for the impact-calculation surface.

Fix code:

```python
print("lca.total_impacts['GWP100'] / lca.system.design_flow  # per m^3 treated")
print("lca.total_impacts['GWP100'] / total_COD_removed       # per kg COD removed")
```

- [ ] **Step 4: Execute and clear outputs.**

```bash
jupyter nbconvert --to notebook --execute --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
jupyter nbconvert --clear-output --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
```

Expected: no errors.

- [ ] **Step 5: Commit.**

```bash
git -C QSDsan add docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
git -C QSDsan commit -m "docs(modeling-notes): add §3 TEA and LCA"
```

### Task 17: Add §4 — Behavior inherited from BioSTEAM

- [ ] **Step 1: §4.1 — Pressure and energy balance idealizations**

Subsection title: `### 4.1. Pressure and energy balance idealizations <a class="anchor" id="s4.1"></a>`

Symptom markdown:

> A pumping unit in the system reports zero power demand, or a digester reports an unphysical T after a "warming" step.

Symptom code:

```python
print("(symptom: utility cost = 0 when the unit clearly needs electricity)")
```

Why markdown:

> Many QSDsan units inherit BioSTEAM's defaults: pressure ambient unless set, temperature unchanged unless set, no spontaneous phase change. For most WW analyses this is fine; for pumping headworks, anaerobic gas handling, or thermal hydrolysis, the idealization matters and must be overridden.

Fix markdown:

> Where energy matters, set inlet/outlet pressures or temperatures explicitly, or subclass the unit to add the energy balance. Audit utility costs after every initial system build — a zero is a smell.

Fix code:

```python
print("u.outs[0].P = 2e5  # 2 bar")
print("# or, in a subclass:")
print("def _design(self):")
print("    super()._design()")
print("    self.power_utility(rho*g*h*Q / eta)")
```

`<div class="alert alert-info"><b>Coming from BioSTEAM?</b> Same defaults apply here; the difference is that WW pumping cases are common enough that you will hit this in regular use.</div>`

- [ ] **Step 2: §4.2 — Degrees of freedom and over-constraint**

Subsection title: `### 4.2. Degrees of freedom and over-constraint <a class="anchor" id="s4.2"></a>`

Symptom markdown:

> A user adds a process specification to fix a stream property; the System refuses to converge, or it converges to a value inconsistent with the spec.

Symptom code:

```python
print("(symptom: System with N+1 specs on N DOF stalls or silently picks one)")
```

Why markdown:

> A System has a finite number of design degrees of freedom. Each `Specification` consumes one. Adding more specs than DOF is mathematically inconsistent: the solver either oscillates or honors one spec and quietly violates the others.

Fix markdown:

> Before adding a spec, count DOF (inlets, outlets, parameters with known values). Adding a spec usually means removing a previously-fixed parameter so it becomes the variable the spec drives. The forthcoming Process Specifications tutorial covers this systematically.

Fix code:

```python
print("# Adding a spec on TSS: also unfix the parameter the spec adjusts.")
print("# Don't: fix HRT and fix TSS together if HRT was the only knob TSS responds to.")
```

- [ ] **Step 3: §4.3 — `simulate()` runs more than `_run`**

Subsection title: `### 4.3. simulate() runs more than _run <a class="anchor" id="s4.3"></a>`

Symptom markdown:

> A parameter sweep over a kinetic constant changes the unit's outlet concentrations cycle to cycle, but `installed_cost` is the same every cycle. The user assumes the cost correlation ignores composition.

Symptom code:

```python
print("(symptom: outlet streams change across sweep iterations, installed_cost doesn't)")
```

Why markdown:

> The sweep is calling `unit._run()` (or assigning into the unit's parameters and then reading outlets) without invoking `_design` and `_cost`. `unit.simulate()` calls `_run` then `_summary`, which itself calls `_design` and `_cost`. The cost results from the previous full simulation stay cached otherwise.

Fix markdown:

> Always go through `unit.simulate()` (or `system.simulate()`) when downstream code reads design or cost results. See [tutorial 5 §1.1](5_SanUnit_advanced.ipynb#1.1.-Fundamental-methods) for the full call graph.

Fix code:

```python
print("for k in k_range:")
print("    unit.k = k")
print("    unit.simulate()  # not just unit._run()")
print("    record(unit.installed_cost)")
```

- [ ] **Step 4: Execute and clear outputs.**

```bash
jupyter nbconvert --to notebook --execute --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
jupyter nbconvert --clear-output --inplace QSDsan/docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
```

Expected: no errors.

- [ ] **Step 5: Commit.**

```bash
git -C QSDsan add docs/source/tutorials/14_Modeling_Notes_and_Pitfalls.ipynb
git -C QSDsan commit -m "docs(modeling-notes): add §4 BioSTEAM-inherited behavior"
```

---

## Phase E — Final build, cross-link verification, changelog

### Task 18: Cross-link verification

- [ ] **Step 1: Full docs rebuild.**

```powershell
cd QSDsan\docs ; Remove-Item -Recurse -Force build ; .\make.bat html *>&1 | Tee-Object -FilePath $env:TEMP\build.log
```

Force a clean build by removing the `build/` directory first.

- [ ] **Step 2: Scan for warnings.**

```bash
grep -E "WARNING|ERROR|undefined label|not_included|broken link" "$TEMP/build.log"
```

Expected: zero output. Any line that appears is a regression to fix before the next step.

- [ ] **Step 3: Click-through verification.**

Open the following pages and follow every cross-link the spec promised, confirming each lands on the right anchor:

```bash
start QSDsan/docs/build/html/index.html                           # home -> FAQ grid card
start QSDsan/docs/build/html/faq/index.html                       # FAQ index -> 4 sub-pages
start QSDsan/docs/build/html/faq/errors.html                      # errors -> notebook callout
start QSDsan/docs/build/html/tutorials/14_Modeling_Notes_and_Pitfalls.html  # notebook -> faq/errors callout, tutorial 5 §1.1 link, tutorial 7/8/11 links
start QSDsan/docs/build/html/tutorials/4_SanUnit_basic.html       # tutorial 4 simulate comment (comment is in code, not a hyperlink, but should read clearly)
start QSDsan/docs/build/html/tutorials/5_SanUnit_advanced.html    # tutorial 5 §1.1 -> notebook §4.3 callout
start QSDsan/docs/build/html/FAQ.html                             # redirect to faq/index
```

For each: confirm the rendered links resolve. Document any broken link in a follow-up commit before proceeding.

### Task 19: Changelog and final commit

**Files:**
- Modify: `QSDsan/CHANGELOG.rst`

- [ ] **Step 1: Add a changelog bullet.**

Under the next unreleased version section in `QSDsan/CHANGELOG.rst`, add (under a Docs subsection, creating it if it doesn't exist):

```rst
- Added ``tutorials/14_Modeling_Notes_and_Pitfalls.ipynb`` covering common modeling surprises across streams/components, unit design, TEA/LCA, and BioSTEAM-inherited behavior. Extended tutorial 5 §1.1 with ``simulate()`` / ``_summary`` call-graph mechanics.
- Restructured the FAQ into a four-page subdirectory (``faq/{errors,tips,styling,ai_assisted_coding}.rst``), absorbing the previously-stranded AI-assisted coding page from the tutorials section. Old ``FAQ.html`` bookmark redirects to the new index.
```

Per [feedback_changelog_convention](../../../memory/feedback_changelog_convention.md), EXPOsan has no changelog of its own; only QSDsan's CHANGELOG.rst gets updated.

- [ ] **Step 2: Final commit.**

```bash
git -C QSDsan add CHANGELOG.rst
git -C QSDsan commit -m "docs: log Modeling Notes & Pitfalls + FAQ restructure"
```

- [ ] **Step 3: Review the full diff before handing back.**

```bash
git -C QSDsan log --oneline main..HEAD
git -C QSDsan diff main..HEAD --stat
```

Surface the commit list and file-change summary to the user. Per [feedback_git_workflow](../../../memory/feedback_git_workflow.md), do not push; the user reviews and pushes.

---

## Self-review notes

**Spec coverage check** (against `2026-05-29-modeling-notes-and-pitfalls-design.md`):

- §"Deliverables" item 1 (notebook) — Tasks 12, 14, 15, 16, 17.
- §"Deliverables" item 2 (tutorial 5 §1.1) — Task 1.
- §"Deliverables" item 3 (faq/ subdirectory) — Tasks 3-7.
- §"Deliverables" item 4 (FAQ.rst deletion) — Task 8 (replaced with stub, not deleted, per spec risk note).
- §"Deliverables" item 5 (tutorials/14_AI_Assisted move) — Task 7.
- §"Deliverables" item 6 (tutorials/index.rst toctree) — Task 10.
- §"Deliverables" item 7 (docs/source/index.rst toctree) — Task 9.
- §"Tutorial 5 §1.1" inline-not-renumbered approach — Task 1.
- §"Cross-links" all three — covered in Tasks 1, 3 (notebook ↔ FAQ), 4 (errors → notebook), 12 (notebook → faq/errors).
- §"Verification" steps — Tasks 11, 13, 18.
- §"Risks/tradeoffs" FAQ link-rot stub — Task 8.
- §"Implementation order" — followed verbatim.

**Placeholder scan:** Several Phase D code cells use `print("(symptom: ...)")` as a stand-in where a full runnable reproducer would require setting up a System or dynamic unit. These are intentional: the spec called for "minimal toy snippets" and demanded total notebook execution under 30 s. They communicate the symptom in prose; the user can decide during execution whether to expand any specific entry to a full reproducer.

**Type/name consistency:** `qs.CEPCI` and `qs.CEPCI_by_year` used consistently. `simulate()`/`_summary`/`_run`/`_design`/`_cost` naming matches BioSTEAM/QSDsan conventions. `WasteStream`/`SanStream`/`Stream` consistently capitalized. Notebook entry IDs `s1.1` through `s4.3` consistent with the spec's anchor scheme. Tutorial-5 §1.1 anchor referenced from §4.3 as `5_SanUnit_advanced.ipynb#1.1.-Fundamental-methods` (the nbsphinx case-preserving / dot-preserving / hyphenated-space slug — verify the actual built anchor matches at Task 18).
