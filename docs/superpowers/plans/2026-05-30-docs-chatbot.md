# QSDsan + EXPOsan Docs Chatbot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a self-hosted RAG chatbot widget embedded on the QSDsan readthedocs site that answers questions about QSDsan (API + tutorials) and EXPOsan (example systems) strictly from indexed docs, with citations and out-of-scope refusal.

**Architecture:** Four loosely-coupled units under the QSDsan repo. An *indexer* (two source adapters → one chunk schema with embeddings → `index.json`), a *FastAPI endpoint* (embed question → cosine top-k → Claude with cached system prompt → cited answer), a *Furo-injected widget* (floating button + chat panel), and a *readthedocs post_build hook* that reindexes on every build. Generation orchestration lives in a pure `engine.py` so it is testable with mocked Claude/Voyage clients.

**Tech Stack:** Python 3.12, FastAPI + uvicorn, Anthropic SDK (Claude Haiku 4.5, config-switchable to Sonnet, prompt caching), Voyage `voyage-4-lite` embeddings, numpy (flat cosine search), `requests`, `docutils`/`BeautifulSoup` for parsing, pytest, vanilla JS/CSS widget, Render (host) + readthedocs (reindex hook).

---

## Conventions for every task

- **TDD:** failing test first, watch it fail, minimal implementation, watch it pass, commit.
- **Test command (from repo root):** `pytest docs/chatbot/tests/<file> -v`
- **Scope:** touch only files under `docs/chatbot/`, `docs/source/_static/`, `docs/source/conf.py`, `.readthedocs.yml`, `render.yaml`, and `.gitignore`. NEVER touch `qsdsan/_tea.py` or `tests/test_tea.py` (parallel `beta` work in the other worktree).
- **No network in tests:** all Voyage/Anthropic/GitHub calls are injected and mocked.
- **Commit message prefix:** `chatbot:` on every commit.
- **Show the diff and ask before committing.** Never push without explicit approval.
- Minimize em-dashes in user-facing prose (widget strings, refusal text, README).

## File structure (what each file owns)

| File | Responsibility |
|------|----------------|
| `docs/chatbot/DESIGN.md` | Spec (already written). |
| `docs/chatbot/__init__.py` | Marks the dir importable in tests; empty. |
| `docs/chatbot/config.py` | All tunables (models, top-k, threshold, URLs) read from env with defaults. Single source of truth. |
| `docs/chatbot/chunking.py` | Pure heading-splitters: `split_rst_by_heading`, `split_html_by_heading`. Shared by both adapters (DRY). |
| `docs/chatbot/index_docs.py` | Two adapters (`build_exposan_chunks`, `build_qsdsan_chunks`) + `build_index` (adds embeddings) + `main()` CLI. Writes `index.json`. |
| `docs/chatbot/embeddings.py` | Voyage wrapper `embed_texts(texts, input_type, client=None)`. Used by indexer (documents) and server (query). |
| `docs/chatbot/retrieval.py` | `load_index(path_or_url)`, `cosine_top_k(query_vec, records, k)`. Pure numpy. |
| `docs/chatbot/prompts.py` | System prompt, `build_user_prompt`, `append_code_disclaimers`, refusal text. Guardrail behavior. |
| `docs/chatbot/engine.py` | `answer_question(...)` orchestration (embed → retrieve → threshold/refuse → Claude → cite). Pure, dependency-injected. |
| `docs/chatbot/server.py` | Thin FastAPI app: loads index at startup, wires real clients into `engine.answer_question`. |
| `docs/chatbot/requirements.txt` | Endpoint runtime deps. |
| `docs/chatbot/README.md` | How to run/deploy/reindex. |
| `docs/chatbot/eval/questions.yaml` | Seeded regression questions + expectations. |
| `docs/chatbot/eval/run_eval.py` | Runs questions against a live `/ask`, asserts expectations. |
| `docs/chatbot/tests/` | All pytest tests + fixtures. Self-contained. |
| `docs/source/_static/js/chatbot.js` | Widget behavior. |
| `docs/source/_static/css/chatbot.css` | Widget styling (Furo light/dark aware). |
| `docs/source/conf.py` | Register widget JS/CSS (modify lines 99-106). |
| `render.yaml` | Render Blueprint (rootDir `docs/chatbot`). |
| `.readthedocs.yml` | Add `build.jobs.post_build` reindex step. |
| `.gitignore` | Already ignores `index.json` (line 29); no change needed. |

## Shared data contracts (used across tasks — keep names exact)

**Chunk record** (one entry in `index.json`):
```python
{
  "text": str,        # the chunk body
  "title": str,       # heading text
  "url": str,         # resolvable citation URL
  "source": str,      # "qsdsan" | "exposan"
  "type": str,        # "tutorial" | "api" | "readme"
  "embedding": list[float],  # added by build_index; absent on raw adapter output
}
```

**Heading split** (output of `chunking` functions): `{"title": str, "text": str, "anchor": str}`.

**`answer_question` return:** `{"answer": str, "citations": list[dict]}` where each citation is `{"n": int, "url": str, "title": str, "source": str}`.

---

### Task 1: Package scaffold + config

**Files:**
- Create: `docs/chatbot/__init__.py`
- Create: `docs/chatbot/config.py`
- Create: `docs/chatbot/tests/__init__.py`
- Create: `docs/chatbot/tests/conftest.py`
- Test: `docs/chatbot/tests/test_config.py`

- [ ] **Step 1: Create the empty package marker and the test path bootstrap**

`docs/chatbot/__init__.py`: empty file.

`docs/chatbot/tests/__init__.py`: empty file.

`docs/chatbot/tests/conftest.py`:
```python
import sys
from pathlib import Path

# Make modules in docs/chatbot importable as top-level (import config, chunking, ...)
CHATBOT_DIR = Path(__file__).resolve().parents[1]
if str(CHATBOT_DIR) not in sys.path:
    sys.path.insert(0, str(CHATBOT_DIR))
```

- [ ] **Step 2: Write the failing test**

`docs/chatbot/tests/test_config.py`:
```python
import importlib


def test_defaults_present():
    config = importlib.import_module("config")
    assert config.GEN_MODEL == "claude-haiku-4-5-20251001"
    assert config.EMBED_MODEL == "voyage-4-lite"
    assert config.TOP_K == 6
    assert 0.0 < config.SIMILARITY_THRESHOLD < 1.0
    assert config.QSDSAN_DOCS_BASE.startswith("https://qsdsan.readthedocs.io")
    assert config.EXPOSAN_RAW_BASE.startswith("https://raw.githubusercontent.com/QSD-Group/EXPOsan")
    assert config.EXPOSAN_BLOB_BASE.startswith("https://github.com/QSD-Group/EXPOsan")


def test_env_override(monkeypatch):
    monkeypatch.setenv("CHATBOT_GEN_MODEL", "claude-sonnet-4-6")
    import config
    importlib.reload(config)
    assert config.GEN_MODEL == "claude-sonnet-4-6"
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_config.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'config'`).

- [ ] **Step 4: Write minimal implementation**

`docs/chatbot/config.py`:
```python
"""Central configuration for the QSDsan/EXPOsan docs chatbot.

Every tunable is read from an environment variable with a safe default so the
same code runs locally, on Render, and in the readthedocs build.
"""
import os

# --- Models ---
GEN_MODEL = os.getenv("CHATBOT_GEN_MODEL", "claude-haiku-4-5-20251001")
EMBED_MODEL = os.getenv("CHATBOT_EMBED_MODEL", "voyage-4-lite")

# --- Retrieval ---
TOP_K = int(os.getenv("CHATBOT_TOP_K", "6"))
SIMILARITY_THRESHOLD = float(os.getenv("CHATBOT_SIMILARITY_THRESHOLD", "0.4"))

# --- Source locations ---
QSDSAN_DOCS_BASE = os.getenv(
    "CHATBOT_QSDSAN_DOCS_BASE", "https://qsdsan.readthedocs.io/en/latest/"
)
EXPOSAN_RAW_BASE = os.getenv(
    "CHATBOT_EXPOSAN_RAW_BASE",
    "https://raw.githubusercontent.com/QSD-Group/EXPOsan/main/exposan/",
)
EXPOSAN_BLOB_BASE = os.getenv(
    "CHATBOT_EXPOSAN_BLOB_BASE",
    "https://github.com/QSD-Group/EXPOsan/blob/main/exposan/",
)

# Where the published index lives; the server loads this on startup.
INDEX_URL = os.getenv(
    "CHATBOT_INDEX_URL",
    "https://qsdsan.readthedocs.io/en/latest/_static/chatbot/index.json",
)

# Docs search fallback used in the refusal message.
DOCS_SEARCH_URL = os.getenv(
    "CHATBOT_DOCS_SEARCH_URL",
    "https://qsdsan.readthedocs.io/en/latest/search.html?q=",
)
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_config.py -v`
Expected: PASS (2 passed).

- [ ] **Step 6: Commit**

```bash
git add docs/chatbot/__init__.py docs/chatbot/config.py docs/chatbot/tests/__init__.py docs/chatbot/tests/conftest.py docs/chatbot/tests/test_config.py
git commit -m "chatbot: add package scaffold and config"
```

---

### Task 2: RST heading splitter

**Files:**
- Create: `docs/chatbot/chunking.py`
- Test: `docs/chatbot/tests/test_chunking_rst.py`

- [ ] **Step 1: Write the failing test**

`docs/chatbot/tests/test_chunking_rst.py`:
```python
import chunking

SAMPLE_RST = """\
=========
BSM1
=========

Summary
-------
This system models BSM1.

Load and run
------------
Use ``create_system`` to build it.
"""


def test_splits_each_heading_into_a_chunk():
    chunks = chunking.split_rst_by_heading(SAMPLE_RST)
    titles = [c["title"] for c in chunks]
    assert "Summary" in titles
    assert "Load and run" in titles


def test_chunk_text_contains_body_under_heading():
    chunks = chunking.split_rst_by_heading(SAMPLE_RST)
    load = next(c for c in chunks if c["title"] == "Load and run")
    assert "create_system" in load["text"]
    assert "models BSM1" not in load["text"]  # body of a different section


def test_anchor_is_github_style_slug():
    chunks = chunking.split_rst_by_heading(SAMPLE_RST)
    load = next(c for c in chunks if c["title"] == "Load and run")
    assert load["anchor"] == "load-and-run"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_chunking_rst.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'chunking'`).

- [ ] **Step 3: Write minimal implementation**

`docs/chatbot/chunking.py`:
```python
"""Pure heading-based chunking for RST source and built HTML.

Both indexer adapters reuse these so chunk boundaries stay consistent.
"""
import re

# RST section adornment characters per the docutils convention.
_ADORNMENT = set("=-~^\"'`#*+.:_")


def slugify(title: str) -> str:
    """GitHub/Sphinx-style anchor slug: lowercase, non-alphanumerics to hyphens."""
    slug = re.sub(r"[^a-z0-9]+", "-", title.strip().lower())
    return slug.strip("-")


def split_rst_by_heading(text: str) -> list[dict]:
    """Split RST into [{title, text, anchor}] by underlined (and overlined) headings.

    A heading is a non-blank title line whose following line is an adornment run
    of one repeated character at least as long as the title.
    """
    lines = text.splitlines()
    sections: list[dict] = []
    current = None
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        nxt = lines[i + 1] if i + 1 < n else ""
        is_underline = (
            line.strip()
            and nxt.strip()
            and len(set(nxt.strip())) == 1
            and nxt.strip()[0] in _ADORNMENT
            and len(nxt.strip()) >= len(line.strip())
        )
        # Skip an over-line directly above the title (e.g. ==== / BSM1 / ====).
        prev_is_overline = (
            i > 0
            and lines[i - 1].strip()
            and len(set(lines[i - 1].strip())) == 1
            and lines[i - 1].strip()[0] in _ADORNMENT
        )
        if is_underline:
            title = line.strip()
            if current is not None:
                sections.append(current)
            current = {"title": title, "text": "", "anchor": slugify(title)}
            i += 2  # consume title + underline
            continue
        if prev_is_overline and not line.strip():
            i += 1
            continue
        if current is not None:
            current["text"] += line + "\n"
        i += 1
    if current is not None:
        sections.append(current)
    for s in sections:
        s["text"] = s["text"].strip()
    return sections
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_chunking_rst.py -v`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add docs/chatbot/chunking.py docs/chatbot/tests/test_chunking_rst.py
git commit -m "chatbot: add RST heading splitter"
```

---

### Task 3: HTML heading splitter

**Files:**
- Modify: `docs/chatbot/chunking.py` (add `split_html_by_heading`)
- Modify: `docs/chatbot/requirements.txt` (created here for the first time)
- Test: `docs/chatbot/tests/test_chunking_html.py`

- [ ] **Step 1: Create requirements.txt (needed for BeautifulSoup)**

`docs/chatbot/requirements.txt`:
```text
fastapi>=0.110
uvicorn[standard]>=0.29
anthropic>=0.40
voyageai>=0.3
numpy>=1.26
requests>=2.31
beautifulsoup4>=4.12
pydantic>=2.6
pyyaml>=6.0
```

Install locally for development:
```bash
pip install -r docs/chatbot/requirements.txt
```

- [ ] **Step 2: Write the failing test**

`docs/chatbot/tests/test_chunking_html.py`:
```python
import chunking

SAMPLE_HTML = """
<html><body>
<section id="create-a-wastestream">
  <h2>Create a WasteStream</h2>
  <p>Use <code>WasteStream</code> to define flow and COD.</p>
</section>
<section id="run-a-simulation">
  <h2>Run a simulation</h2>
  <p>Call <code>sys.simulate()</code> to run it.</p>
</section>
</body></html>
"""


def test_splits_html_sections_by_heading():
    chunks = chunking.split_html_by_heading(SAMPLE_HTML)
    titles = [c["title"] for c in chunks]
    assert "Create a WasteStream" in titles
    assert "Run a simulation" in titles


def test_html_anchor_comes_from_section_id():
    chunks = chunking.split_html_by_heading(SAMPLE_HTML)
    ws = next(c for c in chunks if c["title"] == "Create a WasteStream")
    assert ws["anchor"] == "create-a-wastestream"


def test_html_text_is_plain_no_tags():
    chunks = chunking.split_html_by_heading(SAMPLE_HTML)
    ws = next(c for c in chunks if c["title"] == "Create a WasteStream")
    assert "flow and COD" in ws["text"]
    assert "<code>" not in ws["text"]
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_chunking_html.py -v`
Expected: FAIL (`AttributeError: module 'chunking' has no attribute 'split_html_by_heading'`).

- [ ] **Step 4: Add the implementation**

Append to `docs/chatbot/chunking.py`:
```python
from bs4 import BeautifulSoup

_HEADING_TAGS = ("h1", "h2", "h3", "h4")


def split_html_by_heading(html: str) -> list[dict]:
    """Split built HTML into [{title, text, anchor}] by <section> heading.

    Anchor is the nearest enclosing element id (Sphinx/Furo put the section id on
    the wrapping <section>/<div>), falling back to a slug of the title.
    """
    soup = BeautifulSoup(html, "html.parser")
    chunks: list[dict] = []
    for heading in soup.find_all(_HEADING_TAGS):
        title = heading.get_text(strip=True)
        if not title:
            continue
        # Find the id on the heading or the nearest ancestor section/div.
        anchor = heading.get("id")
        if not anchor:
            container = heading.find_parent(
                lambda tag: tag.name in ("section", "div") and tag.get("id")
            )
            anchor = container.get("id") if container else slugify(title)
        # Body text = the section's text minus the heading itself.
        container = heading.find_parent("section") or heading.parent
        text = container.get_text(" ", strip=True)
        text = text.replace(title, "", 1).strip()
        chunks.append({"title": title, "text": text, "anchor": anchor})
    return chunks
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_chunking_html.py -v`
Expected: PASS (3 passed).

- [ ] **Step 6: Commit**

```bash
git add docs/chatbot/chunking.py docs/chatbot/requirements.txt docs/chatbot/tests/test_chunking_html.py
git commit -m "chatbot: add HTML heading splitter and runtime requirements"
```

---

### Task 4: EXPOsan adapter (GitHub-raw fetch + parse, mocked)

**Files:**
- Create: `docs/chatbot/index_docs.py`
- Create: `docs/chatbot/tests/fixtures/bsm1_README.rst`
- Test: `docs/chatbot/tests/test_adapter_exposan.py`

- [ ] **Step 1: Create the README fixture**

`docs/chatbot/tests/fixtures/bsm1_README.rst`:
```rst
====
bsm1
====

Summary
-------
Benchmark Simulation Model No. 1 for activated sludge.

Load and run
------------
Run ``from exposan.bsm1 import create_system; sys = create_system()`` to build
and then ``sys.simulate(t_span=(0, 10))``.
```

- [ ] **Step 2: Write the failing test**

`docs/chatbot/tests/test_adapter_exposan.py`:
```python
from pathlib import Path

import index_docs

FIXTURE = (Path(__file__).parent / "fixtures" / "bsm1_README.rst").read_text()


def fake_fetch(url):
    # Only the bsm1 README exists in this fake world; everything else 404s (None).
    if url.endswith("exposan/bsm1/README.rst"):
        return FIXTURE
    return None


def test_builds_chunks_only_for_systems_with_readme():
    chunks = index_docs.build_exposan_chunks(
        systems=["bsm1", "ghost"], fetch=fake_fetch
    )
    assert chunks, "expected at least one chunk"
    assert all(c["source"] == "exposan" for c in chunks)
    assert all(c["type"] == "readme" for c in chunks)
    # 'ghost' has no README, so no chunks reference it.
    assert not any("ghost" in c["url"] for c in chunks)


def test_exposan_url_is_resolvable_github_blob():
    chunks = index_docs.build_exposan_chunks(systems=["bsm1"], fetch=fake_fetch)
    run = next(c for c in chunks if c["title"] == "Load and run")
    assert run["url"].startswith(
        "https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bsm1/README.rst"
    )
    assert run["url"].endswith("#load-and-run")


def test_every_chunk_has_text_and_title():
    chunks = index_docs.build_exposan_chunks(systems=["bsm1"], fetch=fake_fetch)
    for c in chunks:
        assert c["title"]
        assert c["text"]
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_adapter_exposan.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'index_docs'`).

- [ ] **Step 4: Write minimal implementation**

`docs/chatbot/index_docs.py`:
```python
"""Indexer: build chunk records from QSDsan built HTML and EXPOsan READMEs.

Source-pluggable. Both adapters emit the same chunk schema; build_index adds
embeddings and writes index.json.
"""
from __future__ import annotations

import config
import chunking


def _http_get(url: str) -> str | None:
    """Default fetcher: return text on 200, None on 404. Network only at runtime."""
    import requests

    resp = requests.get(url, timeout=30)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    return resp.text


def build_exposan_chunks(systems: list[str], fetch=_http_get) -> list[dict]:
    """Fetch each system's README.rst from GitHub raw and chunk by heading.

    Systems whose README is missing (fetch returns None) are skipped, so the
    list can be a superset of what actually exists.
    """
    chunks: list[dict] = []
    for sys_name in systems:
        raw_url = f"{config.EXPOSAN_RAW_BASE}{sys_name}/README.rst"
        text = fetch(raw_url)
        if not text:
            continue
        blob_url = f"{config.EXPOSAN_BLOB_BASE}{sys_name}/README.rst"
        for section in chunking.split_rst_by_heading(text):
            if not section["text"]:
                continue
            chunks.append(
                {
                    "text": section["text"],
                    "title": section["title"],
                    "url": f"{blob_url}#{section['anchor']}",
                    "source": "exposan",
                    "type": "readme",
                }
            )
    return chunks
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_adapter_exposan.py -v`
Expected: PASS (3 passed).

- [ ] **Step 6: Commit**

```bash
git add docs/chatbot/index_docs.py docs/chatbot/tests/fixtures/bsm1_README.rst docs/chatbot/tests/test_adapter_exposan.py
git commit -m "chatbot: add EXPOsan README adapter"
```

---

### Task 5: QSDsan adapter (built-HTML walk + parse)

**Files:**
- Modify: `docs/chatbot/index_docs.py` (add `build_qsdsan_chunks`)
- Create: `docs/chatbot/tests/fixtures/html/tutorials/3_WasteStream.html`
- Create: `docs/chatbot/tests/fixtures/html/api/streams.html`
- Test: `docs/chatbot/tests/test_adapter_qsdsan.py`

- [ ] **Step 1: Create the built-HTML fixtures**

`docs/chatbot/tests/fixtures/html/tutorials/3_WasteStream.html`:
```html
<html><body><section id="wastestream">
  <h1>WasteStream</h1>
  <section id="create-a-wastestream">
    <h2>Create a WasteStream</h2>
    <p>Use WasteStream to set flow and COD.</p>
  </section>
</section></body></html>
```

`docs/chatbot/tests/fixtures/html/api/streams.html`:
```html
<html><body><section id="qsdsan-sanstream-wastestream">
  <h1>WasteStream API</h1>
  <p>Class reference for WasteStream.</p>
</section></body></html>
```

- [ ] **Step 2: Write the failing test**

`docs/chatbot/tests/test_adapter_qsdsan.py`:
```python
from pathlib import Path

import index_docs

HTML_DIR = Path(__file__).parent / "fixtures" / "html"


def test_builds_chunks_from_built_html_tree():
    chunks = index_docs.build_qsdsan_chunks(str(HTML_DIR))
    assert chunks
    assert all(c["source"] == "qsdsan" for c in chunks)


def test_tutorial_pages_tagged_tutorial_api_pages_tagged_api():
    chunks = index_docs.build_qsdsan_chunks(str(HTML_DIR))
    tut = next(c for c in chunks if "flow and COD" in c["text"])
    api = next(c for c in chunks if "Class reference" in c["text"])
    assert tut["type"] == "tutorial"
    assert api["type"] == "api"


def test_qsdsan_url_is_readthedocs_page_plus_anchor():
    chunks = index_docs.build_qsdsan_chunks(str(HTML_DIR))
    ws = next(c for c in chunks if c["title"] == "Create a WasteStream")
    assert ws["url"] == (
        "https://qsdsan.readthedocs.io/en/latest/"
        "tutorials/3_WasteStream.html#create-a-wastestream"
    )
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_adapter_qsdsan.py -v`
Expected: FAIL (`AttributeError: module 'index_docs' has no attribute 'build_qsdsan_chunks'`).

- [ ] **Step 4: Add the implementation**

Append to `docs/chatbot/index_docs.py`:
```python
import os


def _page_type(rel_path: str) -> str:
    """Tag a built page by its location under the docs tree."""
    head = rel_path.replace(os.sep, "/").split("/", 1)[0]
    if head == "tutorials":
        return "tutorial"
    if head == "api":
        return "api"
    return "tutorial"


def build_qsdsan_chunks(html_dir: str, base_url: str | None = None) -> list[dict]:
    """Walk a built Sphinx HTML tree and chunk each page by heading.

    Citation URL = readthedocs base + page path relative to html_dir + #anchor.
    Only files under tutorials/ and api/ are indexed.
    """
    base_url = base_url or config.QSDSAN_DOCS_BASE
    chunks: list[dict] = []
    for root, _dirs, files in os.walk(html_dir):
        for fname in files:
            if not fname.endswith(".html"):
                continue
            abs_path = os.path.join(root, fname)
            rel_path = os.path.relpath(abs_path, html_dir).replace(os.sep, "/")
            head = rel_path.split("/", 1)[0]
            if head not in ("tutorials", "api"):
                continue
            ptype = _page_type(rel_path)
            html = open(abs_path, encoding="utf-8").read()
            page_url = f"{base_url.rstrip('/')}/{rel_path}"
            for section in chunking.split_html_by_heading(html):
                if not section["text"]:
                    continue
                chunks.append(
                    {
                        "text": section["text"],
                        "title": section["title"],
                        "url": f"{page_url}#{section['anchor']}",
                        "source": "qsdsan",
                        "type": ptype,
                    }
                )
    return chunks
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_adapter_qsdsan.py -v`
Expected: PASS (3 passed).

- [ ] **Step 6: Commit**

```bash
git add docs/chatbot/index_docs.py docs/chatbot/tests/fixtures/html docs/chatbot/tests/test_adapter_qsdsan.py
git commit -m "chatbot: add QSDsan built-HTML adapter"
```

---

### Task 6: Voyage embeddings wrapper

**Files:**
- Create: `docs/chatbot/embeddings.py`
- Test: `docs/chatbot/tests/test_embeddings.py`

- [ ] **Step 1: Write the failing test**

`docs/chatbot/tests/test_embeddings.py`:
```python
import embeddings


class FakeResult:
    def __init__(self, vectors):
        self.embeddings = vectors


class FakeClient:
    def __init__(self):
        self.calls = []

    def embed(self, texts, model, input_type):
        self.calls.append((tuple(texts), model, input_type))
        return FakeResult([[0.1, 0.2, 0.3] for _ in texts])


def test_embed_texts_returns_one_vector_per_input():
    client = FakeClient()
    vecs = embeddings.embed_texts(["a", "b"], input_type="document", client=client)
    assert len(vecs) == 2
    assert vecs[0] == [0.1, 0.2, 0.3]


def test_embed_texts_passes_model_and_input_type():
    client = FakeClient()
    embeddings.embed_texts(["q"], input_type="query", client=client)
    (_texts, model, input_type) = client.calls[0]
    assert model == "voyage-4-lite"
    assert input_type == "query"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_embeddings.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'embeddings'`).

- [ ] **Step 3: Write minimal implementation**

`docs/chatbot/embeddings.py`:
```python
"""Thin Voyage embedding wrapper used by the indexer and the server.

input_type is "document" when embedding chunks, "query" when embedding a question.
"""
import config


def _default_client():
    import voyageai

    return voyageai.Client()  # reads VOYAGE_API_KEY from the environment


def embed_texts(texts: list[str], input_type: str, client=None) -> list[list[float]]:
    """Return one embedding vector per input string."""
    if not texts:
        return []
    client = client or _default_client()
    result = client.embed(texts, model=config.EMBED_MODEL, input_type=input_type)
    return list(result.embeddings)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_embeddings.py -v`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add docs/chatbot/embeddings.py docs/chatbot/tests/test_embeddings.py
git commit -m "chatbot: add Voyage embeddings wrapper"
```

---

### Task 7: Indexer orchestration + CLI (`build_index`, `main`)

**Files:**
- Modify: `docs/chatbot/index_docs.py` (add `build_index`, `main`)
- Test: `docs/chatbot/tests/test_build_index.py`

- [ ] **Step 1: Write the failing test**

`docs/chatbot/tests/test_build_index.py`:
```python
import json
from pathlib import Path

import index_docs

HTML_DIR = Path(__file__).parent / "fixtures" / "html"
FIXTURE = (Path(__file__).parent / "fixtures" / "bsm1_README.rst").read_text()


def fake_fetch(url):
    if url.endswith("exposan/bsm1/README.rst"):
        return FIXTURE
    return None


def fake_embed(texts, input_type, client=None):
    assert input_type == "document"
    return [[float(len(t)), 1.0] for t in texts]


def test_build_index_merges_sources_and_adds_embeddings():
    records = index_docs.build_index(
        html_dir=str(HTML_DIR),
        systems=["bsm1"],
        fetch=fake_fetch,
        embed_fn=fake_embed,
    )
    sources = {r["source"] for r in records}
    assert sources == {"qsdsan", "exposan"}
    assert all(len(r["embedding"]) == 2 for r in records)


def test_main_writes_index_json(tmp_path, monkeypatch):
    out = tmp_path / "index.json"
    monkeypatch.setattr(index_docs, "_http_get", fake_fetch)
    monkeypatch.setattr(index_docs, "embed_documents", fake_embed)
    index_docs.main(
        html_dir=str(HTML_DIR),
        systems=["bsm1"],
        out_path=str(out),
    )
    data = json.loads(out.read_text())
    assert isinstance(data, list) and data
    assert {r["source"] for r in data} == {"qsdsan", "exposan"}
    assert all("embedding" in r for r in data)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_build_index.py -v`
Expected: FAIL (`AttributeError: module 'index_docs' has no attribute 'build_index'`).

- [ ] **Step 3: Add the implementation**

Append to `docs/chatbot/index_docs.py`:
```python
import json

import embeddings

# Candidate EXPOsan systems. The adapter skips any without a README, so this can
# safely be a superset; it is refreshed by simply editing this list.
EXPOSAN_SYSTEMS = [
    "adm", "asm", "biobinder", "biogenic_refinery", "bsm1", "bsm2", "bwaise",
    "cas", "eco_san", "hap", "htl", "metab", "metro", "new_generator",
    "pm2_batch", "pm2_ecorecover", "pou_disinfection", "reclaimer", "saf", "werf",
]


def embed_documents(texts, input_type, client=None):
    """Indirection point so tests can monkeypatch embedding without network."""
    return embeddings.embed_texts(texts, input_type=input_type, client=client)


def build_index(html_dir, systems=None, fetch=_http_get, embed_fn=None) -> list[dict]:
    """Build merged chunk records from both adapters and attach embeddings."""
    systems = systems if systems is not None else EXPOSAN_SYSTEMS
    embed_fn = embed_fn or embed_documents
    records = build_qsdsan_chunks(html_dir) + build_exposan_chunks(systems, fetch=fetch)
    if records:
        vectors = embed_fn([r["text"] for r in records], input_type="document")
        for r, vec in zip(records, vectors):
            r["embedding"] = list(vec)
    return records


def main(html_dir=None, systems=None, out_path=None) -> None:
    """CLI entry: build the index and write it as JSON.

    Defaults resolve build-time locations: $READTHEDOCS_OUTPUT/html for the built
    docs, and the in-tree _static path for the output.
    """
    import os

    html_dir = html_dir or os.path.join(
        os.environ.get("READTHEDOCS_OUTPUT", "docs/build"), "html"
    )
    out_path = out_path or "docs/source/_static/chatbot/index.json"
    records = build_index(html_dir, systems=systems, fetch=_http_get, embed_fn=embed_documents)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(records, fh)
    print(f"Wrote {len(records)} chunks to {out_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_build_index.py -v`
Expected: PASS (2 passed).

- [ ] **Step 5: Run the full chatbot suite so far**

Run: `pytest docs/chatbot/tests -v`
Expected: PASS (all green).

- [ ] **Step 6: Commit**

```bash
git add docs/chatbot/index_docs.py docs/chatbot/tests/test_build_index.py
git commit -m "chatbot: add indexer orchestration and CLI"
```

---

### Task 8: Retrieval (load index + cosine top-k)

**Files:**
- Create: `docs/chatbot/retrieval.py`
- Test: `docs/chatbot/tests/test_retrieval.py`

- [ ] **Step 1: Write the failing test**

`docs/chatbot/tests/test_retrieval.py`:
```python
import json

import retrieval

RECORDS = [
    {"text": "a", "title": "A", "url": "u1", "source": "qsdsan", "type": "api",
     "embedding": [1.0, 0.0]},
    {"text": "b", "title": "B", "url": "u2", "source": "exposan", "type": "readme",
     "embedding": [0.0, 1.0]},
    {"text": "c", "title": "C", "url": "u3", "source": "qsdsan", "type": "tutorial",
     "embedding": [0.9, 0.1]},
]


def test_top_k_orders_by_cosine_similarity():
    hits = retrieval.cosine_top_k([1.0, 0.0], RECORDS, k=2)
    assert [r["title"] for r, _score in hits] == ["A", "C"]
    assert hits[0][1] > hits[1][1]


def test_top_score_is_one_for_identical_vector():
    hits = retrieval.cosine_top_k([0.0, 1.0], RECORDS, k=1)
    record, score = hits[0]
    assert record["title"] == "B"
    assert abs(score - 1.0) < 1e-6


def test_load_index_reads_json_file(tmp_path):
    p = tmp_path / "index.json"
    p.write_text(json.dumps(RECORDS))
    loaded = retrieval.load_index(str(p))
    assert len(loaded) == 3
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_retrieval.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'retrieval'`).

- [ ] **Step 3: Write minimal implementation**

`docs/chatbot/retrieval.py`:
```python
"""Load the index and do flat cosine top-k. Corpus is small; no vector DB."""
import json

import numpy as np


def load_index(path_or_url: str) -> list[dict]:
    """Load index.json from a local path or an http(s) URL."""
    if path_or_url.startswith(("http://", "https://")):
        import requests

        resp = requests.get(path_or_url, timeout=30)
        resp.raise_for_status()
        return resp.json()
    with open(path_or_url, encoding="utf-8") as fh:
        return json.load(fh)


def cosine_top_k(query_vec, records: list[dict], k: int) -> list[tuple[dict, float]]:
    """Return the top-k (record, similarity) pairs by cosine similarity."""
    if not records:
        return []
    matrix = np.array([r["embedding"] for r in records], dtype=float)
    query = np.array(query_vec, dtype=float)
    norms = np.linalg.norm(matrix, axis=1) * np.linalg.norm(query)
    norms[norms == 0] = 1e-12
    sims = matrix @ query / norms
    order = np.argsort(-sims)[:k]
    return [(records[i], float(sims[i])) for i in order]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_retrieval.py -v`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add docs/chatbot/retrieval.py docs/chatbot/tests/test_retrieval.py
git commit -m "chatbot: add cosine top-k retrieval"
```

---

### Task 9: Prompts and guardrail text

**Files:**
- Create: `docs/chatbot/prompts.py`
- Test: `docs/chatbot/tests/test_prompts.py`

- [ ] **Step 1: Write the failing test**

`docs/chatbot/tests/test_prompts.py`:
```python
import prompts

RETRIEVED = [
    {"text": "Use WasteStream to set flow and COD.", "title": "Create a WasteStream",
     "url": "https://qsdsan.readthedocs.io/en/latest/tutorials/3_WasteStream.html#x",
     "source": "qsdsan", "type": "tutorial"},
]


def test_system_prompt_states_grounding_and_citation_rules():
    text = prompts.SYSTEM_PROMPT.lower()
    assert "only" in text and "cite" in text
    assert "qsdsan" in text and "exposan" in text


def test_user_prompt_numbers_excerpts_and_includes_urls():
    up = prompts.build_user_prompt("How do I make a WasteStream?", RETRIEVED)
    assert "[1]" in up
    assert RETRIEVED[0]["url"] in up
    assert "How do I make a WasteStream?" in up


def test_append_code_disclaimers_adds_note_after_each_code_block():
    md = "Intro\n\n```python\nx = 1\n```\n\nOutro"
    out = prompts.append_code_disclaimers(md)
    assert out.count("Draft from the docs") == 1
    assert "verify against the linked pages" in out


def test_refusal_message_links_docs_search():
    msg = prompts.refusal_message()
    assert "couldn't find this in the QSDsan/EXPOsan docs" in msg
    assert "search.html" in msg
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_prompts.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'prompts'`).

- [ ] **Step 3: Write minimal implementation**

`docs/chatbot/prompts.py`:
```python
"""System/user prompts and guardrail text. Behavior of the bot lives here."""
import re

import config

SYSTEM_PROMPT = (
    "You are the QSDsan and EXPOsan documentation assistant. "
    "Answer ONLY using the numbered excerpts provided in the user message. "
    "Cite every claim with its excerpt number and URL, like [1]. "
    "If the excerpts do not contain the answer, say you could not find it in the "
    "QSDsan/EXPOsan docs and do not guess. "
    "Refuse questions unrelated to QSDsan or EXPOsan. "
    "When you include code, keep it minimal and faithful to the excerpts."
)

CODE_DISCLAIMER = "Draft from the docs - verify against the linked pages."


def build_user_prompt(question: str, retrieved: list[dict]) -> str:
    """Render numbered excerpts followed by the question."""
    lines = ["Excerpts:"]
    for i, chunk in enumerate(retrieved, start=1):
        lines.append(
            f"[{i}] ({chunk['source']}) {chunk['title']} - {chunk['url']}\n{chunk['text']}"
        )
    lines.append("")
    lines.append(f"Question: {question}")
    return "\n".join(lines)


def append_code_disclaimers(answer_md: str) -> str:
    """Append the verify-against-docs note right after each fenced code block."""
    pattern = re.compile(r"```.*?```", re.DOTALL)

    def repl(match):
        return match.group(0) + f"\n\n*{CODE_DISCLAIMER}*"

    return pattern.sub(repl, answer_md)


def refusal_message() -> str:
    """Out-of-scope / low-similarity refusal with a docs-search link."""
    return (
        "I couldn't find this in the QSDsan/EXPOsan docs. "
        f"Try the docs search: {config.DOCS_SEARCH_URL}"
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_prompts.py -v`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add docs/chatbot/prompts.py docs/chatbot/tests/test_prompts.py
git commit -m "chatbot: add prompts and guardrail text"
```

---

### Task 10: Engine orchestration (grounding / citation / refusal contract)

**Files:**
- Create: `docs/chatbot/engine.py`
- Test: `docs/chatbot/tests/test_engine.py`

- [ ] **Step 1: Write the failing test**

`docs/chatbot/tests/test_engine.py`:
```python
import engine

RECORDS = [
    {"text": "Use WasteStream to set flow and COD.", "title": "Create a WasteStream",
     "url": "https://qsdsan.readthedocs.io/en/latest/tutorials/3_WasteStream.html#x",
     "source": "qsdsan", "type": "tutorial", "embedding": [1.0, 0.0]},
    {"text": "Run create_system then sys.simulate().", "title": "Load and run",
     "url": "https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bsm1/README.rst#load-and-run",
     "source": "exposan", "type": "readme", "embedding": [0.0, 1.0]},
]


def embed_query_high(texts, input_type, client=None):
    # Aligns with the first record (cosine ~1.0).
    return [[1.0, 0.0]]


def embed_query_low(texts, input_type, client=None):
    # Orthogonal to everything (cosine 0 -> below threshold).
    return [[1.0, 1.0]]  # equal sim to both, still below a high threshold? see test


class FakeClaude:
    """Minimal stand-in for anthropic.Anthropic().messages.create."""

    def __init__(self, text):
        self._text = text
        self.received = None

        class Messages:
            def __init__(self, outer):
                self._outer = outer

            def create(self, **kwargs):
                self._outer.received = kwargs

                class Block:
                    type = "text"
                    text = self._outer._text

                class Resp:
                    content = [Block()]

                return Resp()

        self.messages = Messages(self)


def test_grounded_answer_returns_citations():
    claude = FakeClaude("Use WasteStream [1].")
    out = engine.answer_question(
        "How do I make a WasteStream?",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=2,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    assert "WasteStream" in out["answer"]
    assert out["citations"]
    assert out["citations"][0]["url"] == RECORDS[0]["url"]
    assert out["citations"][0]["source"] == "qsdsan"


def test_low_similarity_refuses_without_calling_claude():
    claude = FakeClaude("SHOULD NOT BE CALLED")
    out = engine.answer_question(
        "How do I integrate SCADA/Modbus?",
        records=RECORDS,
        embed_fn=lambda t, input_type, client=None: [[-1.0, 0.0]],  # cosine -1
        claude_client=claude,
        top_k=2,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    assert claude.received is None  # Claude never invoked
    assert "couldn't find this" in out["answer"]
    assert out["citations"] == []


def test_system_prompt_is_cached():
    claude = FakeClaude("Answer [1].")
    engine.answer_question(
        "How do I make a WasteStream?",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=1,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    system = claude.received["system"]
    # System prompt sent as a cacheable block.
    assert isinstance(system, list)
    assert system[0]["cache_control"] == {"type": "ephemeral"}


def test_code_block_gets_disclaimer_appended():
    claude = FakeClaude("Here:\n\n```python\nws = WasteStream('ws')\n```\n")
    out = engine.answer_question(
        "Make a WasteStream",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=1,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    assert "Draft from the docs" in out["answer"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_engine.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'engine'`).

- [ ] **Step 3: Write minimal implementation**

`docs/chatbot/engine.py`:
```python
"""Orchestrate one question: embed -> retrieve -> threshold/refuse -> Claude -> cite.

Pure and dependency-injected: embed_fn and claude_client are passed in so the
whole contract is testable without network.
"""
import prompts
import retrieval


def _used_citation_numbers(answer: str, n: int) -> list[int]:
    """Return excerpt numbers actually referenced as [i] in the answer, in order."""
    import re

    seen = []
    for m in re.finditer(r"\[(\d+)\]", answer):
        i = int(m.group(1))
        if 1 <= i <= n and i not in seen:
            seen.append(i)
    return seen


def answer_question(
    question,
    records,
    embed_fn,
    claude_client,
    *,
    top_k,
    threshold,
    gen_model,
):
    """Return {"answer": str, "citations": list[dict]}."""
    query_vec = embed_fn([question], input_type="query")[0]
    hits = retrieval.cosine_top_k(query_vec, records, k=top_k)

    if not hits or hits[0][1] < threshold:
        return {"answer": prompts.refusal_message(), "citations": []}

    retrieved = [record for record, _score in hits]
    user_prompt = prompts.build_user_prompt(question, retrieved)

    resp = claude_client.messages.create(
        model=gen_model,
        max_tokens=1024,
        system=[
            {
                "type": "text",
                "text": prompts.SYSTEM_PROMPT,
                "cache_control": {"type": "ephemeral"},
            }
        ],
        messages=[{"role": "user", "content": user_prompt}],
    )
    answer = "".join(
        block.text for block in resp.content if getattr(block, "type", None) == "text"
    )
    answer = prompts.append_code_disclaimers(answer)

    used = _used_citation_numbers(answer, len(retrieved))
    if not used:  # answer cited nothing explicitly: surface all retrieved sources
        used = list(range(1, len(retrieved) + 1))
    citations = [
        {
            "n": i,
            "url": retrieved[i - 1]["url"],
            "title": retrieved[i - 1]["title"],
            "source": retrieved[i - 1]["source"],
        }
        for i in used
    ]
    return {"answer": answer, "citations": citations}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_engine.py -v`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add docs/chatbot/engine.py docs/chatbot/tests/test_engine.py
git commit -m "chatbot: add answer engine with grounding, citation, and refusal"
```

---

### Task 11: FastAPI server (thin wiring + contract test)

**Files:**
- Create: `docs/chatbot/server.py`
- Test: `docs/chatbot/tests/test_server.py`

- [ ] **Step 1: Write the failing test**

`docs/chatbot/tests/test_server.py`:
```python
from fastapi.testclient import TestClient

import server


def test_ask_returns_answer_and_citations(monkeypatch):
    def fake_answer(question, **kwargs):
        return {
            "answer": "Use WasteStream [1].",
            "citations": [
                {"n": 1, "url": "https://example/x", "title": "Create a WasteStream",
                 "source": "qsdsan"}
            ],
        }

    monkeypatch.setattr(server, "_answer", fake_answer)
    client = TestClient(server.app)
    resp = client.post("/ask", json={"question": "How do I make a WasteStream?"})
    assert resp.status_code == 200
    body = resp.json()
    assert body["answer"] == "Use WasteStream [1]."
    assert body["citations"][0]["source"] == "qsdsan"


def test_ask_rejects_empty_question(monkeypatch):
    monkeypatch.setattr(server, "_answer", lambda question, **kw: {"answer": "", "citations": []})
    client = TestClient(server.app)
    resp = client.post("/ask", json={"question": "   "})
    assert resp.status_code == 422


def test_healthcheck_ok():
    client = TestClient(server.app)
    assert client.get("/health").status_code == 200
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_server.py -v`
Expected: FAIL (`ModuleNotFoundError: No module named 'server'`).

- [ ] **Step 3: Write minimal implementation**

`docs/chatbot/server.py`:
```python
"""FastAPI endpoint. Loads the index once at startup and wires real clients into
engine.answer_question. Kept thin: all logic is tested in engine/retrieval/prompts.
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, field_validator

import config
import embeddings
import engine
import retrieval

app = FastAPI(title="QSDsan/EXPOsan docs chatbot")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # internal-first; tighten before any public launch
    allow_methods=["POST", "GET"],
    allow_headers=["*"],
)

_STATE = {"records": None, "claude": None}


class AskRequest(BaseModel):
    question: str

    @field_validator("question")
    @classmethod
    def _not_blank(cls, v):
        if not v or not v.strip():
            raise ValueError("question must not be blank")
        return v.strip()


class Citation(BaseModel):
    n: int
    url: str
    title: str
    source: str


class AskResponse(BaseModel):
    answer: str
    citations: list[Citation]


def _records():
    if _STATE["records"] is None:
        _STATE["records"] = retrieval.load_index(config.INDEX_URL)
    return _STATE["records"]


def _claude():
    if _STATE["claude"] is None:
        import anthropic

        _STATE["claude"] = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY
    return _STATE["claude"]


def _answer(question: str, **kwargs):
    """Indirection seam so tests can monkeypatch the whole pipeline."""
    return engine.answer_question(
        question,
        records=_records(),
        embed_fn=embeddings.embed_texts,
        claude_client=_claude(),
        top_k=config.TOP_K,
        threshold=config.SIMILARITY_THRESHOLD,
        gen_model=config.GEN_MODEL,
    )


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/ask", response_model=AskResponse)
def ask(req: AskRequest):
    return _answer(req.question)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_server.py -v`
Expected: PASS (3 passed).

- [ ] **Step 5: Run the full suite**

Run: `pytest docs/chatbot/tests -v`
Expected: PASS (all green).

- [ ] **Step 6: Commit**

```bash
git add docs/chatbot/server.py docs/chatbot/tests/test_server.py
git commit -m "chatbot: add FastAPI endpoint"
```

---

### Task 12: Widget (JS + CSS) and Sphinx registration

**Files:**
- Create: `docs/source/_static/js/chatbot.js`
- Create: `docs/source/_static/css/chatbot.css`
- Modify: `docs/source/conf.py:99-106`

> This task has no unit test (browser UI). Verification is a local docs build + manual smoke check, described in Step 4.

- [ ] **Step 1: Create the widget script**

`docs/source/_static/js/chatbot.js`:
```javascript
(function () {
  // Point this at the deployed Render endpoint. Overridable via a global set in conf.
  const ENDPOINT = window.QSDSAN_CHATBOT_ENDPOINT || "https://qsdsan-chatbot.onrender.com/ask";

  function el(tag, cls, html) {
    const node = document.createElement(tag);
    if (cls) node.className = cls;
    if (html !== undefined) node.innerHTML = html;
    return node;
  }

  // Minimal, safe-ish markdown: code fences, inline code, links, paragraphs.
  function renderMarkdown(md) {
    const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
    let html = esc(md);
    html = html.replace(/```(\w*)\n([\s\S]*?)```/g, (_m, _lang, code) => `<pre><code>${code}</code></pre>`);
    html = html.replace(/`([^`]+)`/g, "<code>$1</code>");
    html = html.replace(/\[([^\]]+)\]\((https?:[^)]+)\)/g, '<a href="$2" target="_blank" rel="noopener">$1</a>');
    html = html.replace(/\*(.+?)\*/g, "<em>$1</em>");
    html = html.split(/\n{2,}/).map((p) => (p.startsWith("<pre>") ? p : `<p>${p.replace(/\n/g, "<br>")}</p>`)).join("");
    return html;
  }

  function citationsHtml(citations) {
    if (!citations || !citations.length) return "";
    const items = citations
      .map((c) => `<li>[${c.n}] <a href="${c.url}" target="_blank" rel="noopener">${c.title}</a> <span class="qsd-src">(${c.source})</span></li>`)
      .join("");
    return `<ul class="qsd-citations">${items}</ul>`;
  }

  function init() {
    const button = el("button", "qsd-chat-button", "Ask the docs");
    button.setAttribute("aria-label", "Open the QSDsan docs assistant");
    const panel = el("div", "qsd-chat-panel qsd-hidden");
    panel.innerHTML = `
      <div class="qsd-chat-header">QSDsan / EXPOsan docs assistant
        <button class="qsd-close" aria-label="Close">&times;</button></div>
      <div class="qsd-chat-log"></div>
      <form class="qsd-chat-form">
        <input class="qsd-chat-input" type="text" placeholder="Ask about QSDsan or EXPOsan..." autocomplete="off"/>
        <button type="submit">Send</button>
      </form>`;
    document.body.appendChild(button);
    document.body.appendChild(panel);

    const log = panel.querySelector(".qsd-chat-log");
    const form = panel.querySelector(".qsd-chat-form");
    const input = panel.querySelector(".qsd-chat-input");

    button.addEventListener("click", () => panel.classList.toggle("qsd-hidden"));
    panel.querySelector(".qsd-close").addEventListener("click", () => panel.classList.add("qsd-hidden"));

    function addMessage(cls, html) {
      const m = el("div", "qsd-msg " + cls, html);
      log.appendChild(m);
      log.scrollTop = log.scrollHeight;
      return m;
    }

    form.addEventListener("submit", async (e) => {
      e.preventDefault();
      const question = input.value.trim();
      if (!question) return;
      addMessage("qsd-user", renderMarkdown(question));
      input.value = "";
      const pending = addMessage("qsd-bot", "<em>Searching the docs...</em>");
      try {
        const resp = await fetch(ENDPOINT, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ question }),
        });
        const data = await resp.json();
        pending.innerHTML = renderMarkdown(data.answer || "") + citationsHtml(data.citations);
      } catch (err) {
        pending.innerHTML = "<em>Sorry, the docs assistant is unavailable right now.</em>";
      }
      log.scrollTop = log.scrollHeight;
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
```

- [ ] **Step 2: Create the widget styles (Furo light/dark aware)**

`docs/source/_static/css/chatbot.css`:
```css
/* Uses Furo CSS variables so it follows the active light/dark theme. */
.qsd-chat-button {
  position: fixed;
  right: 1.25rem;
  bottom: 1.25rem;
  z-index: 1000;
  padding: 0.6rem 1rem;
  border: none;
  border-radius: 999px;
  cursor: pointer;
  background: var(--color-brand-primary, #2b8a3e);
  color: #fff;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.25);
}
.qsd-chat-panel {
  position: fixed;
  right: 1.25rem;
  bottom: 4.5rem;
  z-index: 1000;
  width: min(380px, 90vw);
  height: min(520px, 70vh);
  display: flex;
  flex-direction: column;
  border-radius: 12px;
  overflow: hidden;
  background: var(--color-background-primary, #fff);
  color: var(--color-foreground-primary, #1a1a1a);
  border: 1px solid var(--color-background-border, #ccc);
  box-shadow: 0 6px 24px rgba(0, 0, 0, 0.3);
}
.qsd-hidden { display: none; }
.qsd-chat-header {
  padding: 0.6rem 0.8rem;
  font-weight: 600;
  background: var(--color-background-secondary, #f0f0f0);
  display: flex;
  justify-content: space-between;
  align-items: center;
}
.qsd-close { background: none; border: none; font-size: 1.2rem; cursor: pointer; color: inherit; }
.qsd-chat-log { flex: 1; overflow-y: auto; padding: 0.8rem; }
.qsd-msg { margin-bottom: 0.75rem; line-height: 1.4; }
.qsd-user { text-align: right; }
.qsd-user p { display: inline-block; background: var(--color-background-secondary, #eef); padding: 0.4rem 0.6rem; border-radius: 10px; }
.qsd-bot pre {
  background: var(--color-code-background, #f5f5f5);
  padding: 0.6rem;
  border-radius: 8px;
  overflow-x: auto;
}
.qsd-citations { font-size: 0.85em; margin-top: 0.4rem; padding-left: 1.1rem; }
.qsd-src { opacity: 0.7; }
.qsd-chat-form { display: flex; gap: 0.4rem; padding: 0.6rem; border-top: 1px solid var(--color-background-border, #ccc); }
.qsd-chat-input { flex: 1; padding: 0.45rem; border-radius: 8px; border: 1px solid var(--color-background-border, #ccc); background: var(--color-background-primary, #fff); color: inherit; }
.qsd-chat-form button { padding: 0.45rem 0.8rem; border: none; border-radius: 8px; cursor: pointer; background: var(--color-brand-primary, #2b8a3e); color: #fff; }
```

- [ ] **Step 3: Register the assets in conf.py**

Modify `docs/source/conf.py` lines 99-106 to:
```python
html_css_files = [
	'css/qsdsan.css',
 	'css/copybutton.css',
 	'css/chatbot.css',
# 	'css/theme_overrides.css',
 	]
html_js_files = [
    'js/unit-operation-filters.js',
    'js/chatbot.js',
]
```

- [ ] **Step 4: Build the docs locally and smoke-check the widget**

Run (from `docs/`):
```bash
make html
```
Expected: build succeeds. Open `docs/build/html/index.html` in a browser; confirm an "Ask the docs" button renders bottom-right and the panel opens/closes. (The endpoint will not answer until Render is deployed; that is expected here.)

- [ ] **Step 5: Commit**

```bash
git add docs/source/_static/js/chatbot.js docs/source/_static/css/chatbot.css docs/source/conf.py
git commit -m "chatbot: add Furo widget and register assets"
```

---

### Task 13: Deploy config (Render) + reindex hook (readthedocs) + README

**Files:**
- Create: `render.yaml`
- Modify: `.readthedocs.yml`
- Create: `docs/chatbot/README.md`

> No unit test. Verification is YAML validity + a dry-run of the indexer command.

- [ ] **Step 1: Create the Render Blueprint**

`render.yaml` (repo root; Render requires it there):
```yaml
services:
  - type: web
    name: qsdsan-chatbot
    runtime: python
    rootDir: docs/chatbot
    plan: free
    buildCommand: pip install -r requirements.txt
    startCommand: uvicorn server:app --host 0.0.0.0 --port $PORT
    envVars:
      - key: ANTHROPIC_API_KEY
        sync: false   # set in the Render dashboard; never committed
      - key: VOYAGE_API_KEY
        sync: false
      - key: CHATBOT_INDEX_URL
        value: https://qsdsan.readthedocs.io/en/latest/_static/chatbot/index.json
```

- [ ] **Step 2: Add the post_build reindex hook to .readthedocs.yml**

Modify `.readthedocs.yml` — add a `jobs.post_build` under `build:` (the indexer needs `voyageai`, `requests`, `numpy`, `beautifulsoup4`, which are not in the docs extra, so install them inline). `VOYAGE_API_KEY` is set in the readthedocs dashboard env, not in this file:
```yaml
build:
  os: ubuntu-20.04
  tools:
    python: "3.12"
  jobs:
    post_build:
      - pip install voyageai requests numpy beautifulsoup4
      - python docs/chatbot/index_docs.py
```
Keep the existing `sphinx:` and `python:` sections unchanged. The indexer's `main()` reads `$READTHEDOCS_OUTPUT/html` for the built docs and writes `docs/source/_static/chatbot/index.json`, which readthedocs serves at the `CHATBOT_INDEX_URL` above.

- [ ] **Step 3: Write the README**

`docs/chatbot/README.md`:
```markdown
# QSDsan + EXPOsan docs chatbot

A self-hosted RAG assistant embedded on the QSDsan readthedocs site. It answers
questions about QSDsan (API + tutorials) and EXPOsan (example systems) strictly
from indexed docs, with citations and out-of-scope refusal. See `DESIGN.md` for
the full spec.

## Components

- `index_docs.py` - builds `index.json` from QSDsan built HTML + EXPOsan READMEs.
- `server.py` - FastAPI `/ask` endpoint (deployed on Render).
- `../source/_static/js/chatbot.js` + `css/chatbot.css` - the widget.
- `.readthedocs.yml` post_build - reindexes on every docs build.

## Local development

```bash
pip install -r docs/chatbot/requirements.txt
pytest docs/chatbot/tests -v
```

## Build the index locally

```bash
# After `cd docs && make html` so docs/build/html exists:
export VOYAGE_API_KEY=...   # needed to embed chunks
python docs/chatbot/index_docs.py
```

## Run the endpoint locally

```bash
export ANTHROPIC_API_KEY=...
export VOYAGE_API_KEY=...
export CHATBOT_INDEX_URL=docs/source/_static/chatbot/index.json   # local file is fine
cd docs/chatbot && uvicorn server:app --reload
```

## Configuration

All knobs live in `config.py` and read from env vars: `CHATBOT_GEN_MODEL`
(default Haiku 4.5; set to `claude-sonnet-4-6` to switch), `CHATBOT_TOP_K`,
`CHATBOT_SIMILARITY_THRESHOLD`, `CHATBOT_INDEX_URL`.

## Secrets

`ANTHROPIC_API_KEY` and `VOYAGE_API_KEY` live ONLY in the Render dashboard and the
readthedocs build environment. Never commit them.
```

- [ ] **Step 4: Validate the YAML and dry-run the indexer help**

Run:
```bash
python -c "import yaml; yaml.safe_load(open('render.yaml')); yaml.safe_load(open('.readthedocs.yml')); print('yaml ok')"
python -c "import index_docs; print('systems:', len(index_docs.EXPOSAN_SYSTEMS))" 
```
Expected: `yaml ok` and a system count printed. (Run the second command from `docs/chatbot/` or with that dir on `PYTHONPATH`.)

- [ ] **Step 5: Commit**

```bash
git add render.yaml .readthedocs.yml docs/chatbot/README.md
git commit -m "chatbot: add Render deploy config, readthedocs reindex hook, and README"
```

---

### Task 14: Eval harness (regression guard)

**Files:**
- Create: `docs/chatbot/eval/questions.yaml`
- Create: `docs/chatbot/eval/run_eval.py`
- Test: `docs/chatbot/tests/test_eval_questions.py`

- [ ] **Step 1: Write the seeded eval questions**

`docs/chatbot/eval/questions.yaml`:
```yaml
# Regression guard for the docs chatbot. Each item is run against a live /ask
# endpoint by run_eval.py. `expect` is one of:
#   cite_substring: the answer's citation URLs must include this substring
#   refuse: the answer must be the out-of-scope refusal (no citations)
- question: "create a WasteStream with flow and COD"
  expect: cite_substring
  value: "tutorials/3_WasteStream"
- question: "dynamic simulation pitfalls"
  expect: cite_substring
  value: "tutorials/14_Modeling_Notes_and_Pitfalls"
- question: "how do I load and run the BSM1 system?"
  expect: cite_substring
  value: "exposan/bsm1/README.rst"
- question: "how do I set up SCADA/Modbus integration?"
  expect: refuse
```

- [ ] **Step 2: Write the failing test (validates the eval file shape)**

`docs/chatbot/tests/test_eval_questions.py`:
```python
from pathlib import Path

import yaml

EVAL = Path(__file__).resolve().parents[1] / "eval" / "questions.yaml"


def test_questions_file_is_well_formed():
    items = yaml.safe_load(EVAL.read_text())
    assert len(items) >= 4
    for item in items:
        assert item["question"]
        assert item["expect"] in {"cite_substring", "refuse"}
        if item["expect"] == "cite_substring":
            assert item["value"]


def test_seeded_cases_present():
    items = yaml.safe_load(EVAL.read_text())
    questions = " ".join(i["question"].lower() for i in items)
    assert "wastestream" in questions
    assert "pitfalls" in questions
    assert "bsm1" in questions
    assert "scada" in questions or "modbus" in questions
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest docs/chatbot/tests/test_eval_questions.py -v`
Expected: FAIL (`FileNotFoundError` if run before Step 1, otherwise import error for `run_eval` is not yet relevant). If Step 1 is done, this should pass — in that case write the test first by temporarily renaming the eval file, or accept that the data file precedes its validator. To honor TDD, create the test first, watch it fail with `FileNotFoundError`, then add `questions.yaml`.

> TDD note: do Step 2 before Step 1. Watch `FileNotFoundError`, then create `questions.yaml` and re-run.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest docs/chatbot/tests/test_eval_questions.py -v`
Expected: PASS (2 passed).

- [ ] **Step 5: Write the eval runner**

`docs/chatbot/eval/run_eval.py`:
```python
"""Run the seeded eval questions against a live /ask endpoint.

Usage:
    python docs/chatbot/eval/run_eval.py --endpoint http://localhost:8000/ask

Exits non-zero if any case fails, so it can gate a release manually. Requires a
running server with real keys; not part of the unit test suite.
"""
import argparse
import sys
from pathlib import Path

import requests
import yaml

QUESTIONS = Path(__file__).parent / "questions.yaml"


def run(endpoint: str) -> int:
    items = yaml.safe_load(QUESTIONS.read_text())
    failures = 0
    for item in items:
        resp = requests.post(endpoint, json={"question": item["question"]}, timeout=120)
        data = resp.json()
        answer = data.get("answer", "")
        citations = data.get("citations", [])
        if item["expect"] == "refuse":
            ok = "couldn't find this" in answer and not citations
        else:
            urls = " ".join(c["url"] for c in citations)
            ok = item["value"] in urls
        status = "PASS" if ok else "FAIL"
        if not ok:
            failures += 1
        print(f"[{status}] {item['question']}")
    print(f"\n{len(items) - failures}/{len(items)} passed")
    return 1 if failures else 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint", default="http://localhost:8000/ask")
    args = parser.parse_args()
    sys.exit(run(args.endpoint))
```

- [ ] **Step 6: Final full-suite run**

Run: `pytest docs/chatbot/tests -v`
Expected: PASS (all green across every test file).

- [ ] **Step 7: Commit**

```bash
git add docs/chatbot/eval/questions.yaml docs/chatbot/eval/run_eval.py docs/chatbot/tests/test_eval_questions.py
git commit -m "chatbot: add eval harness as a regression guard"
```

---

## Self-review against the spec

**Spec coverage:**
- Indexer, QSDsan adapter (built HTML, heading split, RTD anchor, source="qsdsan") → Tasks 3, 5, 7. ✓
- Indexer, EXPOsan adapter (GitHub raw fetch, README-only, blob URL, source="exposan") → Tasks 2, 4, 7. ✓
- Same chunk schema `{text,title,url,source,type}` → enforced in Tasks 4/5; embeddings added in Task 7 (documented as the index artifact extension). ✓
- `index.json` gitignored → already at `.gitignore:29`; no task needed (noted in file table). ✓
- Query endpoint: embed → cosine top-k (flat numpy) → cached system prompt + chunks → `{answer, citations[]}` → Tasks 8, 9, 10, 11. ✓
- Loads index from published URL on startup → Task 11 `_records()` + `config.INDEX_URL`. ✓ (Refresh-when-changed is simplified to load-on-startup; a periodic refresh is out of scope for internal-first and can be added later — flagged below.)
- Deploy to Render via render.yaml → Task 13. ✓
- Widget: floating button, chat panel, markdown + clickable citations, Furo light/dark, registered like unit-operation-filters.js → Task 12. ✓
- Reindex hook: one post_build line runs the indexer → Task 13. ✓
- Guardrails: answer only from excerpts + cite; below-threshold refusal without calling Claude + docs-search link; code-block disclaimer; refuse unrelated → Tasks 9, 10 (`test_low_similarity_refuses_without_calling_claude`, disclaimer test). ✓
- Decisions: voyage-4-lite, Haiku 4.5 switchable to Sonnet, prompt caching, Render free tier, secrets only in env → config.py defaults + cache_control in engine (Task 10) + render.yaml `sync: false` (Task 13). ✓
- Testing: chunking tests both adapters with resolvable-url + source assertions (Tasks 2-5); EXPOsan fetch+parse with mocked fixture (Task 4); endpoint contract tests with mocked Claude (Task 10/11); eval harness seeded with the four cases (Task 14). ✓

**Known simplification (flag for review at the checkpoint):** the spec says the server "refreshes when [index.json] changes." This plan loads once at startup (sufficient for internal-first, where Render restarts on deploy/cold-start re-pull the latest index). If live refresh-without-restart is required, add a small TTL/etag check in `_records()`. Raising this rather than silently dropping it.

**Placeholder scan:** no TBD/TODO/"handle edge cases"/"similar to Task N" — every code step contains full code.

**Type consistency:** chunk record keys (`text,title,url,source,type,embedding`), heading-split keys (`title,text,anchor`), citation keys (`n,url,title,source`), and function signatures (`build_exposan_chunks(systems, fetch)`, `build_qsdsan_chunks(html_dir, base_url=None)`, `build_index(html_dir, systems, fetch, embed_fn)`, `embed_texts(texts, input_type, client)`, `cosine_top_k(query_vec, records, k)`, `answer_question(question, records, embed_fn, claude_client, *, top_k, threshold, gen_model)`) are consistent across all tasks. ✓

---

## Review checkpoints during execution

Pause for the user to review after:
- **Task 5** (both adapters done — confirm chunk schema + URL shape before embeddings).
- **Task 11** (endpoint contract complete — confirm guardrail behavior before UI/deploy).
- **Task 14** (feature complete — confirm before any deploy/push).
