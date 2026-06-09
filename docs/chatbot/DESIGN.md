# QSDsan + EXPOsan Documentation Chatbot — Design Spec

An embedded RAG widget on the QSDsan readthedocs site that answers questions
about QSDsan (the package: API + tutorials) and EXPOsan (the runnable example
systems). Self-hosted lightweight RAG: Claude API for generation, a Furo-injected
widget for the UI. This is the single source of truth for the chatbot's design;
it supersedes the earlier separate design and "next steps" notes.

## Status and scope

- **Internal-first.** Released to the team and close contributors to validate
  value before any public launch.
- A NotebookLM pilot validated that QSDsan's docs ground well and that strict
  grounding makes out-of-scope refusal reliable. The risk to guard against is
  confidently-wrong API in generated code, so the index favors structured,
  high-signal sources (real signatures, runnable doctests, known-correct tests)
  over raw prose.
- The bot covers QSDsan (API + tutorials) and the EXPOsan example systems
  (bsm1, bsm2, adm, metab, htl, the sanitation systems, etc.), the latter via
  code-extracted system wiring (see the indexer below).
- Code lives in the QSD-Group/QSDsan repo, on the `chatbot` worktree
  branched off `beta`. All chatbot code lives under `docs/chatbot/` and
  `docs/source/_static/`. Keep ALL commits scoped to those paths; do NOT touch
  the parallel `beta` work (`qsdsan/_tea.py`, `tests/test_tea.py`).

## Architecture (Path B — self-hosted lightweight RAG)

Four loosely-coupled, independently-testable units, all under the QSDsan repo.

### 1. Indexer — `docs/chatbot/index_docs.py` (+ `code_adapter.py`)

Source-PLUGGABLE. Several adapters write the SAME chunk schema into
`docs/source/_static/chatbot/index.json`:

```
{text, title, url, source, type}
```

- **QSDsan HTML adapter** (`index_docs.build_qsdsan_chunks`) — runs at docs-build
  time over the BUILT HTML (tutorials + API + FAQ + homepage). Splits by heading
  (`chunking.split_html_by_heading`). Citation url = readthedocs page + section
  anchor. `source = "qsdsan"`; `type` = `tutorial` / `api` / `faq` / `guide`.
- **Code adapter** (`code_adapter.build_code_chunks`) — AST-walks the installed
  `qsdsan` package source (located via `find_spec`, never imported) and emits
  structured, high-signal chunks. `source = "code"`, citation url = GitHub blob
  URL + line anchor (`.../QSDsan/blob/main/<path>#L<start>-L<end>`):
  - **Doctests** (`type = "doctest"`) — one chunk per docstring that contains a
    `>>>` block (parsed with stdlib `doctest`), the full runnable interaction,
    titled with the qualified name (e.g. `qsdsan.WasteStream.from_composite`).
  - **API surface** (`type = "api"`) — one chunk per public class/function:
    signature with arg defaults, decorators, base classes, return annotation,
    and the one-line summary. Restores the accurate signatures the autodoc HTML
    buries. Undocumented functions/methods are dropped as low-signal noise;
    classes are always kept.
- **Example adapter** (`code_adapter.build_example_chunks`) — `source =
  "example"`:
  - **Tests** (`type = "test"`) — each `test_*` function in the repo `tests/`
    tree as a self-contained, known-correct usage snippet.
  - **EXPOsan systems** (`type = "system"`) — each module-level `create_system`
    function (the system-assembly wiring), when EXPOsan is installed. EXPOsan is
    optional: if it is not importable the adapter yields nothing rather than
    failing. Citation url = `.../EXPOsan/blob/main/exposan/<sys>/...`.

Robustness: the walkers skip files that fail to parse/decode (`SyntaxError` /
`UnicodeDecodeError`) instead of failing the build, suppress `SyntaxWarning` at
the parse point (many contributor sources have non-raw strings), and cap snippet
chunks (`_MAX_SNIPPET_LINES`) so one oversized `create_system` cannot dominate
the retrieved context — the citation still spans the whole function.

`index.json` is gitignored. Each chunk's `source` tag lets answers label
provenance and lets the prompt prefer prose for explanation and code for exact
usage.

> **History:** An EXPOsan README adapter (fetching per-system `README.rst` from
> GitHub raw) was tried and removed 2026-05-31 — the README chunks produced
> low-quality answers. EXPOsan is now covered by structured extraction of
> `create_system` wiring (above), NOT by re-dumping README or raw `.py` text.

### 2. Query endpoint — `docs/chatbot/server.py` (FastAPI)

`POST /ask {question}`:

1. Embed the question.
2. Cosine top-k over `index.json` (flat numpy search; corpus is small, no
   vector DB).
3. Call Claude with a cached system prompt + retrieved chunks.
4. Return `{answer, citations[]}`.

Loads `index.json` from the published docs URL on startup. Deploys to Render via
`render.yaml`. The answer pipeline lives in `engine.answer_question` (pure,
dependency-injected: `embed_fn` and `claude_client` are passed in).

### 3. Widget — `docs/source/_static/js/chatbot.js` + `css/chatbot.css`

Floating button opening a chat panel. Renders answer markdown with clickable
citation links. Furo light/dark aware. Registered via `html_js_files` /
`html_css_files` in `docs/source/conf.py` (same mechanism as
`unit-operation-filters.js`).

### 4. Reindex hook — `.readthedocs.yml`

One `build.jobs.post_build` line runs the indexer (all adapters), so the index
always matches published docs + current source. Reindex on every build/release;
this is what fixes the stale-page problem.

## Guardrails (balanced + disclaimers)

Behavior lives in `docs/chatbot/prompts.py` and `engine.py`.

- Answer ONLY from retrieved excerpts; cite each claim with its number + URL.
- Prefer `(qsdsan)` prose excerpts to EXPLAIN concepts and `(code)` excerpts
  (doctests/signatures) for EXACT, runnable usage and API.
- If top retrieval similarity is below a threshold, DON'T call Claude — return
  "I couldn't find this in the QSDsan/EXPOsan docs" + a docs-search link
  (out-of-scope refusal).
- Greetings / capability questions get a friendly welcome (`is_smalltalk`).
- **EXPOsan routing** (`is_exposan_question`): catalog/discovery questions
  ("what systems exist", "what have you built") are pointed to the Systems page
  + the EXPOsan repo. A specific system asked about in a build/run/code way
  ("how do I run BSM1?") instead falls through to retrieval, so it is answered
  from the indexed `create_system` wiring.
- Auto-append after every code block: "Draft from the docs — verify against the
  linked pages."
- Refuse questions unrelated to QSDsan/EXPOsan.

## Decisions

- **Structured extraction, not raw text.** Code/example chunks are extracted via
  AST + `doctest` into high-signal units; raw `.py` and README dumps were
  rejected as too noisy.
- **Embeddings:** Voyage `voyage-4-lite` (first 200M tokens free, effectively $0
  at this corpus size). Needs `VOYAGE_API_KEY` in the readthedocs build env AND
  on the Render endpoint. (Optional future: `voyage-code-3` for code chunks —
  better code retrieval — gated behind `config.EMBED_MODEL`; weigh the
  complexity of mixing models within one index first.)
- **Generation model:** start on Claude Haiku 4.5 (cheap/fast) with prompt
  caching; config-switchable to Sonnet.
- **Endpoint host:** Render (free tier fine; ~30-60s cold start after idle is
  acceptable for internal use).
- **Secrets** (`ANTHROPIC_API_KEY`, `VOYAGE_API_KEY`) live only in Render +
  readthedocs env. NEVER commit them.

## Testing (TDD — write tests first)

- **Indexer:** chunking unit tests for every adapter; assert each chunk has a
  resolvable url (readthedocs anchor for `qsdsan`; GitHub blob URL for
  `code`/`example`) and a correct `source`/`type` tag.
- **Code/example adapters:** test against small recorded module strings and a
  fixture package tree — no installed package and no live network. Cover
  doctest grouping, signature rendering, the undocumented-noise filter,
  module-level `create_system` only, snippet truncation, and graceful skips on
  unparseable/odd-encoding sources and on a missing EXPOsan.
- **Endpoint/engine:** contract tests with a MOCKED Claude for grounding /
  citation / refusal / smalltalk / EXPOsan-routing behavior.
- **Eval harness** at `docs/chatbot/eval/questions.yaml` as a regression guard,
  seeded with at least:
  - "create a WasteStream with flow+COD" → cites the QSDsan streams page
  - "dynamic simulation pitfalls" → cites QSDsan tutorial 14
  - "how do I load and run the BSM1 system?" → cites the EXPOsan `system.py`
    (`create_system` wiring)
  - "create a Component from an existing chemical" → cites a QSDsan source
    doctest (only the code adapter can answer this)
  - "SCADA/Modbus integration" → refuses (out of scope)

## Project conventions

- Use TDD (failing test first, then implement); recorded fixtures, no live
  network in tests.
- Update docstrings/READMEs alongside code; stale docs are bugs.
- Keep ALL commits scoped to chatbot files only. Do NOT touch or commit the
  parallel `beta` work (`qsdsan/_tea.py`, `tests/test_tea.py`).
- Show the diff and ask before committing. NEVER push without explicit approval.
- Minimize em-dashes in user-facing prose.

## Remaining / optional

- `voyage-code-3` embeddings for code chunks (see Decisions) — the last
  not-yet-done item; everything else above is implemented.
