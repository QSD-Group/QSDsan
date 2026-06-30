# QSDsan + EXPOsan Documentation Chatbot — Design Spec

An embedded RAG widget on the QSDsan readthedocs site that answers questions
about both QSDsan (the package: API + tutorials) and EXPOsan (the runnable
example systems). Self-hosted lightweight RAG: Claude API for generation, a
Furo-injected widget for the UI.

> **Update (2026-05-31):** The EXPOsan README adapter was removed after internal
> testing showed the README chunks produced low-quality answers. EXPOsan is no
> longer indexed. Instead, the query engine detects EXPOsan / "what systems"
> questions and returns a pointer to the QSDsan Systems page and the EXPOsan
> repository (`prompts.is_exposan_question` / `exposan_pointer_message`). The
> indexer now covers QSDsan built HTML only. The rest of this spec is unchanged.

## Status and scope

- **Internal-first.** Released to the team and close contributors to validate
  value before any public launch.
- A NotebookLM pilot validated that QSDsan's docs ground well and that strict
  grounding makes out-of-scope refusal reliable. The risk to guard against is
  confidently-wrong API in generated code.
- The bot covers BOTH QSDsan (API + tutorials) and EXPOsan (the example
  systems: bsm1, bsm2, adm, metab, htl, the sanitation systems, etc.).
- Code lives in the QSD-Group/QSDsan repo, on the `docs-chatbot` worktree
  branched off `beta`. All chatbot code lives under `docs/chatbot/` and
  `docs/source/_static/`.

## Architecture (Path B — self-hosted lightweight RAG)

Four loosely-coupled, independently-testable units, all under the QSDsan repo.

### 1. Indexer — `docs/chatbot/index_docs.py`

Source-PLUGGABLE. Two adapters write the SAME chunk schema into
`docs/source/_static/chatbot/index.json`:

```
{text, title, url, source, type}
```

- **QSDsan adapter** — runs at docs-build time over the BUILT HTML (tutorials +
  API). Splits by heading into chunks. Citation url = readthedocs page +
  section anchor. `source = "qsdsan"`.
- **EXPOsan adapter** — FETCHES the 19 per-system `README.rst` (plus the
  top-level README) from GitHub raw. EXPOsan has no docs build, so this must NOT
  depend on EXPOsan being installed or checked out. Chunk by heading. Citation
  url = `github.com/QSD-Group/EXPOsan/blob/main/exposan/<sys>/README.rst`.
  `source = "exposan"`. Index READMEs ONLY, not the `.py` source (too noisy).

`index.json` is gitignored. Each chunk's `source` tag lets answers label
provenance.

### 2. Query endpoint — `docs/chatbot/server.py` (FastAPI)

`POST /ask {question}`:

1. Embed the question.
2. Cosine top-k over `index.json` (flat numpy search; corpus is small, no
   vector DB).
3. Call Claude with a cached system prompt + retrieved chunks.
4. Return `{answer, citations[]}`.

Loads `index.json` from the published docs URL on startup, refreshes when it
changes. Deploys to Render via `render.yaml`. The endpoint, widget, and
guardrails are identical regardless of source — the adapters only differ inside
the indexer.

### 3. Widget — `docs/source/_static/js/chatbot.js` + `css/chatbot.css`

Floating button opening a chat panel. Renders answer markdown with clickable
citation links. Furo light/dark aware. Registered via `html_js_files` /
`html_css_files` in `docs/source/conf.py` (same mechanism as
`unit-operation-filters.js`).

### 4. Reindex hook — `.readthedocs.yml`

One `build.jobs.post_build` line runs the indexer (both adapters), so the index
always matches published docs + current EXPOsan READMEs. Reindex on every
build/release; this is what fixes the stale-page problem.

## Guardrails (balanced + disclaimers)

Behavior lives in `docs/chatbot/prompts.py`.

- Answer ONLY from retrieved excerpts; cite each claim with its number + URL.
- If top retrieval similarity is below a threshold, DON'T call Claude — return
  "I couldn't find this in the QSDsan/EXPOsan docs" + a docs-search link
  (out-of-scope refusal).
- Auto-append to every code block: "Draft from the docs — verify against the
  linked pages."
- Refuse questions unrelated to QSDsan/EXPOsan.

## Decisions

- **Embeddings:** Voyage `voyage-4-lite` (first 200M tokens free, effectively $0
  at this corpus size). Needs `VOYAGE_API_KEY` in the readthedocs build env AND
  on the Render endpoint.
- **Generation model:** start on Claude Haiku 4.5 (cheap/fast) with prompt
  caching; config-switchable to Sonnet.
- **Endpoint host:** Render (free tier fine; ~30-60s cold start after idle is
  acceptable for internal use).
- **Secrets** (`ANTHROPIC_API_KEY`, `VOYAGE_API_KEY`) live only in Render +
  readthedocs env. NEVER commit them.

## Testing (TDD — write tests first)

- **Indexer:** chunking unit tests for BOTH adapters; assert every chunk has a
  resolvable url (readthedocs anchor for qsdsan, GitHub blob URL for exposan)
  and a correct `source` tag.
- **EXPOsan adapter:** test the GitHub-raw fetch + parse with a recorded/mocked
  README fixture (no live network in tests).
- **Endpoint:** contract tests with a MOCKED Claude for grounding/citation/
  refusal behavior.
- **Eval harness** at `docs/chatbot/eval/questions.yaml` as a regression guard,
  seeded with:
  - "create a WasteStream with flow+COD" → cites the QSDsan streams page
  - "dynamic simulation pitfalls" → cites QSDsan tutorial 14
  - "how do I load and run the BSM1 system?" → cites the EXPOsan bsm1 README
  - "SCADA/Modbus integration" → refuses (out of scope)

## Project conventions

- Use TDD (failing test first, then implement).
- Update docstrings/READMEs alongside code; stale docs are bugs.
- Keep ALL commits scoped to chatbot files only. Do NOT touch or commit the
  parallel `beta` work (`qsdsan/_tea.py`, `tests/test_tea.py`) — that's a
  separate effort in the other worktree.
- Show the diff and ask before committing. NEVER push without explicit approval.
- Minimize em-dashes in user-facing prose.
