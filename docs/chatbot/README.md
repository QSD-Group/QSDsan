# QSDsan + EXPOsan docs chatbot

A self-hosted RAG assistant embedded on the QSDsan readthedocs site. It answers
questions about QSDsan (API + tutorials) strictly from indexed docs, with citations
and out-of-scope refusal. EXPOsan and "what systems" questions are routed to a
pointer (the Systems page + the EXPOsan repo) rather than answered from the index.
See `DESIGN.md` for the original spec and `../superpowers/plans/2026-05-30-docs-chatbot.md`
for the build plan.

## Components

- `index_docs.py` - builds `index.json` from the QSDsan built HTML (tutorials + API).
- `embeddings.py` - Voyage embedding wrapper (documents at index time, queries at ask time).
- `retrieval.py` - flat numpy cosine top-k over the index.
- `prompts.py` - system/user prompts, guardrail text, and the EXPOsan pointer.
- `engine.py` - the answer pipeline (EXPOsan pointer, embed, retrieve, refuse-below-threshold, Claude, cite).
- `server.py` - FastAPI `/ask` endpoint (deployed on Render).
- `../source/_static/js/chatbot.js` + `css/chatbot.css` - the Furo widget.
- `.readthedocs.yml` post_build - reindexes on every docs build.

## Local development

```bash
pip install -r docs/chatbot/requirements.txt
python -m pytest docs/chatbot/tests -v -o addopts=""
```

## Preview the widget locally

The widget renders on any local docs build (it calls the deployed Render endpoint
for answers, so the panel works even before you set up the index):

```bash
cd docs
make html
# open docs/build/html/index.html and click "Ask the docs" (bottom right)
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

## Deploy (Render)

`render.yaml` (repo root) defines the `qsdsan-chatbot` web service:
`https://qsdsan-chatbot.onrender.com` (endpoint `.../ask`). Create it via Render's
Blueprint flow, set `ANTHROPIC_API_KEY` and `VOYAGE_API_KEY` in the dashboard
(they are declared `sync: false` so they are never committed), and deploy from the
docs branch. Free tier sleeps after idle with a ~30-60s cold start, which is
acceptable for internal use.

## Configuration

All knobs live in `config.py` and read from env vars: `CHATBOT_GEN_MODEL`
(default Haiku 4.5; set to `claude-sonnet-4-6` to switch), `CHATBOT_TOP_K`,
`CHATBOT_SIMILARITY_THRESHOLD`, `CHATBOT_INDEX_URL`.

## Secrets

`ANTHROPIC_API_KEY` and `VOYAGE_API_KEY` live ONLY in the Render dashboard and the
readthedocs build environment. Never commit them.
