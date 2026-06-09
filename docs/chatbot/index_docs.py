"""Indexer: build chunk records from the QSDsan built HTML.

Walks the built Sphinx HTML (tutorials + API), splits by heading into chunks,
attaches embeddings, and writes index.json. EXPOsan is intentionally not indexed;
EXPOsan questions are routed to a pointer response by the query engine instead.
"""
from __future__ import annotations

import json
import os

import config
import chunking
import code_adapter
import embeddings


# Pages indexed from the built docs: the tutorials / API / FAQ subtrees plus the
# homepage (index.html), which carries the install and overview content. Read the
# Docs build scaffolding is skipped.
_INDEXED_DIRS = ("tutorials", "api", "faq")
_INDEXED_ROOT_PAGES = ("index.html",)
_SKIP_PAGES = {"genindex.html", "search.html", "py-modindex.html", "404.html"}


def _page_type(rel_path: str) -> str:
    """Tag a built page by its location under the docs tree."""
    head = rel_path.replace(os.sep, "/").split("/", 1)[0]
    if head == "tutorials":
        return "tutorial"
    if head == "api":
        return "api"
    if head == "faq":
        return "faq"
    return "guide"


def build_qsdsan_chunks(html_dir: str, base_url: str | None = None) -> list[dict]:
    """Walk a built Sphinx HTML tree and chunk each page by heading.

    Citation URL = readthedocs base + page path relative to html_dir + #anchor.
    Indexes the tutorials/, api/, and faq/ subtrees plus the homepage; skips Read
    the Docs build scaffolding (search, genindex, ...).
    """
    base_url = base_url or config.QSDSAN_DOCS_BASE
    chunks: list[dict] = []
    for root, _dirs, files in os.walk(html_dir):
        for fname in files:
            if not fname.endswith(".html") or fname in _SKIP_PAGES:
                continue
            abs_path = os.path.join(root, fname)
            rel_path = os.path.relpath(abs_path, html_dir).replace(os.sep, "/")
            parts = rel_path.split("/")
            is_root_page = len(parts) == 1
            if not (
                parts[0] in _INDEXED_DIRS
                or (is_root_page and fname in _INDEXED_ROOT_PAGES)
            ):
                continue
            ptype = _page_type(rel_path)
            with open(abs_path, encoding="utf-8") as fh:
                html = fh.read()
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


def build_all_chunks(
    html_dir, base_url=None, code_chunks_fn=None, example_chunks_fn=None
) -> list[dict]:
    """Merge the HTML adapter's chunks with the code and example adapters'.

    All write the same schema. ``code_chunks_fn`` (API + doctests from the
    package) and ``example_chunks_fn`` (test snippets + EXPOsan systems) are
    indirection points (like ``embed_fn``) so tests can supply deterministic
    chunks instead of walking the installed packages.
    """
    code_chunks_fn = code_chunks_fn or code_adapter.build_code_chunks
    example_chunks_fn = example_chunks_fn or code_adapter.build_example_chunks
    chunks = build_qsdsan_chunks(html_dir, base_url)
    chunks.extend(code_chunks_fn())
    chunks.extend(example_chunks_fn())
    return chunks


def embed_documents(texts, input_type, client=None):
    """Indirection point so tests can monkeypatch embedding without network."""
    return embeddings.embed_texts(texts, input_type=input_type, client=client)


def build_index(
    html_dir, embed_fn=None, code_chunks_fn=None, example_chunks_fn=None
) -> list[dict]:
    """Build QSDsan HTML + code + example chunk records and attach embeddings."""
    embed_fn = embed_fn or embed_documents
    records = build_all_chunks(
        html_dir, code_chunks_fn=code_chunks_fn, example_chunks_fn=example_chunks_fn
    )
    if records:
        vectors = embed_fn([r["text"] for r in records], input_type="document")
        if len(vectors) != len(records):
            raise ValueError(
                f"embed_fn returned {len(vectors)} vectors for {len(records)} records"
            )
        for r, vec in zip(records, vectors):
            r["embedding"] = list(vec)
    return records


def main(html_dir=None, out_path=None) -> None:
    """CLI entry: build the index and write it as JSON.

    On readthedocs the index is written into the build OUTPUT static dir
    ($READTHEDOCS_OUTPUT/html/_static/chatbot/index.json) so it is published;
    locally it defaults to the in-tree source _static path.
    """
    rtd_output = os.environ.get("READTHEDOCS_OUTPUT")
    html_dir = html_dir or os.path.join(rtd_output or "docs/build", "html")
    if out_path is None:
        if rtd_output:
            out_path = os.path.join(
                rtd_output, "html", "_static", "chatbot", "index.json"
            )
        else:
            out_path = "docs/source/_static/chatbot/index.json"
    records = build_index(html_dir, embed_fn=embed_documents)
    if not records:
        raise SystemExit(
            f"Indexer produced no chunks (html_dir={html_dir!r}); aborting so a "
            "broken build does not publish an empty index."
        )
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(records, fh)
    print(f"Wrote {len(records)} chunks to {out_path}")


if __name__ == "__main__":
    main()
