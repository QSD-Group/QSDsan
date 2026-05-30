"""Indexer: build chunk records from QSDsan built HTML and EXPOsan READMEs.

Source-pluggable. Both adapters emit the same chunk schema; build_index adds
embeddings and writes index.json.
"""
from __future__ import annotations

import os

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
