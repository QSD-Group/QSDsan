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
