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
