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
