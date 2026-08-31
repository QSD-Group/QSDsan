import json
from pathlib import Path

import pytest

import index_docs

HTML_DIR = Path(__file__).parent / "fixtures" / "html"


def fake_embed(texts, input_type, client=None):
    assert input_type == "document"
    return [[float(len(t)), 1.0] for t in texts]


def fake_code_chunks():
    return [
        {
            "text": "ws = WasteStream()",
            "title": "qsdsan.WasteStream",
            "url": "https://github.com/QSD-Group/QSDsan/blob/main/"
            "qsdsan/_waste_stream.py#L1-L2",
            "source": "code",
            "type": "doctest",
        }
    ]


def fake_example_chunks():
    return [
        {
            "text": "def test_make_stream(): ...",
            "title": "tests.test_streams.test_make_stream",
            "url": "https://github.com/QSD-Group/QSDsan/blob/main/"
            "tests/test_streams.py#L1-L2",
            "source": "example",
            "type": "test",
        }
    ]


def test_build_index_indexes_qsdsan_and_adds_embeddings():
    records = index_docs.build_index(
        html_dir=str(HTML_DIR),
        embed_fn=fake_embed,
        code_chunks_fn=lambda: [],
        example_chunks_fn=lambda: [],
    )
    assert records
    assert {r["source"] for r in records} == {"qsdsan"}
    assert all(len(r["embedding"]) == 2 for r in records)


def test_build_index_merges_code_and_example_chunks_with_embeddings():
    records = index_docs.build_index(
        html_dir=str(HTML_DIR),
        embed_fn=fake_embed,
        code_chunks_fn=fake_code_chunks,
        example_chunks_fn=fake_example_chunks,
    )
    assert {r["source"] for r in records} == {"qsdsan", "code", "example"}
    extra = [r for r in records if r["source"] in ("code", "example")]
    assert len(extra) == 2
    assert all(len(r["embedding"]) == 2 for r in extra)


def test_main_writes_index_json(tmp_path, monkeypatch):
    out = tmp_path / "index.json"
    monkeypatch.setattr(index_docs, "embed_documents", fake_embed)
    index_docs.main(html_dir=str(HTML_DIR), out_path=str(out))
    data = json.loads(out.read_text())
    assert isinstance(data, list) and data
    # main() runs the real adapters over the installed qsdsan (HTML + code +
    # the test-suite example snippets), so all three sources are published.
    sources = {r["source"] for r in data}
    assert {"qsdsan", "code", "example"} <= sources
    assert all("embedding" in r for r in data)


def test_build_index_raises_on_embedding_count_mismatch():
    def bad_embed(texts, input_type, client=None):
        return [[0.0, 1.0]]  # one vector regardless of how many texts

    with pytest.raises(ValueError):
        index_docs.build_index(
            html_dir=str(HTML_DIR),
            embed_fn=bad_embed,
            code_chunks_fn=lambda: [],
            example_chunks_fn=lambda: [],
        )


def test_main_writes_bare_filename(tmp_path, monkeypatch):
    # out_path with no directory component must not crash.
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(index_docs, "embed_documents", fake_embed)
    index_docs.main(html_dir=str(HTML_DIR), out_path="index.json")
    assert (tmp_path / "index.json").exists()


def test_main_writes_to_rtd_output_when_set(tmp_path, monkeypatch):
    # On readthedocs, the index must land in the BUILD OUTPUT _static so it is
    # published, not in the source tree (which was already copied pre-post_build).
    monkeypatch.setenv("READTHEDOCS_OUTPUT", str(tmp_path))
    monkeypatch.setattr(index_docs, "embed_documents", fake_embed)
    # html_dir passed explicitly to our fixtures; out_path left to default.
    index_docs.main(html_dir=str(HTML_DIR))
    published = tmp_path / "html" / "_static" / "chatbot" / "index.json"
    assert published.exists(), "index.json must be written into the RTD html output"


def test_main_defaults_to_source_static_without_rtd(tmp_path, monkeypatch):
    # Locally (no READTHEDOCS_OUTPUT), default to the in-tree source _static path.
    monkeypatch.delenv("READTHEDOCS_OUTPUT", raising=False)
    monkeypatch.chdir(tmp_path)
    (tmp_path / "docs" / "source" / "_static" / "chatbot").mkdir(parents=True)
    monkeypatch.setattr(index_docs, "embed_documents", fake_embed)
    index_docs.main(html_dir=str(HTML_DIR))
    assert (tmp_path / "docs" / "source" / "_static" / "chatbot" / "index.json").exists()
