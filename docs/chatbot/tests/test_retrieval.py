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


def test_cosine_top_k_empty_records_returns_empty():
    assert retrieval.cosine_top_k([1.0, 0.0], [], k=3) == []


def test_load_index_reads_json_file(tmp_path):
    p = tmp_path / "index.json"
    p.write_text(json.dumps(RECORDS))
    loaded = retrieval.load_index(str(p))
    assert len(loaded) == 3
