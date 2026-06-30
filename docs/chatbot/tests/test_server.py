from fastapi.testclient import TestClient

import server


def test_ask_returns_answer_and_citations(monkeypatch):
    def fake_answer(question, **kwargs):
        return {
            "answer": "Use WasteStream [1].",
            "citations": [
                {"n": 1, "url": "https://example/x", "title": "Create a WasteStream",
                 "source": "qsdsan"}
            ],
        }

    monkeypatch.setattr(server, "_answer", fake_answer)
    client = TestClient(server.app)
    resp = client.post("/ask", json={"question": "How do I make a WasteStream?"})
    assert resp.status_code == 200
    body = resp.json()
    assert body["answer"] == "Use WasteStream [1]."
    assert body["citations"][0]["source"] == "qsdsan"


def test_ask_rejects_empty_question(monkeypatch):
    monkeypatch.setattr(server, "_answer", lambda question, **kw: {"answer": "", "citations": []})
    client = TestClient(server.app)
    resp = client.post("/ask", json={"question": "   "})
    assert resp.status_code == 422


def test_healthcheck_ok():
    client = TestClient(server.app)
    assert client.get("/health").status_code == 200
