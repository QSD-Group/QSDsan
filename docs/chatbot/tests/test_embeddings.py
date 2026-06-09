import embeddings


class FakeResult:
    def __init__(self, vectors):
        self.embeddings = vectors


class FakeClient:
    def __init__(self):
        self.calls = []

    def embed(self, texts, model, input_type):
        self.calls.append((tuple(texts), model, input_type))
        return FakeResult([[0.1, 0.2, 0.3] for _ in texts])


def test_embed_texts_returns_one_vector_per_input():
    client = FakeClient()
    vecs = embeddings.embed_texts(["a", "b"], input_type="document", client=client)
    assert len(vecs) == 2
    assert vecs[0] == [0.1, 0.2, 0.3]


def test_embed_texts_passes_model_and_input_type():
    client = FakeClient()
    embeddings.embed_texts(["q"], input_type="query", client=client)
    (_texts, model, input_type) = client.calls[0]
    assert model == "voyage-4-lite"
    assert input_type == "query"


def test_embed_texts_empty_returns_empty():
    client = FakeClient()
    assert embeddings.embed_texts([], input_type="document", client=client) == []


def test_embed_texts_batches_requests_under_provider_limit():
    # Voyage rejects batches over 1000 inputs; the wrapper must split large
    # corpora into capped requests and still return one vector per input in order.
    client = FakeClient()
    n = embeddings._MAX_BATCH * 2 + 1
    texts = [f"t{i}" for i in range(n)]
    vecs = embeddings.embed_texts(texts, input_type="document", client=client)
    assert len(vecs) == n
    assert all(len(call_texts) <= embeddings._MAX_BATCH for call_texts, *_ in client.calls)
    assert len(client.calls) == 3  # 2 full batches + 1 remainder
