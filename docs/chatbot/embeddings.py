"""Thin Voyage embedding wrapper used by the indexer and the server.

input_type is "document" when embedding chunks, "query" when embedding a question.
"""
import config

# Voyage rejects batches over 1000 inputs per request; stay well under to also
# keep each request's total token count comfortable for the whole corpus.
_MAX_BATCH = 128


def _default_client():
    import voyageai

    return voyageai.Client()  # reads VOYAGE_API_KEY from the environment


def embed_texts(texts: list[str], input_type: str, client=None) -> list[list[float]]:
    """Return one embedding vector per input string.

    Splits large corpora into capped requests so the provider's per-request batch
    limit is never exceeded; vectors are returned in input order.
    """
    if not texts:
        return []
    client = client or _default_client()
    vectors: list[list[float]] = []
    for start in range(0, len(texts), _MAX_BATCH):
        batch = texts[start : start + _MAX_BATCH]
        result = client.embed(batch, model=config.EMBED_MODEL, input_type=input_type)
        vectors.extend(result.embeddings)
    return list(vectors)
