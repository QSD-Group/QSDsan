"""Thin Voyage embedding wrapper used by the indexer and the server.

input_type is "document" when embedding chunks, "query" when embedding a question.
"""
import config


def _default_client():
    import voyageai

    return voyageai.Client()  # reads VOYAGE_API_KEY from the environment


def embed_texts(texts: list[str], input_type: str, client=None) -> list[list[float]]:
    """Return one embedding vector per input string."""
    if not texts:
        return []
    client = client or _default_client()
    result = client.embed(texts, model=config.EMBED_MODEL, input_type=input_type)
    return list(result.embeddings)
