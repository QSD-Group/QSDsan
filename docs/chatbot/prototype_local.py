"""Local interactive prototype for the docs chatbot.

Builds an in-memory index from the code + example adapters (and the built HTML
if ``docs/build/html`` exists), then runs the real RAG loop
(embed -> retrieve -> threshold/refuse -> Claude -> cite) in a REPL so you can
see how retrieval over qsdsan source / doctests / tests / EXPOsan systems behaves
before deploying.

This is a dev harness, not part of the served app (the endpoint is ``server.py``).
It only orchestrates already-tested functions (build_* / engine.answer_question).

Requires VOYAGE_API_KEY and ANTHROPIC_API_KEY in the environment.
Run from docs/chatbot:

    python prototype_local.py

Optional: build the HTML docs first (``cd docs && make html``) to also index the
tutorials/API pages; otherwise the prototype covers code + examples only.
"""
import os
import sys

import code_adapter
import config
import embeddings
import engine
import index_docs

_EMBED_BATCH = 128  # Voyage caps inputs per request; embed documents in batches.


def _require_keys() -> None:
    missing = [k for k in ("VOYAGE_API_KEY", "ANTHROPIC_API_KEY") if not os.getenv(k)]
    if missing:
        sys.exit(f"Set these environment variables first: {', '.join(missing)}")


def _build_records() -> list[dict]:
    """Gather chunks from every available adapter and attach embeddings."""
    chunks: list[dict] = []
    html_dir = os.path.join("..", "build", "html")  # docs/build/html
    if os.path.isdir(html_dir):
        html_chunks = index_docs.build_qsdsan_chunks(html_dir)
        chunks.extend(html_chunks)
        print(f"  HTML chunks:         {len(html_chunks)}")
    else:
        print("  HTML chunks:         0 (run `cd docs && make html` to include them)")

    code_chunks = code_adapter.build_code_chunks("qsdsan")
    example_chunks = code_adapter.build_example_chunks("qsdsan")
    chunks.extend(code_chunks)
    chunks.extend(example_chunks)
    print(f"  code chunks:         {len(code_chunks)}")
    print(f"  example chunks:      {len(example_chunks)}")

    texts = [c["text"] for c in chunks]
    vectors: list[list[float]] = []
    for i in range(0, len(texts), _EMBED_BATCH):
        batch = texts[i : i + _EMBED_BATCH]
        vectors.extend(embeddings.embed_texts(batch, input_type="document"))
        print(f"  embedding {min(i + _EMBED_BATCH, len(texts))}/{len(texts)}", end="\r")
    print()
    for chunk, vec in zip(chunks, vectors):
        chunk["embedding"] = list(vec)
    return chunks


def main() -> None:
    _require_keys()
    print("Building in-memory index...")
    records = _build_records()
    if not records:
        sys.exit("No chunks built; nothing to query.")

    import anthropic

    client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY
    print(
        f"\nReady: {len(records)} chunks "
        f"(model={config.GEN_MODEL}, top_k={config.TOP_K}, "
        f"threshold={config.SIMILARITY_THRESHOLD}).\n"
        "Ask a question (blank line or Ctrl-C to quit).\n"
    )

    while True:
        try:
            question = input("ask> ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break
        if not question:
            break
        result = engine.answer_question(
            question,
            records=records,
            embed_fn=embeddings.embed_texts,
            claude_client=client,
            top_k=config.TOP_K,
            threshold=config.SIMILARITY_THRESHOLD,
            gen_model=config.GEN_MODEL,
        )
        print("\n" + result["answer"] + "\n")
        if result["citations"]:
            print("Citations:")
            for c in result["citations"]:
                print(f"  [{c['n']}] ({c['source']}) {c['title']}\n       {c['url']}")
        print()


if __name__ == "__main__":
    main()
