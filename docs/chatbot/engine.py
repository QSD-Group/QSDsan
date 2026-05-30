"""Orchestrate one question: embed -> retrieve -> threshold/refuse -> Claude -> cite.

Pure and dependency-injected: embed_fn and claude_client are passed in so the
whole contract is testable without network.
"""
import re

import prompts
import retrieval


def _used_citation_numbers(answer: str, n: int) -> list[int]:
    """Return excerpt numbers actually referenced as [i] in the answer, in order."""
    seen = []
    for m in re.finditer(r"\[(\d+)\]", answer):
        i = int(m.group(1))
        if 1 <= i <= n and i not in seen:
            seen.append(i)
    return seen


def answer_question(
    question,
    records,
    embed_fn,
    claude_client,
    *,
    top_k,
    threshold,
    gen_model,
):
    """Return {"answer": str, "citations": list[dict]}."""
    query_vec = embed_fn([question], input_type="query")[0]
    hits = retrieval.cosine_top_k(query_vec, records, k=top_k)

    if not hits or hits[0][1] < threshold:
        return {"answer": prompts.refusal_message(), "citations": []}

    retrieved = [record for record, _score in hits]
    user_prompt = prompts.build_user_prompt(question, retrieved)

    resp = claude_client.messages.create(
        model=gen_model,
        max_tokens=1024,
        system=[
            {
                "type": "text",
                "text": prompts.SYSTEM_PROMPT,
                "cache_control": {"type": "ephemeral"},
            }
        ],
        messages=[{"role": "user", "content": user_prompt}],
    )
    answer = "".join(
        block.text for block in resp.content if getattr(block, "type", None) == "text"
    )
    answer = prompts.append_code_disclaimers(answer)

    used = _used_citation_numbers(answer, len(retrieved))
    if not used:  # answer cited nothing explicitly: surface all retrieved sources
        used = list(range(1, len(retrieved) + 1))
    citations = [
        {
            "n": i,
            "url": retrieved[i - 1]["url"],
            "title": retrieved[i - 1]["title"],
            "source": retrieved[i - 1]["source"],
        }
        for i in used
    ]
    return {"answer": answer, "citations": citations}
