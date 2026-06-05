import engine

RECORDS = [
    {"text": "Use WasteStream to set flow and COD.", "title": "Create a WasteStream",
     "url": "https://qsdsan.readthedocs.io/en/latest/tutorials/3_WasteStream.html#x",
     "source": "qsdsan", "type": "tutorial", "embedding": [1.0, 0.0]},
    {"text": "Run create_system then sys.simulate().", "title": "Load and run",
     "url": "https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bsm1/README.rst#load-and-run",
     "source": "exposan", "type": "readme", "embedding": [0.0, 1.0]},
]


def embed_query_high(texts, input_type, client=None):
    # Aligns with the first record (cosine ~1.0).
    return [[1.0, 0.0]]


class FakeClaude:
    """Minimal stand-in for anthropic.Anthropic().messages.create."""

    def __init__(self, text):
        self._text = text
        self.received = None

        class Messages:
            def __init__(self, outer):
                self._outer = outer

            def create(self, **kwargs):
                self._outer.received = kwargs

                class Block:
                    type = "text"
                    text = self._outer._text

                class Resp:
                    content = [Block()]

                return Resp()

        self.messages = Messages(self)


def test_grounded_answer_returns_citations():
    claude = FakeClaude("Use WasteStream [1].")
    out = engine.answer_question(
        "How do I make a WasteStream?",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=2,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    assert "WasteStream" in out["answer"]
    assert out["citations"]
    assert out["citations"][0]["url"] == RECORDS[0]["url"]
    assert out["citations"][0]["source"] == "qsdsan"


def test_low_similarity_refuses_without_calling_claude():
    claude = FakeClaude("SHOULD NOT BE CALLED")
    out = engine.answer_question(
        "How do I integrate SCADA/Modbus?",
        records=RECORDS,
        embed_fn=lambda t, input_type, client=None: [[-1.0, 0.0]],  # cosine -1
        claude_client=claude,
        top_k=2,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    assert claude.received is None  # Claude never invoked
    assert "couldn't find this" in out["answer"]
    assert out["citations"] == []


def test_system_prompt_is_cached():
    claude = FakeClaude("Answer [1].")
    engine.answer_question(
        "How do I make a WasteStream?",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=1,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    system = claude.received["system"]
    assert isinstance(system, list)
    assert system[0]["cache_control"] == {"type": "ephemeral"}


def test_code_block_gets_disclaimer_appended():
    claude = FakeClaude("Here:\n\n```python\nws = WasteStream('ws')\n```\n")
    out = engine.answer_question(
        "Make a WasteStream",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=1,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    assert "Draft from the docs" in out["answer"]


def test_answer_without_explicit_citation_falls_back_to_all_retrieved():
    claude = FakeClaude("A general answer with no bracket markers.")
    out = engine.answer_question(
        "How do I make a WasteStream?",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=2,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    # No [n] markers in the answer -> surface all retrieved as citations.
    assert len(out["citations"]) == 2


def test_citations_are_bounded_and_subset_only():
    claude = FakeClaude("Only the second excerpt is relevant; see [2]. Ignore [5].")
    out = engine.answer_question(
        "How do I set properties on a WasteStream?",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=2,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    # Only the cited, in-range excerpt [2] is surfaced; [5] (out of range) is dropped.
    assert [c["n"] for c in out["citations"]] == [2]
    assert out["citations"][0]["url"] == RECORDS[1]["url"]


def test_exposan_catalog_question_returns_pointer_without_calling_claude():
    # Catalog/discovery questions are pointed to the Systems page + repo. (A code
    # how-to like "how do I run BSM1?" instead falls through to retrieval now.)
    claude = FakeClaude("SHOULD NOT BE CALLED")
    out = engine.answer_question(
        "What systems does QSDsan include?",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=2,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    assert claude.received is None  # Claude never invoked
    assert "github.com/QSD-Group/EXPOsan" in out["answer"]
    assert "systems/index.html" in out["answer"]
    assert out["citations"] == []


def test_smalltalk_returns_greeting_without_calling_claude():
    claude = FakeClaude("SHOULD NOT BE CALLED")
    out = engine.answer_question(
        "hello",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=2,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    assert claude.received is None  # Claude never invoked
    assert "QSDsan documentation assistant" in out["answer"]
    assert "couldn't find" not in out["answer"]
    assert out["citations"] == []


def test_empty_model_answer_refuses():
    claude = FakeClaude("")  # model returned no usable text
    out = engine.answer_question(
        "How do I make a WasteStream?",
        records=RECORDS,
        embed_fn=embed_query_high,
        claude_client=claude,
        top_k=2,
        threshold=0.5,
        gen_model="claude-haiku-4-5-20251001",
    )
    assert "couldn't find this" in out["answer"]
    assert out["citations"] == []
