import prompts

RETRIEVED = [
    {"text": "Use WasteStream to set flow and COD.", "title": "Create a WasteStream",
     "url": "https://qsdsan.readthedocs.io/en/latest/tutorials/3_WasteStream.html#x",
     "source": "qsdsan", "type": "tutorial"},
]


def test_system_prompt_states_grounding_and_citation_rules():
    text = prompts.SYSTEM_PROMPT.lower()
    assert "only" in text and "cite" in text
    assert "qsdsan" in text and "exposan" in text


def test_user_prompt_numbers_excerpts_and_includes_urls():
    up = prompts.build_user_prompt("How do I make a WasteStream?", RETRIEVED)
    assert "[1]" in up
    assert RETRIEVED[0]["url"] in up
    assert "How do I make a WasteStream?" in up


def test_append_code_disclaimers_adds_note_after_each_code_block():
    md = "Intro\n\n```python\nx = 1\n```\n\nOutro"
    out = prompts.append_code_disclaimers(md)
    assert out.count("Draft from the docs") == 1
    assert "verify against the linked pages" in out


def test_append_code_disclaimers_handles_multiple_blocks():
    md = "```python\na = 1\n```\ntext\n```python\nb = 2\n```"
    out = prompts.append_code_disclaimers(md)
    assert out.count("Draft from the docs") == 2


def test_append_code_disclaimers_noop_without_code():
    md = "Just prose, no code."
    assert prompts.append_code_disclaimers(md) == md


def test_refusal_message_points_to_nav_search():
    msg = prompts.refusal_message()
    assert "couldn't find this in the QSDsan/EXPOsan docs" in msg
    assert "search function at the top of the navigation bar" in msg
