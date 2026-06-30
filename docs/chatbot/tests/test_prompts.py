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


def test_refusal_message_points_to_search_and_github():
    msg = prompts.refusal_message()
    assert "couldn't find this" in msg
    assert "search bar" in msg
    assert "github.com/QSD-Group/QSDsan" in msg


def test_is_smalltalk_detects_greetings_and_meta():
    for q in ["hello", "Hi!", "hey", "thanks", "good morning",
              "help", "what can you do?", "what can you help me with"]:
        assert prompts.is_smalltalk(q), q


def test_is_smalltalk_ignores_technical_questions():
    assert not prompts.is_smalltalk("How do I create a WasteStream?")
    assert not prompts.is_smalltalk("help me create a WasteStream")
    assert not prompts.is_smalltalk("what systems are available")


def test_greeting_message_is_a_welcome():
    msg = prompts.greeting_message()
    assert "QSDsan documentation assistant" in msg
    assert "couldn't find" not in msg


def test_is_exposan_question_detects_systems_and_exposan():
    assert prompts.is_exposan_question("How do I run BSM1?")
    assert prompts.is_exposan_question("Tell me about EXPOsan")
    assert prompts.is_exposan_question("what systems does QSDsan include?")
    assert prompts.is_exposan_question("what have you built?")


def test_is_exposan_question_ignores_plain_qsdsan_questions():
    assert not prompts.is_exposan_question("How do I create a WasteStream?")
    assert not prompts.is_exposan_question("How do I run a dynamic simulation?")


def test_exposan_pointer_message_has_both_links():
    msg = prompts.exposan_pointer_message()
    assert "systems/index.html" in msg
    assert "github.com/QSD-Group/EXPOsan" in msg
