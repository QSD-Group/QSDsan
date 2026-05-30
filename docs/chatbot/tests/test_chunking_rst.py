import chunking

SAMPLE_RST = """\
=========
BSM1
=========

Summary
-------
This system models BSM1.

Load and run
------------
Use ``create_system`` to build it.
"""


def test_splits_each_heading_into_a_chunk():
    chunks = chunking.split_rst_by_heading(SAMPLE_RST)
    titles = [c["title"] for c in chunks]
    assert "Summary" in titles
    assert "Load and run" in titles


def test_chunk_text_contains_body_under_heading():
    chunks = chunking.split_rst_by_heading(SAMPLE_RST)
    load = next(c for c in chunks if c["title"] == "Load and run")
    assert "create_system" in load["text"]
    assert "models BSM1" not in load["text"]  # body of a different section


def test_anchor_is_github_style_slug():
    chunks = chunking.split_rst_by_heading(SAMPLE_RST)
    load = next(c for c in chunks if c["title"] == "Load and run")
    assert load["anchor"] == "load-and-run"
