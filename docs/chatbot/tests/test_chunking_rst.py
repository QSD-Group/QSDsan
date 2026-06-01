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


LITERAL_BLOCK_RST = """\
System Simulation
-----------------
Run the system::

    >>> biobinder.simulate_and_print(sys)
    biobinder
    ---------
    Received flowsheet: None

More prose after the block.
"""


def test_indented_literal_block_is_not_a_heading():
    chunks = chunking.split_rst_by_heading(LITERAL_BLOCK_RST)
    titles = [c["title"] for c in chunks]
    # The indented "biobinder / ---------" inside the literal block must NOT
    # become its own section.
    assert titles == ["System Simulation"]
    sim = chunks[0]
    # The real section keeps the whole literal block, including the trailing prose.
    assert "simulate_and_print" in sim["text"]
    assert "More prose after the block" in sim["text"]
