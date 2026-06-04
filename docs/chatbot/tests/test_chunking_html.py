import chunking

SAMPLE_HTML = """
<html><body>
<section id="create-a-wastestream">
  <h2>Create a WasteStream</h2>
  <p>Use <code>WasteStream</code> to define flow and COD.</p>
</section>
<section id="run-a-simulation">
  <h2>Run a simulation</h2>
  <p>Call <code>sys.simulate()</code> to run it.</p>
</section>
</body></html>
"""


def test_splits_html_sections_by_heading():
    chunks = chunking.split_html_by_heading(SAMPLE_HTML)
    titles = [c["title"] for c in chunks]
    assert "Create a WasteStream" in titles
    assert "Run a simulation" in titles


def test_html_anchor_comes_from_section_id():
    chunks = chunking.split_html_by_heading(SAMPLE_HTML)
    ws = next(c for c in chunks if c["title"] == "Create a WasteStream")
    assert ws["anchor"] == "create-a-wastestream"


def test_html_text_is_plain_no_tags():
    chunks = chunking.split_html_by_heading(SAMPLE_HTML)
    ws = next(c for c in chunks if c["title"] == "Create a WasteStream")
    assert "flow and COD" in ws["text"]
    assert "<code>" not in ws["text"]


NESTED_FURO_HTML = """
<html><body>
<section id="wastestream">
  <h1>WasteStream<a class="headerlink" href="#wastestream" title="Link to this heading">¶</a></h1>
  <p>A WasteStream tracks flows and properties.</p>
  <section id="create-a-wastestream">
    <h2>Create a WasteStream<a class="headerlink" href="#create-a-wastestream">¶</a></h2>
    <p>Body A about flow and COD.</p>
  </section>
  <section id="run-a-simulation">
    <h2>Run a simulation<a class="headerlink" href="#run-a-simulation">¶</a></h2>
    <p>Body B about simulate.</p>
  </section>
</section>
</body></html>
"""


def test_headerlink_pilcrow_not_in_title_or_text():
    chunks = chunking.split_html_by_heading(NESTED_FURO_HTML)
    for c in chunks:
        assert "¶" not in c["title"]
        assert "¶" not in c["text"]


def test_heading_text_not_leaked_into_body():
    chunks = chunking.split_html_by_heading(NESTED_FURO_HTML)
    create = next(c for c in chunks if c["title"] == "Create a WasteStream")
    assert "Body A about flow and COD" in create["text"]
    assert not create["text"].startswith("Create a WasteStream")


def test_h1_chunk_excludes_nested_subsection_bodies():
    chunks = chunking.split_html_by_heading(NESTED_FURO_HTML)
    top = next(c for c in chunks if c["title"] == "WasteStream")
    assert "WasteStream tracks flows" in top["text"]
    # The h1 chunk must NOT contain the child subsections' bodies.
    assert "Body A about flow and COD" not in top["text"]
    assert "Body B about simulate" not in top["text"]
    assert top["anchor"] == "wastestream"
