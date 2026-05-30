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
