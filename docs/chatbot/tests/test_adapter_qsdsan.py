from pathlib import Path

import index_docs

HTML_DIR = Path(__file__).parent / "fixtures" / "html"


def test_builds_chunks_from_built_html_tree():
    chunks = index_docs.build_qsdsan_chunks(str(HTML_DIR))
    assert chunks
    assert all(c["source"] == "qsdsan" for c in chunks)


def test_tutorial_pages_tagged_tutorial_api_pages_tagged_api():
    chunks = index_docs.build_qsdsan_chunks(str(HTML_DIR))
    tut = next(c for c in chunks if "flow and COD" in c["text"])
    api = next(c for c in chunks if "Class reference" in c["text"])
    assert tut["type"] == "tutorial"
    assert api["type"] == "api"


def test_qsdsan_url_is_readthedocs_page_plus_anchor():
    chunks = index_docs.build_qsdsan_chunks(str(HTML_DIR))
    ws = next(c for c in chunks if c["title"] == "Create a WasteStream")
    assert ws["url"] == (
        "https://qsdsan.readthedocs.io/en/latest/"
        "tutorials/3_WasteStream.html#create-a-wastestream"
    )


def test_homepage_install_section_is_indexed():
    # Installation lives on the docs homepage (index.html), not under tutorials/api,
    # so the indexer must cover the homepage for the bot to answer install questions.
    chunks = index_docs.build_qsdsan_chunks(str(HTML_DIR))
    inst = next(c for c in chunks if c["title"] == "Installation")
    assert inst["type"] == "guide"
    assert "pip install qsdsan" in inst["text"]
    assert inst["url"] == (
        "https://qsdsan.readthedocs.io/en/latest/index.html#installation"
    )


def test_contributing_dev_setup_is_indexed():
    # CONTRIBUTING.html is a top-level guide page (not under tutorials/api), so it
    # must be whitelisted for the bot to answer dev-setup questions like the uv
    # install workflow.
    chunks = index_docs.build_qsdsan_chunks(str(HTML_DIR))
    setup = next(c for c in chunks if c["title"] == "Setting up a development environment")
    assert setup["type"] == "guide"
    assert "uv" in setup["text"]
    assert setup["url"] == (
        "https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html#setting-up"
    )
