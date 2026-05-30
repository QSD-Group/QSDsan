from pathlib import Path

import index_docs

FIXTURE = (Path(__file__).parent / "fixtures" / "bsm1_README.rst").read_text()


def fake_fetch(url):
    # Only the bsm1 README exists in this fake world; everything else 404s (None).
    if url.endswith("exposan/bsm1/README.rst"):
        return FIXTURE
    return None


def test_builds_chunks_only_for_systems_with_readme():
    chunks = index_docs.build_exposan_chunks(
        systems=["bsm1", "ghost"], fetch=fake_fetch
    )
    assert chunks, "expected at least one chunk"
    assert all(c["source"] == "exposan" for c in chunks)
    assert all(c["type"] == "readme" for c in chunks)
    # 'ghost' has no README, so no chunks reference it.
    assert not any("ghost" in c["url"] for c in chunks)


def test_exposan_url_is_resolvable_github_blob():
    chunks = index_docs.build_exposan_chunks(systems=["bsm1"], fetch=fake_fetch)
    run = next(c for c in chunks if c["title"] == "Load and run")
    assert run["url"].startswith(
        "https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bsm1/README.rst"
    )
    assert run["url"].endswith("#load-and-run")


def test_every_chunk_has_text_and_title():
    chunks = index_docs.build_exposan_chunks(systems=["bsm1"], fetch=fake_fetch)
    for c in chunks:
        assert c["title"]
        assert c["text"]
