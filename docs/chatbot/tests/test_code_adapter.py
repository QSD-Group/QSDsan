from pathlib import Path

import code_adapter

BLOB_BASE = "https://github.com/QSD-Group/QSDsan/blob/main"
EXPOSAN_BLOB = "https://github.com/QSD-Group/EXPOsan/blob/main"
CODEPKG_ROOT = Path(__file__).parent / "fixtures" / "codepkg"
REPO_ROOT = Path(__file__).resolve().parents[3]

# A tiny module with exactly one doctest, laid out so line numbers are known:
#   line  8: the method docstring opens ('"""Build ...')
#   line 10: >>> ws = WasteStream()      (first doctest source line)
#   line 11: >>> ws.from_composite(2)
#   line 12: 4                           (last doctest output line)
# so the doctest block spans file lines 10-12.
ONE_DOCTEST_MODULE = '''\
"""Top-level module docstring, no doctest here."""


class WasteStream:
    """A stream of waste."""

    def from_composite(self, x):
        """Build a WasteStream from a composite.

        >>> ws = WasteStream()
        >>> ws.from_composite(2)
        4
        """
        return x * 2
'''


def test_one_doctest_yields_one_chunk_with_right_schema():
    chunks = code_adapter.extract_doctest_chunks(
        ONE_DOCTEST_MODULE,
        module_name="qsdsan",
        rel_path="qsdsan/_waste_stream.py",
        blob_base=BLOB_BASE,
    )
    assert len(chunks) == 1
    chunk = chunks[0]
    assert chunk["title"] == "qsdsan.WasteStream.from_composite"
    assert chunk["source"] == "code"
    assert chunk["type"] == "doctest"
    assert chunk["url"] == (
        "https://github.com/QSD-Group/QSDsan/blob/main/"
        "qsdsan/_waste_stream.py#L10-L12"
    )
    assert ">>> ws.from_composite(2)" in chunk["text"]


# Two methods with doctests, plus prose-only docstrings that must NOT produce
# chunks (module docstring and a summary-only method).
TWO_DOCTEST_MODULE = '''\
"""Module docstring with no doctest."""


class Stream:
    """A stream; this prose-only docstring has no doctest."""

    def cod(self):
        """Return the COD.

        >>> Stream().cod()
        42
        """
        return 42

    def describe(self):
        """Summary line only, no example here."""
        return "a stream"

    def tss(self):
        """Return the TSS.

        >>> Stream().tss()
        7
        """
        return 7
'''


def test_one_chunk_per_docstring_with_a_doctest():
    chunks = code_adapter.extract_doctest_chunks(
        TWO_DOCTEST_MODULE,
        module_name="qsdsan",
        rel_path="qsdsan/_stream.py",
        blob_base=BLOB_BASE,
    )
    titles = sorted(c["title"] for c in chunks)
    assert titles == ["qsdsan.Stream.cod", "qsdsan.Stream.tss"]


# A doctest whose statement spans multiple lines, exercising "..." continuation
# reconstruction and multi-line span counting.
CONTINUATION_MODULE = '''\
def build():
    """Build a dict.

    >>> d = {
    ...     'a': 1,
    ... }
    >>> d['a']
    1
    """
    return {}
'''


def test_continuation_lines_are_reconstructed_with_ellipsis():
    chunks = code_adapter.extract_doctest_chunks(
        CONTINUATION_MODULE,
        module_name="qsdsan",
        rel_path="qsdsan/_build.py",
        blob_base=BLOB_BASE,
    )
    assert len(chunks) == 1
    text = chunks[0]["text"]
    assert ">>> d = {" in text
    assert "...     'a': 1," in text
    assert "... }" in text
    # Block spans the first '>>>' (line 4) through the last output '1' (line 8).
    assert chunks[0]["url"].endswith("qsdsan/_build.py#L4-L8")


# Real qsdsan docstrings contain malformed doctest directives (e.g. "+NUMBER",
# which is not a registered doctest option in a plain-Python build). doctest
# raises ValueError on these; one bad docstring must not crash the whole index
# build or lose the good doctests. (A name pytest never registers is used so the
# test reproduces the production failure even though pytest registers +NUMBER.)
INVALID_OPTION_MODULE = '''\
def good():
    """Good example.

    >>> good()
    1
    """
    return 1


def bad():
    """Has a bogus doctest option.

    >>> bad()  # doctest: +QSDSAN_NOT_AN_OPTION
    2
    """
    return 2
'''


def test_invalid_doctest_option_is_skipped_not_fatal():
    chunks = code_adapter.extract_doctest_chunks(
        INVALID_OPTION_MODULE,
        module_name="qsdsan",
        rel_path="qsdsan/_x.py",
        blob_base=BLOB_BASE,
    )
    titles = [c["title"] for c in chunks]
    assert titles == ["qsdsan.good"]


def test_walks_package_dir_for_doctests_with_dotted_names_and_blob_paths():
    pkg_dir = CODEPKG_ROOT / "mypkg"
    chunks = code_adapter.build_code_chunks_from_dir(
        str(pkg_dir), str(CODEPKG_ROOT), blob_base=BLOB_BASE
    )
    by_title = {c["title"]: c for c in chunks}
    # Dotted module names reflect the package layout, including subpackages.
    assert "mypkg.core.Widget.size" in by_title
    assert "mypkg.sub.helpers.doubled" in by_title
    # Blob URL path is relative to the repo root (the package's parent dir).
    assert by_title["mypkg.core.Widget.size"]["url"].startswith(
        f"{BLOB_BASE}/mypkg/core.py#L"
    )
    assert by_title["mypkg.sub.helpers.doubled"]["url"].startswith(
        f"{BLOB_BASE}/mypkg/sub/helpers.py#L"
    )
    assert all(c["source"] == "code" for c in chunks)


# --- API-surface extractor ---------------------------------------------------

# A module exercising defaults, a private function (must be skipped), and a
# multi-line docstring whose summary is only the first line. Line numbers:
#   line  9: 'def from_composite(' opens the public function
#   line 10: its docstring opens, so the signature span is L9-L9.
API_FUNCTION_MODULE = '''\
"""Module docstring."""


def _private(x):
    """Hidden helper."""
    return x


def from_composite(x, unit='mg/L'):
    """Build from a composite.

    A longer explanation that must NOT appear in the one-line summary.
    """
    return x
'''


def test_api_chunk_renders_signature_and_summary_skips_private():
    chunks = code_adapter.extract_api_chunks(
        API_FUNCTION_MODULE,
        module_name="qsdsan",
        rel_path="qsdsan/_waste_stream.py",
        blob_base=BLOB_BASE,
    )
    by_title = {c["title"]: c for c in chunks}
    # Private (underscore-prefixed) names are not part of the public API surface.
    assert "qsdsan._private" not in by_title
    chunk = by_title["qsdsan.from_composite"]
    assert chunk["source"] == "code"
    assert chunk["type"] == "api"
    assert "def from_composite(x, unit='mg/L')" in chunk["text"]
    # One-line summary only, not the longer explanation.
    assert "Build from a composite." in chunk["text"]
    assert "longer explanation" not in chunk["text"]
    assert chunk["url"] == (
        f"{BLOB_BASE}/qsdsan/_waste_stream.py#L9-L9"
    )


# A class with a base, a decorated property method, and a return annotation.
API_CLASS_MODULE = '''\
"""Module docstring."""


class WasteStream(SanStream):
    """A stream of waste."""

    @property
    def COD(self) -> float:
        """Chemical oxygen demand."""
        return 1.0
'''


def test_api_chunk_for_class_bases_and_decorated_method():
    chunks = code_adapter.extract_api_chunks(
        API_CLASS_MODULE,
        module_name="qsdsan",
        rel_path="qsdsan/_waste_stream.py",
        blob_base=BLOB_BASE,
    )
    by_title = {c["title"]: c for c in chunks}
    cls = by_title["qsdsan.WasteStream"]
    assert "class WasteStream(SanStream):" in cls["text"]
    assert "A stream of waste." in cls["text"]
    assert cls["type"] == "api"
    meth = by_title["qsdsan.WasteStream.COD"]
    assert "@property" in meth["text"]
    assert "def COD(self) -> float" in meth["text"]
    assert "Chemical oxygen demand." in meth["text"]


# An undocumented public method is low signal (a bare signature with no summary)
# and dilutes retrieval; classes are kept for their structural value.
API_UNDOCUMENTED_MODULE = '''\
class Tool:
    """A tool."""

    def documented(self):
        """Does a thing."""
        return 1

    def undocumented(self, x):
        return x
'''


def test_api_skips_undocumented_functions_but_keeps_classes():
    chunks = code_adapter.extract_api_chunks(
        API_UNDOCUMENTED_MODULE,
        module_name="qsdsan",
        rel_path="qsdsan/_t.py",
        blob_base=BLOB_BASE,
    )
    titles = {c["title"] for c in chunks}
    assert "qsdsan.Tool" in titles
    assert "qsdsan.Tool.documented" in titles
    assert "qsdsan.Tool.undocumented" not in titles


def test_package_walker_emits_both_api_and_doctest_chunks():
    chunks = code_adapter.build_code_chunks_from_dir(
        str(CODEPKG_ROOT / "mypkg"), str(CODEPKG_ROOT), blob_base=BLOB_BASE
    )
    types = {c["type"] for c in chunks}
    assert {"api", "doctest"} <= types
    pairs = {(c["title"], c["type"]) for c in chunks}
    # The class itself yields an API chunk; the method yields both API + doctest.
    assert ("mypkg.core.Widget", "api") in pairs
    assert ("mypkg.core.Widget.size", "api") in pairs
    assert ("mypkg.core.Widget.size", "doctest") in pairs


# --- test-function snippets (source="example") -------------------------------

# Test functions are known-correct usage; helpers and non-test code are not.
# Line numbers: 'def test_creates_widget' opens at line 8, ends at line 10.
TEST_MODULE = '''\
import pytest


def _helper():
    return 1


def test_creates_widget():
    w = Widget()
    assert w.size() == 3


class TestStream:
    def test_cod(self):
        assert Stream().cod() == 42
'''


def test_test_snippets_capture_test_functions_only():
    chunks = code_adapter.extract_test_chunks(
        TEST_MODULE,
        module_name="tests.test_widget",
        rel_path="tests/test_widget.py",
        blob_base=BLOB_BASE,
    )
    by_title = {c["title"]: c for c in chunks}
    assert "tests.test_widget._helper" not in by_title  # non-test skipped
    snippet = by_title["tests.test_widget.test_creates_widget"]
    assert snippet["source"] == "example"
    assert snippet["type"] == "test"
    assert "def test_creates_widget():" in snippet["text"]
    assert "assert w.size() == 3" in snippet["text"]
    assert snippet["url"] == f"{BLOB_BASE}/tests/test_widget.py#L8-L10"
    # Test methods on a class are captured and dedented to a runnable snippet.
    method = by_title["tests.test_widget.TestStream.test_cod"]
    assert method["text"].startswith("def test_cod(self):")


# --- EXPOsan create_system wiring (source="example") -------------------------

# Only the system-assembly function is high-signal wiring; module helpers are not.
SYSTEM_MODULE = '''\
"""bsm1 system."""

import qsdsan as qs


def create_system(flowsheet=None):
    """Create the BSM1 system."""
    flowsheet = flowsheet or qs.Flowsheet('bsm1')
    A1 = qs.sanunits.CSTR('A1')
    sys = qs.System('bsm1_sys', path=(A1,))
    return sys


def _helper():
    return 1
'''


def test_system_chunk_captures_create_system_wiring():
    chunks = code_adapter.extract_system_chunks(
        SYSTEM_MODULE,
        module_name="exposan.bsm1.system",
        rel_path="exposan/bsm1/system.py",
        blob_base=EXPOSAN_BLOB,
    )
    assert len(chunks) == 1
    chunk = chunks[0]
    assert chunk["title"] == "exposan.bsm1.system.create_system"
    assert chunk["source"] == "example"
    assert chunk["type"] == "system"
    assert "def create_system(flowsheet=None):" in chunk["text"]
    assert "qs.System('bsm1_sys', path=(A1,))" in chunk["text"]
    assert chunk["url"] == f"{EXPOSAN_BLOB}/exposan/bsm1/system.py#L6-L11"
    assert "truncated" not in chunk["text"]  # short function: no truncation


# Only module-level create_system is true system wiring; a same-named method on
# an analysis-model class is not.
SYSTEM_WITH_METHOD_MODULE = '''\
def create_system():
    """Module-level system."""
    return 1


class Model:
    def create_system(self):
        """A model helper method, not system wiring."""
        return 2
'''


def test_system_predicate_is_module_level_only():
    chunks = code_adapter.extract_system_chunks(
        SYSTEM_WITH_METHOD_MODULE,
        module_name="exposan.x.systems",
        rel_path="exposan/x/systems.py",
        blob_base=EXPOSAN_BLOB,
    )
    titles = [c["title"] for c in chunks]
    assert titles == ["exposan.x.systems.create_system"]


def test_oversized_snippet_is_truncated_but_citation_spans_full_function():
    body = "\n".join(f"    x{i} = {i}" for i in range(200))
    src = f'def create_system():\n    """Big system."""\n{body}\n    return 1\n'
    chunks = code_adapter.extract_system_chunks(
        src, "exposan.x.systems", "exposan/x/systems.py", blob_base=EXPOSAN_BLOB
    )
    text = chunks[0]["text"]
    assert "truncated" in text
    assert len(text.splitlines()) <= code_adapter._MAX_SNIPPET_LINES + 1
    # The citation still points to the WHOLE function, not the capped extract.
    start_s, end_s = chunks[0]["url"].split("#L")[1].split("-L")
    assert int(start_s) == 1
    assert int(end_s) > code_adapter._MAX_SNIPPET_LINES


def test_build_test_chunks_over_the_real_suite():
    chunks = code_adapter.build_test_chunks_from_dir(
        str(REPO_ROOT / "tests"), str(REPO_ROOT), blob_base=BLOB_BASE
    )
    assert chunks
    assert all(c["source"] == "example" and c["type"] == "test" for c in chunks)
    assert all(c["title"].split(".")[-1].startswith("test") for c in chunks)
    assert any(c["url"].startswith(f"{BLOB_BASE}/tests/") for c in chunks)


def test_build_exposan_chunks_is_graceful_when_package_absent():
    # An unimportable package must yield no chunks, not raise.
    assert code_adapter.build_exposan_chunks("not_a_real_pkg_qsdsan_xyz") == []


def test_parsing_invalid_escape_sequences_does_not_emit_syntaxwarning():
    # Real EXPOsan/qsdsan sources contain non-raw strings with invalid escapes
    # (e.g. '\d', '\$'); ast.parse warns on those. The adapter must not spew
    # SyntaxWarnings into the build log when walking such files.
    import warnings

    src = "def f():\n    \"\"\"Doc.\"\"\"\n    return '\\d+'\n"
    with warnings.catch_warnings():
        warnings.simplefilter("error", SyntaxWarning)  # any SyntaxWarning -> error
        code_adapter.extract_api_chunks(
            src, "qsdsan", "qsdsan/_x.py", blob_base=BLOB_BASE
        )
