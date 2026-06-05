"""Code adapter: emit structured, high-signal chunks from Python source.

Unlike the HTML/RST adapters, this walks Python source with ``ast`` and the
stdlib ``doctest`` parser to surface curated usage examples that the autodoc
HTML buries. Chunks use the same schema as every other adapter
(``{text, title, url, source, type}``) so retrieval treats them uniformly.

Pure functions over source strings: no installed package and no network, so
tests run against small recorded module strings.
"""
from __future__ import annotations

import ast
import doctest
import importlib.util
import os
import textwrap
import warnings

import config

_FUNC_TYPES = (ast.FunctionDef, ast.AsyncFunctionDef)

# Some EXPOsan create_system functions run to hundreds of lines (the systems are
# written by many contributors with varied style). Cap snippet chunks so one huge
# excerpt cannot dominate the top-k prompt; the citation still spans the full def.
_MAX_SNIPPET_LINES = 80
_TRUNCATION_NOTE = "# ... (truncated - see linked source)"


def _reconstruct_doctest(examples) -> str:
    """Rebuild the interactive doctest block (>>> / ... / output) as text."""
    lines: list[str] = []
    for ex in examples:
        src_lines = ex.source.rstrip("\n").split("\n")
        lines.append(">>> " + src_lines[0])
        for cont in src_lines[1:]:
            lines.append("... " + cont)
        if ex.want:
            lines.extend(ex.want.rstrip("\n").split("\n"))
    return "\n".join(lines)


def _doctest_span(docstring_lineno: int, examples) -> tuple[int, int]:
    """File line range [start, end] covered by a docstring's doctest examples.

    ``example.lineno`` is the 0-based offset of the example's first source line
    within the (uncleaned) docstring, whose first line sits on the same physical
    line as the opening quote, so file line = ``docstring_lineno + offset``.
    """
    first = examples[0]
    last = examples[-1]
    start = docstring_lineno + first.lineno
    source_lines = last.source.count("\n")
    want_lines = last.want.count("\n")
    end = docstring_lineno + last.lineno + source_lines + want_lines - 1
    return start, end


_DEF_TYPES = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)


def _iter_nodes(source: str):
    """Yield (qualname_parts, node) for every class/function, nested included."""
    with warnings.catch_warnings():
        # Many real sources have non-raw strings with invalid escapes (regex,
        # LaTeX); they parse fine but ast.parse warns. Don't pollute build logs.
        warnings.simplefilter("ignore", SyntaxWarning)
        tree = ast.parse(source)

    def walk(node, prefix):
        for child in node.body:
            if isinstance(child, _DEF_TYPES):
                parts = prefix + [child.name]
                yield parts, child
                yield from walk(child, parts)

    yield from walk(tree, [])


def extract_doctest_chunks(
    source: str,
    module_name: str,
    rel_path: str,
    blob_base: str,
) -> list[dict]:
    """Extract one chunk per docstring that contains a doctest.

    The chunk text is the full runnable doctest block (its interactions build on
    each other, so the block is the natural retrieval unit). Title is the
    qualified name (``module_name`` + class/function path); the citation URL is
    the GitHub blob URL with a line anchor spanning the doctest block.
    """
    parser = doctest.DocTestParser()
    base = blob_base.rstrip("/")
    chunks: list[dict] = []
    for parts, node in _iter_nodes(source):
        docstring = ast.get_docstring(node, clean=False)
        if not docstring:
            continue
        try:
            examples = parser.get_examples(docstring)
        except ValueError:
            # Malformed doctest directive (e.g. a bogus "# doctest: +OPTION").
            # Skip this one docstring rather than failing the whole build.
            continue
        if not examples:
            continue
        docstring_lineno = node.body[0].lineno
        start, end = _doctest_span(docstring_lineno, examples)
        chunks.append(
            {
                "text": _reconstruct_doctest(examples),
                "title": ".".join([module_name, *parts]),
                "url": f"{base}/{rel_path}#L{start}-L{end}",
                "source": "code",
                "type": "doctest",
            }
        )
    return chunks


def _is_public(parts) -> bool:
    """Public API = no name component starts with an underscore (private/dunder)."""
    return all(not p.startswith("_") for p in parts)


def _render_signature(node) -> str:
    """Render a class or def header (decorators, args/defaults, bases, returns)."""
    decorators = [f"@{ast.unparse(d)}" for d in node.decorator_list]
    if isinstance(node, ast.ClassDef):
        bases = [ast.unparse(b) for b in node.bases]
        bases += [
            f"{kw.arg}={ast.unparse(kw.value)}" if kw.arg else f"**{ast.unparse(kw.value)}"
            for kw in node.keywords
        ]
        base_str = f"({', '.join(bases)})" if bases else ""
        header = f"class {node.name}{base_str}:"
    else:
        kw = "async def" if isinstance(node, ast.AsyncFunctionDef) else "def"
        returns = f" -> {ast.unparse(node.returns)}" if node.returns is not None else ""
        header = f"{kw} {node.name}({ast.unparse(node.args)}){returns}:"
    return "\n".join([*decorators, header])


def _signature_span(node) -> tuple[int, int]:
    """File line range [start, end] for a def/class header, decorators included."""
    start = node.decorator_list[0].lineno if node.decorator_list else node.lineno
    if node.body:
        end = node.body[0].lineno - 1
        if end < start:
            end = node.lineno
    else:
        end = node.end_lineno or node.lineno
    return start, end


def _one_line_summary(node) -> str:
    """First line of the (cleaned) docstring, or '' when there is no docstring."""
    doc = ast.get_docstring(node, clean=True)
    if not doc or not doc.strip():
        return ""
    return doc.strip().splitlines()[0].strip()


def extract_api_chunks(
    source: str,
    module_name: str,
    rel_path: str,
    blob_base: str | None = None,
) -> list[dict]:
    """Extract one chunk per public class/function: signature + one-line summary.

    Restores the accurate signatures (defaults, decorators, bases, return type)
    that the autodoc HTML buries, titled with the qualified name and cited to the
    GitHub blob URL for the header lines.
    """
    blob_base = blob_base or config.QSDSAN_BLOB_BASE
    base = blob_base.rstrip("/")
    chunks: list[dict] = []
    for parts, node in _iter_nodes(source):
        if not _is_public(parts):
            continue
        summary = _one_line_summary(node)
        # A bare, undocumented function/method signature is low-signal noise that
        # dilutes retrieval; keep classes for their structural value regardless.
        if not summary and not isinstance(node, ast.ClassDef):
            continue
        text = _render_signature(node)
        if summary:
            text += "\n\n" + summary
        start, end = _signature_span(node)
        chunks.append(
            {
                "text": text,
                "title": ".".join([module_name, *parts]),
                "url": f"{base}/{rel_path}#L{start}-L{end}",
                "source": "code",
                "type": "api",
            }
        )
    return chunks


def _node_source_span(source: str, node) -> tuple[str, int, int]:
    """Return (dedented full source, start_line, end_line) for a def, decorators included."""
    lines = source.splitlines()
    start = node.decorator_list[0].lineno if node.decorator_list else node.lineno
    end = node.end_lineno or node.lineno
    snippet = "\n".join(lines[start - 1 : end])
    return textwrap.dedent(snippet), start, end


def _cap_snippet(text: str) -> str:
    """Bound a snippet to _MAX_SNIPPET_LINES, appending a truncation note."""
    lines = text.split("\n")
    if len(lines) <= _MAX_SNIPPET_LINES:
        return text
    return "\n".join(lines[:_MAX_SNIPPET_LINES]) + "\n" + _TRUNCATION_NOTE


def _snippet_chunks(source, module_name, rel_path, blob_base, predicate, type_):
    """Emit one ``source="example"`` chunk per function the predicate selects.

    The chunk text is the function's FULL source (a self-contained, known-correct
    usage snippet), unlike the API extractor which keeps only the signature.
    """
    base = blob_base.rstrip("/")
    chunks: list[dict] = []
    for parts, node in _iter_nodes(source):
        if not isinstance(node, _FUNC_TYPES) or not predicate(parts, node):
            continue
        text, start, end = _node_source_span(source, node)
        chunks.append(
            {
                "text": _cap_snippet(text),
                "title": ".".join([module_name, *parts]),
                "url": f"{base}/{rel_path}#L{start}-L{end}",
                "source": "example",
                "type": type_,
            }
        )
    return chunks


def extract_test_chunks(
    source: str,
    module_name: str,
    rel_path: str,
    blob_base: str | None = None,
) -> list[dict]:
    """Each ``test_*`` function as a self-contained, known-correct usage snippet."""
    blob_base = blob_base or config.QSDSAN_BLOB_BASE
    return _snippet_chunks(
        source,
        module_name,
        rel_path,
        blob_base,
        predicate=lambda parts, node: parts[-1].startswith("test"),
        type_="test",
    )


def extract_system_chunks(
    source: str,
    module_name: str,
    rel_path: str,
    blob_base: str | None = None,
) -> list[dict]:
    """Each module-level ``create_system`` function: system-assembly wiring only.

    Restricted to module level so a same-named method on an analysis-model class
    (not actual system wiring) is not mistaken for a system entry point.
    """
    blob_base = blob_base or config.EXPOSAN_BLOB_BASE
    return _snippet_chunks(
        source,
        module_name,
        rel_path,
        blob_base,
        predicate=lambda parts, node: len(parts) == 1 and node.name == "create_system",
        type_="system",
    )


def _module_name_and_rel_path(abs_path: str, pkg_root: str) -> tuple[str, str]:
    """Map a source file to its dotted module name and repo-relative path.

    ``pkg_root`` is the directory that CONTAINS the top package, so the path
    relative to it is both the import path (dots) and the blob path (slashes).
    ``__init__.py`` maps to its package, not a ``.__init__`` submodule.
    """
    rel_path = os.path.relpath(abs_path, pkg_root).replace(os.sep, "/")
    mod_path = rel_path[: -len(".py")]
    if mod_path.endswith("/__init__"):
        mod_path = mod_path[: -len("/__init__")]
    module_name = mod_path.replace("/", ".")
    return module_name, rel_path


def build_code_chunks_from_dir(
    pkg_dir: str,
    pkg_root: str,
    blob_base: str | None = None,
) -> list[dict]:
    """Walk a package directory and extract API + doctest chunks from each module.

    Filesystem only (reads source, no import), so it is safe to run over a
    fixture tree and never executes package side effects.
    """
    blob_base = blob_base or config.QSDSAN_BLOB_BASE
    chunks: list[dict] = []
    for root, _dirs, files in os.walk(pkg_dir):
        for fname in sorted(files):
            if not fname.endswith(".py"):
                continue
            abs_path = os.path.join(root, fname)
            module_name, rel_path = _module_name_and_rel_path(abs_path, pkg_root)
            try:
                with open(abs_path, encoding="utf-8") as fh:
                    source = fh.read()
                chunks.extend(
                    extract_api_chunks(source, module_name, rel_path, blob_base)
                )
                chunks.extend(
                    extract_doctest_chunks(source, module_name, rel_path, blob_base)
                )
            except (SyntaxError, UnicodeDecodeError):
                continue  # skip files Python can't parse/decode, don't fail build
    return chunks


def build_code_chunks(
    package: str = "qsdsan",
    blob_base: str | None = None,
) -> list[dict]:
    """Locate an installed package and extract API + doctest chunks from source.

    Uses ``find_spec`` to find the package directory WITHOUT importing it, so the
    package's import-time side effects never run during indexing.
    """
    spec = importlib.util.find_spec(package)
    if spec is None or not spec.submodule_search_locations:
        raise ModuleNotFoundError(f"package not found or not a package: {package!r}")
    pkg_dir = list(spec.submodule_search_locations)[0]
    pkg_root = os.path.dirname(pkg_dir)
    return build_code_chunks_from_dir(pkg_dir, pkg_root, blob_base=blob_base)


def build_test_chunks_from_dir(
    tests_dir: str,
    repo_root: str,
    blob_base: str | None = None,
) -> list[dict]:
    """Walk a ``tests/`` tree and emit one example chunk per ``test_*`` function."""
    blob_base = blob_base or config.QSDSAN_BLOB_BASE
    chunks: list[dict] = []
    for root, _dirs, files in os.walk(tests_dir):
        for fname in sorted(files):
            if not (fname.startswith("test_") and fname.endswith(".py")):
                continue
            abs_path = os.path.join(root, fname)
            module_name, rel_path = _module_name_and_rel_path(abs_path, repo_root)
            try:
                with open(abs_path, encoding="utf-8") as fh:
                    source = fh.read()
                chunks.extend(
                    extract_test_chunks(source, module_name, rel_path, blob_base)
                )
            except (SyntaxError, UnicodeDecodeError):
                continue
    return chunks


def build_exposan_chunks(
    package: str = "exposan",
    blob_base: str | None = None,
) -> list[dict]:
    """EXPOsan ``create_system`` wiring, when EXPOsan is installed.

    EXPOsan is optional: if it is not importable this returns ``[]`` rather than
    failing, so the indexer works whether or not EXPOsan is present.
    """
    blob_base = blob_base or config.EXPOSAN_BLOB_BASE
    base = blob_base.rstrip("/")
    try:
        spec = importlib.util.find_spec(package)
    except (ImportError, ValueError):
        return []
    if spec is None or not spec.submodule_search_locations:
        return []
    pkg_dir = list(spec.submodule_search_locations)[0]
    chunks: list[dict] = []
    for root, _dirs, files in os.walk(pkg_dir):
        for fname in sorted(files):
            if not fname.endswith(".py"):
                continue
            abs_path = os.path.join(root, fname)
            # Path/module relative to the EXPOsan package dir; base already ends
            # in ".../exposan", so the final blob URL is .../exposan/<sys>/...
            rel = os.path.relpath(abs_path, pkg_dir).replace(os.sep, "/")
            mod = rel[: -len(".py")]
            if mod.endswith("/__init__"):
                mod = mod[: -len("/__init__")]
            module_name = f"{package}.{mod.replace('/', '.')}"
            try:
                with open(abs_path, encoding="utf-8") as fh:
                    source = fh.read()
                chunks.extend(
                    extract_system_chunks(source, module_name, rel, blob_base=base)
                )
            except (SyntaxError, UnicodeDecodeError):
                continue
    return chunks


def build_example_chunks(
    package: str = "qsdsan",
    exposan_package: str = "exposan",
) -> list[dict]:
    """Aggregate example chunks: the QSDsan test suite plus EXPOsan systems.

    The test suite lives beside the package (repo root); EXPOsan is optional.
    """
    chunks: list[dict] = []
    spec = importlib.util.find_spec(package)
    if spec and spec.submodule_search_locations:
        repo_root = os.path.dirname(list(spec.submodule_search_locations)[0])
        tests_dir = os.path.join(repo_root, "tests")
        if os.path.isdir(tests_dir):
            chunks.extend(build_test_chunks_from_dir(tests_dir, repo_root))
    chunks.extend(build_exposan_chunks(exposan_package))
    return chunks


def _preview(package: str = "qsdsan") -> None:
    """Dev helper: walk a package and print/dump all adapter chunks for inspection.

    Usage: ``python code_adapter.py [package] [out.json]``. No embeddings or
    network; just the structured code + example chunks the indexer will merge.
    """
    import collections
    import json
    import sys

    out_path = sys.argv[2] if len(sys.argv) > 2 else None
    chunks = build_code_chunks(package) + build_example_chunks(package)
    files = {c["url"].split("#")[0] for c in chunks}
    by_type = dict(collections.Counter(c["type"] for c in chunks))
    print(f"{len(chunks)} chunks from {len(files)} files in {package!r}: {by_type}")
    for c in chunks[:5]:
        print(f"\n# ({c['source']}/{c['type']}) {c['title']}\n{c['url']}\n{c['text'][:200]}")
    if out_path:
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(chunks, fh, ensure_ascii=False, indent=2)
        print(f"\nWrote {len(chunks)} chunks to {out_path}")


if __name__ == "__main__":
    import sys

    _preview(sys.argv[1] if len(sys.argv) > 1 else "qsdsan")
