"""Pure heading-based chunking for RST source and built HTML.

Both indexer adapters reuse these so chunk boundaries stay consistent.
"""
import re

from bs4 import BeautifulSoup

# RST section adornment characters per the docutils convention.
_ADORNMENT = set("=-~^\"'`#*+.:_")


def slugify(title: str) -> str:
    """GitHub/Sphinx-style anchor slug: lowercase, non-alphanumerics to hyphens."""
    slug = re.sub(r"[^a-z0-9]+", "-", title.strip().lower())
    return slug.strip("-")


def split_rst_by_heading(text: str) -> list[dict]:
    """Split RST into [{title, text, anchor}] by underlined (and overlined) headings.

    A heading is a non-blank title line whose following line is an adornment run
    of one repeated character at least as long as the title.
    """
    lines = text.splitlines()
    sections: list[dict] = []
    current = None
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        nxt = lines[i + 1] if i + 1 < n else ""
        is_underline = (
            line.strip()
            and line == line.lstrip()  # title must be at column 0 (not indented)
            and nxt.strip()
            and len(set(nxt.strip())) == 1
            and nxt.strip()[0] in _ADORNMENT
            and len(nxt.strip()) >= len(line.strip())
        )
        # Skip an over-line directly above the title (e.g. ==== / BSM1 / ====).
        prev_is_overline = (
            i > 0
            and lines[i - 1].strip()
            and len(set(lines[i - 1].strip())) == 1
            and lines[i - 1].strip()[0] in _ADORNMENT
        )
        if is_underline:
            title = line.strip()
            if current is not None:
                sections.append(current)
            current = {"title": title, "text": "", "anchor": slugify(title)}
            i += 2  # consume title + underline
            continue
        if prev_is_overline and not line.strip():
            i += 1
            continue
        if current is not None:
            current["text"] += line + "\n"
        i += 1
    if current is not None:
        sections.append(current)
    for s in sections:
        s["text"] = s["text"].strip()
    return sections


_HEADING_TAGS = ("h1", "h2", "h3", "h4")


def split_html_by_heading(html: str) -> list[dict]:
    """Split built HTML into [{title, text, anchor}] by section heading.

    Furo/Sphinx wraps each section in <section id="..."> with the heading as the
    first child and nested subsections as later siblings. For each heading the
    chunk text is the heading's following siblings up to (but excluding) the next
    nested <section>, so a parent heading does not absorb its child subsections.
    The headerlink pilcrow ("¶") is stripped from titles and bodies.
    """
    soup = BeautifulSoup(html, "html.parser")
    chunks: list[dict] = []
    for heading in soup.find_all(_HEADING_TAGS):
        # Drop the Furo headerlink (the "¶" pilcrow) before reading the title.
        for link in heading.select("a.headerlink"):
            link.decompose()
        title = heading.get_text(strip=True)
        if not title:
            continue
        # Anchor: id on the heading, else nearest ancestor section/div with an id.
        anchor = heading.get("id")
        if not anchor:
            container = heading.find_parent(
                lambda tag: tag.name in ("section", "div") and tag.get("id")
            )
            anchor = container.get("id") if container else slugify(title)
        # Body = siblings after the heading, stopping at the first nested section.
        parts = []
        for sib in heading.find_next_siblings():
            name = getattr(sib, "name", None)
            if name == "section":
                break
            if name in _HEADING_TAGS:
                break
            text = sib.get_text(" ", strip=True) if name else str(sib).strip()
            if text:
                parts.append(text)
        chunks.append(
            {"title": title, "text": " ".join(parts).strip(), "anchor": anchor}
        )
    return chunks
