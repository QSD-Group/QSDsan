"""Pure heading-based chunking for RST source and built HTML.

Both indexer adapters reuse these so chunk boundaries stay consistent.
"""
import re

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
