# Consolidate QSDsan-website into the QSDsan docs

**Date:** 2026-05-21
**Status:** Design — pending review
**Author:** Yalin Li (with Claude)

## Problem

QSDsan currently maintains **two** public-facing properties:

1. **`QSDsan-website`** — a Jekyll (minimal-mistakes) marketing site at `qsdsan.com`
   (hosted via GitHub Pages; apex domain set through a `CNAME` file).
2. **`QSDsan/docs`** — Sphinx/Furo documentation at `qsdsan.readthedocs.io`
   (built by ReadTheDocs from `QSD-Group/QSDsan`).

The Jekyll site carries little unique content beyond the docs. Its only standalone
material is a curated **Research** page (6 papers), a **Learning** page (2 workshops),
an empty **App** "Coming Soon" placeholder, and a marketing splash home page. Everything
else (Docs, YouTube) is an outbound link, and "Docs" already points to ReadTheDocs.

Maintaining two sites is redundant. We want a single canonical site and a proper
publications page.

## Goal

- Make the **Sphinx docs the single canonical site**, served at **`qsdsan.com`** via a
  ReadTheDocs custom domain.
- Add a **Publications** page (hybrid: featured showcase + complete list).
- Migrate the still-useful content (Learning/workshops, a welcoming landing page,
  footer links) and retire the Jekyll site.

## Non-goals

- The **App** is *not* migrated into the docs — it will live separately at
  `https://qsdsan.app/`. The docs only link to it.
- No redesign of existing docs content (Tutorials, API, Systems, FAQ).
- No new publication-tracking automation (Scholar/Zotero sync) in this pass.

## Decisions (from brainstorming)

| Question | Decision |
|----------|----------|
| Publications page scope | **Hybrid** — featured papers (rich) + complete list (below) |
| Complete-list source of truth | **BibTeX file in repo** rendered by `sphinxcontrib-bibtex` |
| Which website content to keep | **Learning/workshops** + a **marketing-style landing** |
| App page | External link to `https://qsdsan.app/` (not migrated) |
| Newsletter + YouTube | **Footer links only** |
| DNS control | User controls the `qsdsan.com` registrar/DNS |
| Landing page approach | **Enhance the current RST landing** (hero + feature highlights + cards via `sphinx-design`) |
| Initial `publications.bib` | **User will provide the full list**; ship featured 6 + bib structure meanwhile |

## Architecture

The Sphinx docs (`QSDsan/docs/source`) become canonical. ReadTheDocs serves them at
`qsdsan.com`; `qsdsan.readthedocs.io` continues to work and redirects to the canonical
domain. The `QSDsan-website` repo is archived.

### Site structure

**Sidebar nav (toctree in `index.rst`):**

```
Home (landing)
Tutorials
API
Systems
Publications      <- new
Learning          <- migrated from website
FAQ
Contributing
Changelog
App  ↗            <- external link to https://qsdsan.app/
```

**Footer (Furo `footer_icons` in `html_theme_options`):**
GitHub · PyPI/Docs · YouTube · Email · Newsletter signup.

### Components

#### 1. Publications page (`docs/source/publications/index.rst` + `publications.bib`)

- **Featured section (top):** rich panels for selected papers — image + synopsis +
  citation + action buttons (Read Paper / Source Code) — built with `sphinx-design`
  cards/grids. Seeded from the current 6 on `research.md`:
  QSD, QSDsan, DMsan, Biogenic Refinery, NEWgenerator, HTL. The featured set is
  hand-curated and easy to edit later.
- **Complete list (below):** a `.. bibliography::` directive from `publications.bib`,
  sorted **newest-first**. New papers are added via PR to the `.bib` file.
- **Images:** copy the local `images/research/*.png` from the website into
  `docs/source/images/publications/`. Replace publisher hot-links (e.g. `pubs.acs.org`)
  with **local copies** — more robust and avoids hotlink/copyright fragility.

**Interface:** a contributor adds a paper by appending one BibTeX entry to
`publications.bib` (and, optionally, adding it to the featured set). No code changes.

#### 2. Learning page (`docs/source/learning/index.rst`)

Migrate the AEESP 2022 + EES workshop entries verbatim (materials + video links),
using the same `sphinx-design` panel style as Publications for visual consistency.

#### 3. Landing page (`docs/source/index.rst`)

Keep the existing content (What is QSDsan / Installation / Join the Community /
References) and **enhance** the top with:
- a hero image (reuse `images/index/splash.png`),
- the 3 feature highlights (integrated / under uncertainty / with visualization) as a
  `sphinx-design` grid,
- cards linking to Publications and Learning.

The existing `.. grid::` cards (Tutorials/API/Systems/FAQ) stay.

#### 4. Footer + external nav

- Furo `footer_icons` for GitHub, PyPI, YouTube, Email, Newsletter.
- External `App` link added as a toctree entry: ``App <https://qsdsan.app/>``.

#### 5. Dependencies / config

- Add `sphinxcontrib-bibtex` to the docs requirements (the `docs` extra used by
  `.readthedocs.yml`), and confirm `sphinx-design` is present.
- Register `sphinxcontrib.bibtex` in `conf.py` `extensions` and set
  `bibtex_bibfiles = ['publications/publications.bib']` (relative to the Sphinx source
  dir `docs/source/`, where `conf.py` lives).
- Configure the bibliography style for newest-first ordering (custom pybtex style or a
  sort option), confirmed during implementation.

### Domain / RTD migration

1. **RTD admin → Domains:** add `qsdsan.com`, set it **canonical**; RTD provisions a
   Let's Encrypt certificate.
2. **DNS (user-controlled):** `qsdsan.com` is an **apex** domain. Two viable setups,
   chosen at implementation time based on registrar capabilities:
   - **Apex canonical (preferred if supported):** ALIAS/ANAME record on `qsdsan.com`
     pointing to the target shown in RTD admin.
   - **www canonical + apex redirect:** CNAME `www.qsdsan.com` → RTD target, and
     redirect apex `qsdsan.com` → `www` at the registrar/DNS layer.
   The exact CNAME/ALIAS **target value is read from the RTD Domains panel** during
   implementation (not guessed here).
3. Remove the GitHub-Pages `CNAME` mechanism (it's superseded by the DNS records above).
4. Verify `qsdsan.readthedocs.io` still resolves and redirects to `qsdsan.com`.

### Retiring the Jekyll site

- Archive the `QSD-Group/QSDsan-website` GitHub repo (read-only) with a README note
  pointing to the new location.
- Disable its GitHub Pages deployment so it no longer competes for the domain.

## Data flow

```
Contributor edits docs/source/*  +  publications.bib
        │  (git push / PR to QSD-Group/QSDsan)
        ▼
ReadTheDocs build (Sphinx + furo + sphinx-design + sphinxcontrib-bibtex)
        │
        ▼
Served at qsdsan.com (canonical)  ←─ qsdsan.readthedocs.io redirects here
```

## Error handling / edge cases

- **Apex DNS:** if the registrar lacks ALIAS/ANAME, fall back to www-canonical + redirect
  (see above) so `qsdsan.com` never breaks for existing visitors.
- **SSL provisioning lag:** RTD's Let's Encrypt cert can take time after DNS propagates;
  expect a short window where HTTPS is not yet valid.
- **Old inbound links:** papers and the footer reference `qsdsan.readthedocs.io`; RTD's
  canonical-domain redirect preserves them.
- **Broken/hot-linked images:** localize publisher images to avoid future 404s.
- **bibtex build failures:** a malformed `.bib` entry fails the RTD build; PRs are
  validated by the normal docs build check.

## Testing / verification

- `docs` build passes locally (`make html`) and on RTD with the new extensions.
- Publications page renders: featured panels + complete bibliography, newest-first.
- Learning page renders with working material/video links.
- Landing page renders hero + feature grid + cards.
- Footer shows GitHub / PyPI / YouTube / Email / Newsletter; App link works.
- After DNS cutover: `https://qsdsan.com` serves the docs with valid SSL;
  `https://qsdsan.readthedocs.io` redirects to it.

## Open inputs (from user)

- **`publications.bib` content** — the complete list of QSDsan papers (BibTeX/Zotero
  export). Until provided, the page ships with the featured 6 and a clearly-marked
  placeholder section.

## Out of scope / future

- Automated publication sync from Google Scholar / Zotero.
- Any redesign of the App at `qsdsan.app`.
