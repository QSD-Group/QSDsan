# Agent Context

QSDsan is the reusable modeling engine for sanitation and resource-recovery systems. Put generic components, streams, units, processes, system behavior, TEA/LCA utilities, and public APIs here.

For cross-repo structure, dependency, import, or `sanunits` namespace work, use the architecture skill for the agent runtime in use:

```text
.codex/skills/qsdsan-exposan-architecture/SKILL.md
.claude/skills/qsdsan-exposan-architecture/SKILL.md
```

Keep EXPOsan out of QSDsan runtime dependencies. EXPOsan-dependent checks belong in explicit integration coverage.

## Changelog & releases

- The changelog is `CHANGELOG.rst` at the repo root (surfaced in the docs via `docs/source/CHANGELOG.rst`, which just `.. include::`s it). EXPOsan has no changelog of its own, so record **major** EXPOsan changes here; minor ones such as typo fixes need no entry.
- Versions are static in each repo's `pyproject.toml` (`version = "X.Y.Z"`). EXPOsan also pins an exact `qsdsan==X.Y.Z`, which must be bumped to match for coordinated releases.
- For a release, add a `` `X.Y.Z`_ `` section at the top of `CHANGELOG.rst` with a matching footer link target (`.. _X.Y.Z: https://github.com/QSD-Group/QSDsan/releases/tag/vX.Y.Z`), then bump the `pyproject` version(s).

## Packaging data files

- Any non-`.py` file QSDsan reads at runtime (e.g., `units_of_measure.txt`, files under `data/`) **must** be declared in `[tool.setuptools.package-data]` in `pyproject.toml`. Editable installs read the source tree, so a missing entry still works locally — but the built wheel/sdist omits the file, which breaks `import qsdsan` for non-editable installs (and fails every EXPOsan test at collection, since importing qsdsan is the first thing they do).
- After adding or relocating such a file, verify it ships: build a wheel (`uv build --wheel`) and confirm the file is inside it. Local (editable) tests passing is **not** sufficient to catch this class of bug.
