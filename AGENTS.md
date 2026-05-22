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
