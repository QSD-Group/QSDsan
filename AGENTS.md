# Agent Context

QSDsan is the reusable modeling engine for sanitation and resource-recovery systems. Put generic components, streams, units, processes, system behavior, TEA/LCA utilities, and public APIs here.

For cross-repo structure, dependency, import, or `sanunits` namespace work, use the architecture skill for the agent runtime in use:

```text
.codex/skills/qsdsan-exposan-architecture/SKILL.md
.claude/skills/qsdsan-exposan-architecture/SKILL.md
```

Keep EXPOsan out of QSDsan runtime dependencies. EXPOsan-dependent checks belong in explicit integration coverage.
