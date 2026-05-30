# flexaidds — Grok Skill

Project-scoped skill that turns Grok into an expert FlexAIDdS developer.

## What It Gives You

- Deep knowledge of the C++26 core engine (`LIB/`)
- Python package (`python/flexaidds/`) and pybind11 bindings
- Full CMake build system + all important targets and options
- Testing discipline (`ctest` + pytest with `@requires_core`)
- Architecture (genetic algorithm → Voronoi scoring → StatMech → BindingMode clustering → vibrational + Shannon entropy → cavity detection)
- Strict workflow rules the project demands:
  - Always verify with actual build/test runs before claiming anything is done
  - Use `todo_write` for tasks with 3+ steps
  - Commit + push immediately after changes (conventional prefixes)
  - Fresh builds after touching CMake or adding sources

## How to Use

| Method                    | Command                          |
|---------------------------|----------------------------------|
| Direct slash command      | `/flexaidds`                     |
| Skills menu               | `/skills flexaidds`              |
| Automatic activation      | Just mention "FlexAIDdS", "tENCoM", "statmech", "BindingMode", etc. |

It combines beautifully with other skills (`/review`, `/implement`, `/check`, etc.).

## Why This Is Committed

This lives in `.grok/skills/` (not `~/.grok/skills/`) on purpose.  
Anyone who clones the repo gets the skill automatically. No extra setup required for teammates.

## Maintenance

- Edit `SKILL.md` when project structure, build commands, or core rules change.
- The authoritative workflow rules live in the repo root `AGENTS.md`. Keep this skill aligned with it.
- Keep the frontmatter description accurate (it drives auto-invocation).
- Full technical depth lives in [CLAUDE.md](../../CLAUDE.md).

Created: 2026-05-27
Updated: 2026-05-27 (added cross-tool AGENTS.md system)
