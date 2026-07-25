---
name: flexaid-docking
description: >
  Thin alias for FlexAID / FlexAIDδS docking workflows. Canonical policy and
  scripts live in the repo skill `.grok/skills/flexaidds/` (`/flexaidds`).
  Use for the same triggers as flexaidds: molecular docking, redocking, ensemble
  analysis, thermodynamic ledger, DatasetRunner, CF/contact-function scoring proxy.
  Natural language: "FlexAid docking", "run docking", "redock", "binding mode analysis".
user_invocable: true
metadata:
  short-description: "Alias → /flexaidds (canonical docking skill)"
---

# flexaid-docking (thin alias)

**Do not use this directory as a second source of scientific policy.**

| Authority | Path |
|-----------|------|
| **Source of truth** | `AGENTS.md` (repo root) |
| **Canonical skill** | `.grok/skills/flexaidds/SKILL.md` (`/flexaidds`) |
| **Benchmark ops** | `.agents/skills/flexaidds-benchmarking/SKILL.md` |
| **Methodology** | `METHODOLOGY.md` (cite by section; never fork numbers) |

## Agent contract

1. Immediately load and follow **`.grok/skills/flexaidds/SKILL.md`**.
2. Run scripts only from the canonical tree:
   ```bash
   python3 .grok/skills/flexaidds/scripts/validate_skill.py
   python3 .grok/skills/flexaidds/scripts/resolve_build.py --check
   python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py --check
   ```
3. Prefer `export FLEXAIDDS_REQUIRE_BUILD=1` in claim/agent sessions (hard-fail stale builds).
4. Obey the **deception-proof claim contract** in the flexaidds skill: no success language without real execution + on-disk `result.csv` / `RUN_RECEIPT` + RMSD≤2.0 Å and PoseBusters for modern claims.
5. **Local-first** live OUT (`$FLEXAIDDS_LOCAL_ROOT`); iCloud is thin durable mirror only.

## Why this alias exists

Hosts and users historically invoke `/flexaid-docking` or `~/.grok/skills/flexaid-docking`. A full fork of SKILL.md there went stale (wrong paths, missing resolve_build, pre-local-first iCloud rituals). This in-repo file is the **redirect**, so agents never invent science from a duplicate tree.

If a home install under `~/.grok/skills/flexaid-docking` still has full scripts/data, treat it as a **cache only** — re-ensure from the git checkout before any dock.
