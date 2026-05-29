# AGENTS.md — FlexAIDdS Agent Instructions

> **Primary agent instruction file for this repo is `CLAUDE.md`.** This AGENTS.md exists for compatibility with Codex, Claude Code, and other tools that look for it, and to document the packaged Grok skill.

## Skill Location & Invocation

The canonical FlexAID / FlexAIDdS / FlexAID∆S skill lives at:

```
.grok/skills/flexaid-docking/SKILL.md
```

**User-facing trigger phrases / slash commands:**
- `/flexaid-docking`
- `/FlexAid docking`
- `/FlexAidDS`
- `FlexAIDdS`, `FlexAID∆S`, "FlexAID docking", "ensemble analysis", "thermodynamic ledger", "CF/contact-function scoring proxy"

When any of the above appear, the agent **must** load and follow `.grok/skills/flexaid-docking/SKILL.md`.

## Mandatory Startup Ritual (Every Session)

```bash
git status
find . -maxdepth 4 -iname '*skill*' -o -iname 'SKILL.md' -o -iname '*.xml' -o -iname 'AGENTS.md'
python3 .grok/skills/flexaid-docking/scripts/validate_skill.py
```

## Guardrails (Copied from Skill)

- Inspect repo state first; validate every claim against files/commits/tests/logs.
- **avoid unsafe git** operations; **never merge branches or rewrite history** without explicit user confirmation.
- Separate "CF/contact-function scoring proxy" language from real thermodynamic ledger claims. Do not overclaim true ∆G.
- **preserve current ranking** behavior unless explicitly told to change the thermodynamic integration.
- Thermodynamic/ensemble work only behind tests and feature flags.
- Produce **chunked implementation plans** only.
- Never change scientific docking, scoring, or ranking code without a packaging/test requirement for a tiny non-behavioral import/path fix.

## Validation Commands (Run Before "Done")

```bash
python3 .grok/skills/flexaid-docking/scripts/validate_skill.py
python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
```

## What Agents Must Not Do

- Modify LIB/ or python/flexaidds/ scientific kernels (statmech, Vcontacts, BindingMode, tENCoM, etc.) except for the smallest non-behavioral packaging/import fix required by a test.
- Claim slash commands are natively registered beyond what the host TUI actually supports.
- Invent content from inaccessible Grok share links (only the title "Grok Fixes FlexAID Skill XML" was visible; body was not fetched).
- Skip the validator or test runs.

## References

- Full skill + guardrails: [.grok/skills/flexaid-docking/SKILL.md](.grok/skills/flexaid-docking/SKILL.md)
- Terminology & distinctions: [.grok/skills/flexaid-docking/references/flexaid-docking-guidance.md](.grok/skills/flexaid-docking/references/flexaid-docking-guidance.md)
- Primary detailed instructions: [CLAUDE.md](CLAUDE.md)
- Project README skill section (added for this packaging fix)

This file + the skill were added as part of the FlexAIDdS skill packaging fix (branch: feature/docs-validation-boundary).
