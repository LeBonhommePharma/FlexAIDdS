---
name: flexaid-docking
description: >
  Use this skill for FlexAID, FlexAIDdS, and FlexAID∆S docking workflows, including
  safe repo review, implementation planning, XML/package validation, and
  docking/thermodynamic-roadmap task decomposition. Triggered by: /flexaid-docking,
  /FlexAid docking, /FlexAidDS, FlexAIDdS, FlexAID∆S, "FlexAID docking", or any
  mention of ensemble analysis, thermodynamic ledger, or CF/contact-function scoring
  proxy work on the FlexAIDdS codebase.
user_invocable: true
metadata:
  short-description: "FlexAID / FlexAIDdS / FlexAID∆S docking, validation, safe planning"
---

# FlexAID / FlexAIDdS / FlexAID∆S Skill

**Primary invocations (documented aliases):**
- `/flexaid-docking`
- `/FlexAid docking`
- `/FlexAidDS`
- Direct phrases: `FlexAIDdS`, `FlexAID∆S`, `FlexAID docking`, `ensemble analysis`, `thermodynamic ledger`

This skill activates for any task involving the FlexAID or FlexAIDdS (FlexAID∆S) molecular docking engine, its Python package `flexaidds`, benchmarks, thermodynamics layer, or related packaging.

## Mandatory First Actions (ALWAYS)

1. **Always inspect repo state first** before any other action. Run exactly these discovery commands:
   ```bash
   git status
   find . -maxdepth 4 -iname '*skill*' -o -iname 'SKILL.md' -o -iname '*.xml' -o -iname 'AGENTS.md'
   ```
2. **validate claims against files**, commits, tests, and logs — never trust memory or prior summaries.
3. Then run the project-specific skill validator:
   ```bash
   python3 .grok/skills/flexaid-docking/scripts/validate_skill.py
   ```
3. Inspect repo structure with `list_dir`, `read_file` on README.md, CLAUDE.md, docs/, python/flexaidds/, LIB/ key headers only as needed. Never assume layout.

## Core Guardrails (Non-Negotiable)

- **Inspect first, claim never**: Every factual statement about code, behavior, or history must be validated against actual files, `git log`, test output, or build logs in the current session. Do not trust prior conversation summaries.
- **Git safety** and **avoid unsafe git** operations: Never run `git push`, `git merge`, `git rebase`, `git reset --hard`, or any history-rewriting command without explicit user confirmation. **never merge branches or rewrite history** without explicit confirmation. Prefer read-only inspection.
- **No unsafe operations**: Do not force-push, delete branches, or edit `.git/` directly.
- **Separate scoring proxy from thermodynamics**:
  - The core engine uses **CF/contact-function scoring proxy** (VoronoiCF, Vcontacts) for pose ranking during GA search.
  - True thermodynamic quantities (Helmholtz F, entropy S, Cv, Boltzmann weights) come from the StatMechEngine / BindingMode layer on top of the ensemble.
  - Never claim "computed true binding free energy ΔG" unless the full partition function + vibrational corrections (tENCoM) + explicit solvent/conc terms are active and validated against experimental ITC or known benchmarks.
  - Use precise language: "CF/contact-function scoring proxy", "ensemble-derived free energy estimate", "thermodynamic ledger (F, H, -TS, Cv)".
- **Preserve ranking behavior** and **preserve current ranking**: Do not alter pose ranking, clustering, or final output order unless the user explicitly requests a change to the thermodynamic integration or WHAM procedure. Any such change requires new tests + feature flag.
- **Thermodynamic / ensemble work gated** and **thermodynamic/ensemble work only behind tests** and feature flags: All new ensemble analysis, ΔS contributions, or free-energy ledger features must be implemented behind tests (`ctest`, `pytest`) and optional feature flags. Never enable in default paths without passing validation.
- **Chunked plans only** and **produce chunked implementation plans**: When asked for implementation work (Codex, Claude Code, Grok Build, or human), always produce small, reviewable chunks with explicit test gates between chunks. Never deliver monolithic diffs.
- **Terminology preservation** (do not rename or dilute):
  - FlexAID (legacy)
  - FlexAIDdS / FlexAID∆S (entropy-augmented)
  - docking, ensemble analysis, thermodynamic ledger, CF/contact-function scoring proxy, Voronoi contact function.

## What This Skill Must NOT Do

- Change scientific formulas, docking ranking, or scoring behavior without explicit request + tests.
- Overclaim thermodynamic accuracy (e.g., "exact ΔG" vs. "ensemble estimate from partition function").
- Delete or overwrite existing skill content (preserve in references/ or git history).
- Invent or assume the content of inaccessible external links (e.g., Grok share pages) — only use what is locally verifiable or explicitly provided in the current prompt.
- Assume slash commands beyond what the host TUI actually supports (document as user-facing trigger phrases + `/flexaid-docking` shorthand).

## Validation & Packaging

The skill itself is packaged under:
```
.grok/skills/flexaid-docking/
├── SKILL.md
├── scripts/
│   └── validate_skill.py
├── references/
│   └── flexaid-docking-guidance.md
└── assets/ (optional)
```

**Local validation commands (run these before any claim of "done"):**
```bash
python3 .grok/skills/flexaid-docking/scripts/validate_skill.py
python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
```

The validator enforces:
- Valid SKILL.md YAML frontmatter (`name`, `description`)
- Zero malformed XML anywhere (well-formedness, single root element, escaped ampersands, UTF-8, no illegal nesting/IDs)
- No broken relative links in SKILL.md
- All required aliases and guardrail phrases present

## References

See [references/flexaid-docking-guidance.md](references/flexaid-docking-guidance.md) for preserved scientific terminology, scoring proxy vs. thermodynamic ledger distinctions, and historical context from the FlexAIDdS implementation roadmap.

## Workflow for Typical Tasks

1. Discovery (git status + find + validator) — mandatory.
2. Read relevant source (never edit LIB/ or python/flexaidds/ scientific kernels without tests).
3. If implementation requested: produce chunked plan with per-chunk test commands.
4. Validate claims with `git diff`, build, and test runs — never skip.
5. Update this skill or its validator if packaging or guardrails evolve.
6. Commit only after validator + tests pass (see README for commit rules).

This skill exists to keep all FlexAID / FlexAIDdS / FlexAID∆S work safe, reproducible, and correctly scoped between scoring proxies and real statistical mechanics.
