---
name: flexaid-docking
description: >
  Use this skill for FlexAID and FlexAIDδS docking workflows, including
  safe repo review, implementation planning, XML/package validation, and
  docking/thermodynamic-roadmap task decomposition.

  Natural language triggers include:
  - Any mention of FlexAID, FlexAIDδS, FlexAIDdS, "molecular docking", "perform docking",
    "run docking", "docking simulation", "redock", "redocking", "binding mode analysis",
    "thermodynamic analysis", "ensemble docking", "pose ranking", or "vibrational entropy".
  - Skill maintenance: "update the flexaid-docking skill", "update the docking skill",
    "refresh flexaid skill", "pull latest flexaid-docking", "update the skill".
  - Any request involving the FlexAIDδS binary, flexaidds Python package, tENCoM,
    StatMechEngine, BindingMode, or thermodynamic ledger work.

  When a docking-related request is detected, the skill should ask clarifying questions
  about organism/species, biological target (protein/RNA/DNA + chains), ligand source
  (PDB ID, MOL2, SMILES, SDF, residue name), intent (self-docking/redocking vs cross-docking
  vs screening), thermodynamic requirements, and any special constraints before proceeding.
user_invocable: true
metadata:
  short-description: "FlexAID / FlexAIDδS docking, validation, safe planning"
---

# FlexAID / FlexAIDδS Skill

**Primary invocations (documented aliases):**
- `/flexaid-docking`
- `/FlexAid docking`
- `/FlexAIDδS`, `/FlexAIDdS`
- Natural language (strongly supported):
  - "update the flexaid-docking skill", "update the docking skill", "refresh the flexaid skill"
  - "dock this ligand", "perform molecular docking", "redock the co-crystallized ligand",
    "run FlexAIDδS on this target", "analyze the thermodynamic ledger", "binding mode prediction with entropy"

This skill activates for any task involving the FlexAID or FlexAIDδS molecular docking engine, its Python package `flexaidds`, benchmarks, thermodynamics layer, or related packaging.

**Conversational behavior (important):**  
When activated by any docking-related natural language request, the skill MUST ask clarifying questions before taking action. Key dimensions to establish:
- Biological context (organism / species)
- Target macromolecule (protein, RNA, DNA; specific chain(s); PDB ID or local file)
- Ligand(s) (name, SMILES, MOL2/SDF/PDB residue, or "extract from the PDB co-crystal")
- Docking intent (self-docking / redocking of known ligand vs. cross-docking vs. virtual screening)
- Thermodynamic depth required (full ensemble free energy / partition function, tENCoM vibrational entropy, temperature, etc.)
- Special constraints (covalent attachment, modified residues, NMR multi-model, bio-unit .pdb1 preference, user-specified receptor/ligand chains)
- Input/output preferences (local paths vs. automatic RCSB download + splitting via redock_from_pdb.py)

Never guess these details. Ask focused, numbered questions and wait for the user to provide the missing information.

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
  - FlexAIDδS (entropy-augmented)
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
│   ├── validate_skill.py
│   ├── ensure_docking_data.py                  # unified data matrix management + --source
│   └── update_skill.py                         # built-in autoupdate for the skill + all sub-components
                                                #   (dry-run by default, --source, --yes, auto-validator)
├── data/
│   └── README.md                  # Documents the required MC_*.dat files
├── references/
│   └── flexaid-docking-guidance.md
└── assets/ (optional)
```

**Local validation commands (run these before any claim of "done"):**
```bash
python3 .grok/skills/flexaid-docking/scripts/validate_skill.py
python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
```

Before any real docking run, run the unified data ensure script:
```bash
python3 .grok/skills/flexaid-docking/scripts/ensure_docking_data.py
```

If you have a known-good FlexAIDδS installation elsewhere, use the deeply integrated `--source` flag:
```bash
python3 .grok/skills/flexaid-docking/scripts/ensure_docking_data.py \
    --source /path/to/your/working/flexaidds/install
```

You can also combine it with an explicit binary:
```bash
python3 .grok/skills/flexaid-docking/scripts/ensure_docking_data.py \
    --source /path/to/good/install \
    --binary /path/to/current/build/FlexAIDδS
```

### Keeping the Skill Up to Date (New in 2026-05)

The skill now includes a first-class, safe autoupdate tool:

```bash
# Always start here (completely safe)
python3 .grok/skills/flexaid-docking/scripts/update_skill.py --dry-run -v

# When you are ready (requires a full FlexAIDδS checkout as source)
python3 .grok/skills/flexaid-docking/scripts/update_skill.py --yes

# Using an explicit source (works great for portable copies too)
python3 .grok/skills/flexaid-docking/scripts/update_skill.py \
    --source ~/FlexAIDdS \
    --yes \
    --data          # optional: also refresh bundled matrices
```

The updater:
- Is **dry-run by default**
- Detects full checkouts automatically (or via `--source` / `FLEXAIDDS_ROOT`)
- Refreshes scripts, references, docs, bin/ shortcuts, and (optionally) data
- Always runs the validator at the end
- Never modifies anything without explicit `--yes`

See the script header and `--help` for all options.

The validator enforces:
- Valid SKILL.md YAML frontmatter (`name`, `description`)
- Zero malformed XML anywhere (well-formedness, single root element, escaped ampersands, UTF-8, no illegal nesting/IDs)
- No broken relative links in SKILL.md
- All required aliases and guardrail phrases present

## Critical Runtime Data Management (Interaction Matrices)

The FlexAIDδS binary depends on precomputed atom-type interaction matrices (`MC_*.dat` files) for the Voronoi CF/contact-function scoring proxy. These files are **not** part of the main source tree.

This skill treats them as first-class managed assets:

- The `data/` directory inside the skill ships with the required matrices, making the skill self-contained and portable.
- `scripts/ensure_docking_data.py` is the primary, production-grade tool. It supports both automatic discovery and explicit sourcing from another installation (`--source`).
- A convenience wrapper `copy_docking_data_from_install.py` exists for the common “I have a known-good install” workflow.
- Running the ensure script is considered part of the standard pre-docking workflow.

See `data/README.md` for the full rationale, file list, and maintenance guidance.

**Recommended before any docking task:**
```bash
python3 .grok/skills/flexaid-docking/scripts/ensure_docking_data.py
```

## References

See [references/flexaid-docking-guidance.md](references/flexaid-docking-guidance.md) for preserved scientific terminology, scoring proxy vs. thermodynamic ledger distinctions, and historical context from the FlexAIDδS implementation roadmap.

## Workflow for Typical Tasks

1. Discovery (git status + find + validator) — mandatory.
2. Read relevant source (never edit LIB/ or python/flexaidds/ scientific kernels without tests).
3. If implementation requested: produce chunked plan with per-chunk test commands.
4. Validate claims with `git diff`, build, and test runs — never skip.
5. Update this skill or its validator if packaging or guardrails evolve.
   Use the built-in updater: `scripts/update_skill.py --dry-run` then `--yes`.
6. Commit only after validator + tests pass (see README for commit rules).

### Convenience Shortcuts (`bin/` directory)

For ergonomics, the skill provides short commands in `bin/`:

```bash
.grok/skills/flexaid-docking/bin/ensure-docking-data
.grok/skills/flexaid-docking/bin/validate-skill
.grok/skills/flexaid-docking/bin/copy-docking-data
.grok/skills/flexaid-docking/bin/update-skill          # built-in autoupdate (dry-run by default)
```

**These are pure symlinks.** Running them executes the exact same code as the real scripts. They change nothing about behavior or verification requirements.

**Important:** These shortcuts are for convenience only. They never replace running the actual FlexAIDδS binary, the full validator, or any scientific analysis. No scientific claim is ever valid without executing the real code.

## Quickstart for Actual Docking + Thermodynamics

For users who want to run real FlexAIDδS jobs (not just review code), start here:

→ **[QUICKSTART.md](QUICKSTART.md)** — End-to-end guide for preparing inputs, running docking, and computing the thermodynamic ledger.

This skill exists to keep all FlexAID / FlexAIDδS work safe, reproducible, and correctly scoped between scoring proxies and real statistical mechanics.
