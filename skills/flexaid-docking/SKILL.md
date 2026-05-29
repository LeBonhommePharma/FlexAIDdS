---
name: flexaid-docking
description: Use this skill for FlexAID, FlexAIDdS, and FlexAID∆S docking workflows, including safe repo review, implementation planning, XML/package validation, and docking/thermodynamic-roadmap task decomposition.
---

# FlexAID Docking Skill

This skill guides an agent (Claude Code, Codex/OpenAI, Grok Build, or similar)
through tasks against the **FlexAIDdS** (aka **FlexAID∆S**) repository — an
entropy-aware molecular docking engine combining a genetic algorithm with a
Voronoi contact-function (CF) scoring proxy and a statistical-mechanics
"thermodynamic ledger" over pose ensembles.

The skill is intentionally narrow: it makes the agent *cautious* about the
docking scientific stack, and *thorough* about packaging, validation, and
planning.

---

## Invocation aliases

The following user-facing trigger phrases all map to this skill. They are
**trigger phrases**, not guaranteed native slash-command registrations — the
host environment (Claude Code, Codex, Grok) may or may not register them as
real `/`-commands. Treat any of these as a request to load this skill:

- `/FlexAid docking`
- `/FlexAidDS`
- `FlexAIDdS`
- `FlexAID∆S`
- `FlexAID delta-S` (plain-text fallback for the ∆ character)

If the host does register slash commands natively, the canonical command name
is `flexaid-docking`. Otherwise the agent should recognize the phrases above in
user prompts and switch into this skill's behavior.

---

## When to use this skill

Load this skill when the user is asking the agent to:

1. Inspect, plan changes to, or review code in the FlexAIDdS repository.
2. Run a docking workflow against a target/ligand pair using the existing
   `FlexAID` binary, Python bindings, or `python -m flexaidds`.
3. Generate an implementation roadmap (e.g. for Codex, Claude Code, or Grok
   Build) covering scoring, ensemble analysis, or the thermodynamic ledger.
4. Validate or repair skill/package XML, `SKILL.md`, or related agent
   metadata in this repo.
5. Triage a docking-result directory and produce a structured report.

If the user is asking for *new science* (a new scoring term, a different
thermodynamic estimator, a different ranking rule), this skill must surface
guardrails (see "What this skill must not do" below) before any code edit.

---

## What this skill must not do

1. **Do not modify scientific terminology.** Preserve all of:
   *FlexAID*, *FlexAIDdS*, *FlexAID∆S*, *docking*, *ensemble analysis*,
   *thermodynamic ledger*, *CF / contact-function scoring proxy*.
2. **Do not overclaim free-energy semantics.** The CF score is a *proxy*. Do
   not relabel it as a true thermodynamic binding free energy in
   documentation, output, or commit messages.
3. **Do not change docking ranking behavior** without an explicit user
   request that names the file and the change. Ranking changes must land
   behind a feature flag and tests.
4. **Do not merge branches or rewrite git history** without explicit user
   confirmation in the same conversation.
5. **Do not introduce GPL/AGPL dependencies.** See `THIRD_PARTY_LICENSES.md`
   and `docs/licensing/clean-room-policy.md`.
6. **Do not delete skill content** — preserve old versions in
   `references/` or rely on git history.
7. **Do not claim slash commands are natively registered** unless the host
   environment has explicit support; describe them as user-facing trigger
   phrases.

---

## Required behavior

When this skill is invoked, the agent should:

1. **Inspect repo state first.** Run `git status`, look at the current
   branch, and confirm there is no in-flight work that the user has forgotten
   about.
2. **Avoid destructive git ops.** No `--force`, no `reset --hard`, no
   `branch -D`, no force-push unless the user typed the exact command.
3. **Validate claims against the repo.** When the user (or memory) says "X
   exists at path Y", `Read` or `Grep` for it before acting.
4. **Separate scoring-proxy language from real thermodynamic claims.**
   Anywhere outputs say "free energy", "∆G", or "binding affinity", check the
   computation that produced the number. If it is a CF-proxy score, the
   surrounding text must say so.
5. **Preserve current ranking behavior** unless explicitly requested. A
   refactor that touches `LIB/BindingMode.cpp`, `LIB/Vcontacts.cpp`,
   `LIB/statmech.cpp`, or `LIB/gaboom.cpp` must include a before/after
   ranking comparison or a test that pins ranking.
6. **Implement thermodynamic / ensemble work behind tests and flags.** New
   estimators land in a guarded code path with a unit test and (when
   relevant) a regression test against a known pose set.
7. **Produce chunked plans on request.** When the user asks for a Codex /
   Claude Code / Grok Build roadmap, return a sequence of small,
   independently-mergeable steps with a clear test for each.
8. **Never auto-merge or rewrite history.** Branch creation is fine; merge
   and rebase require an explicit user instruction.

---

## Files in this skill

| Path | Purpose |
|------|---------|
| `SKILL.md` | This file. Metadata + behavior contract. |
| `references/flexaid-docking-guidance.md` | Long-form agent guidance (scoring vs thermodynamics, ledger, planning template). |
| `references/skill-manifest.xml` | Well-formed reference XML manifest of the skill (for hosts that prefer XML). |
| `scripts/validate_skill.py` | Standalone validator: frontmatter, alias text, XML well-formedness, file references. |
| `assets/` | Reserved for future static assets (icons, sample inputs). Contains a `.gitkeep` placeholder. |

---

## Validation

Run from the repo root:

```bash
python3 skills/flexaid-docking/scripts/validate_skill.py
python3 -m pytest tests/test_flexaid_skill.py
```

The validator exits non-zero on any of:

- Missing required skill/package files or directories.
- Missing `SKILL.md` frontmatter `name` or `description`.
- `name` field that is not the literal string `flexaid-docking`.
- `description` field that does not match the required OpenAI/Codex skill
  metadata.
- Missing alias documentation (`/FlexAid docking`, `/FlexAidDS`,
  `FlexAIDdS`, or `FlexAID∆S`).
- Malformed XML in `references/` or `assets/`.
- XML manifest metadata, aliases, IDs, or `<file path="..."/>` entries that
  drift from the canonical `SKILL.md` contract.
- A local Markdown or file-table reference that does not exist on disk.

---

## Relationship to FlexAIDdS code

This skill **wraps** the existing FlexAIDdS codebase; it does not replace any
of it. The actual docking engine lives in `LIB/` and is built via
`CMakeLists.txt`. The Python entry points are documented in
[CLAUDE.md](../../CLAUDE.md) under "Usage Modes":

- `./FlexAID config.inp ga.inp` — legacy native CLI
- `./flexaids dock receptor.pdb ligand.mol2` — modern native CLI
- `import flexaidds` / `python -m flexaidds` — Python API and CLI
- PyMOL plugin under `pymol_plugin/`

This skill is **packaging and guidance only**. It must not be confused with
the docking engine itself.
