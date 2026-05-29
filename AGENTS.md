# AGENTS.md — FlexAIDdS

This file tells coding agents (Codex/OpenAI, Claude Code, Grok Build, and
similar) how to behave inside this repository.

The repo also ships a skill package — see `skills/flexaid-docking/` — that
encodes the same rules in a host-loadable form. This file is the
plain-English mirror for hosts that read `AGENTS.md` directly.

For deeper development guidance (build commands, test commands, repo
structure, code conventions) the authoritative source is [CLAUDE.md](CLAUDE.md).
This file does not duplicate that content; it focuses on **agent behavior**.

---

## Where the skill lives

```
skills/flexaid-docking/
├── SKILL.md                              # metadata + behavior contract
├── references/
│   ├── flexaid-docking-guidance.md       # long-form agent guidance
│   └── skill-manifest.xml                # XML mirror of the metadata
├── scripts/
│   └── validate_skill.py                 # packaging validator
└── assets/                               # reserved (.gitkeep placeholder)
```

`SKILL.md` is the canonical source of truth. `skill-manifest.xml` mirrors
the same metadata for hosts that consume XML. Both are validated by
`scripts/validate_skill.py`.

---

## How to invoke the skill

The skill responds to any of these trigger phrases:

- `/FlexAid docking`
- `/FlexAidDS`
- `FlexAIDdS`
- `FlexAID∆S`
- `FlexAID delta-S`

These are **trigger phrases**, not guaranteed native slash-command
registrations. The host environment decides which (if any) are wired up as
real `/`-commands. The canonical skill ID is `flexaid-docking`.

---

## How to validate the skill

From the repo root:

```bash
python3 skills/flexaid-docking/scripts/validate_skill.py
python3 -m pytest tests/test_flexaid_skill.py
```

The validator exits 0 on success and 1 on any failure. It is dependency-free
(Python 3 stdlib only) and safe to run in CI. The pytest test is a lightweight
wrapper that fails CI if the packaged skill, XML manifest, metadata, aliases,
or local references drift.

---

## What agents must not do in this repo

These are hard rules. They apply whether the skill is loaded or not.

1. **Do not modify docking scoring or ranking behavior** without an explicit
   user instruction that names the file. The scientific stack lives under
   `LIB/` — in particular `gaboom.cpp`, `Vcontacts.cpp`, `vcfunction.cpp`,
   `statmech.cpp`, `BindingMode.cpp`, `encom.cpp`, and everything under
   `LIB/ShannonThermoStack/`, `LIB/tENCoM/`, `LIB/LigandRingFlex/`.
2. **Do not relabel CF proxy scores as "binding free energy" or "∆G."**
   The Voronoi contact-function output is a *proxy*. Use the vocabulary
   contract in `skills/flexaid-docking/references/flexaid-docking-guidance.md`.
3. **Do not introduce GPL or AGPL dependencies.** See
   `THIRD_PARTY_LICENSES.md` and `docs/licensing/clean-room-policy.md`.
4. **Do not merge branches or rewrite git history** without explicit user
   confirmation in the same conversation.
5. **Do not run destructive git operations** (`push --force`,
   `reset --hard` on a dirty tree, `branch -D`, `clean -fd` outside
   `build*/` and `WRK/`) without explicit user confirmation.
6. **Do not delete skill content.** Old versions stay in git history or in
   `skills/flexaid-docking/references/`.
7. **Do not claim slash commands are natively registered** unless the host
   has explicitly enabled them.

---

## What agents should do

1. **Inspect repo state first** (`git status`, current branch) before any
   non-trivial edit.
2. **Verify claims against files.** If memory or the user says "X exists at
   path Y", read or grep before acting.
3. **Validate before pushing.** Build with CMake, run `ctest`, run
   `pytest` for the Python package, and the skill validator above.
4. **Plan in chunks** when asked for a roadmap. See the planning template in
   `references/flexaid-docking-guidance.md` section 5.
5. **Gate new science behind tests and flags.** Any new estimator or
   ranking term lands behind a CMake/Python flag that defaults OFF, with a
   regression test pinning the off-path behavior.
6. **Use the right vocabulary.** "CF score", "ensemble ∆F estimate",
   "vibrational entropy" — not vague "binding affinity" claims.

---

## Relationship to CLAUDE.md

`CLAUDE.md` is the development guide: build commands, test commands, file
layout, conventions. Read it before doing any work that touches code.

`AGENTS.md` (this file) is the behavior contract: what an agent must and
must not do, regardless of how it discovered the repository.

`SKILL.md` (under `skills/flexaid-docking/`) is the packaged, host-loadable
form of the same contract, with structured metadata and references.

All three should stay consistent. If they drift, `CLAUDE.md` wins on
development mechanics, `AGENTS.md` and `SKILL.md` win on agent behavior.
