# FlexAID / FlexAIDdS / FlexAID∆S Docking Guidance (Reference)

This document preserves exact scientific terminology and guardrails for any agent or human working on the FlexAIDdS codebase via the `flexaid-docking` skill.

## Canonical Names (Never Alter)

- **FlexAID** — legacy single-structure docking engine (original 2015 codebase).
- **FlexAIDdS** / **FlexAID∆S** — the entropy-augmented fork that adds full canonical (and grand-canonical) ensemble thermodynamics on top of the GA + Voronoi CF pipeline.
- **∆S** — explicit entropy contributions (configurational via ShannonThermoStack, vibrational via tENCoM/ENCoM, and higher-order terms).

## Scoring vs. Thermodynamics — Critical Distinction

| Layer                        | What It Computes                          | Language You MUST Use                          | What It Is NOT                  |
|------------------------------|-------------------------------------------|------------------------------------------------|---------------------------------|
| Voronoi CF / Vcontacts       | Contact-function score (proxy)            | "CF/contact-function scoring proxy", "Voronoi CF score" | Not a free energy               |
| StatMechEngine + BindingMode | Partition function Z, F = -kT ln Z, H_eff, -T∆S, Cv, Boltzmann weights | "ensemble-derived free energy", "thermodynamic ledger (F, H, -TS, Cv)" | Not experimental ∆G_bind unless validated vs ITC |
| tENCoM / ENCoM               | Vibrational entropy differential ∆S_vib   | "vibrational entropy correction (tENCoM)"      | Not full anharmonic solvent entropy |
| Full TI / WHAM / Grand PF    | Concentration-dependent occupancy, selectivity | "grand-partition-function analysis"            | Requires explicit experimental anchors |

**Rule**: Never claim "we computed the true binding free energy" from a single docking run. The ledger is an *estimate* from the sampled ensemble under the chosen force field and solvent model.

## Historical Context (for Roadmap Tasks)

- Phase 1: StatMechEngine core (partition function, WHAM, TI, parallel tempering).
- Phase 2: Python bindings, `flexaidds` package, `load_results()`, `dock()` API, ENCoM wrapper.
- Phase 3: tENCoM backbone flexibility + vibrational ∆S integration.
- Current boundary (see docs/VALIDATED_CAPABILITIES.md): Core 1.0 supports the GA + CF + StatMech + tENCoM path on Linux/macOS/Windows. Many Swift/TS/PWA/Fleet layers remain experimental.

## Implementation Rules (When Producing Plans)

1. **Chunked only** — each chunk < 300 LOC changed, includes its own test command (`ctest --output-on-failure -R <name>`, `pytest python/tests/test_*.py -q -k <keyword>`).
2. **Feature flags** — any new thermodynamic term, ensemble method, or ranking change must be off-by-default until the validator + full benchmark suite passes.
3. **No behavior change to ranking** — the order of poses/modes returned by the legacy `./FlexAID` or `flexaidds.dock()` must remain identical unless the user writes "change the ranking to use full vibrational ledger".
4. **XML / packaging hygiene** — any new skill, prompt pack, or Codex/OpenAI agent definition must be SKILL.md + YAML frontmatter (never raw XML with unescaped `&`, `∆`, or multi-root documents). The validator script enforces this.
5. **Git discipline** — branch, commit, push only after validator + tests. Never force, never rewrite shared history without explicit `/confirm` from the human.

## Common Pitfalls to Flag

- Assuming a high Voronoi CF score == tight binder (ignores entropy).
- Reporting "∆G = X kcal/mol" from a 300 K single-temperature GA run without the full ledger + concentration term.
- Editing `LIB/statmech.cpp`, `LIB/BindingMode.cpp`, `LIB/Vcontacts.cpp`, or `python/flexaidds/thermodynamics.py` without adding or updating the corresponding GoogleTest / pytest.
- Treating the Grok share link body as authoritative when the fetch only returned the title "Grok Fixes FlexAID Skill XML" — always fall back to local files and the current prompt.
- Using a legacy `AMINO8/12/26.def` with modern `MC_*.dat` matrices (different atom type numbering → wrong typing and scoring).
- Forgetting that `FLEDIH` lines in `AMINO.def` are what actually enable side-chain sampling in the GA. Residues without FLEDIH entries (ALA, GLY, PRO) get no side-chain flexibility from this mechanism.

## Definition Files (AMINO*.def / NUCLEOTIDES*.def) — Practical Notes

These files (bundled in the skill's `data/`) are required alongside the matrices:

- `AMINO.def` (2011.12.08) is the current standard. It defines the 20 amino acids with:
  - `ATMTYP` lines (serial, type code, name, r/m flag, parents)
  - `CONECT` (covalent bonds)
  - `FLEDIH` (which dihedrals are rotatable for GA sampling)
- The `AMINO8/12/26` variants are legacy and use incompatible atom type numbers.
- `FLEDIH` entries directly determine the side-chain flexibility that will be explored. See the full per-residue counts via `inspect_definition_files.py` or the source `AMINO.def`.

Always run the ensure script (and optionally the inspector) before production docking jobs.

## References Back to Source

- `docs/thermodynamics.md` — full ledger equations and component breakdown.
- `docs/VALIDATED_CAPABILITIES.md` + `docs/EXPERIMENTAL_CAPABILITIES.md`
- `CLAUDE.md` (this repo's primary agent instruction file)
- `python/flexaidds/` and `LIB/` for implementation.

Keep this file in sync with any evolution of the thermodynamic boundary or skill packaging rules.
