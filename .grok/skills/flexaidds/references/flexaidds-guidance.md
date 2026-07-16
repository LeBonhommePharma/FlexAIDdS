# FlexAID / FlexAIDdS / FlexAID∆S Docking Guidance (Reference)

This document preserves exact scientific terminology and guardrails for any agent or human working on the FlexAIDdS codebase via the `flexaidds` skill.

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

### DatasetRunner / campaign CSV columns (compatibility names)

Live campaign CSVs keep historical column names. Agents **must** use the semantics below in prose, reports, and claims:

| CSV column        | Actual meaning                                                                 | Forbidden language                          |
|-------------------|---------------------------------------------------------------------------------|---------------------------------------------|
| `best_score`      | CF/contact-function scoring proxy of the elected pose (same concept as `elected_cf` / REMARK CF) | "free energy", "ΔG", "binding affinity"     |
| `predicted_dG`    | Ensemble free-energy *estimate* F when StatMech ledger is present; may fall back to CF | "experimental ΔG", "true binding free energy" without full Z+vib+solvent validation |
| `predicted_dH` / `predicted_TdS` | Configurational ledger proxies when available                         | Calorimetric ΔH / TΔS without calibration   |
| `elected_cf` / `cf_native` / `cf_best_cluster` | Explicit CF fields (correct names)                        | Equating CF to ΔG                           |

Prefer saying **cf_score** / **CF proxy** in new documentation while leaving CSV headers unchanged.

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

## Softβ / FO / TEMPER (election vs sampling)

| Name | What it is | Default / note |
|------|------------|----------------|
| Softβ \(\tilde G=\tilde H-T\tilde S\) | Soft free energy over **mode members** on **CF** (≡ ACF) | DatasetRunner S1: **OFF** (`FLEXAIDDS_SOFTBETA_ELECTION=0`) |
| Engine TEMPER + CLUSTA FO | Arm B: density clustering + soft-T emission order | **Not** the same as DatasetRunner Softβ flag |
| TEMPER value (e.g. 21) | Soft temperature β=1/T on CF a.u. | **Not** physical K or \(k_B T\) in kcal unless calibrated |

Softβ **reorders** modes; it does **not** sample. If BCR=0 (no head ≤2 Å), Softβ cannot invent S1 success. Policy: `docs/implementation/softbeta_election_policy.md`.

## Common Pitfalls to Flag

- Assuming a high Voronoi CF score == tight binder (ignores entropy).
- Reporting "∆G = X kcal/mol" from a 300 K single-temperature GA run without the full ledger + concentration term.
- Equating Softβ \(\tilde G\) with experimental ΔG or with FO@TEMPER21 pilot arm B.
- Re-ranking BCR=0 ensembles with Softβ and claiming docking success.
- Emitting `SHARESCL 0.20` in `ga.inp` (typo; production is **10**).
- Using a binary that drops the last LIG.inp HETTYP atom (`latm = atm_cnt-1` bug) — integrity gate must fail closed.
- Claiming ranking science when native CF oracle fails (crystal CF hundreds of units worse than best GA CF).
- Editing `LIB/statmech.cpp`, `LIB/BindingMode.cpp`, `LIB/Vcontacts.cpp`, or `python/flexaidds/thermodynamics.py` without adding or updating the corresponding GoogleTest / pytest.
- Treating the Grok share link body as authoritative when the fetch only returned the title "Grok Fixes FlexAID Skill XML" — always fall back to local files and the current prompt.
- Using a legacy `AMINO8/12/26.def` with modern `MC_*.dat` matrices (different atom type numbering → wrong typing and scoring). `AMINO26.def` on disk is unused unless `DEFTYP` is set.
- Forgetting that `FLEDIH` lines in `AMINO.def` are what actually enable side-chain sampling in the GA.
- Missing Lovell_LIB.dat, rotobs.lst, or SYBYL_emat.dat — these are part of the complete runtime data pack required by the binary.
- Baking `/Users/<name>/...` absolute paths into skills or scripts (forbidden; use `FLEXAIDDS_*` or repo-relative resolution).

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

## Reproducibility & Provenance Tooling (2026-05+)

For any work that will be shared, published, or audited, the skill now provides a general, reusable reproducibility layer:

- `dataset_runner.py --package` → full `REPRODUCIBILITY_MANIFEST.json` + professional one-pager `VALIDATION_SUMMARY.md` (includes every critical data file hash, binary SHA256, git state, rich environment capture).
- `inspect_definition_files.py --reproducibility` → compact ready-to-paste JSON block + hash table (perfect for redock reports and lab notebooks).

These tools are deliberately **not** tied to a single workflow. They work for DatasetRunner campaigns, `redock_from_pdb`, manual `dock()` calls, and future tooling. Use them consistently when precision and defensibility matter.
