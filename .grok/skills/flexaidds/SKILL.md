---
name: flexaidds
description: >
  Use this skill for FlexAID and FlexAIDδS docking workflows, including
  safe repo review, implementation planning, XML/package validation, and
  docking/thermodynamic-roadmap task decomposition. Includes auto-generation of
  high-end publication figures + 6s animated cover art (Grok Imagine / imagine-tools)
  for the best-scoring binding mode after Gate 6 success, with NRDD-quality aesthetics,
  thermodynamic equation, FlexAID∆S branding, reproducibility metadata overlay, and
  blue→red entropy heatmap.

  Natural language triggers include:
  - Any mention of FlexAID, FlexAIDδS, FlexAIDdS, "molecular docking", "perform docking",
    "run docking", "docking simulation", "redock", "redocking", "binding mode analysis",
    "thermodynamic analysis", "ensemble docking", "pose ranking", or "vibrational entropy".
  - Figure / viz triggers: "generate figure", "cover art", "animated cover", "imagine figure",
    "publication figure", "Nature Reviews cover", "NRDD figure", "best mode visualization",
    "entropy heatmap figure", "promotional animation".
  - Skill maintenance: "update the flexaidds skill", "update the docking skill",
    "refresh flexaid skill", "pull latest flexaidds", "update the skill".
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

**Source of truth:** `AGENTS.md` (repo root). This skill derives from `AGENTS.md` and defers to it when rules conflict. For Astex entropy benchmark launch/monitor/resume work, also read `.agents/skills/flexaidds-benchmarking/SKILL.md`. For DatasetRunner campaigns, use `.grok/skills/flexaidds-dataset-runner/SKILL.md` (in-repo thin launcher; home `~/.grok/skills/flexaidds-dataset-runner` is a mirror). `/flexaid-docking` is a **thin alias** (`.grok/skills/flexaid-docking/`) that redirects here — never treat a stale home fork as policy.

**Repository hygiene:** Never commit `.env` / secret files. Never add machine-specific absolute paths (`/Users/...`) to committed skills or shared scripts — use repo-relative paths or `FLEXAIDDS_*` environment variables. Run `python3 scripts/check_repo_hygiene.py` before pushing skill changes.

---

## Science / ops contract (post pilot8 — DO NOT MISLEAD)

These rules override any older “Softβ ON by default”, “entropy will lift S1”, or “SHARESCL 0.20” text elsewhere.

### CF proxy vs Softβ vs true ΔG

| Layer | Role | Default |
|-------|------|---------|
| **GA search** | Samples **CF/contact-function proxy** (Voronoi VCT) | Always CF search |
| **Engine TEMPER + CLUSTA FO** (arm **B**) | Density modes + soft-T **ACF** emission when TEMPER>0 | Arm B protocol (`TEMPER 21`, **not** kcal \(k_BT\)) |
| **DatasetRunner Softβ S1** | Optional election: \(\tilde G=\tilde H-T\tilde S\) over **already clustered** modes | **`FLEXAIDDS_SOFTBETA_ELECTION=0` (OFF)** |
| **StatMech / tENCoM / solvent ledger** | True-ish thermo only when full path validated | Not implied by Softβ |

- **Softβ ≠ FO@TEMPER21.** Pilot arm B was engine FO + TEMPER21, **not** DatasetRunner Softβ rescoring of CF heads. See `docs/implementation/softbeta_election_policy.md`.
- **Softβ cannot create ≤2 Å poses if BCR=0** (no near-native among emitted heads). Never re-rank BCR=0 pilots expecting S1 success.
- Prefer language: “CF soft-β ranking proxy \(\tilde G\)” — **never** “true binding free energy ΔG” unless full ledger is active and labeled.

### Classic three-engine red-pair (A / B0 / B)

| Arm | Engine | TEMPER | CLUSTA | Ranking story |
|-----|--------|--------|--------|---------------|
| **A** | FlexAID 2015-era pin | 0 | CF | CF red bar |
| **B0** | master FlexAID | 0 | CF | CF control |
| **B** | master FlexAID | **21** | **FO** (single literature MinPts) | Entropy arm = engine soft free energy on modes |
| **C0** | FlexAIDdS DatasetRunner | separate | separate | Out of band until FO dual-suffix election verified |

- **No dual-launch** of heavy GA on one Mac. Serial **A → B0 → B**. Local-first I/O (`scripts/use_local_first_benchmark_storage.sh`); sync iCloud later.
- **Matrix pin:** `MC_st0r5.2_6.dat` MD5 **`72d7c7396702331d96ff12d18f831796`**.
- **PSHARE:** production **`SHARESCL 10`**, **`SHAREPEK 5`**, **`SHAREALF 4`** (`scripts/generate_flexaid_inp.py`). **Never** ship `SHARESCL 0.20` (pilot typo; ~50× niche radius). Override only via `FLEXAIDDS_GA_SHARESCL` with receipt.
- **AMINO.def** is the live type file from DEPSPA unless `DEFTYP` is set. `AMINO26.def` on disk ≠ used.
- **Ligand emission:** `LIB/read_lig.cpp` must set **inclusive** `latm = atm_cnt` so last HETTYP atom is emitted (fix for missing 90017/90027). Rebuild binary after that fix; integrity gate catches regressions.

### Fail-closed prep / science gates (before claims)

```bash
# Prep (wired into generate_flexaid_inp by default)
python3 scripts/clean_target_apo.py          # strip HOH/metals for redock TARGET
python3 scripts/validate_ligand_integrity.py --work <work> --max-bond 3.0
# After FlexAID emits INI / poses:
python3 scripts/validate_ligand_integrity.py --work <work> --require-ini
python3 scripts/native_cf_oracle_gate.py --work <work> --results <out>/<pdb>
# Canary driver
bash scripts/run_pilot8_canary_gates.sh --arm B0 --pdb 1P62,1T40 ...
```

- **Native CF oracle:** FAIL (exit 1) when `CF_native > best_ga_cf + tol` → **ranking / Softβ / entropy claims forbidden**. Softβ does not repair a CF landscape that rejects the crystal.
- **3Dsig success metric (red bars):** **S_top10** = any of ranks 0..9 RMSD ≤ 2.0 Å; median over 10k bootstrap; 10 sims × 2e6 evals. Deck targets ~**0.66 / 0.69** (FlexAID / FlexAIDdS) on Astex Diverse N=85 — **not** pilot8 rates.
- Modern claim packages still need PoseBusters (+ tENCoM where required by benchmarking skill). RMSD-only is not full claim success.

### Deception-proof claim contract (normative — refuse without evidence)

**Refuse** all “docking success”, “recognition success”, “claim-ready”, “benchmark pass”, or numerical success-rate language unless **every** applicable gate below is satisfied from **on-disk artifacts** in the current session. Memory, chat history, and log fragments alone are never enough.

| Gate | Required evidence | Applies to |
|------|-------------------|------------|
| **Real execution** | Engine process ran; binary SHA256 matches `resolve_build.py --check` / pin | All docks & claims |
| **Runtime data** | `ensure_docking_data.py --check` OK (matrix + defs next to binary) | All docks |
| **Durable receipt** | `result.csv` present for DatasetRunner / modern packages; `RUN_RECEIPT` (or `.json`) with binary SHA256 + matrix MD5 for classic red-pair | Campaigns / claims |
| **Geometry + physics (modern)** | Rank-0 (or elected) **RMSD ≤ 2.0 Å** **and** PoseBusters pass on that same pose | Modern / DatasetRunner / PB claim tables |
| **STRICT / claim_ready** | Plus official `bust_cli`, tENCoM/Eigen on exact pose SHA-256, protocol eligibility, score–pose consistency — see `benchmarks/protocols/admission_metrics_contract.md` | STRICT packages |
| **Classic 3Dsig primary** | Report **S_top10** (any of ranks 0..9 RMSD ≤ 2.0 Å); S1/BCR diagnostic only | Classic A/B0/B arms |
| **Science gate** | If S_top10=0/N and BCR=0/N → **DOCKING COMPLETE — SCIENCE GATE FAIL**; no Softβ/ranking science claims | All |
| **Native CF oracle** | `native_cf_oracle_gate.py` must not fail when ranking/Softβ/entropy claims are made | Ranking claims |
| **Terminology** | CF/`best_score` = scoring proxy only; ledger F/H/−TS/Cv = ensemble estimates; never “true experimental ΔG” without full validated path | All prose |
| **DoF budget** | Claim runs: **fixed generations**, scale **population** via `FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1` (optional `FLEXAIDDS_BUDGET_SCALE`); never freeze base 1000×6000 as if it were the effective budget — read `[EVAL-BUDGET]` logs | Claim campaigns |
| **Docking semantics** | Self-docking vs cross-docking confirmed (orchestrator `native` vs `non_native`); never mix rates without labeling | Benchmarks |
| **Methodology** | Cite `METHODOLOGY.md` §N for parity / determinism / Astex-85 / ctest — **do not restate or fork numbers** in skills | Validation handoffs |

**Hard refuse phrases when evidence is missing:**
- “Success rate = …” without reading `result.csv` / success_pb / admission metrics for that exact OUT.
- “PoseBusters passed” without a PB receipt on the elected pose.
- “tENCoM validated” without Eigen/diff output tied to the pose SHA-256.
- “Local-first complete” when live GA wrote only to CloudDocs without local staging.

**Build pin discipline:** Before any real dock or claim, run `python3 .grok/skills/flexaidds/scripts/resolve_build.py --check`. For fail-closed CI/agent sessions set `FLEXAIDDS_REQUIRE_BUILD=1` so missing/stale builds are hard errors, not WARN.

**Storage:** Live GA/OUT/logs/binaries → **local** `$FLEXAIDDS_LOCAL_ROOT` (default `~/flexaidds_results`). iCloud is a **thin durable mirror** (`result.csv`, RUN_RECEIPT, thin OPS) only — see `docs/ICLOUD_BENCHMARK_STORAGE.md` and `AGENTS.md` § Benchmark storage. Never claim from iCloud-only live GA trees that hang FileProvider.

### Ops monitor scope

`scripts/run_benchmark_ops_monitor.sh` / `benchmark_ops_monitor.py` track **three_engine red-pair only** (`A|B0|B` / `3dsig_r10`). Do **not** treat C0_claim/C0_legacy as the live red-pair science path unless explicitly re-scoped.

**Primary invocations (documented aliases):**
- `/flexaidds`
- `/FlexAid docking`
- `/FlexAidDS`
- `/FlexAIDδS`, `/FlexAIDdS`
- `FlexAIDdS`, `FlexAID∆S`
- Natural language (strongly supported):
  - "update the flexaidds skill", "update the docking skill", "refresh the flexaid skill"
  - "dock this ligand", "perform molecular docking", "redock the co-crystallized ligand",
    "run FlexAIDδS on this target", "analyze the thermodynamic ledger", "binding mode prediction with entropy"
  - "run DatasetRunner", "benchmark on Astex", "run casf2016 benchmark", "distributed docking campaign", "dataset benchmarking"
  - "generate the cover figure", "create NRDD animation for the best mode", "imagine figure after docking", "add promotional cover art + 6s animation"

This skill activates for any task involving the FlexAID or FlexAIDδS molecular docking engine, its Python package `flexaidds`, **DatasetRunner** benchmarking campaigns, thermodynamics layer, or related packaging.

**Why leading researchers and pharma teams use this skill**
- **Pharma-grade reproducibility out of the box**: Every run (via DatasetRunner or manual) captures git SHA, binary SHA256, *complete* hashes of every critical runtime file (all matrices + 16 definition files + Lovell_LIB.dat + rotobs.lst + SYBYL_emat + scoring support), rich conda/pip + system environment, and produces a professional validation package on demand (`--package`).
- Beautiful one-pager `VALIDATION_SUMMARY.md` + `REPRODUCIBILITY_MANIFEST.json` — ready for papers, internal audits, collaboration, or regulatory packages.
- `inspect_definition_files --reproducibility` gives the same high-quality snapshot for one-off redocking and manual work.
- Production-grade DatasetRunner for systematic benchmarking on public and proprietary sets with professional reports.
- Self-contained critical data (no more "missing MC_*.dat or AMINO.def" surprises).
- Strong scientific guardrails and precise terminology (never confuses CF proxy with thermodynamic ledger).
- Extremely low-friction for both quick experiments and large distributed campaigns.

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
   python3 .grok/skills/flexaidds/scripts/validate_skill.py
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
- Assume slash commands beyond what the host TUI actually supports (document as user-facing trigger phrases + `/flexaidds` shorthand).

## Validation & Packaging

The skill itself is packaged under:
```
.grok/skills/flexaidds/
├── SKILL.md
├── scripts/
│   ├── validate_skill.py
│   ├── ensure_docking_data.py                  # unified runtime data (matrices + *.def files) + --source
│   ├── dataset_runner.py                       # high-quality wrapper for FlexAIDδS DatasetRunner (benchmarks, distributed runs, reports)
│   └── update_skill.py                         # built-in autoupdate for the skill + all sub-components
                                                #   (dry-run by default, --source, --yes, auto-validator)
├── data/
│   └── README.md                  # Documents MC_*.dat matrices + all AMINO*.def / NUCLEOTIDES*.def files
├── references/
│   └── flexaidds-guidance.md
└── assets/ (optional)
```

**Local validation commands (run these before any claim of "done"):**
```bash
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python3 .grok/skills/flexaidds/scripts/resolve_build.py --check
python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
```

**Production build resolution (autonomous, SHA-pinned):**
```bash
# After any C++ rebuild — refreshes ~/.flexaidds_env and ~/.flexaidds/active_build.json
python3 .grok/skills/flexaidds/scripts/resolve_build.py --sync-env

# Resume a campaign on the exact same engine binary:
export FLEXAIDDS_ENGINE_SHA256="<sha-from-prior-run-manifest>"
python3 .grok/skills/flexaidds/scripts/resolve_build.py --check
```

Before any real docking run, run the unified data ensure script:
```bash
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py
```

If you have a known-good FlexAIDδS installation elsewhere, use the deeply integrated `--source` flag:
```bash
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py \
    --source /path/to/your/working/flexaidds/install
```

You can also combine it with an explicit binary:
```bash
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py \
    --source /path/to/good/install \
    --binary /path/to/current/build/FlexAIDδS
```

### Keeping the Skill Up to Date (New in 2026-05)

The skill now includes a first-class, safe autoupdate tool:

```bash
# Always start here (completely safe)
python3 .grok/skills/flexaidds/scripts/update_skill.py --dry-run -v

# When you are ready (requires a full FlexAIDδS checkout as source)
python3 .grok/skills/flexaidds/scripts/update_skill.py --yes

# Using an explicit source (works great for portable copies too)
python3 .grok/skills/flexaidds/scripts/update_skill.py \
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

## Critical Runtime Data Management (Interaction Matrices + Definition Files)

The FlexAIDδS binary depends on two categories of runtime data files that are **not** part of the main source tree:

1. **Interaction matrices** (`MC_*.dat`) — used for the Voronoi contact-function (CF) scoring proxy during genetic algorithm search.
2. **Definition files** (`*.def`) — used for atom typing, covalent connectivity, and side-chain flexibility sampling.

### Definition Files (`*.def`) and Additional Runtime Data

The skill also bundles:
- `AMINO*.def` + `NUCLEOTIDES*.def` (atom typing, connectivity, and side-chain flexibility via `FLEDIH` entries)
- Supporting files (`Lovell_LIB.dat`, `rotobs.lst`, `SYBYL_emat.dat`, scoring matrices, etc.)

**Key practical points:**
- `AMINO.def` (2011 version) is the current standard. Legacy variants (AMINO8/12/26) use different atom type numbering and should be avoided with modern matrices.
- `FLEDIH` lines in `AMINO.def` directly control which side-chain torsions the GA will sample.
- All these files must live next to the binary at runtime.

See `data/README.md` for the full file list, format details, and per-residue FLEDIH mapping. Use `ensure_docking_data.py --info` or `inspect-definition-files` for diagnostics.

### Management in This Skill

This skill treats all these files as first-class managed assets:

- The `data/` directory ships with the complete runtime set (matrices + all `*.def` + Lovell_LIB, rotobs.lst, SYBYL_emat, scoring support files, etc.), making the skill fully self-contained.
- `scripts/ensure_docking_data.py` automatically discovers and places both matrices **and** definition files next to the binary (supports `--source`, `--dry-run`, `--check`, etc.).
- Use the ensure script **before every real docking task**.

See `data/README.md` for the complete file list and deeper format details (including full FLEDIH mappings per residue).

**Recommended before any docking task (now covers matrices + all definition files + extra runtime data):**
```bash
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py
```

The tools now automatically choose the right balance:
- In normal interactive use → rich diagnostics (what `--info` used to require).
- In CI or resource-constrained environments → automatic lightweight behavior.

You can still force modes with `--info` or `--quick` if needed. The `inspect-definition-files` helper follows the same smart logic.

## DatasetRunner — Distributed Benchmarking

The skill provides first-class support for the FlexAIDδS `DatasetRunner`, a powerful orchestrator for running systematic benchmarking campaigns across standard datasets.

**What it does:**
- Discovers and runs docking on curated datasets (Astex Diverse, CASF-2016, ITC-187, DUD-E subsets, psychopharmacology sets, etc.)
- Supports tiered execution (Tier 1 = fast sanity, Tier 2 = full comprehensive)
- Computes docking power, scoring power, and thermodynamic/entropy-related metrics
- Produces structured JSON + beautiful Markdown reports
- Supports local parallel, thread-pool, and MPI-distributed execution

**Typical usage via the skill:**

```bash
# Ensure all runtime data is present first (critical)
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py

# Run a single well-known dataset (Tier 1 for speed)
python3 -m flexaidds.dataset_runner --dataset astex_diverse --tier 1

# Full campaign with reports (prefer local results under $FLEXAIDDS_LOCAL_ROOT; thin-mirror later)
python3 -m flexaidds.dataset_runner --all --tier 2 --results-dir "${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}/benchmarks_$(date +%Y%m%d)"

# Distributed run (launch with mpirun)
mpirun -n 8 python -m flexaidds.dataset_runner --all --tier 2 --distributed

# Dry-run to validate pipeline without actual docking
python3 -m flexaidds.dataset_runner --dataset casf2016 --tier 1 --dry-run
```

**Important guardrails when using DatasetRunner through this skill:**
- Always run `ensure_docking_data.py` first (or the inspector) — missing matrices or definition files will cause silent or noisy failures.
- Use `--dry-run` liberally before committing large compute resources.
- Respect the distinction between CF/contact-function scoring proxy (used during search) and the full thermodynamic ledger (computed afterward).
- **Softβ S1 election defaults OFF** (`FLEXAIDDS_SOFTBETA_ELECTION=0` / `FLEXAIDDS_ELECTION_SHANNON_F=0`). Opt in only with explicit intent + log `[SOFTBETA-ELECT] Softβ S1 ON`. Never equate Softβ with true ΔG or with arm-B FO@TEMPER21.
- Do **not** enable Softβ / ranking experiments when `native_cf_oracle_gate.py` fails on canaries (CF rejects native).
- For any published benchmark results, **always** pass `--package` (or run the inspector with `--reproducibility`). The resulting `VALIDATION_SUMMARY.md` + manifest gives you complete, auditable provenance (binary + every data file hash + environment).
- Classic FlexAID A/B0/B red-pair uses `scripts/run_flexaid_arm_pilot8.sh` + `generate_flexaid_inp.py` (not DatasetRunner alone). Pin matrix MD5 + binary SHA256 in `RUN_RECEIPT.json`.

**Per-entry processing & Master Manager (new automation)**
The DatasetRunner now automatically saves and resumes *individual entries* (one target + structural state = one work item). A `EntryTaskManager` master coordinator allocates these fine-grained tasks to workers.

- Use `--resume` on long or expensive campaigns. It skips any target that already has a complete per-entry JSON result.
- Results layout: `results/<slug>/tierN/<target>_<state>.json` + `_entry_manifest.json` (with full per-entry wall time + cost in CPU-seconds)
- **Hybrid MPI**: Non-root ranks respect `--workers` locally (true MPI + threading).
- **Cost-aware scheduling**: On resume, previous costs are auto-loaded (with EMA history) to schedule cheaper entries first.
- All of the above appears in the final Markdown reports and in the skill's reproducibility validation package.
- The manager controls the worker pool size (`--workers`) and makes crash recovery + resource tracking first-class.

**CI Validation**: A dedicated GitHub Actions job (`.github/workflows/ci.yml`) now runs on every PR/push:
```bash
python3 .grok/skills/flexaidds/scripts/dataset_runner.py --dataset astex_diverse --tier 1 --dry-run --resume --package
```
It verifies the full reproducibility package + per-entry artifacts are produced correctly.

See `examples/small_real_benchmark_1stp.sh` for a minimal real-world-style example using a single complex.

**Reproducibility & Audit Packages (new in 2026-05)**
```bash
# Recommended for anything you intend to share or publish
python3 .grok/skills/flexaidds/scripts/dataset_runner.py \
    --all --tier 2 --package
# Note: prefer local OUT ($FLEXAIDDS_LOCAL_ROOT); thin-mirror result.csv/RUN_RECEIPT to iCloud after success.
# Legacy iCloud-only active trees are deprecated (FileProvider hang risk).

# For manual redocking or one-off work, capture a snapshot at inspection time
python3 .grok/skills/flexaidds/scripts/inspect_definition_files.py --reproducibility
```
The generated package contains:
- `REPRODUCIBILITY_MANIFEST.json` (machine-readable, full hashes + conda/pip capture)
- `VALIDATION_SUMMARY.md` (beautiful one-pager with tables, instructions, precise terminology, and regulatory notes)
- Your results/ directory

This is the general, reusable solution that works for DatasetRunner campaigns, redock_from_pdb workflows, and future tooling.

See the full CLI and library interface via:
```bash
python -m flexaidds.dataset_runner --help
```

Detailed dataset configurations live in `python/flexaidds/dataset_runner/datasets/`.

## Canonical full BindingMode protocol (298 K / 310 K) — **local-first**

**Goal:** best BindingMode (lowest ensemble free_energy after full thermo ledger + entropy corrections) for target+ligand molecular recognition at a stated temperature.

**Storage policy (overrides any older “iCloud-only live results” wording):**
| Layer | Where |
|-------|--------|
| Live GA / OUT / logs / binaries | **Local** `$FLEXAIDDS_LOCAL_ROOT` (default `~/flexaidds_results`) via `scripts/ensure_local_first_layout.sh` + `scripts/claim_local_staging_paths.sh` |
| Durable thin mirror | iCloud `$FLEXAIDDS_ICLOUD` / `$FLEXAIDDS_RESULTS` — **only** `result.csv`, RUN_RECEIPT, thin OPS (sync after local success) |
| CloudDocs I/O | **Must** use `scripts/icloud_safe_io.py`; never raw `find`/`rglob` under Mobile Documents |

See `AGENTS.md` § Benchmark storage and `docs/ICLOUD_BENCHMARK_STORAGE.md`. Agents that write live GA traffic only to CloudDocs are **out of contract**.

### Mandatory ritual (every time)
```bash
git status
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python3 .grok/skills/flexaidds/scripts/resolve_build.py --check
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py --check
# Prefer hard-fail in agent sessions:
# export FLEXAIDDS_REQUIRE_BUILD=1
```

### 1. Local layout + free resources
```bash
bash scripts/ensure_local_first_layout.sh
source scripts/use_local_first_benchmark_storage.sh 2>/dev/null || true
# Close competing heavy agents if needed; clear only stale /tmp/flexaidds* (not active OUT).
```

### 2. Re-ensure runtime data (full, no --quick)
```bash
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py
```

### 3. Launch the 4 canonical campaigns (local OUT; Metal pre-flight still required)
Confirm self-docking (`astex_diverse` / native) vs cross-docking (`astex_nonnative` / non_native) **before** launch. Softβ remains **OFF** unless explicitly opted in.

```bash
# 298 K
bash .grok/skills/flexaidds/scripts/launch_full_benchmark.sh astex_diverse 298 astex_diverse_298K
bash .grok/skills/flexaidds/scripts/launch_full_benchmark.sh astex_nonnative 298 astex_nonnative_298K

# 310 K
bash .grok/skills/flexaidds/scripts/launch_full_benchmark.sh astex_diverse 310 astex_diverse_310K
bash .grok/skills/flexaidds/scripts/launch_full_benchmark.sh astex_nonnative 310 astex_nonnative_310K
```

Launcher expectations:
- Sources `~/.flexaidds_env` when present; pin engine with `resolve_build.py --sync-env` after rebuilds.
- `validate_skill` + `ensure_docking_data` + Metal pre-flight (`.metallib`, Metal framework link).
- Early `run_status.json` (temperature, pids, **local** `output_dir`, binary path, command).
- Detach with nohup/disown; keep `binary.log` + `stderr.log` on **local** OUT.
- After success: thin-sync receipts to iCloud (`scripts/sync_claim_local_to_icloud.sh` / campaign archive helpers).

### 4. Analyze only valid local results
```bash
tail -f "$OUT_DIR/binary.log"
cat "$OUT_DIR/run_status.json"
python3 .grok/skills/flexaidds/scripts/summarize_campaign.py "$OUT_DIR" --verbose --extract-best-mode
# Accept only: real RMSD (not 999 placeholders), modes/poses > 0, temp == requested,
# returncode 0 (or still running), free_energy from ledger — not CF proxy alone.
```

Post-finish verify:
```bash
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python3 .grok/skills/flexaidds/scripts/summarize_campaign.py "$OUT_DIR" --extract-best-mode
# Require deception-proof gates above before any success language.
# Quarantine 999 placeholders, 0 modes, temp drift, or missing result.csv / RUN_RECEIPT.
```

### 5. The exact requested answer
Only after summarize + receipt gates pass: report the best BindingMode at the requested T (ensemble free_energy sort from the full ledger). Label docking mode (self vs cross) and never mix CF scores into “ΔG” language.

Historical note: older text described exclusive iCloud live trees; that path is **deprecated** because CloudDocs FileProvider hangs. Local-first + thin mirror is mandatory.


## Workflow for Typical Tasks

1. Discovery (git status + find + validator) — mandatory.
2. Read relevant source (never edit LIB/ or python/flexaidds/ scientific kernels without tests).
3. If implementation requested: produce chunked plan with per-chunk test commands.
4. Validate claims with `git diff`, build, and test runs — never skip.
5. Update this skill or its validator if packaging or guardrails evolve.
   Use the built-in updater: `scripts/update_skill.py --dry-run` then `--yes`.
6. Commit only after validator + tests pass (see README for commit rules).

### Convenience Shortcuts (`bin/` directory)

For ergonomics, the skill provides executable shell wrappers in `bin/` (never symlinks into `scripts/` — editing a symlink would corrupt the underlying Python tools):

```bash
.grok/skills/flexaidds/bin/ensure-docking-data
.grok/skills/flexaidds/bin/validate-skill
.grok/skills/flexaidds/bin/copy-docking-data
.grok/skills/flexaidds/bin/update-skill          # built-in autoupdate (dry-run by default)
.grok/skills/flexaidds/bin/dataset-runner        # DatasetRunner campaigns with safety + diagnostics integration
```

**These are pure symlinks.** Running them executes the exact same code as the real scripts. They change nothing about behavior or verification requirements.

**Important:** These shortcuts are for convenience only. They never replace running the actual FlexAIDδS binary, the full validator, or any scientific analysis. No scientific claim is ever valid without executing the real code.

## Quickstart for Actual Docking + Thermodynamics

For users who want to run real FlexAIDδS jobs (not just review code), start here:

→ **[QUICKSTART.md](QUICKSTART.md)** — End-to-end guide for preparing inputs, running docking, and computing the thermodynamic ledger.

## Publication-Quality Figure & Animation Generation (Imagine Integration + Gate 6)

The skill now includes first-class support for automatically producing a **static high-end cover figure + 6-second promotional animation** focused on the **best-scoring binding mode** (top-ranked by the ensemble free energy / thermodynamic ledger) after a successful docking.

**This is a deliberate post-processing / packaging feature** (purely additive; never affects scoring, ranking, poses, or any scientific computation). It is **gated on Gate 6 success** (the F/S cross-check from the thermodynamic audit / provenance.gate_results already produced by full ledger runs).

### When the gate fires (skill agent execution contract)
After any docking workflow reaches successful completion **and** Gate 6 (F/S cross-check) has passed (detected via audit JSON, reproducibility.json, or explicit user confirmation), the agent **MUST**:

1. Ensure the run produced a results directory (with PDBs + REMARKs or the JSON repro).
2. Run the preparation step (writes prompts + metadata derived from real values):
   ```bash
   # Preferred (via the run wrapper that the skill uses)
   VISUALIZE=1 bash scripts/run_flexaidds.sh ... --visualize
   # or
   python -c '
   from pathlib import Path
   from flexaidds.figures import prepare_publication_figures
   prepare_publication_figures(Path("results/my_run"), visualize=True, require_gate6=True)
   '
   ```
3. **Then invoke the imagine tools** available to the host (image_gen / imagine_text_to_image for the cover, video_gen / imagine_image_to_video for the 6s animation, image_edit for refinement/overlay polish if text is soft). Use the dynamically generated `prompt_cover.txt` and `prompt_animation.txt` (they already contain the exact real ΔG/ΔH/-TΔS numbers, ligand/target, git sha, gate status, etc.).
4. Save the materialized assets inside `<results_dir>/figures/` with canonical names:
   - `cover_best_mode.png` (or .jpg)
   - `animation_6s.mp4`
5. Confirm the assets contain the required elements (the prompts are engineered to request them):
   - Sleek bottom banner: `/flexaids-docking • FlexAID∆S`
   - Thermodynamic equation `ΔG=ΔH−TΔS` with the actual run values calligraphed.
   - Reproducibility metadata overlay (gate6:PASS, short git, date, run id).
   - Cyan/teal accents + deep navy gradients matching thebonhomme.com + LeBonhommePharma/FlexAIDdS identity.
   - Entropy heatmap (blue→red), induced-fit side chains, PyMOL-style publication base + promotional styling.
   - "Proudly suitable for the cover of Nature Reviews Drug Discovery" + high-end X scientific post aesthetic (cinematic, clean, SwitchCraft-inspired elegant MD viz quality).
6. The prompts are plain-text and AI-tool compatible (Grok, ChatGPT, Claude, etc.).

### Quick usage example (the one requested)
```bash
FLEXAIDDS_SOURCE=/path/to/FlexAIDdS \
SKIP_REBUILD=1 \
bash run_flexaidds.sh 1stp biotin.mol2 --temperature 298.15 -o results/test_run --visualize
```
This produces `results/test_run/figures/` containing the prompts + metadata (and later the rendered cover + animation) alongside the usual reproducibility artifacts.

### Aesthetics & prompt contract (redesigned)
The prompts are built in `python/flexaidds/figures.py` from real docking output. They enforce the NRDD-cover + reference-video aesthetic (deep navy #0a0e14 gradients, #22D3EE teal/cyan, gold for ΔG, terra for entropy, hybrid clean scientific rendering with subtle entropy wash, exact banner/equation/footer baked in, JetBrains Mono / thebonhomme.com typography for all labels). 

**PLIP integration for interactions**: Prompts emulate the clean, professional 3D interaction diagrams from PLIP (Protein-Ligand Interaction Profiler, https://github.com/pharmai/plip) — color-coded per its standard legend (blue for H-bonds, grey dashed for hydrophobics, etc.), with emphasis on the *most favourable contacts and those contributing most to the CF/Voronoi score*. If you have PLIP installed, run `plip -f <best_pose.pdb> -p -y` in the results dir before/after the skill step; the generated PNG/.pse makes an excellent base image for image_to_image or manual refinement (the prepare step will auto-detect and note a `base_plip_interactions.png` when possible). This gives pixel-accurate, publication-grade interaction viz baked into the promotional cover/anim.

See the module for the canonical TEMPLATE_COVER / TEMPLATE_ANIMATION.

**Guardrail**: Figure generation is post-hoc promotional only. Use precise language in any sharing ("best-scoring binding mode by the ensemble-derived thermodynamic ledger", "visualization generated from run outputs").

All new visualization work lives behind the existing "chunked plans + tests + validator" discipline.

## References (updated)

See [references/flexaidds-guidance.md](references/flexaidds-guidance.md) for preserved scientific terminology, scoring proxy vs. thermodynamic ledger distinctions, and historical context from the FlexAIDδS implementation roadmap.

This skill exists to keep all FlexAID / FlexAIDδS work safe, reproducible, and correctly scoped between scoring proxies and real statistical mechanics.

## Agent Instruction Maintenance

When workflow rules, build commands, or constraints change:

1. Update `AGENTS.md` first.
2. Propagate the delta into this file, `CLAUDE.md`, `.agents/skills/flexaidds-benchmarking/SKILL.md`, and `docs/custom-instructions/` (Claude, Codex/Cursor, Grok Build, ChatGPT).
3. Run `python3 scripts/check_repo_hygiene.py` and `python3 .grok/skills/flexaidds/scripts/validate_skill.py`.

Platform-specific packs live under `docs/custom-instructions/` — see the table in `AGENTS.md` → "Agent Instruction Files".
