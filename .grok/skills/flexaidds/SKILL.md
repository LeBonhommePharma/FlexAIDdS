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

**Source of truth:** `AGENTS.md` (repo root). This skill derives from `AGENTS.md` and defers to it when rules conflict. For Astex entropy benchmark launch/monitor/resume work, also read `.agents/skills/flexaidds-benchmarking/SKILL.md`.

**Repository hygiene:** Never commit `.env` / secret files. Never add machine-specific absolute paths (`/Users/...`) to committed skills or shared scripts — use repo-relative paths or `FLEXAIDDS_*` environment variables. Run `python3 scripts/check_repo_hygiene.py` before pushing skill changes.

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
python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
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

# Full campaign with reports
python3 -m flexaidds.dataset_runner --all --tier 2 --results-dir results/benchmarks_2026

# Distributed run (launch with mpirun)
mpirun -n 8 python -m flexaidds.dataset_runner --all --tier 2 --distributed

# Dry-run to validate pipeline without actual docking
python3 -m flexaidds.dataset_runner --dataset casf2016 --tier 1 --dry-run
```

**Important guardrails when using DatasetRunner through this skill:**
- Always run `ensure_docking_data.py` first (or the inspector) — missing matrices or definition files will cause silent or noisy failures.
- Use `--dry-run` liberally before committing large compute resources.
- Respect the distinction between CF/contact-function scoring proxy (used during search) and the full thermodynamic ledger (computed afterward).
- For any published benchmark results, **always** pass `--package` (or run the inspector with `--reproducibility`). The resulting `VALIDATION_SUMMARY.md` + manifest gives you complete, auditable provenance (binary + every data file hash + environment).

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

## M3 Pro iCloud Canonical Best-BindingMode Protocol (exact 298 K / 310 K full runs)

**This is the production, zero-friction path to the exact requested answer: the best BindingMode (lowest free_energy after full thermo ledger + entropy corrections) from FlexAIDdS for a target+ligand molecular recognition event at precise temperature.**

Only this path is supported for the 4 canonical full first runs (Astex Diverse self-docking + Astex Non-Native cross-docking at 298 K and 310 K).

### Mandatory ritual (every time)
```bash
git status
find . -maxdepth 4 -iname '*skill*' -o -iname 'SKILL.md' -o -iname 'AGENTS.md'
python3 .grok/skills/flexaidds/scripts/validate_skill.py   # must PASS
```

### 1. Close hogs + space (iCloud green-lit)
User-confirmed: only ~325 GB occupied out of 2 TB iCloud. No aggressive cleanup in results/ needed.
- Close/pause Codex + openclaw etc (frees CPU/GPU for GA + Metal kernels in Shannon/tENCoM etc).
- Local tmp only: rm stale /private/tmp/flexaidds* and /tmp/flexaidds* (protect any active benchmark's omp-build/campaign dirs).

### 2. Re-ensure (full, no --quick)
```bash
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py
# Heavy Non-Native will require+verify Lovell_LIB.dat rotobs.lst SYBYL_emat.dat etc.
```

### 3. Launch the 4 canonical (exact command; iCloud-only results + Metal pre-flight enforced)
All output exclusively under $FLEXAIDDS_RESULTS (iCloud /Mobile Documents/.../FlexAIDdS/results).

```bash
# 298 K
bash .grok/skills/flexaidds/scripts/launch_full_benchmark.sh astex_diverse 298 astex_diverse_298K
bash .grok/skills/flexaidds/scripts/launch_full_benchmark.sh astex_nonnative 298 astex_nonnative_298K

# 310 K
bash .grok/skills/flexaidds/scripts/launch_full_benchmark.sh astex_diverse 310 astex_diverse_310K
bash .grok/skills/flexaidds/scripts/launch_full_benchmark.sh astex_nonnative 310 astex_nonnative_310K
```

The launcher (inside):
- Sources ~/.flexaidds_env (iCloud BUILD + RESULTS + BINARY).
- export PATH=... ; validate_skill + ensure (full).
- Metal pre-flight (counts .metallib (7+), system_profiler Metal 4 on M3 Pro, otool links Metal.framework, prints the exact grep for logs).
- Heavy dataset guard + symlinks for Non-Native.
- Early run_status.json (status, exact temperature, pids, iCloud output_dir/binary, command).
- nohup + disown portable detach (macOS-safe, separate stderr.log).
- Prints the 4 OUT_DIRs (full-298K-...-TS etc on iCloud), tail cmds, post-verify cmds, and "safely log out".

Example created (all iCloud):
- full-298K-astex_diverse_298K-1780383769
- full-310K-astex_diverse_310K-1780383775
- full-298K-astex_nonnative_298K-1780383775
- full-310K-astex_nonnative_310K-1780383775

### 4. Analyze along the way + valid results only
Live (monitors or):
```bash
# The 4 binary.logs (iCloud)
tail -f $OUT_DIR/binary.log
tail -f $OUT_DIR/stderr.log
cat $OUT_DIR/run_status.json   # temperature exact, status, returncode when done
```

Use the skill helper (Metal + health + validity + extract):
```bash
python3 .grok/skills/flexaidds/scripts/summarize_campaign.py $OUT_DIR --verbose --extract-best-mode
# Strict heuristic: real RMSD > placeholder, modes/prepared signals >0, temp==requested, returncode 0 or running, has subdir.
# --extract-best-mode: scans for rank-1 / lowest free_energy BindingMode (thermo ledger) + prints REMARKs + pointers to the PDB/JSON.
# (Full power: python -c 'from flexaidds.results import load_results; r=load_results(str(p)); print(r.top_mode())' — sorts by free_energy.)
```

Grep for the science (in any of the 4 logs):
```bash
grep -iE 'metal|backend|dispatch|shannon|using metal|success rate|RMSD|Binding Mode|Prepared|GA:|temperature|free energy|entropy' $OUT_DIR/binary.log | tail -30
```

iCloud FS green-light (no fuckup):
```bash
python benchmarks/re-dock/icloud_fs_check.py --path $OUT_DIR
# All tests (json, churn, executable, nested, bit persist) must PASS.
```

Post-finish full verify (ritual + only keep valid):
```bash
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python benchmarks/re-dock/icloud_fs_check.py --path $OUT_DIR
python3 .grok/skills/flexaidds/scripts/summarize_campaign.py $OUT_DIR --extract-best-mode
# Inspect the surfaced best BindingMode: exact T, lowest F from full ledger (partition + vib + Shannon config), real RMSD <<2 where success, >0 modes/poses, Metal dispatch lines for entropy kernels, returncode 0.
# Quarantine or discard anything with 999 placeholders, 0 modes, temp drift, high bogus energies, or no Metal when build had metallibs.
```

### 5. The "exact requested answer"
Once a dir passes the above (✅ from summarize, fs green, best mode extract shows sane thermo at exact T), that is the best BindingMode for that target+ligand at that T from the full FlexAIDdS entropy-driven simulation.

The 4 canonicals (launched 2026-06-02) are the reference full runs.

Monitors (persistent filtered tails on the 4 iCloud binary.logs) stream events for metal/dispatch/success/RMSD/temp/prepared/binding/entropy as they progress.

All large results, logs, status, configs, outputs live exclusively on iCloud (325 GB / 2 TB per user; fs_check + launcher + summarize enforce the paths).

### Why this eliminates friction for the best BindingMode
- Single canonical launcher (PATH, ritual, ensure, Metal preflight, temp fidelity in config+status, iCloud-only, detachable, early status).
- Early diagnosis (run_status, monitors, summarize health).
- Validity gate before accepting "this is the answer" (no more 999-as-success, 0-mode silent, temp 300, no Metal, premature tables).
- Best surfaced explicitly (free_energy sort from full thermo ledger, not CF proxy).
- Reproducible + auditable (skill validate, fs_check, manifests, git/binary/data hashes via ensure/package).

See also the launched dirs' run_status.json + binary.log for live proof.

(Protocol added/expanded after the 4 canonical 298/310 launches succeeded with all iCloud/Metal/T pre-flights green.)


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
