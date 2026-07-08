# AGENTS.md — FlexAIDdS AI Coding Agents

**This is the single source of truth for all AI coding agents working in this repository (Claude, Codex, Grok, GPT, etc.).**

Every other instruction file — `CLAUDE.md`, `.grok/skills/flexaidds/SKILL.md`, the ChatGPT/custom-agent instructions under `docs/custom-instructions/` — must derive from or explicitly defer to this document. When rules conflict, **this file takes precedence**.

All agents are expected to internalize and strictly follow the rules below on every task.

---

## Core Workflow Rules (Non-Negotiable)

These rules exist because this is a complex, performance-critical scientific codebase. Sloppy execution creates expensive bugs and broken builds.

1. **Verify with actual execution before claiming success.**
   Never say "done", "fixed", "implemented", "should work", or "complete" until you have run the relevant build and/or test commands and shown clean, passing output.
   - C++ changes → full test build + `ctest --output-on-failure`.
   - Python changes → the affected pytest suite.
   - CMake or source-list changes → a fresh configure + build.

2. **Use a todo list for any task with 3 or more distinct actions.**
   Open the list at the start of the task. Keep **exactly one** item in progress at a time. Mark items completed immediately when finished — never batch. Re-read the list before ending any turn that still has pending or in-progress work.

3. **After any code change, commit and push immediately.**
   Do not batch multiple logical changes. Use conventional commit prefixes (`Fix:`, `Add:`, `Update:`, `Refactor:`, etc.). If a git operation appears to hang, kill stale git processes (`kill $(pgrep -f git)`) and retry; check `git config core.fsmonitor` if problems persist.

4. **Fresh builds after touching the build system or sources.**
   Never assume a target still builds after editing `CMakeLists.txt` or adding/removing `.cpp`/`.h` files under `LIB/`. Run a clean configure + build and confirm linking succeeds. Watch for disk-space issues before long builds.

5. **Zero test failures before any push.**
   Run the full relevant suite (`ctest --output-on-failure` for C++, `pytest` for Python) after changes that could affect behavior. Fix all failures in the same session. Never push with red tests.

6. **Prioritized lists are completed in full.**
   When given P0/P1/P2 (or similar), finish every item before stopping.

7. **Run the thing when asked.**
   If the user explicitly asks you to *run* a command, benchmark, or test, do it — do not spend 20+ tool calls exploring first.

---

## Scientific Guardrails (Non-Negotiable)

- **Inspect first, claim never.** Every factual statement about code, behavior, or history must be validated against actual files, `git log`, test output, or build logs in the current session. Do not trust prior conversation summaries.
- **Avoid unsafe git.** Never `git push --force`, `git merge`, `git rebase`, `git reset --hard`, delete branches, or rewrite history without explicit user confirmation. Prefer read-only inspection.
- **Separate the scoring proxy from thermodynamics.** The GA search ranks poses with the **CF/contact-function scoring proxy** (VoronoiCF, Vcontacts). True thermodynamic quantities (Helmholtz F, entropy S, Cv, Boltzmann weights) come from the StatMechEngine / BindingMode layer on top of the ensemble. Never claim "computed true binding free energy ΔG" unless the full partition function + vibrational corrections (tENCoM) + explicit solvent/concentration terms are active and validated. Use precise language: "CF/contact-function scoring proxy", "ensemble-derived free energy estimate", "thermodynamic ledger (F, H, -TS, Cv)".
- **Preserve current ranking.** Do not alter pose ranking, clustering, or output order unless the user explicitly requests a change to the thermodynamic integration or WHAM procedure. Any such change requires new tests + a feature flag.
- **Thermodynamic/ensemble work only behind tests** and feature flags. Never enable new ΔS/free-energy features in default paths without passing validation.
- **Produce chunked implementation plans** with explicit test gates between chunks. Never deliver monolithic diffs.
- **Strict licensing:** Apache-2.0 only. Never introduce GPL/AGPL dependencies or use GPL code as inspiration (see `docs/licensing/clean-room-policy.md` and `THIRD_PARTY_LICENSES.md`).

---

## Repository Hygiene (Non-Negotiable)

- **Never commit secrets or local environment files.** Do not add `.env`, `.env.*`, `.envrc`, or any file containing API keys, tokens, or machine-local credentials. Use `.env.example` (placeholders only) when documenting required variables. `.gitignore` already excludes these — do not force-add them with `git add -f`.
- **Never commit machine-specific absolute paths in agent skills or shared scripts.** Committed automation under `.agents/`, `.grok/skills/`, `docs/custom-instructions/`, and new `scripts/` helpers must resolve paths from the repo root (`Path(__file__).resolve().parents[...]`, `git rev-parse --show-toplevel`) or from documented environment variables (`FLEXAIDDS_ROOT`, `FLEXAIDDS_ICLOUD`, `FLEXAIDDS_POSEBUSTERS_BIN`, `FLEXAIDDS_TENCOM_BIN`). Do not bake in `/Users/<username>/...` paths.
- **Run the hygiene check before pushing skill or agent-instruction changes:**
  ```bash
  python3 scripts/check_repo_hygiene.py
  ```

---

## Essential Commands

**C++ (recommended starting point for most work)**
```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
ctest --test-dir build --output-on-failure
```

**Python package only** (many tasks do not need the full C++ build)
```bash
pip install -e ./python
pytest tests/ -q
```

**Python bindings smoke**
```bash
cmake -B build -DBUILD_PYTHON_BINDINGS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
```

**Skill / runtime-data maintenance**
```bash
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python3 .grok/skills/flexaidds/scripts/resolve_build.py --check    # reject stale builds; pin engine SHA256
python3 .grok/skills/flexaidds/scripts/resolve_build.py --sync-env  # refresh ~/.flexaidds_env after rebuilds
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py --check
python3 .grok/skills/flexaidds/scripts/update_skill.py --dry-run -v
```

**Benchmark workflows**
```bash
python3 .grok/skills/flexaidds/scripts/dataset_runner.py --dataset astex_diverse --tier 1 --dry-run --resume --package
bash .grok/skills/flexaidds/scripts/launch_full_benchmark.sh astex_diverse 298 astex_diverse_298K
bash scripts/run_benchmark_production.sh --dry-run
python3 scripts/validate_benchmark_results.py <results-dir>/summary.csv --manifest benchmarks/datasets/astex_diverse.yaml --shannon-log-dir <results-dir>/astex_diverse --out-dir <results-dir>/figures
```

See `CLAUDE.md` → "Build System" for the full table of CMake options and targets.

---

## The FlexAIDdS Skill

The canonical agent skill lives at `.grok/skills/flexaidds/SKILL.md` (single directory, single slash command).

**User-facing trigger phrases / slash commands:** `/flexaidds`, `/FlexAid docking`, `/FlexAidDS`, `FlexAIDdS`, `FlexAID∆S`, "FlexAID docking", "ensemble analysis", "thermodynamic ledger", "CF/contact-function scoring proxy".

**Validation commands (run before claiming "done" on any skill change):**
```bash
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py --check
```

The skill is self-contained: `scripts/` (validator, data-ensure, dataset runner, updater), `data/` (matrices + `*.def` runtime files), `references/flexaidds-guidance.md` (terminology contract), `examples/`, and `bin/` convenience symlinks.

### Astex Entropy / Benchmark Orchestration Skill

For Astex entropy/FlexAIDdS benchmark launch, monitor, resume, restart, audit, or handoff work, use `.agents/skills/flexaidds-benchmarking/SKILL.md` before touching live runs. It is the shared contract for Grok Build, Claude Code, Claude Science, Claude Cowork Dispatch, Codex, and other agents.

PoseBusters and tENCoM/Eigen are mandatory validators for benchmark claims. A pose is successful only when `RMSD <= 2.0 A` and PoseBusters passes; never report RMSD-only success as benchmark success.

---

## Project Snapshot

FlexAIDdS is an entropy-driven molecular docking engine that combines genetic algorithms with statistical mechanics (partition functions, free energy, vibrational + configurational entropy).

**Primary languages:** C++26 (core engine in `LIB/`), Python (`python/flexaidds/` + bindings), Objective-C++ (Metal GPU), optional CUDA.

**Key architectural layers** (in execution order):
1. Genetic Algorithm (`LIB/gaboom.cpp`)
2. Voronoi Contact Scoring (`LIB/Vcontacts.cpp`)
3. Statistical Mechanics engine (`LIB/statmech.cpp`)
4. BindingMode clustering + thermodynamic integration (`LIB/BindingMode.cpp`)
5. Vibrational entropy (ENCoM + tENCoM)
6. Shannon configurational entropy with hardware acceleration (`LIB/ShannonThermoStack/`)
7. Cavity detection (`LIB/CavityDetect/`)

Full details live in `CLAUDE.md`.

---

## Agent Instruction Files (Per Platform)

| Agent / platform | Primary file(s) |
|------------------|-----------------|
| **All agents** | `AGENTS.md` (source of truth) |
| **Claude** | `CLAUDE.md` + `docs/custom-instructions/claude-instructions.md` |
| **Codex / Cursor** | `docs/custom-instructions/codex-cursor-instructions.md` (or root `.cursorrules` pointing here) |
| **Grok Build CLI** | `.grok/skills/flexaidds/SKILL.md` + `docs/custom-instructions/grok-build-cli-instructions.md` |
| **ChatGPT** | `docs/custom-instructions/chatgpt-instructions.md` |
| **Benchmark orchestration** (all agents) | `.agents/skills/flexaidds-benchmarking/SKILL.md` |

## Maintaining These Instructions (Production Discipline)

- `AGENTS.md` is the **authoritative source** for workflow rules and high-level constraints.
- When you change a core rule or command, update `AGENTS.md` first, then propagate the delta into the derived files:
  - `CLAUDE.md` (Claude's richer, more detailed reference)
  - `.grok/skills/flexaidds/SKILL.md` (Grok `/flexaidds` skill)
  - `.agents/skills/flexaidds-benchmarking/SKILL.md` (shared benchmark contract)
  - `docs/custom-instructions/` (Claude, Codex/Cursor, Grok, ChatGPT platform packs)
- Every derived file must state that `AGENTS.md` is the source of truth and defer to it, rather than restating rules that can drift. Prefer a short pointer over a duplicated paragraph.
- Keep the "Core Workflow Rules" wording identical across files if you must repeat it, to prevent drift.
- After any agent-instruction edit, run `python3 scripts/check_repo_hygiene.py` and `python3 .grok/skills/flexaidds/scripts/validate_skill.py`.

---

## When in Doubt

1. Read `AGENTS.md` first (this file).
2. Then read the tool-specific file for the agent you are currently using.
3. When the task involves running builds or tests, actually run them and show the output.

These rules exist to protect the quality and velocity of a real scientific codebase. Follow them without exception.
