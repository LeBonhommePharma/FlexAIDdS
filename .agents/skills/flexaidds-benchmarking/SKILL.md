---
name: flexaidds-benchmarking
description: Coordinate FlexAIDdS and Astex entropy benchmarks across agents. Use when any agent, including Grok Build, Claude Code, Claude Science, Claude Cowork Dispatch, or Codex, needs to launch, monitor, resume, restart, audit, or explain FlexAIDdS benchmark campaigns, especially benchmarks/astex_entropy, Astex Diverse native/non-native runs, iCloud result archival, PoseBusters validation, tENCoM/Eigen thermodynamic rescoring, no-duplicate launch hygiene, and benchmark handoff/status reporting.
---

# FlexAIDdS Benchmarking

**Source of truth:** `AGENTS.md` (repo root). This skill is the shared benchmark contract for Grok Build, Claude Code, Claude Science, Claude Cowork Dispatch, Codex, and other agents.

## Operating Rule

Treat this skill as the shared benchmark contract. Any agent can work the benchmark, but it must leave artifacts and status clear enough that another agent can continue without asking LP which AI handled the prior step.

**Primary live goal:** three-engine **classic FlexAID** Astex Diverse red-pair — **A** (2015 CF) vs **B0** (master CF) vs **B** (master FO + TEMPER21). FlexAIDdS **C0/C** DatasetRunner is a **separate** packaging path until FO dual-suffix + Softβ policy are verified. Prefer `docs/implementation/3dsig_red_pair_protocol.md` over stale `LIVE_QUEUE.md` when they conflict.

### Hard bans (misgiving prevention)

| Do | Don't |
|----|--------|
| Serial **one** heavy arm at a time (A→B0→B) | Dual-launch A with B0/B or C0 on one Mac |
| **SHARESCL 10** / SHAREPEK 5 in `ga.inp` | Reintroduce **SHARESCL 0.20** (pilot typo) |
| Softβ DatasetRunner election **OFF** by default | Claim Softβ will fix **BCR=0** / re-rank pilot heads for ≤2 Å |
| Call arm B **FO@TEMPER21** | Call arm B “Softβ S1 rescoring of CF ensembles” |
| Fail-closed ligand integrity + native CF oracle | Claim ranking science when oracle fails (native CF ≫ decoy CF) |
| Local-first OUT/work for classic arms | Write live GA traffic only to hanging iCloud paths |
| Rebuild binary after `read_lig` latm fix | Run “fixed science” on old Mach-O that drops last HETTYP atom |

### Metrics (be precise)

- **3Dsig red-bar success:** **S_top10** (any of ranks 0..9 RMSD ≤ 2.0 Å), 10 sims × 2e6 evals, 10k bootstrap median. Deck Astex Diverse: FlexAID **~0.66**, FlexAIDdS **~0.69**.
- **S1** = rank-0 only; **BCR** = min RMSD over cluster heads (diagnostic sampling ceiling). Softβ/FO election cannot raise S1 if BCR>2.
- Modern packages: success for claims may also require PoseBusters — RMSD-only is not enough for PB claim tables.

Before touching a live run:

1. Read `AGENTS.md` (repo root).
2. Read `docs/implementation/3dsig_red_pair_protocol.md` + `docs/implementation/softbeta_election_policy.md`.
3. Resolve and pin the active C++ build (rejects stale `FLEXAIDDS_BUILD` paths):

```bash
python3 .grok/skills/flexaidds/scripts/resolve_build.py --check
```

4. Ops snapshot (three_engine red-pair only):

```bash
bash scripts/run_benchmark_ops_monitor.sh
# optional older path:
python3 .agents/skills/flexaidds-benchmarking/scripts/astex_entropy_status.py
```

5. Preflight canaries when prep changed:

```bash
bash scripts/run_pilot8_canary_gates.sh --arm B0 --pdb 1P62,1T40 \
  --work-root "$FLEXAID_WORK_ROOT" --results-root "<OUT>"
```

Do not launch duplicate benchmark work if a matching live process or current orchestrator run already exists. Monitor or resume the existing namespace instead.

## Canonical Paths (environment variables — never hardcode machine paths)

Resolve everything from the git checkout root or documented env vars. Never commit `/Users/<username>/...` paths in skills or shared scripts.

| Purpose | Variable | Typical default |
|---------|----------|-----------------|
| Repo / workspace root | `FLEXAIDDS_ROOT` | `git rev-parse --show-toplevel` |
| Astex entropy module | (under root) | `benchmarks/astex_entropy` (local workspace; gitignored) |
| Long-term iCloud results | `FLEXAIDDS_ICLOUD` | `~/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy` |
| PoseBusters `bust` binary | `FLEXAIDDS_POSEBUSTERS_BIN` | `$FLEXAIDDS_ROOT/.venv-posebusters/bin/bust` |
| tENCoM/Eigen diff binary | `FLEXAIDDS_TENCOM_BIN` | `$FLEXAIDDS_ROOT/build_lto/tencom_entropy_diff` |
| Astex entropy venv | `FLEXAIDDS_ASTEX_VENV` | `$FLEXAIDDS_ROOT/.venv-astex-entropy` |

Orchestrator summaries: `<FLEXAIDDS_ICLOUD>/orchestrator_runs/<run_id>/orchestrator_summary.{json,md}`

## Required Validators

PoseBusters and tENCoM/Eigen are mandatory. No PoseBusters, no tENCoM/Eigen, no benchmark claim.

The orchestrator preflights both through `benchmarks/astex_entropy/orchestrate.py`. Do not disable this check. If either validator is missing, stop and report the exact missing path plus the rebuild/install step needed.

## One-command Workflow

Use the orchestrator for normal work. It prepares the manifest, runs selected pose generators, rescoring each generated pose set with Shannon collapse, tENCoM/Eigen, thermodynamic `G_bind`, RMSD, and PoseBusters.

All examples assume `cd` to `$FLEXAIDDS_ROOT` (repo root) and `.venv-astex-entropy` exists or use `$FLEXAIDDS_ASTEX_VENV`.

Native one-target wiring smoke:

```bash
cd "${FLEXAIDDS_ROOT:-.}"
.venv-astex-entropy/bin/python -m benchmarks.astex_entropy orchestrate --mode native --tools flexaidds --max-targets 1 --skip-rescore
```

Native full head-to-head:

```bash
cd "${FLEXAIDDS_ROOT:-.}"
.venv-astex-entropy/bin/python -m benchmarks.astex_entropy orchestrate --mode native --tools flexaidds,vina,rdock,boltz
```

Non-native smoke:

```bash
cd "${FLEXAIDDS_ROOT:-.}"
.venv-astex-entropy/bin/python -m benchmarks.astex_entropy orchestrate --mode non_native --tools flexaidds --max-targets 3 --download-missing --skip-rescore
```

Native plus non-native:

```bash
cd "${FLEXAIDDS_ROOT:-.}"
.venv-astex-entropy/bin/python -m benchmarks.astex_entropy orchestrate --mode all --tools flexaidds,vina,rdock,boltz --download-missing
```

Use `--dry-run` when validating command generation only. Use `--continue-on-error` only when LP explicitly wants partial progress over fail-fast rigor.

## Monitoring And Resume

Use `scripts/astex_entropy_status.py` first (optionally `--work-dir "$FLEXAIDDS_ICLOUD"`). It reports active benchmark processes, latest orchestrator summary, pose CSV counts, rescored CSV counts, and the generated FlexAIDdS command.

Resume rules:

- If a matching process is active, monitor it. Do not relaunch.
- If no process is active and a summary reports an error, read the summary and the tool log before restarting.
- If rerunning a partially completed FlexAIDdS output namespace, do not use `--force` unless LP explicitly asks for a destructive restart. The runner skips completed targets by `result.csv`.
- If preserving smoke and full results matters, set a new `work_dir` in a copied config instead of overwriting the same iCloud namespace.
- For external tools, rerunning `run` rewrites the mode/tool pose CSV. Keep this in mind before relaunching a full tool batch.

Restart rules:

- Prefer a new orchestrator run with the same `work_dir` when resuming.
- Prefer a new timestamped `work_dir` when comparing methods, changing config, or rerunning after a bad setup.
- Record the exact command in the final answer or handoff.

## Success Definition

**3Dsig red-pair (classic A/B0/B) primary statistic:** S_top10 — any of top-10 ranked modes has RMSD ≤ 2.0 Å (deck contract). Report S1 and BCR as diagnostics; do not replace S_top10 with Softβ-only tables.

**Modern / DatasetRunner claim packages:** a pose is claim-successful only when **both** are true:

```text
RMSD <= 2.0 A
PoseBusters passes
```

RMSD alone is not full claim success. PoseBusters failure means failure even if RMSD is good.

**Science gate after docking:** if S_top10 = 0/N and BCR = 0/N, report **DOCKING COMPLETE — SCIENCE GATE FAIL**. Do not start Softβ election experiments; fix prep/emission/sampling first (`read_lig` latm, SHARESCL 10, clean apo, native CF oracle).

## Reporting Contract

Every closeout or handoff must include:

- Exact command launched or resumed.
- Modes and tools.
- iCloud work dir and orchestrator run ID.
- PoseBusters path and tENCoM/Eigen path verified by preflight.
- Active PIDs if any process remains running.
- Pose CSV and rescored CSV paths.
- Success rates from `success_pb`, not RMSD-only counts.
- Any missing validators, failed tools, interrupted jobs, or partial outputs.

Never report benchmark numbers from memory or logs alone when CSV artifacts exist. Read the CSV or summary JSON.

## Repository Hygiene

- Never commit `.env` files or API keys.
- Never add machine-specific absolute paths to this skill or its scripts. Use env vars above.
- After edits, run `python3 scripts/check_repo_hygiene.py` from repo root.