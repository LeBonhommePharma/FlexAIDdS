---
name: flexaidds-benchmarking
description: Coordinate FlexAIDdS and Astex entropy benchmarks across agents. Use when any agent, including Grok Build, Claude Code, Claude Science, Claude Cowork Dispatch, or Codex, needs to launch, monitor, resume, restart, audit, or explain FlexAIDdS benchmark campaigns, especially benchmarks/astex_entropy, Astex Diverse native/non-native runs, iCloud result archival, PoseBusters validation, tENCoM/Eigen thermodynamic rescoring, no-duplicate launch hygiene, and benchmark handoff/status reporting.
---

# FlexAIDdS Benchmarking

## Operating Rule

Treat this skill as the shared benchmark contract. Any agent can work the benchmark, but it must leave artifacts and status clear enough that another agent can continue without asking LP which AI handled the prior step.

Before touching a live run:

1. Read `/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/AGENTS.md`.
2. Read `/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/CLAUDE_BENCHMARK_HANDOFF.md` if the task could overlap older v90-v94 campaign work.
3. Run the status script:

```bash
python3 /Users/lp.more/Documents/PhD/Programs/FlexAIDdS/.agents/skills/flexaidds-benchmarking/scripts/astex_entropy_status.py
```

Do not launch duplicate benchmark work if a matching live process or current orchestrator run already exists. Monitor or resume the existing namespace instead.

## Canonical Paths

- Benchmark workspace: `/Users/lp.more/Documents/PhD/Programs/FlexAIDdS`
- Live source/data checkout: `/Users/lp.more/Projects/FlexAIDdS`
- Astex entropy module: `/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/benchmarks/astex_entropy`
- Long-term iCloud results: `/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy`
- Orchestrator summaries: `<icloud>/orchestrator_runs/<run_id>/orchestrator_summary.{json,md}`

## Required Validators

PoseBusters and tENCoM/Eigen are mandatory. No PoseBusters, no tENCoM/Eigen, no benchmark claim.

Required paths by default:

```text
/Users/lp.more/Projects/FlexAIDdS/.venv-posebusters/bin/bust
/Users/lp.more/Projects/FlexAIDdS/build_lto/tencom_entropy_diff
```

The orchestrator preflights both through `benchmarks/astex_entropy/orchestrate.py`. Do not disable this check. If either validator is missing, stop and report the exact missing path plus the rebuild/install step needed.

## One-command Workflow

Use the orchestrator for normal work. It prepares the manifest, runs selected pose generators, rescoring each generated pose set with Shannon collapse, tENCoM/Eigen, thermodynamic `G_bind`, RMSD, and PoseBusters.

Native one-target wiring smoke:

```bash
cd /Users/lp.more/Documents/PhD/Programs/FlexAIDdS
.venv-astex-entropy/bin/python -m benchmarks.astex_entropy orchestrate --mode native --tools flexaidds --max-targets 1 --skip-rescore
```

Native full head-to-head:

```bash
cd /Users/lp.more/Documents/PhD/Programs/FlexAIDdS
.venv-astex-entropy/bin/python -m benchmarks.astex_entropy orchestrate --mode native --tools flexaidds,vina,rdock,boltz
```

Non-native smoke:

```bash
cd /Users/lp.more/Documents/PhD/Programs/FlexAIDdS
.venv-astex-entropy/bin/python -m benchmarks.astex_entropy orchestrate --mode non_native --tools flexaidds --max-targets 3 --download-missing --skip-rescore
```

Native plus non-native:

```bash
cd /Users/lp.more/Documents/PhD/Programs/FlexAIDdS
.venv-astex-entropy/bin/python -m benchmarks.astex_entropy orchestrate --mode all --tools flexaidds,vina,rdock,boltz --download-missing
```

Use `--dry-run` when validating command generation only. Use `--continue-on-error` only when LP explicitly wants partial progress over fail-fast rigor.

## Monitoring And Resume

Use `scripts/astex_entropy_status.py` first. It reports active benchmark processes, latest orchestrator summary, pose CSV counts, rescored CSV counts, and the generated FlexAIDdS command.

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

A pose is successful only when both are true:

```text
RMSD <= 2.0 A
PoseBusters passes
```

RMSD alone is not success. PoseBusters failure means failure even if RMSD is good.

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
