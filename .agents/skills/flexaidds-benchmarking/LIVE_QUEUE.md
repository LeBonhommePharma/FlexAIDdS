# Benchmark handoff — live queue & launch (all agents)

**Audience:** Codex, Claude, Grok, Cursor — anyone who might launch, monitor, or resume benchmarks.  
**Author:** Grok Build session · **Updated (UTC):** 2026-07-15T02:55Z  
**Source of truth for rules:** `AGENTS.md` · **Skill:** `.agents/skills/flexaidds-benchmarking/SKILL.md`  

This file is the **tracked** live queue handoff. A local mirror may also exist as root `CLAUDE_BENCHMARK_HANDOFF.md` (gitignored pattern `/*_HANDOFF*.md`).

---

## 0. Storage policy (non-negotiable)

**All new benchmark outputs go to iCloud Drive.**

```bash
source ~/.flexaidds_env   # or: source scripts/use_icloud_benchmark_storage.sh
# FLEXAIDDS_ICLOUD  = $HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks
# FLEXAIDDS_RESULTS = $FLEXAIDDS_ICLOUD/results
# FLEXAIDDS_QUEUE_ROOT = $FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1
```

| Role | Path (via env) |
|------|----------------|
| iCloud root | `$FLEXAIDDS_ICLOUD` |
| Campaign results | `$FLEXAIDDS_RESULTS/campaigns/` |
| Three-engine queue | `$FLEXAIDDS_QUEUE_ROOT` |
| Queue docs | `$FLEXAIDDS_ICLOUD/README_STORAGE.md` |
| **Legacy local** | `~/flexaidds_results/` — archive / binary staging only; **do not start new production campaigns here** |

Staged **binaries** may remain on local disk (queue `bin/` is a symlink to local staging) so iCloud sync does not corrupt Mach-O.

---

## 1. LIVE NOW — C0 full85 (arm C / FlexAIDdS)

| Field | Value |
|-------|--------|
| **Status** | **RUNNING** |
| **PID** | `10963` (`benchmark_datasets`; child `FlexAIDdS` on first targets) |
| **Lock / pid file** | `$FLEXAIDDS_QUEUE_ROOT/logs/C0_full85.{lock,pid}` |
| **Log** | `$FLEXAIDDS_QUEUE_ROOT/logs/C0_full85.log` |
| **Output** | `$FLEXAIDDS_RESULTS/campaigns/C0_full85_defined_cleft_nativeseed_forbidden/` |
| **Receipt** | `$…/C0_full85_…/RUN_RECEIPT.json` |
| **Mode** | `defined-cleft-redock` (native seed **OFF**) |
| **N** | 85 (Astex Diverse, cognate site, pose-blind) |
| **GA** | pop=1000 · gen=6000 · restarts=5 · **T=298 K** |
| **Env knobs** | `FLEXAIDDS_VCT_R0=4` `SHARING_ALPHA=4` `EVAL_SCALE_DIHEDRAL=-1` `FLEXAIDDS_NATIVE_SEED_FRAC=0` `FLEXAIDDS_SEED_ELITISM=0` |
| **Matrix MD5** | `72d7c7396702331d96ff12d18f831796` |
| **Binary SHA256 (C)** | `c7166f0291e6a59e0fbab40d8a44a67d28b1d1a355ebf0429eb087c4da2d37a6` |
| **Runner SHA256** | `32fa19af0dd66c19d74c7a21525bd60ef66b3bcce157905393d77cace7d9f05d` |
| **Started (UTC)** | 2026-07-15T02:54:52Z |
| **First target** | 1G9V (pose-blind + defined-cleft confirmed in log) |

### Do NOT

- Dual-launch into the same `OUT` dir.
- Kill a healthy run to “restart.”
- Report success as RMSD-only; primary S1 = elected RMSD ≤ 2 Å; S2 = S1 ∧ PoseBusters; do **not** sell BCR/S3 as abstract success.
- Mix with **seeded** oracle-ceiling numbers (see §3).

### Monitor

```bash
source ~/.flexaidds_env
PID=$(cat "$FLEXAIDDS_QUEUE_ROOT/logs/C0_full85.pid")
kill -0 "$PID" && echo LIVE || echo DEAD
tail -f "$FLEXAIDDS_QUEUE_ROOT/logs/C0_full85.log"
find "$FLEXAIDDS_RESULTS/campaigns/C0_full85_defined_cleft_nativeseed_forbidden" -name result.csv | wc -l
```

### Resume (only if dead mid-run)

Runner skips existing `result.csv` unless `--force`. Re-run the same launch script (preflight + dual-launch lock):

```bash
source ~/.flexaidds_env
"$FLEXAIDDS_QUEUE_ROOT/scripts/preflight_strict.sh"
"$FLEXAIDDS_QUEUE_ROOT/scripts/run_C0_full85.sh"
```

---

## 2. Three-engine entropy queue (protocol)

**Protocol:** `benchmarks/protocols/three_engine_entropy_comparison.md`  
**Queue status:** `$FLEXAIDDS_QUEUE_ROOT/STATUS.md`

| Stage | Status | Notes |
|-------|--------|--------|
| P0 preflight | **OK** | `scripts/preflight_strict.sh` |
| C0 smoke (2) | **DONE** | legacy local ref under queue `results/*_local_ref` |
| C0 pilot8 | **DONE** | S1 **2/8 (25%)**; all `seed_echo=0` / `native_pose_seeded=0`; election gap visible (S3 often better than S1) |
| **C0 full85** | **RUNNING** | §1 |
| FlexAID A / B0 / B | **Not wired** | binaries staged; ProcessLigand `.inp` harness missing |
| Full multi-arm P3 | **Blocked** | needs A/B harness + pilot multi-arm gate |

**Pilot8 targets:** 1G9V, 1GPK, 1MEH, 1P62, 1Q4G, 1R9O, 1T40, 2BYS  
**S1 hits (pilot):** 1P62, 1Q4G only.

Launch scripts (iCloud queue):

```text
$FLEXAIDDS_QUEUE_ROOT/scripts/preflight_strict.sh
$FLEXAIDDS_QUEUE_ROOT/scripts/run_C0_full85.sh [--dry-run]
```

---

## 3. Oracle-ceiling restore (seeded, separate science track)

| Field | Value |
|-------|--------|
| Path | `~/flexaidds_results/oracle_ceiling_restore_v43proto_r3` (local; completed) |
| Status | **COMPLETE** N=85 |
| BCR ≤ 2 Å ceiling | **83/85 = 97.65%** (`exceeds_90_ceiling: true`) |
| BCR fails | **1HNN**, **1HQ2** |
| success_pb | **40/85 ≈ 47%** (chemistry/PB gap ≠ placement ceiling) |
| Seed | **native seed ON** — **not** comparable to C0 full85 claim table |

Aggregate: use `scripts/aggregate_oracle_ceiling.py`.  
Do **not** dual-launch this dir; do **not** mix into three-engine headline rates.

---

## 4. Scientific guardrails (short)

- CF / contact-function scoring proxy ≠ experimental ΔG.
- TIER-1 claim path: no native seed, single matrix pin, S1 primary.
- PoseBusters + tENCoM/Eigen for claim-ready / modern secondary metrics on FlexAIDdS arms.
- Local free disk was tight (~14 GiB true free at launch); watch mid-run disk even though destination is iCloud (local APFS cache).

---

## 5. Codex / agent checklist when you open this repo

1. `source ~/.flexaidds_env` (or `scripts/use_icloud_benchmark_storage.sh`).
2. Read this file + `AGENTS.md` + `.agents/skills/flexaidds-benchmarking/SKILL.md`.
3. If PID live → **monitor only**; never relaunch same OUT.
4. Report rates from CSV/JSON on disk, not memory.
5. After any code change affecting the runner: fresh build + pin SHAs in queue `provenance.json` before claiming a new campaign.

---

## 6. Exact launch command (for forensics)

```text
caffeinate -i -s $QUEUE/bin/C/benchmark_datasets \
  --benchmark crossdock_json:$QUEUE/inputs/astex_native_85.json \
  --mode defined-cleft-redock \
  --output $FLEXAIDDS_RESULTS/campaigns/C0_full85_defined_cleft_nativeseed_forbidden/ \
  --threads 1 --omp-threads 4 \
  --ga-population 1000 --ga-generations 6000 \
  --temperature 298 --job-timeout-seconds 10800
```

(`$QUEUE` = `$FLEXAIDDS_QUEUE_ROOT`)
