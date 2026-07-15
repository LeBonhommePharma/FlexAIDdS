# Astex-85 benchmark run — HANDOFF (Bonhomme → Grok agent)

**Handoff time:** Jul 11 2026 ~19:38 EDT
**Owning agent going forward:** Grok Build CLI session `019f1b82-2a99-7c93-8dce-c26656432aab`
**Reason:** user delegated run ownership to Grok; Bonhomme is standing down monitoring to avoid two agents touching the same processes.

## What is running
Detached benchmark launched via `benchmarks/astex_repro/run_full.sh`, logging to `benchmarks/astex_repro/full_FIXED.log`.
- Launch: Jul 10 17:05 EDT (launcher PID 40998 — now cross-session, shows EPERM to this sandbox but process pool is live)
- Protocol: `--mode autonomous`, GetCleft cavity-confined (FLEXAIDDS_CLEFT_SPHERE_DIR), NO oracle leakage (0 ORACLE/COGNATE injects verified)
- GA: pop=1000, gen=2000 base (scaled per flexibility), FLEXAIDDS_RESTARTS=10, FLEXAIDDS_NO_SEC=1
- Per-job timeout: 10800 s (3 h)
- Energy matrix PINNED: FLEXAIDDS_DATA_DIR=engine/ → MC_st0r5.2_6.dat md5 9dc93717dfed0698006d88dd6a9627bc (dock + provenance both record it)
- Engine pinned: FLEXAIDDS_BINARY=benchmarks/astex_repro/engine/FlexAIDdS

## Progress at handoff
- Completed: 8/85 (1G9V 1GM8 1GPK 1HP0 1HNN 1HQ2 1IA1 1IGJ). Currently docking 1J3J (#9, live, poses already emitted r6).
- ~3 h/target because 10 concurrent restarts saturate 18 GB RAM (swap) → every target burns the full wall.
- ETA at this rate: ~10 days for all 85.

## KNOWN ISSUE — zero-pose targets (recoverable)
- **1HNN, 1HQ2** completed with num_poses=0 (0/9 restarts finished GA loop before the 3h wall under heavy overnight load).
- MECHANISM (verified): each restart that finishes its GA loop (prints "TIMING SUMMARY") emits exactly 10 poses. total_poses = completed_restarts × 10. Zero-pose = 0 restarts crossed the wall. This is RESOURCE STARVATION, not an engine/scoring defect.
- RESUME BEHAVIOR: skip predicate (LIB/DatasetRunner.cpp:5375) fires only on num_poses>0 OR success. 1HNN/1HQ2 have num_poses=0 AND success=0 → they WILL redock on any resume. Safe.
- RECOMMENDATION: after main run finishes/crashes, do a targeted resume for any target with <~3 completed restarts (re-run run_full.sh; --skip-done keeps the good ones), ideally when the co-running benchmark session is idle so RAM contention is lower.

## Resume command (idempotent — skips completed targets)
    bash benchmarks/astex_repro/run_full.sh
(reads FLEXAIDDS_* from the script; completed targets with poses are skipped)

## RAM watch
18 GB cap, been at 0.8–1.2 GB free the whole run. Crash risk real when the other benchmarking session is also active. A crash costs only the in-progress target's 3 h; everything else is on disk.
