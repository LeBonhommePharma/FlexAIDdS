# Offline Benchmarks queue — no docking while another session owns the box

**Schema:** multi-session coordination (Sol #9 / Codex)  
**Updated:** 2026-07-26  
**Mechanism:** `scripts/benchmark_coord.py` — `BENCHMARK_HOLD.json` + atomic `BENCHMARK_DOCK_LOCK/`

## Role split

| Role | May dock? | WORKERS | Notes |
|------|:---------:|--------:|-------|
| **Dock owner** (holds lock) | **yes** | ≤4 | Finishes current benchmarks; stamps binary into OUT |
| **Benchmarks / analysis** | **no** | n/a | Offline queue only until owner releases lock / hold |

Two independent agents share **no memory**. Convention alone cannot prevent dual docks.  
**Refuse to launch** if hold exists or lock dir is present.

## Before any analysis — read the other session’s commits

```bash
git -C "$FLEXAIDDS_ROOT" log -8 --oneline
```

Do **not** re-derive: HEM, population election, six baselines, matrix debates already on main.

## Offline queue (do these first)

1. **Full-population ceiling** on SEARCH-MISS `1J3J 1K3U 1L7F 1N1M 1M2Z` — frozen poses only; may reframe Phase 4.  
2. **Per-term CF decomposition** of SCORING-LOCKED gaps (+17.9 / +28.8 / +70.2 on 1OQ5/1SQ5/1YGC) — frozen poses, `probe_cf --config`.  
3. **Matrix pin resolve** — `72d7c739…` vs `9dc93717…` (`md5` / git history only).  
4. **Crystal-reference PoseBusters ceiling** — no docking; bounds every strict claim.

## Dock owner checklist

```bash
python3 scripts/benchmark_coord.py status
# free ≥ 20 GiB (or explicit FLEXAIDDS_DISK_FLOOR_OVERRIDE=1 emergency only)
python3 scripts/benchmark_coord.py preflight \
  --out ~/flexaidds_results/<run> --workers 2 \
  --binary build/FlexAIDdS --owner "<session-id>"
# use stamped binary path from JSON; never rebuild while peer run is live
# when finished:
python3 scripts/benchmark_coord.py release --token <token>
```

## Operator hold (manual freeze)

```bash
python3 scripts/benchmark_coord.py hold --owner ops --reason "disk emergency / dual-session freeze"
python3 scripts/benchmark_coord.py unhold
```

## Live G4.1 note

If `g4_1_boom_frac_*` (or any dock) is already running **without** this lock, do **not** steal ownership.  
New launches must call preflight; the in-flight owner finishes, then releases / unhold.
