# Grok Build Worker Dispatch — Oracle Record Recovery (2026-06-29)

> **Mandatory for every Grok Build worker session** touching benchmarks, `LIB/Vcontacts.cpp`,
> launch scripts, or full-85 campaigns. Machine-readable mirror:
> `docs/dev/grok_build_worker_orders.json`

## Context

| Run | Success | vs v109 (80/85) |
|-----|---------|-----------------|
| v109 record | 80/85 (94.1%) | — |
| v127 | 78/85 | −2 (Wave 1: parameter/selector) |
| v130 | 73/85 | −7 (+ Wave 2: Vcontacts/SoA path) |
| v131 HEAD | ~72/85 | −8 (broken binary — stop burning compute) |

**Persistent structural failures:** `1HNN`, `1N2V`, `1TW6` (until expB site + sulfo + holo fixes land).

## Active Phase DAG (do not skip)

```
NOW     Bisect Vcontacts (smoke-12 × 3)
        └─ scripts/queue_bisect_vcontacts.py --daemon  [PID in queue_bisect_vcontacts.state.json]
        └─ variants: safe → head_soa_off → head_soa_on

NEXT    v131_safe_full85 (if bisect passes ≥8/12, 0/3 guard fail)
        └─ scripts/launch_v131_safe_full85.py --skip-build --ignore-smoke-gate
        └─ binary: v131_safe @ 82ad51f4+sulfo+holo (pre-27e68e51 escape hatch)

THEN    Cherry-pick only proven Vcontacts fix onto HEAD
        └─ guilty commits since 82ad51f4: 27e68e51, f9c80fe5, d4d68592
        └─ rebuild HEAD; re-run smoke-12 before any full-85

LATER   r0=7 + expB/sulfo bundle (1G9V already recovering in v130/v131)
        └─ only after 78/85 recovery confirmed on safe binary

CEILING 83–84/85 after structural fixes on 1HNN/1N2V
```

## BLOCKED until bisect names guilty commit

| Campaign | Script | Reason |
|----------|--------|--------|
| v132 ablation ladder (all steps) | `queue_v132_ablation_ladder.py`, `launch_v132_ablation.py` | Combined knob turns — paused |
| v131 HEAD full-85 | `launch_v131_full85.py` | Wave-2 broken binary |
| v132 consensus/hbond/logsumexp ablations | `launch_v132_ablation.py *` | Deferred until Vcontacts fixed |
| New combined knob turns | any ad-hoc launcher | Audit violation |

## ALLOWED now

| Action | Script |
|--------|--------|
| Vcontacts bisect watcher | `queue_bisect_vcontacts.py --daemon` |
| Single bisect smoke variant | `launch_vcontacts_bisect_smoke.py {safe,head_soa_off,head_soa_on}` |
| Build bisect binaries | `build_vcontacts_bisect.sh`, `build_v131_safe.sh` |
| v131 safe smoke-12 | `launch_v131_smoke12.py` |
| Scalar perf (quiet queue) | `queue_after_v130_scalar_perf.py` (does not compete with bisect) |

## Worker startup checklist

```bash
cat docs/dev/grok_build_worker_orders.json
cat ~/Documents/PhD/Programs/FlexAIDdS/results/queue_bisect_vcontacts.state.json
cat ~/Documents/PhD/Programs/FlexAIDdS/results/queue_v132_ablation_ladder.state.json
# v132 status MUST be paused_for_vcontacts_bisect
```

## Do NOT

- Launch `launch_v131_full85.py` or any HEAD `build_lto` full-85 until bisect completes.
- Resume `queue_v132_ablation_ladder.py` until `vcontacts_bisect_summary.json` exists and names fix.
- Claim a new oracle record from v131 HEAD results (incomplete + regressed).
- Start competing full-85 campaigns while `queue_bisect_vcontacts` is `watching` or `launching`.

## Success gates

| Gate | Criterion |
|------|-----------|
| Bisect pass | `vcontacts_bisect_summary.json` best ≥8/12, regression guards (1HQ2,1S3V,1T40) clean |
| Recovery | `v131_safe_full85` ≥78/85 on per-target `result.csv` |
| HEAD unblocked | smoke-12 on patched HEAD matches safe ±1 after cherry-pick |