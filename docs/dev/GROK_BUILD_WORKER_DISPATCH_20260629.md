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

**Recent progress:** Vcontacts fix landed (revert 27e68e51 + f9c80fe5 identified; `FLEXAIDS_USE_SOA_DISTANCES=OFF` default; PR4 parity). v132_ablation_ladder active. Guard-bisect experiment (1HQ2/1T40) → best 1/2, "deeper selector bisect needed". New commit-bisect tools added for remaining suspects.

## Active Phase DAG (do not skip)

```
NOW     v132 ablation ladder
        └─ scripts/queue_v132_ablation_ladder.py --daemon
        steps: consensus_on → safe_binary → logsumexp_only → hbond_zero
        (Recent guard_bisect on 1HQ2/1T40; deeper commit-bisect tools for remaining Vcontacts suspects)

NEXT    v131_safe_full85_recovery (deferred until after v132 consensus_on baseline)
        └─ scripts/launch_v131_safe_full85.py --skip-build

LATER   r0=7 + expB/sulfo bundle (1G9V already recovering)
        └─ only after v132 baseline + 78/85 recovery confirmed on safe/ patched binary

CEILING 83–84/85 after structural fixes on 1HNN/1N2V + selector work
```

**Vcontacts status (completed):** bisect identified 27e68e51 (catastrophic) + f9c80fe5. Reverted 27e68e51; default `FLEXAIDS_USE_SOA_DISTANCES=OFF`; PR4 scalar parity landed. Guard/commit-bisect ongoing for remaining suspects.

## BLOCKED

| Campaign | Script | Reason |
|----------|--------|--------|
| v131 HEAD full-85 | `launch_v131_full85.py` | Under test via v132 ladder + selector work |
| combined knob turns | any ad-hoc | Audit violation until v132 ladder baselines land |

## ALLOWED now

| Action | Script |
|--------|--------|
| v132 ablation ladder | `queue_v132_ablation_ladder.py --daemon`, `launch_v132_ablation.py` |
| Vcontacts commit bisect (remaining suspects) | `queue_bisect_vcontacts_commits.py`, `launch_vcontacts_commit_bisect_smoke.py`, `build_vcontacts_commit_bisect.sh` |
| v132 guard / isolation | `launch_v132_guard_bisect.py`, `launch_v132_isolation4.py` |
| v131 safe smoke-12 / recovery | `launch_v131_smoke12.py`, `launch_v131_safe_full85.py`, `build_v131_safe.sh` |
| Scalar perf (quiet queue) | `queue_after_v130_scalar_perf.py`, `launch_perf_scalar_quiet.py` |
| Legacy bisect tools | `queue_bisect_vcontacts.py`, `launch_vcontacts_bisect_smoke.py` |

## Worker startup checklist

```bash
cat docs/dev/grok_build_worker_orders.json
cat ~/Documents/PhD/Programs/FlexAIDdS/results/queue_v132_ablation_ladder.state.json
cat ~/Documents/PhD/Programs/FlexAIDdS/results/v132_*guard*/v132_guard_bisect_report.txt 2>/dev/null || true
# v132 status should be "active" or "launching"
```

## Do NOT

- Launch `launch_v131_full85.py` or HEAD `build_lto` full-85 until v132 ladder completes a consensus_on baseline.
- Resume old v132 ablation until current queue state shows unblocked.
- Claim new oracle records from incomplete/HEAD runs.
- Start combined knob turns while ladder or commit-bisect is active.

## Success gates

| Gate | Criterion |
|------|-----------|
| v132 ladder step | consensus_on baseline established; compare vs v130/v131; guards tracked |
| Deeper Vcontacts | commit-bisect summary names next action for f9c80fe5 / d2295cf0 |
| Recovery | `v131_safe_full85` or patched equivalent ≥78/85 (v109 record matched at 80/85 on safe) |
| Unblock further | ladder passes + selector fixes land → re-enable other campaigns |