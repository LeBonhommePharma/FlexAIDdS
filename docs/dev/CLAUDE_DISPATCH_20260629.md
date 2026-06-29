# Claude Dispatch — FlexAIDdS v124–v127 Campaign (2026-06-29)

> **Read first**: `AGENTS.md` (repo root). This is the operational dispatch from the Grok session that synced `origin/master`, resumed benchmarks, queued v127, and spawned monitoring + literature audit.

---

## Dispatch summary

**Mission**: Steward the in-flight Astex 85 oracle-ceiling campaign (v126 → v127), finalize v124 metrics, fix consensus-default drift, and produce a defensible v50b/v124/v126/v127 comparison — without disturbing running daemons or overclaiming thermodynamics.

---

## Live state (refresh from logs)

| Run | Dir | Status |
|-----|-----|--------|
| v124 | `results/v124_full85_20260626_0413_consensus_guard` | **DONE** 85/85 dirs · canary 6/6 · CSV may need regen (was 59 rows) |
| v126 | `results/v126_20260628_2347_optB_smoke` | **RUNNING** PID 13781 · `/tmp/FlexAIDdS_v126` |
| v127 | auto via `scripts/queue_after_v124_v126.py` | **QUEUED** after v126 exits |
| Audit | `results/scientific_robustness_audit_20260629.md` | pending (Grok background agent) |

**Daemons — do not kill**: queue watcher ~51275, monitor ~52272.

```bash
tail -f ~/Documents/PhD/Programs/FlexAIDdS/results/campaign_monitor.log
tail -f ~/Documents/PhD/Programs/FlexAIDdS/results/queue_after_v124_v126.log
```

---

## P0 tasks (execute in order)

1. **Monitor only** — let v126 (PID 13781) finish. No competing full-85 launches.

2. **Confirm v127 launch** — when v126 exits, check `queue_after_v124_v126.log` for `v127 launched`. Fallback:
   ```bash
   python3 scripts/launch_v127_full85.py
   ```

3. **Regenerate v124 summary CSV** if `astex_diverse_results.csv` still has empty rows — aggregate 85 `*/result.csv` files. Report success ≤2Å, sentinels, canary. **Flag mixed binary**: first 59 @ SHA `8a00dfdd`, last 26 @ HEAD.

4. **Read literature audit** when it lands; fold P0/P1 risks into v127 write-up.

5. **After v127 completes** — comparison table: v50b (81.2%) vs v124 vs v126 vs v127, oracle-ceiling only, ≤2Å.

---

## P1 fix

**Consensus default drift** — `DatasetRunner.cpp` defaults `FLEXAIDDS_CONSENSUS_SCORER=0` (v125); docs/scripts say `1`. Restore default `1` post-logsumexp; sync `BENCHMARK_STANDARD.md` + `REPRODUCIBILITY.md`.

---

## Scientific guardrails

- CF scoring proxy language — not true ΔG.
- Label `ORACLE_CEILING` vs `AUTONOMOUS` always.
- Do not revert SMFREE→PSHARE without benchmark (v117: 81%→37%).
- No `LIB/` scoring changes without tests; v127 validates `ba5364d3` H-bond/VCT patch.

---

## Version commits

`15b536f8` v124 consensus-guard · `ce8f3368` v125 diagnostic · `a4056163` v126 logsumexp · `ba5364d3` H-bond/VCT · `d0be195b` queue/v127 launcher · `bfd56f45` monitor

---

## Do NOT

- Kill PIDs 13781 / 51275 / 52272 without user OK.
- Force-push `master` without explicit user OK.
- Claim final rates before 85 populated CSV rows.

---

*Update when v126/v127 complete.*