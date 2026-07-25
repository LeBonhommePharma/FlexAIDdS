# Genuine Astex baseline — autonomous blind (2026-07-24 run)

**Hub:** [`COMPARATIVE_SCIENCE_README.md`](COMPARATIVE_SCIENCE_README.md)  
**Status:** recorded reference (pre-merge vs later PRs; not a 3Dsig S_top10 claim table).  
**Recorded from:** Claude Science relay (2026-07-25); operator confirmation of completed 80/85 archive.  
**Campaign id:** `v_autonomous_20260724_160919`  
**Metric contract:** **genuine** = rank-0 RMSD ≤ 2.0 Å **and** `seed_echo=0` (sentinel-guarded).  
**Not:** JCIM top-10 / 3Dsig S_top10 bootstrap median (different contract).

---

## Headline

| Metric | Value |
|--------|------:|
| **GENUINE top-1 ≤2 Å** | **20 / 79 = 25.3%** |
| Best-cluster ≤2 Å (BCR / sampling ceiling) | 22 / 79 = **27.8%** |
| Election gap (BCR success − genuine) | **2 targets** |
| Seed-echo contamination | **0** |
| Median RMSD (reported) | 3.99 Å |
| Targets with scores in denominator | 79 (80/85 finished; 5 never scored) |

**Success PDBs (20):**  
1HNN, 1HQ2, 1OPK, 1P62, 1Q1G, 1Q41, 1R1H, 1T46, 1TZ8, 1U4D, 1UML, 1V4S, 1W1P, 1Y6B, 1Y6R, 1YQY, 1YWR, 2BM2, 2BSM, 2HB1

---

## Scientific interpretation (Science + forward plan)

1. **Clean multi-target number.** Zero seed-echo; usable as a **pre-merge** genuine baseline under autonomous blind.  
2. **Election gap collapsed.** Only ~2 targets sample sub-2 Å but elect worse (BCR 27.8% vs genuine 25.3%). Consistent with `free_energy_strict` / ACF-strict election work removing largest-cluster popularity bias.  
3. **Bottleneck is sampling.** Ceiling **27.8%** — GA finds a sub-2 Å head on 22/79. Election/scoring alone cannot exceed BCR.  
4. **Route:** FORWARD plan **Wave 3 (sampling / BCR raisers)**; wall oracle / wall-before-memetic remains pacing for refine levers.  
5. **Caveats:** 80/85 not full 85; run **predates** PR #300/#301 merges — re-baseline after claim binary pin if citing post-merge rates.

Near-miss BCR heads worth inspection (election-gap class, not the main 20-point gap): 1OQ5 (~1.06 Å), 1SQ5 (~1.12), 1YGC (~1.24), 1YVF (~1.70) — small headroom.

---

## Relation to comparative P0–P5 pipeline

| Item | Status |
|------|--------|
| Pipeline gates | Wired; 17 unit tests; `--pipeline-dry` |
| Live blocker | **P2** — needs real `native_cf_oracle_gate` JSON (`ok` / `exit_code` / `ranking_forbidden`), not empty/deferred |
| Matrix for claim | **9dc9** (`9dc93717dfed0698006d88dd6a9627bc`) — same as this baseline era pin; do not confuse with 72d7 packing fork |
| This 25.3% figure | **Genuine S1-style**, autonomous — **not** a substitute for arm A/B **S_top10** comparative table |

---

## Operator env (Science)

```bash
export PYTHONPATH=$PWD/python
# CLIs need numpy:
~/.claude-science/conda/envs/python/bin/python scripts/run_comparative_phases.py --pipeline-dry
~/.claude-science/conda/envs/python/bin/python scripts/comparative_phase_gate.py --dry-run
# pytest needs an env with pytest:
PYTHONPATH=$PWD/python ~/.claude-science/conda/envs/cpp-python-core/bin/python -m pytest python/tests/test_comparative_phases.py -q
```

Flag is **`--pipeline-dry`**, not `--dry-run`, on `run_comparative_phases.py`.
