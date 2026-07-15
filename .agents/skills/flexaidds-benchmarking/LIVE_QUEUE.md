# LIVE GOAL — Three-engine Astex Diverse 85 (PRIORITY)

**Updated (UTC):** 2026-07-15T05:04Z  
**Primary science goal (highest priority):**  
**FlexAID 2015-era (A) vs FlexAID master entropy off (B0) vs FlexAID master entropy on (B) vs FlexAIDdS (C0/C)**  
on **Astex Diverse N=85**, TIER-1 cognate, **no native seed**, **identical matrix**.

Oracle-ceiling restore is **complete** (diagnostic / seeded ceiling only — **not** the claim path).

---

## Scientific question

Under a matched cognate-pocket, **no native-seed** redocking protocol and a **single pinned energy matrix**, how do docking power and ranking change across:

| Arm | Engine | Entropy ranking |
|-----|--------|-----------------|
| **A** | FlexAID 2015-era / pre-polish lineage | Off (CF default) |
| **B0** | FlexAID current master | **Off** (`TEMPER 0` → CF clustering) |
| **B** | FlexAID current master | **On** (`TEMPER 298`, BindingMode / FO path) |
| **C0** | FlexAIDdS (pinned binary) | CF / standard DatasetRunner election (claim path) |
| **C** | Same SHA as C0 | Full stack; thermo ledger reported separately |

**Primary KPIs:** S1 (elected RMSD ≤ 2 Å), S2 (S1 ∧ PoseBusters).  
**S3 / BCR:** diagnostic only — never headline success.  
**Admission:** `seed_echo=0`, `native_pose_seeded=0`, matrix MD5 pin.

**Normative docs:**
- `benchmarks/protocols/three_engine_entropy_comparison.md`
- `benchmarks/protocols/admission_metrics_contract.md`
- Aggregator: `python3 scripts/aggregate_claim_metrics.py`

---

## Hard invariants

| Item | Value |
|------|--------|
| Dataset | Astex Diverse **85** |
| Matrix | `MC_st0r5.2_6.dat` MD5 **`72d7c7396702331d96ff12d18f831796`** |
| Search | pop **1000** · gen **6000** · restarts **5** · T **298** (thermo arms) |
| Seed | **Forbidden** for claim rows |
| Storage | **iCloud only** (`$FLEXAIDDS_ICLOUD` / `$FLEXAIDDS_RESULTS`) |
| RAM (this Mac) | **~18 GiB max** — **one heavy campaign at a time** |

---

## Queue status (this machine)

| Job | Status | Notes |
|-----|--------|--------|
| Oracle-ceiling restore (seeded) | **DONE** | 83/85 BCR (97.65%); **not** three-engine claim |
| **C0 full85** FlexAIDdS | **LIVE** | `$FLEXAIDDS_QUEUE_ROOT/logs/C0_full85.pid` |
| FlexAID **A→B0→B pilot8** | **PAUSED** (SIGSTOP, RAM guard) | Resume only when C0 done + free RAM ≥ 3 GiB |
| Full A/B0/B N=85 | **Blocked** | After pilot8 gate + C0 complete |
| Entropy isolation C re-rank | Later | Prefer re-rank frozen C0 ensembles when I/O allows |

**Resume A/B pilot (after C0 finishes):**
```bash
source ~/.flexaidds_env
bash "$FLEXAIDDS_QUEUE_ROOT/scripts/ram_guard_resume_pilot.sh"
export OMP_NUM_THREADS=1 FLEXAIDDS_PARALLEL_RESTARTS=0
bash "$FLEXAIDDS_ROOT/scripts/run_A_pilot8.sh"
bash "$FLEXAIDDS_ROOT/scripts/run_B0_pilot8.sh"
bash "$FLEXAIDDS_ROOT/scripts/run_B_pilot8.sh"
```

**Aggregate claim metrics:**
```bash
python3 scripts/aggregate_claim_metrics.py --c0-full85
python3 scripts/aggregate_claim_metrics.py "$FLEXAIDDS_RESULTS/campaigns/three_engine/A/pilot8"
```

---

## Deprioritized (do not steal cycles from three-engine)

- Oracle-ceiling re-runs / seeded ceiling campaigns  
- Packaging/Homebrew unless blocking science binaries  
- Repo hygiene refactors that do not unblock A/B/C docks  
- Fleet production orchestration until three-engine pilot + full85 land  

---

## Agent operating rule

When in doubt: **advance three-engine Astex 85 (A / B0 / B / C0)** with correct matrix pin, no seed, S1/S2 reporting, and RAM-safe serialization. Everything else waits.
