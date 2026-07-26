# Next campaign step — score-only Native/Elected CF inversion map

**Date:** 2026-07-25  
**Authoring tip:** main @ `b19704bc` lineage (methodology gates landed)  
**Status:** **RECOMMENDED — not launched** (analysis goal; no dock in this workorder)  
**Hub:** [`CAMPAIGN_GATE_SUMMARY.md`](CAMPAIGN_GATE_SUMMARY.md) · audit [`DOCKING_BUG_AUDIT_2026-07-25.md`](DOCKING_BUG_AUDIT_2026-07-25.md)

---

## Primary recommendation (one procedure, zero GA knobs)

### Experiment name
**Native–Elected–BCR CF inversion map** (score-only, production LOCCLF)

### Exactly one controlled choice
**Pose role under a fixed scoring surface** — score three pose classes with the **same** production dock config and `probe_cf` (no GA, no BOOM, no WAL, no memetic):

| Pose role | Source |
|-----------|--------|
| **native** | crystal ligand SDF (`*_ligand.sdf`) |
| **elected** | `elected_pose.pdb` / `elected_pose_path` from a frozen local run leaf |
| **bcr_head** (optional when local PDB exists) | cluster head achieving `best_cluster_rmsd` if path recoverable |

**Fixed (not varied):** matrix unused (score-only); production `ops/gates/configs/{PDB}_dock_config.json`; `FLEXAIDDS_WAL_COERCIVE=0`; `FLEXAIDDS_PB_CLASH_WEIGHT=0` (keep burial off so this map is pure CF.com/wal landscape); workers N/A.

If a later **dock** is needed after this map, it is a **separate** one-variable run (see §Follow-on). Do not bundle.

### Panel (cheap, discriminating)
Reuse the methodology clean probes + election-gap near-misses already used in STEP 3:

`1J3J 1K3U 1L7F 1N1M 1M2Z 1OQ5 1SQ5 1YGC`

Prefer poses from a **local** leaf that still has PDBs, e.g.  
`~/flexaidds_results/pilot_w1_boom_interval_20260725_134740/{PDB}/elected_pose.pdb`  
(not a raw iCloud tree walk). Baseline archive thin `result.csv` alone is insufficient without pose files — materialize only the 8× elected PDBs if needed.

### Protocol (score-only)
```bash
# For each PDB, each pose in {native.sdf, elected_pose.pdb}:
build/probe_cf \
  --receptor ~/.flexaidds/benchmarks/astex_diverse/$PDB/${PDB}_apo.pdb \
  --pose <pose> \
  --ligand ~/.flexaidds/benchmarks/astex_diverse/$PDB/${PDB}_ligand.sdf \
  --config ops/gates/configs/${PDB}_dock_config.json \
  --binary build/FlexAIDdS --data-dir .
# Require --config (methodology: without it CF inflates ~200×).
```

Record: `cf_total`, `cf_com`, `cf_wal`, `cf_clash`, paths, binary sha256, git tip, full `FLEXAIDDS_*` env (should be empty of BOOM/WAL/PB for this map).

### Classification (the inventive product)

For each target with both scores:

| Class | Definition (ε ≈ 0.5 CF units) | Implication for **next** dock lever |
|-------|--------------------------------|--------------------------------------|
| **SCORING-LOCKED** | `cf_elected + ε < cf_native` | Landscape prefers the false min; **do not** expect BOOM/coarse alone to elect native. Prefer scoring/burial opponent work (strong pb_clash decoys) or accept honest fail on that target. |
| **SEARCH-MISS** | `cf_native + ε < cf_elected` | Near-native is **better** on CF but was not found/elected → **sampling** lever. |
| **TIED / NOISE** | \|cf_native − cf_elected\| ≤ ε | Uninformative; need more poses or BCR head. |

Cross with E10/BCR labels from `workorders/E10_election_vs_scoring.md` and pilot `result.csv`:

- SEARCH-MISS ∧ BCR ≤ 2 Å → **election/pool** issue (near-native in pool, not rank-0).  
- SEARCH-MISS ∧ BCR ≫ 2 → pure **sampling ceiling**.  
- SCORING-LOCKED ∧ multi-restart identical CF (e.g. 1N1M elect CF ≈ −99.314 on multiple restarts per audit) → **false-min attractor**, not “need more BOOM”.

### ACCEPT / FAIL (this step only)

| Layer | PASS | FAIL |
|-------|------|------|
| **Instrumentation / completeness** | ≥6/8 targets have both native + elected CF with production config; 1M2Z native spot-check still ~−117.7 (sanity) | Missing configs/poses; accidental no-`--config` inflation |
| **Science product** | Every scored target labeled SCORING-LOCKED / SEARCH-MISS / TIED with dCF table written to a workorder | Ambiguous labels without ε; mixing pb_clash-on into the map |
| **Docking quality** | **Not claimed** — no genuine/BCR rate from this step | Citing inversion map as Astex success-rate change |

### Why this is inventive (what prior gates did not unlock)

1. **Stops lever roulette.** The campaign already burned compute on **unwired** BOOM interval ([`STEP3_PILOT_GATE.md`](STEP3_PILOT_GATE.md), B1) and **structurally unpassable** WAL ([`WALL_ORACLE.md`](WALL_ORACLE.md), B3). BOOM_FRAC=0.1 is **live** ([`BOOM_FRAC_AB.md`](BOOM_FRAC_AB.md)) but left **1N1M on the same elect CF as control** (−99.314) — diversity inject without knowing whether native is even CF-better is cargo-cult sampling.  
2. **Splits the bottleneck that E10 only half-answers.** E10 ([`E10_election_vs_scoring.md`](E10_election_vs_scoring.md)) shows sampling primary **globally** (election-gap ~18.8%). It does not, per target, answer: *is the elected decoy winning on CF against the crystal?* That single bit decides whether the next dollar of compute buys **search** or **scoring**.  
3. **Uses evidence already in hand.** Multi-restart CF consensus on 1N1M ([`DOCKING_BUG_AUDIT_2026-07-25.md`](DOCKING_BUG_AUDIT_2026-07-25.md) §B4/election) is a fingerprint of SCORING-LOCKED attractors; the map **tests** that hypothesis with production LOCCLF rather than another GA.  
4. **Cheaper than another W1 dock** and does not violate R7 (no full-85), does not enable memetic, does not re-run WAL.

### Explicit non-claims

- Not a genuine top-1 or BCR rate.  
- Not permission to set `WALL_PILOT_PASS` or memetic (pb_clash formal PASS remains micro-effect only — [`PB_CLASH_ORACLE.md`](PB_CLASH_ORACLE.md)).  
- Not a product change to claim `boom_inject_fraction`.  
- Not proof that `free_energy_strict` fixed election (still unmeasured on a post-fix full-85).

---

## Follow-on (only after map; still one variable each)

**If SEARCH-MISS dominates the panel:**  
- **One variable:** `FLEXAIDDS_COARSE_ORIENTATIONS=256` vs claim default 64 (requires `coarse_init.enabled` in dock JSON — confirm DatasetRunner emit or env path before launch).  
- Panel: same 8; WORKERS≤2; R≤2; matrix **9dc9** (`md5 9dc93717dfed0698006d88dd6a9627bc`); seed OFF.  
- ACCEPT (docking): BCR count ↑ or mean BCR ↓ on SEARCH-MISS subset; no clean-probe genuine regression; liveness: log/config shows n_orientations=256.  
- Rationale: BOOM re-injects **mid-run** then reconverges; coarse orientations expand **t=0 basin entry** (Wave 3 K2 in `docs/implementation/WAVE3_SAMPLING_BCR_PLAN.md`).

**If SCORING-LOCKED dominates (esp. 1N1M-class):**  
- **One variable:** stronger deep-interpenetration decoy construction for pb_clash (not weight ladder yet) until `cf_clash` is non-trivial — then re-run pb_clash oracle.  
- Do **not** burn another BOOM_FRAC pilot on those targets first.

**If mixed:** dock coarse-orient **only** on SEARCH-MISS codes (`--only-codes`); leave SCORING-LOCKED out of BCR claims.

---

## Rejected alternatives (constrained novelty)

| Rejected | Why not now |
|----------|-------------|
| Re-run **WAL_COERCIVE** | B3 structural no-op; OFF≡ON forever ([`CAMPAIGN_GATE_SUMMARY.md`](CAMPAIGN_GATE_SUMMARY.md)) |
| **BOOM_INTERVAL** alone | B1 invalid; claim frac=0.0 |
| **BOOM_FRAC efficacy panel** as next primary | Liveness already PASS; 1N1M same false-min elect — without inversion map this thrash repeats |
| **Full-85** / dual full-85 | R7; gates incomplete; WORKERS≤4 still not enough reason to launch |
| **Memetic / WALL_PILOT_PASS** | Blocked; micro pb_clash not unlock |
| **NO_SEC full 2000 gens** first | Valid later one-var honesty check, but confounded if target is SCORING-LOCKED (burns wall-clock proving budget on unwinnable landscape) |
| **Niche Cartesian rewrite** | Real Wave 3 item but multi-file product change; needs flag + tests; not the cheapest discriminator |

---

## Campaign constraints (carry forward)

- One variable per **dock** run; WORKERS≤4; OMP=1/worker; no dual full-85.  
- Genuine metric only for rates (seed_echo=0 ∧ rank-0 RMSD&lt;2); never raw `success_rmsd` alone.  
- Matrix **9dc9** when docking; provenance: binary sha256, git tip, full `FLEXAIDDS_*`.  
- Baseline 25.3% remains **pre-`free_energy_strict`** reference, not election-fix proof.

---

## Cadence line (for the next owner)

**Phase:** W0.5 diagnostic (score-only)  
**One variable:** pose role under fixed production CF (native vs elected)  
**Genuine/BCR:** not applicable this step  
**PASS/FAIL:** completeness + classification table (see ACCEPT)  
**Blocks:** full-85, memetic, WAL re-panel, interval-only BOOM  

---

*Analysis-only deliverable. Experiment not executed in the session that authored this file.*
