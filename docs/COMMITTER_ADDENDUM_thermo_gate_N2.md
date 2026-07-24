# Committer addendum — thermo gate claim + N.2 remap validation

**Date:** 2026-07-22  
**HEAD baseline reviewed:** `d623c45ea` / later `a5a1a82c3`  
**Context:** Arm A CF.com blow-up deep-dive + independent source re-verification  
**Not a code change** — claims to correct, work items to schedule.

---

## 1. Thermo impossibility gate (`7a44c7035`) — claim overreach

### What shipped
- Header + unit tests for `thermo_gate::is_impossible` / `apply_gate` (`LIB/ThermodynamicEngine.h`).
- Gate forces `dG_eff → +1000` when `dH>0 ∧ dS<0`, behind `FLEXAIDDS_THERMO_SCORE`.
- Comments and docs (README, BENCHMARK.md, gaboom.cpp) state that the gate **affects selection** / demotes impossible poses so clustering cannot elect them rank-0.

### What the code actually does (verified)
| Claim | Reality |
|---|---|
| "promotes dG_eff to the ranking criterion" | **False for pose QuickSort.** Ranking is `QuickSort` on CF energy at `gaboom.cpp:705`, `:1147`, `:1238` — finalized **before** `thermo_engine->compute()` at `:1308`. |
| "clustering never selects impossible poses rank-0" | **Not wired.** Repo consumers of `dG_eff` / `thermo_impossible` outside the engine are **printf only** (`gaboom.cpp:1338–1358`). |
| Gate runs by default | **No.** Requires `thermo_engine_enabled` (default false) **and** `FLEXAIDDS_THERMO_SCORE` truthy. |

`ProtocolConfig.h` still documents:

> promote ΔG_eff = ⟨CF⟩ − T·H …

`gaboom.cpp:1336–1337` still says:

> Reporting-only unless FLEXAIDDS_THERMO_SCORE=1, which promotes dG_eff to the ranking criterion in place of min(CF).

That promotion path is **not implemented** in the selection pipeline. Arm B (`ops/run_astex85_twoarm.sh`) sets `FLEXAIDDS_THERMO_SCORE=1` and Softβ election — Softβ reorders on CF-derived soft free energy, **not** on the gated `dG_eff` sentinel.

### Required action (pick one)
**Option A — wire it (feature work):**  
After `compute()`, if `thermo_score_enabled()`, re-rank / re-elect using `dG_eff` (or attach sentinel to cluster ACF). Add a unit/integration test that a synthetic impossible population loses election.

**Option B — retract the claim (docs-only, preferred if A is out of scope):**  
- Rewrite commit-message-level language in `ThermodynamicEngine.h`, `ProtocolConfig.h`, `gaboom.cpp` comments, `README.md`, `docs/BENCHMARK.md`.  
- State clearly: **diagnostic / reporting only; does not change pose or cluster selection.**  
- Leave unit tests as pure gate-math tests.

**Do not leave the middle ground** — operators will trust the comment and misread Arm B results.

### Secondary physics caveats (if Option A)
- Gate uses CF as ΔH proxy (`dH_i ≈ CF_i`); sign of CF ≠ calibrated enthalpy.
- A single complex-level `TdS_vib` is broadcast to all poses — "per-pose impossibility" is really one dS gate across the population.

---

## 2. N.2 → N.ar (`3f56c34b1`) — needs A/B, not settlement

### What shipped (unconditional default path)
`LIB/top.cpp:78`:

```cpp
if (!strcmp(s, "N.2"))   return 10;  // N.ar — sp2 imine is an acceptor; N.am (donor) reversed the H-bond sign
```

Also remaps C.1→C.2, I→Br in the same commit. **Not env-gated.**

### Two live, defensible destinations

| Mapping | Row | Entries (live matrix) | Argument |
|---|---|---|---|
| **N.2 → 7** (type-exact) | N.2 | 13 nonzero | Identity: keep imine on its own SYBYL row; auditor preference |
| **N.2 → 10** (shipped) | N.ar | 20–22 nonzero | Chemistry: sp2 imine is acceptor; prior N.am(11) reversed H-bond donor/acceptor sign |

Both rows are **live** (unlike N.3 row 8, which is all-zero — that fix in `3bdde9932` is **sound** and should stay).

### Required action
1. **Do not treat N.2→10 as settled** until pose-quality evidence exists.
2. A/B on a polar / imine-containing subset (or full Astex with typed N.2 ligands):
   - Arm X: `N.2 → 7` (type-exact)
   - Arm Y: `N.2 → 10` (current)
   - Metrics: success@2Å, native CF rank, H-bond REMARK counts on elected poses
3. Gate experimental remaps behind `FLEXAIDDS_ATOM_REMAP_EXPERIMENTAL=1` if shipping further typing changes without canary.

### Keep
- **`3bdde9932` N.3 → N.am (11)** — correct; dead row 8 was zeroing live chemistry.

---

## 3. Unrelated one-liner still open
`BustCli.cpp` mandatory key `no_protein_clashes` vs upstream bust 0.6.5 — 1/79 ctest fail. Independent of Arm A regression; one-line schema pin.

---

## 4. Relation to Arm A CF.com blow-up
These two items are **not** the primary cause of Arm A CF~−3000 failures. That cause is **raw extensive CF.com with `FLEXAIDDS_VCT_NORM` and `FLEXAIDDS_COM_FLOOR` left off** (see `ops/canary_com_tame.env`).  

This addendum is independent hygiene so future operators do not:
- believe the thermo gate ranks poses when it only prints, or
- treat N.2→N.ar as validated chemistry without a canary.
