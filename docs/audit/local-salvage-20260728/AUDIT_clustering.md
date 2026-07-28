# AUDIT — Theme C: Clustering & Guards

**Auditor:** adversarial code auditor (OPS/CI discipline)
**Repo:** /Users/lp.more/Projects/FlexAIDdS @ main `d623c45ea`
**Commits:** `89d8dcd3e` (P2 cluster-rep election), `024ba8068` (revert of spread guard `d7ef67380`), `fb8a31225` (two-gate spread guard replacement)
**Verdict:** MINOR ISSUES — all three theme commits are code-sound; this is a *corrective* theme that cleans up two prior DEFAULT-OFF violations. The only findings are (1) a low-severity reference-point nuance in P2 and (2) a process red flag that is out-of-theme and already remediated by the revert.

---

## 89d8dcd3e — P2 cluster-rep election (FLEXAIDDS_CLUSTER_REP)

**What it does (verified from diff + live code):** Replaces the compile-time `#define OUTPUT_CLUSTER_CENTER` and the default-ON Boltzmann medoid (both from out-of-theme `3e674479c`) with ONE runtime gate `flexaids::cluster_rep_mode()` (LIB/ClusterRepMode.h). Modes: `unset|lowcf` (default), `medoid` (pure geometric), `bmedoid` (Boltzmann-CF-weighted, ablation), `center` (density-peak center). Adds opt-in `.pop.tsv` audit dump gated on `refstructure==1 && FLEXAIDDS_DUMP_POP=1`.

**default_behavior_changed:** **NO** — with no env vars set the mode is `LOWCF`.
- `ClusterRepMode.h:40-58`: unset `FLEXAIDDS_CLUSTER_REP` → falls through; legacy `FLEXAIDDS_MEDOID_REFINE` honored **only when explicitly non-zero** (default-ON removed) → returns `LOWCF`.
- `cluster.cpp:215-216`: medoid block entered only for `MEDOID||BMEDOID`; `LOWCF` skips it entirely (no table mutation, no stdout).
- `cluster.cpp:456`: provenance REMARKs emitted only `if (rep_mode != LOWCF)` → the default PDB is byte-unchanged.
- `DensityPeak_Cluster.cpp`: `output_cluster_center = (mode==CENTER)` = `false` by default → same as the prior `OUTPUT_CLUSTER_CENTER false`; REMARK string resolves identically ("the lowest CF").

**Reference-point nuance (the low finding):** "lowcf = current behavior" is byte-identical to the **pre-medoid baseline** (a HEAD build with `FLEXAIDDS_MEDOID_REFINE=0`), NOT to `3e674479c`-as-shipped, which had the medoid DEFAULT-ON. P2 deliberately turns that default-ON medoid OFF. This is the *intended correction* of a DEFAULT-OFF violation, so it satisfies the invariant (default = prior CORRECT behavior) — but a reviewer diffing only against the immediately-preceding commit will see the effective default change (medoid off). Documented correctly in the commit body; flagged low only so the reference point is not misread.

**Work-order sub-requirements — all satisfied:**
- (a) defaults to `lowcf` = pre-medoid behavior — **YES** (verified above).
- (b) medoid mode is PURE-geometric, not CF-weighted — **YES**: `cluster.cpp:236` initializes `weights(members.size(), 1.0)`; Boltzmann weights computed only `if (boltzmann)`. Pure medoid = `argmin_m Σ_{n≠m} ‖x_m−x_n‖²`, CF plays no role.
- (c) restores rather than removes the density-peak center — **YES**: `3e674479c` hardcoded `OUTPUT_CLUSTER_CENTER true→false` (removed center election); P2 converts it to runtime `output_cluster_center=(mode==center)`, so `center` restores the density-peak center as an option instead of deleting it.

**Determinism (invariant #3):** medoid tie-break is `if (cost < best_cost)` (strict `<`), `best_member` seeded to `old_head`, members collected in ascending chromosome index — lowest index wins ties, seed-independent (`cluster.cpp:264-280`). `rep_shift_src` keyed by chromosome index survives `QuickSort_Clusters` (which permutes cluster tables, not chromosome indices).

**Perf (invariant #6):** `cluster_rep_mode()` = one `getenv` at function entry (`cluster.cpp:210`), not per-eval. Medoid reuses `coord_cache` built once at `cluster.cpp:108` — no geometry rebuild. The `.pop.tsv` per-chromosome `ic2cf`+2×`calc_rmsd` re-score is O(num_chrom), runs **once at the very end** (only `free()`s follow, verified to EOF), and only when `FLEXAIDDS_DUMP_POP=1` — never on the benchmark hot path.

**Correctness:** SOUND · **Severity:** LOW · **Verdict:** makes_sense

---

## 024ba8068 — Revert of spread guard d7ef67380

**What it does (verified):** Reverts `d7ef67380` in full — removes the single-gate spread guard from `DatasetRunner.cpp` and the `cluster_spread_max{15.0f}` field + env/JSON plumbing from `ProtocolConfig.{h,cpp}`.

**Clean-inverse check:** `d7ef67380` diffstat = `+49/−2` across DatasetRunner.cpp/ProtocolConfig.{cpp,h}; revert `024ba8068` = `+2/−49` across the same three files — exact mirror image. The `-S(=O)2-N-` sulfonamide remap (`9450761d4`/PR #297) is correctly KEPT (context-scoped, unrelated).

**default_behavior_changed:** **YES, corrective** — the reverted guard's `cluster_spread_max{15.0f}` default fired on EVERY run (`proto.cluster_spread_max > 0.0f` was always true), demoting correct near-native rank-0 heads (evidence in commit: 1T46 0.16→15.78 Å, 2GBP 0.74→11.50 Å, 64 pass→fail / 0 fail→pass, median 6.4 Å global collapse; provenance.json confirmed `cluster_spread_max:15.0` active for the v132 run). The revert restores the canonical no-guard default → bit-identity to the pre-`d7ef67380` baseline. This change to the default is the fix, not a new risk.

**Correctness:** SOUND · **Severity:** NONE · **Verdict:** makes_sense

---

## fb8a31225 — Two-gate cluster spread guard (default-off replacement)

**What it does (verified from diff + live code):** Reintroduces a spread guard in `select_pose_freq_gated_pooled` requiring BOTH gates to fire before demotion: Gate 1 (rank-0 ≥ θ Å from all top-4 peers AND holds < `cluster_pop_min_fraction` of merged population; θ = 0.70·pocket_radius or population Q75 pairwise-RMSD spread), Gate 2 (< `cluster_consensus_k` restarts within `cluster_consensus_tau` Å). Skipped in oracle mode. On demotion it `std::swap`s rank-0↔rank-1 (the old guard `erase`d rank-0 outright).

**default_behavior_changed:** **NO — genuinely default-off / dead code.**
- `ProtocolConfig.h:87`: `cluster_spread_max{0.0f}`; entry gate `DatasetRunner.cpp:1263` is `if (proto.cluster_spread_max > 0.0f && ...)` → whole block skipped unless `FLEXAIDDS_CLUSTER_SPREAD_MAX` is set.
- All 5 knobs default no-op (`0.0f / 0.35f / 2.0f / 3 / 0.0f`); the only writers are env (`from_env`, guarded by `env_opt_*`) and JSON (`from_json`, guarded by `is_null`). Grep confirms **no DatasetRunner or default config sets it nonzero**.
- Regression test `test_protocol_config.cpp:169` (`DefaultsMatchHistoricalFallbacks`) pins `cluster_spread_max == 0.0f` with a comment citing the d7ef67380 collapse — the class of bug that broke it before is now CI-pinned.

**Cannot repeat the collapse:** the failure mode of `d7ef67380` was a nonzero *default*. Here the default is `0.0f`; even when enabled, demotion requires two independent gates plus a population-minority condition, and it swaps (keeps rank-0 in pool) rather than erasing. The v134 O(N²) `population_q75_spread` is cached (`std::optional q75_cache`, computed at most once per election) and only runs when the guard is on.

**Compile:** `<set>` (line 51) and `<optional>` (line 46) present; in-scope `PoseInfo` (line 1048) carries `restart`/`freq`/`cf`/`path`. Consistent with the claimed LTO build (ctest 77/79; the 2 failures — PoseBust parity, Mol2SdfReader — compile neither changed TU).

**Correctness:** SOUND · **Severity:** NONE · **Verdict:** makes_sense

---

## Theme summary

All three theme commits are code-sound and this is a **corrective/hardening theme**: `89d8dcd3e` fixes the three constraint violations `3e674479c` shipped (default-ON medoid, CF-reinjection, wrong-way DensityPeak `#define`), while `024ba8068`+`fb8a31225` undo a merged 15 Å-default spread guard that caused a 64/85 Astex collapse and replace it with a genuinely default-off two-gate version. Default-OFF/bit-identity, determinism, physics (pure-geometric medoid), and perf (no hot-path work, no per-eval rebuild) all hold with no env vars set. Live findings are limited to one low-severity reference-point nuance in P2.

## Highest-severity finding

**PROCESS red flag (out-of-theme, remediated) — the single most important thing a reviewer needs to know:** the reverted `d7ef67380` was a *merged* "fix" (via PR #298 / `5c136efe9`) that shipped a `cluster_spread_max{15.0f}` DEFAULT — a direct DEFAULT-OFF invariant violation that fired on every run and demoted correct near-native rank-0 heads, collapsing Astex to 18/85 (21.2%) from a 96.5% baseline. It passed review and merged. The theme's job was to undo it, and it does so cleanly (exact inverse diff) and replaces it with a default-off, CI-pinned two-gate guard. No live code in this theme repeats the violation. The lesson: "was this merged" ≠ "this is correct" — the theme demonstrates the invariant working *as remediation*, and the new guard's `0.0f` default is now regression-tested so the same class of bug cannot re-merge silently.
