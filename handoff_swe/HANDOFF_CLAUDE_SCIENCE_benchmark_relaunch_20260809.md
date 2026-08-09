# Handoff → Claude Science: Astex-84 relaunch is CLEARED

**From:** Claude Code (Fable 5) · **Date:** 2026-08-09 · **Main at:** `9035cf34`
**Shannon task:** `flexaidds_astex84_dG_relaunch_20260809`
**Recommended executor:** Claude Science (`benchmark_execute`)

---

## MISSION AND STOP CONDITION

Run the Astex-84 relaunch and produce **one number that is attributable** — a
success rate whose engine state, election rule and inputs are all pinned and
recorded, so it can be compared to a later run rather than being a draw from an
unknown distribution.

**Stop condition:** 84/84 complete, `INPUT INTEGRITY: OK`, `provenance.txt`
written, and the top-1 <2 Å rate computed **by you** from `rmsd_to_crystal`
(symmetry-corrected) — never from the `success` column, which means "docking
ran".

## METHODOLOGY MANDATE — reproducible-first, verify before you claim

This campaign exists because the previous two were not reproducible. Operate
accordingly:

1. **Measure, don't infer.** Every claim in this document is backed by a
   command that was run. Hold your own output to that standard: if you report a
   rate, show the query that produced it.
2. **Hash coordinates, never whole files.** `^(ATOM|HETATM)` lines only.
   `BindingMode.cpp:739` stamps `REMARK FLEXAID.commit=... dirty=... seed=...`
   into every pose on every build, so a whole-file md5 has a ~100% false-positive
   rate. I made exactly this mistake on #405 and wrongly reported a science
   change; the coordinates were byte-identical. Do not repeat it.
3. **Always run a same-arm control.** Before attributing a difference to a
   change, run the *unchanged* arm twice. When I saw the `.rrd` differ between
   two conditions, the control showed two identical-condition runs differed on
   MORE lines — the effect was pre-existing nondeterminism, not the change.
4. **One variable per comparison.** If you want the representative rule
   measured, hold grid, GA trajectory and mode ranking fixed and emit both
   representatives from the SAME run (§5). Two campaigns cannot give that
   cleanly even post-#403.
5. **Fail closed on absent fields.** `MAX_REMARK` is 5000 and `safe_remark_cat`
   truncates silently; `soft_beta_G` / `free_energy` / `proxy_free_energy` are
   appended late (`cluster.cpp:556-561`), so they are what disappears first.
   **Spot-check these REMARKs on the first completed target, not after 84.**
   "Field absent" must never be read as "no change".
6. **Quote no number across an engine-behaviour boundary.** Pre-#403
   multi-thread results are void. If in doubt, re-run rather than pool.

---

## VERDICT: clear to relaunch. Three config changes required.

The variance that made the last two campaigns incomparable is fixed and verified.
Everything below is measured, not asserted.

---

## 1. What was wrong, and what fixed it

**Root cause (not the GA).** `CleftDetector::generate_probes` merged per-thread SURFNET
probe vectors under `#pragma omp critical` in **thread-arrival order**. Probe order fixes
`generate_grid()`'s cleftgrid index assignment (`cleftgrid[FA->num_grd++]`, first-touch),
and that index **is GA gene 0** (`ic_bounds.cpp:15` sets `index_max = num_grd-1` under
`locclf`, set unconditionally at `top.cpp:1742`; `ic2cf.cpp:100-105` does
`grd_idx = (uint)icv[i]` straight into `cleftgrid[grd_idx]`). Same seed → different 3D
anchors → multi-ångström swings.

**Fixed by #403** (`dfb134f5`): per-iteration buckets concatenated ascending after the
parallel region, reproducing the serial order bit-exactly at any thread count and schedule.

**Why your probes looked clean and the campaign didn't:** the divergence needs CPU
contention to expose. Idle box → thread arrival order repeats → looks deterministic.
`R=10` restarts supply exactly that contention (restarts are separate forked processes,
`DatasetRunner.cpp:6437`, not nested OpenMP).

### Evidence

| check | result |
|---|---|
| pre-fix, 4 threads ×6 **under load**, 1GPK | **2 distinct grids** |
| pre-fix, 4 threads ×6 **under load**, 1G9V | **4 distinct grids** |
| pre-fix, 1 thread vs 4 threads | different grids (`aa1a261f` vs `2ca980fc`) |
| post-fix, all 5 canonical targets | 1-thread coords **bit-exact vs pre-merge main**; 4-thread ≡ 1-thread |
| post-fix, under load ×6 | **1 grid, 1 pose set** |
| `test_cleft_cavity` | 16/16 PASS |

Targets verified: 1GPK, 1G9V, 1HNN, 1GM8, 1HP0.

**Single-thread numbers are unchanged, bit-exact. Multi-thread numbers now converge onto
the single-thread answer.** So `--omp-threads 3` is finally comparable to the 1-thread canon.

### REFUTED — do not spend time here
- Static scheduling on the four `gaboom.cpp` pragmas (the original fix request). Loop bodies
  write per-index → order-independent. With the grid held fixed, the GA is bit-reproducible.
- Nested OpenMP across restarts (they're separate processes).
- Restart-timeout pose-pool truncation.
- Unstable MIF energy sort.
- `FLEXAIDDS_CLEFT_SORT` as *the* fix — it imposes a **third** order, neither serial nor
  thread-arrival, and changes 1-thread poses. **Leave it OFF.**

---

## 2. Launch config — exact diff from the killed `astex84_cleftsort_20260809_124322/run.sh`

```diff
- export FLEXAIDDS_CLEFT_SORT=1        # REMOVE: third canonical order, changes 1-thread poses
+ export FLEXAIDDS_SOFTBETA_ELECTION=1 # ADD: ΔG̃ = ΔH̃ − TΔS̃ mode ranking (LP's decision)
```
Plus: **restage `$OUT/bin/`** from a build of main `f5303093` or later — the harness runs its
own copied binaries, so rebuilding in the repo does not reach it.

Everything else in that script stays: `FLEXAID_SEED=12345`, `OMP_NUM_THREADS=3`,
`FLEXAIDDS_RESTARTS=10`, `--threads 2 --omp-threads 3`, 1000×1000, the input-integrity
shasum bracket (keep it — it is the check the first campaign lacked).

### On `FLEXAIDDS_SOFTBETA_ELECTION=1`
Default election is `min finite head CF` — **purely enthalpic**, the pre-2017 baseline.
The entropy-aware ranking is implemented exactly per LP's ISMB/ECCB 2017 (3Dsig) slide 12
in `LIB/SoftBetaFreeEnergy.h:5-14`:
`Z = Σ e^(−CF_i/T)`, `P_i = e^(−CF_i/T)/Z`, `H̃ = Σ P_i·CF_i`, `S̃ = −Σ P_i ln P_i`,
`G̃ = H̃ − T·S̃ = E_min − T·ln Z`. Member CFs come from the `.mcf` sidecar
(`cluster.cpp:590-610`). `T_soft` is **dimensionless CF arbitrary units — not Kelvin,
not 1/k_BT**.

---

## 3. What is invalidated

**Every pre-#403 multi-threaded number**, including the 39.3% (v2) and the 43.5% baseline.
Those were draws from a permuted-grid distribution, not measurements. Do not quote either.
The relaunch establishes a **new baseline from scratch** — there is nothing to preserve
continuity with.

**#405 is NOT a break in the chain.** Initially I reported it as a silent science change;
that was wrong and is retracted. Docked **coordinates are byte-identical** across
`ea044869` → `aff69d2f` → main on both targets (1GPK `f0bfe158`, 1G9V `b99f51b5`). The whole
difference is six added `REMARK thermo_*` lines. I had hashed whole `.pdb` files including
headers. All seven statmech commits are **pose-inert**: five touch no hot-path file, the two
that do change only print labels and REMARKs.

> **Lesson for your own comparisons: hash `^(ATOM|HETATM)` lines only.** Never whole pose
> files. Also exclude `wall_time_s`, any `*_pose_path` column, and `posebusters_input_sha256`
> (it embeds an absolute path on line 1).

---

## 4. Known caveats — read before interpreting results

1. **`1cb0cf9a` is UNMERGED** (on `fix/1hp0-cf-divergence`): *"1HP0 CF=-1700 phantom —
   reconcile emitted CF.app + .mcf with recomputed score."* The ΔG̃ ranking **consumes
   `.mcf` member CFs**. If emitted `.mcf` values do not reconcile with recomputed scores,
   the mode ranking is built on unreliable energies. The large-magnitude case is fixed in
   main (`8e7d17fd`, kills 73M CF.app), but this reconciliation is not.
   **Recommend landing it before drawing conclusions about ΔG̃ ranking.**
2. **Residual nondeterminism downstream.** Post-#403 the `.rrd` cluster-label column still
   differs run-to-run (two identical runs differed on 186 lines). **Elected poses were stable**
   across all such runs. Do not build any metric on `.rrd` column 2.
3. **Statmech fixes are unvalidated on Astex.** The Shannon fence fires on ~100% of
   heavy-tailed populations. Pose-inert, so it cannot move RMSD — but any *thermodynamic*
   number is on newly-changed estimators. First campaign **establishes**, does not confirm.
4. Wall cap **is** in main (`WAL_CONTACT_CAP=50`, `vcfunction.cpp:585`) — election keys are
   not corrupted by uncapped wall energy.
5. `success` column means "docking ran", **not** RMSD<2. Compute the rate yourself from
   `rmsd_to_crystal`, symmetry-corrected.

---

## 5. Proposed next experiment (LP's design — NOT implemented)

Two-level election: **rank modes by ΔG̃**, then compare **two representative rules within
the winning mode**.

⚠️ **As literally stated this is a null experiment.** `P_i = e^(−CF_i/T)/Z` is monotonically
decreasing in `CF_i`, so **argmax P_i ≡ argmin CF_i** — "most probable pose" and "min-CF pose"
are the same pose, by construction, on every target.

For a meaningful contrast the ensemble representative must use the *distribution*, not its
argmax. Recommended: **Boltzmann-weighted medoid** = the member minimising `Σ_j P_j·RMSD(i,j)`.
It is the geometric analogue of `H̃ = Σ P_i·CF_i`, always a real sampled pose, and is exactly
the "centre of the wide well" from slide 11. Ingredients already exist in `cluster.cpp`:
`coord_cache` (member geometries) + `.mcf` (member CFs → `P_i`).

**Run both arms in ONE campaign via dual emission.** The CSV already pairs
`cf_top1_*` with `entropy_top1_*` — extend that pattern. Identical grid, GA trajectory, mode
ranking and elected mode; the only difference is which member is written. Any RMSD delta is
attributable to the representative rule alone. Two separate campaigns cannot give that cleanly.

Precedent: `BindingMode::elect_Representative(bool useOPTICSordering)` (BindingMode.h:110)
already offers lowest-CF vs density-centre, and BindingMode.cpp:582-586 records that emitting
lowest-CF *"made Fast OPTICS output indistinguishable from the CF algorithm"* — the null result
above, already hit once. That path is FOPTICS-only; the default path (`cluster.cpp`,
`Clus_TOP[j] = j` over a CF-sorted population) has **no representative knob at all**.

---

## 6. Reproduce / verify anything here

```bash
determinism_cleft/repro_cleft_determinism.sh 1GPK     # or 1G9V (diverges harder)
```
Prints distinct grid/pose counts for gate OFF vs ON, under load, at 1 and 4 threads.
Full write-up: `determinism_cleft/FINDING_cleft_grid_nondeterminism.md`.

## 7. Open PRs
- **#406** — per-generation any-pose RMSD trace (`FLEXAIDDS_GENTRACE`), default off, verified
  non-perturbing (poses bit-identical on/off). Pending tsan. Answers "was a near-native pose
  ever sampled and then lost" — relevant to the 14 election-loss targets. Not a blocker.
