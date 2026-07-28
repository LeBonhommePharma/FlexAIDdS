# AUDIT — Theme A: SCORING PHYSICS
FlexAIDdS main @ d623c45ea · auditor: adversarial code review (scoring physics)
Patch: theme_A_scoring.patch (8 commits) · verified against live merged tree.

## Bottom line
**minor_issues.** Every genuine scoring change in this theme is env-gated and
default-OFF; with no env vars set, the default CF/ranking path is byte-identical
to before. The DEFAULT-OFF / bit-identity invariant HOLDS across all 8 commits.

The one substantive finding is a **correctness-of-claim defect in the thermo
impossibility gate (7a44c7035)**: the `+1000` sentinel is advertised (commit
msg + header comments) as making "downstream clustering never select impossible
poses rank-0", but **nothing reads it** — `dG_eff` / `thermo_impossible` /
`n_impossible_poses` are consumed ONLY by `printf` in gaboom.cpp. Ranking is
`QuickSort(CF)` at gaboom.cpp:1238, which runs BEFORE `compute()` at 1308; and
`compute()` runs only under `thermo_engine_enabled` (default false). The gate is
inert w.r.t. pose selection. This is DEFAULT-SAFE (no bit-identity breach) but
the feature does not do what it claims. Medium.

---

## KEY SCRUTINY answers (as posed in the task)

**Q1: Is the thermo `dG_eff=+1000` sentinel in the DEFAULT scoring path or gated?**
GATED, three times over, and even when fully enabled it never touches scoring:
- `compute()` (ThermodynamicEngine) is called from gaboom.cpp:1308 ONLY inside
  `if (FA->thermo_engine_enabled && FA->thermo_engine != nullptr)` (1241).
  `thermo_engine_enabled` defaults to **false** (config_parser.cpp:355;
  flexaid.h:678 "default false; zero cost when off"). Default run never even
  computes `dG_eff`.
- The gate body inside compute() is further behind `if
  (flexaids::thermo_score_enabled())` (ThermodynamicEngine.cpp:164), i.e.
  `FLEXAIDDS_THERMO_SCORE` truthy (ProtocolConfig.cpp:309, default OFF).
- **Even with both flags on, the sentinel changes only a printed number.**
  Repo-wide grep for consumers of `dG_eff`/`thermo_impossible`/`apply_gate`
  outside the engine internals returns printf lines only. QuickSort(CF) at
  gaboom.cpp:1238 has already fixed the ranking before compute() runs at 1308.
  => The commit-message claim "clustering can never rank it 0" is FALSE in the
  merged code. No selection path reads the gated field.

**Q2: Is the P1 pocket-presence penalty default-on or gated?**
GATED, default-OFF. `FA->pb_pocket_weight = 0.0` hard-set in top.cpp:539 and
defaulted to 0.0 in config_parser.cpp; `pb_pocket_enabled(w) = (w>0.0)`
(ProtocolConfig.h). The penalty block (vcfunction.cpp:1204) is entered only when
`pb_pocket_enabled` is true. The enclosing PB block (981) is entered only when
`pb_clash_weight>0 || pb_pocket_weight>0` — both zero by default => whole block
skipped => bit-identical.

---

## Per-commit

### 7a44c7035 — thermo impossibility gate (dH>0 & dS<0 -> dG_eff=+1000)
- **What:** adds `thermo_gate::is_impossible/apply_gate` (header-only), a
  Boltzmann-population `dG_eff = <CF> - T*H` diagnostic, and a gate that forces
  `dG_eff` to +1000 when `mean_CF>0 & TdS_vib<0`. All under
  `thermo_engine_enabled` + `FLEXAIDDS_THERMO_SCORE`.
- **default_behavior_changed:** NO. compute() not called in default (engine off);
  even when on, only reporting.
- **correctness:** SUSPECT. The primitive is fine and unit-tested (7/7), and the
  ΔS-source reasoning (use TdS_vib, not Shannon H≥0, else dead code) is sound.
  BUT the advertised integration does not exist: (a) no ranking/clustering
  consumer reads `dG_eff`/`thermo_impossible` — printf-only (gaboom.cpp:1341,
  1355); (b) selection is QuickSort(CF) at gaboom.cpp:1238, BEFORE compute() at
  1308. Secondary physics concerns even if it WERE wired in: the per-pose test
  uses `dH_i = CF_i` (contact-fitness proxy, not a calibrated enthalpy — sign of
  CF is not sign of ΔH), and a SINGLE complex-level `TdS_vib` is applied to every
  pose (r.gate_dS_used is one scalar), so "per-pose impossibility" is really
  "one dS gate broadcast across poses".
- **severity:** MEDIUM (false capability claim + misleading in-code comments an
  operator would trust; not high — defaults byte-identical, physics primitive
  correct).
- **verdict:** NEEDS_CHANGE. Either wire the gate/dG_eff into an actual
  selection consumer, or correct the commit message + ThermodynamicEngine.h
  comments to state "diagnostic only; does NOT affect pose selection."

### d5c33b349 — FLEXAIDDS_COM_FLOOR soft lower clamp on CF.com
- **What:** softplus soft-floor on `FA->optres[j].cf.com` at -F, applied in
  vcfunction.cpp:884 after VCT_NORM.
- **default_behavior_changed:** NO. Pure `getenv("FLEXAIDDS_COM_FLOOR")` guard
  with inner `if(F>0.0)`; unset/≤0 => block skipped => bit-identical.
- **correctness:** SOUND (math). Verified the form: softfloor(x)=-F+F*softplus((x+F)/F)
  with numerically-stable softplus `(z>0?z:0)+log1p(exp(-|z|))`. x»-F -> x
  (near-identity); x->-inf -> -F (bounded); derivative in (0,1] (monotone in the
  com channel). Implementation matches. Caveat: monotone in com ONLY, not in the
  summed CF — it rebalances channels (that is the intended "enabler" effect), and
  the commit itself flags the functional form as RECONSTRUCTED from a summary
  ("Confirm F and the exact squashing against the original work order").
- **severity:** LOW (default-off + author-flagged unverified form).
- **verdict:** MAKES_SENSE (gated); confirm form vs work order before any canary.

### dd292d6dd — P1 pocket-presence soft penalty
- **What:** `pb_pocket_weight*(d-radius)^2` on ligand-centroid-to-nearest-heavy-
  receptor-atom distance, into cf.pb_clash; also fixes config_parser env-override
  ordering so FLEXAIDDS_PB_CLASH_* aren't clobbered by JSON defaults.
- **default_behavior_changed:** NO. weight defaults 0.0 (top.cpp:539); gated by
  pb_pocket_enabled; enclosing PB block skipped when both weights 0. The
  config_parser env-override addition only overrides when the env var is set.
- **correctness:** SOUND. Quadratic ramp is zero inside radius (correct sign,
  dimensionally fine as a CF-unit penalty). Chebyshev-shell nearest-atom search
  with `(R-1)*cell)^2 >= best_d2` early-out is a correct branch-and-bound over the
  hoisted cell list; rec_heavy precomputed and rebuilt atomically with the other
  parallel arrays (vcfunction.cpp:1002-1029). Attribution to first ligand optres
  matches pb_clash/GIST convention.
- **severity:** NONE (for defaults).
- **verdict:** MAKES_SENSE.

### bb8e166e5 — carve out coordinating-metal pairs from pb_clash
- **What:** skip clash pair when either partner is Zn/Fe/Mg/Ca/Mn/Co/Ni/Cu
  (Na/K excluded); receptor flags precomputed in cache.
- **default_behavior_changed:** NO. `pb_metal_carveout` = getenv truthy, default
  OFF; when off both metal flags are hard-false, skip branch never taken, visited
  pair set unchanged. Also only reachable inside pb_clash_weight>0 (non-default).
- **correctness:** SOUND. Physics is right (1.9-2.3 Å coordination bonds sit
  inside 0.75*(vdw_i+vdw_j) and were being scored as clashes vs the metal_coord
  Morse reward). Na/K exclusion is defensible (spectator ions).
- **severity:** NONE.
- **verdict:** MAKES_SENSE.

### 88df9bf96 — precomputed pb_vdw_radius instead of per-eval lookup
- **What:** optionally read atoms[].pb_vdw_radius instead of
  posebusters_vdw_radius(get_element(type)).
- **default_behavior_changed:** NO. `pb_vdw_cached` getenv, default OFF; both
  ligand+receptor sides switched together so a pair never mixes sources.
- **correctness:** SOUND, and honestly gated: authors correctly identified this
  is NOT a pure perf swap (the two sources diverge for some atoms) and gated it +
  wrote a parity test rather than shipping it as default. (The commit's *initial*
  divergence attribution — "type 39/H diverges" — was WRONG, and is corrected in
  074612fba; the CODE was always correct, only the explanatory comment/test was.)
- **severity:** NONE (default-off).
- **verdict:** MAKES_SENSE.

### c856dcdf6 — auto-couple grid-hoist to pb_clash_weight>0 (remove opt-in flag)
- **What:** removes FLEXAIDDS_PB_CLASH_GRID_HOIST; hoist now unconditional when
  the PB block runs. `use_cache = pb_cache.valid`.
- **default_behavior_changed:** NO — THIS IS THE ONE TO SCRUTINIZE, and it is
  clean. The hoist lives ENTIRELY inside the `pb_clash_weight>0 ||
  pb_pocket_weight>0` block (vcfunction.cpp:981). With pb_clash OFF (default) the
  block is never entered, so removing the flag cannot change the pb_clash-OFF
  path. When pb_clash IS on, the hoist only moves grid *construction* out of the
  eval loop; CF math is identical (commit shows single-eval bit-identity across
  legacy/env-hoist/auto-hoist). This is the correct fix for the v134 lesson
  (loop-invariant receptor grid must not be rebuilt per eval).
- **correctness:** SOUND. Cache keyed on atm_cnt+ratio, thread_local, invalidated
  via matches().
- **severity:** NONE.
- **verdict:** MAKES_SENSE.

### 83c8f045b — never diverge silently on GPU/accelerated path
- **What:** warn-once when CUDA/METAL backend + pb_clash_weight>0 (GPU zeroes
  cf.pb_clash); adds FLEXAIDDS_FORCE_CPU to pin CPU backend.
- **default_behavior_changed:** NO. pb_clash_weight default 0 and FORCE_CPU unset
  => neither warning nor backend override fires; select_backend() unchanged.
- **correctness:** SOUND. This is a correctness *guard* — it surfaces a real
  physics divergence that was previously silent. Good.
- **severity:** NONE.
- **verdict:** MAKES_SENSE.

### 074612fba — correct the pb_vdw divergent set to I/Na/K + posebusters CLI autodetect
- **What:** rewrites test_pb_vdw_parity.cpp against real reader tables (true
  divergent set = I, Na, K; H does NOT diverge); adds find_program(POSEBUSTERS_BIN)
  baked-in fallback.
- **default_behavior_changed:** NO. FLEXAIDDS_PB_VDW_CACHED stays default-OFF; CLI
  autodetect only affects the optional bust bridge, gated FLEXAIDDS_POSEBUSTERS_BIN
  define, native NativePoseQC path unchanged when absent.
- **correctness:** SOUND. This is a test/comment correction (the earlier
  "H diverges" reading was wrong); the scoring code is untouched. CLI autodetect
  is provenance-clean (optional external, no new build dep).
- **severity:** NONE.
- **verdict:** MAKES_SENSE.

---

## Cross-cutting notes (not per-commit findings)
- **Reproducibility (invariant #3):** c856dcdf6 documents that GA end-to-end
  output is NOT run-to-run reproducible in this build (best_CF -77.1284 vs
  -77.1467 same binary). This predates the theme, but it is why every bit-identity
  claim here rests on SINGLE-EVAL comparison, not GA-level. Latent invariant-#3
  risk for the engine; flag for the reproducibility auditor.
- **cf.pb_clash overloading:** both pb_clash and pb_pocket accumulate into
  cf.pb_clash. Correct for attribution, but the [PB_CLASH] name now covers two
  physically different terms — a reporting-clarity nit, not a bug.

## Highest-severity finding
7a44c7035 thermo impossibility gate: the `+1000` sentinel and `dG_eff` are
advertised as deranking thermodynamically-impossible poses in clustering, but no
selection/clustering path reads them (printf-only), and ranking (QuickSort on CF,
gaboom.cpp:1238) is finalized before compute() runs (1308). The feature is inert
w.r.t. pose selection and its in-code comments misdescribe the data flow.
DEFAULT-SAFE (byte-identical defaults) but functionally a no-op that claims
otherwise. Medium — needs the wiring done or the claim retracted.
