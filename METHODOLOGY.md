# METHODOLOGY.md — FlexAIDdS canonical benchmarking, determinism & CI methodology

**Single source of truth for HOW work is measured, validated, and gated — for every agent
(Claude, Claude Code, Codex, Grok, GPT) and every human.** `AGENTS.md` governs coding conduct;
this file governs *methodology*. When a skill, script, or agent describes a benchmark, a
determinism check, or a merge gate, it must **defer to this file**, not restate it. If a
procedure here changes, it changes here first and everything else re-reads it. The goal is
**agent-independence**: the same task run by any agent produces the same result because all of
them execute the identical procedure below.

**Comparative / baseline science entry point (A/B/C, genuine rates, pipeline):**  
`docs/implementation/COMPARATIVE_SCIENCE_README.md` (does not replace this file’s gates).

Status of this document is authoritative as of the commit that carries it. Do not fork the
numbers into other files — reference `METHODOLOGY.md §N`.

---

## 0. Environment invariants (MUST hold on every run)

- **cmake:** `/opt/homebrew/bin/cmake` — NOT on PATH. Either `export PATH="/opt/homebrew/bin:$PATH"`
  or call the full path. Same for `ctest`.
- **Build:** `cd build && /opt/homebrew/bin/cmake --build . --target <FlexAIDdS|benchmark_datasets> -j4`
- **Energy matrix:** ALWAYS use `MC_st0r5.2_6.dat`. Set `FLEXAIDDS_DATA_DIR=$PWD/build` (or the dir
  that contains the matrix). Never let provenance fall back to the empty `<bin>/../WRK` path — an
  empty matrix silently corrupts scoring. Verify the run log shows the matrix loaded.
- **Deterministic seed:** `FLEXAID_SEED=12345`. This is the ONLY determinism seed. Do **not** use
  `FLEXAIDDS_SEED_BASE` for determinism (it offsets per-restart seeds).
- **Benchmark GA budget:** **2000 generations, population 1000** (the FlexAID/FlexAIDdS norm =
  2,000,000 evals/restart). Do not change for accuracy runs.
- **Restarts:** `FLEXAIDDS_RESTARTS=<n>`. Published Astex protocol = 10 restarts; a fast A/B may
  use 1–3 but MUST state it.
- **Security channel off for benchmarks:** `FLEXAIDDS_NO_SEC=1`.
- **Thread rule:** workers × omp-threads ≤ physical P-cores. On this host use `--threads 1
  --omp-threads 4` (or 6). Never oversubscribe.
- **RMSD reporting:** report **rank-0 elected pose RMSD** (the elected `_0.pdb`), heavy-atom,
  2.0 Å cutoff, **in-place in the receptor frame — never superposed**. NEVER report seed-elitism
  / `_INI.pdb` RMSD as the result. In-place is what the engine does (`LIB/calc_rmsd.cpp:92-102`,
  and the same in the original FlexAID); a superposed value measures shape, not placement, and
  is not the quantity the 2.0 Å criterion is defined on.
- **RMSD engine — read this before quoting a number.** There are **FIVE** RMSD implementations
  in this tree and they are not interchangeable. An unlabelled RMSD is not reportable.
  **Two of the five are C++ and both are live on the campaign path** — see the warning after
  the list:
  1. **In-repo metric** — `python/flexaidds/benchmark.py::compute_rmsd`, used by
     `dataset_runner` and by **every CI tier**. This is the gate. In-place since #354. Ligand
     selected by residue against the reference atom count since #363/#366 — *not* by unioning
     every non-water HETATM, which merges cofactors into the ligand. **NOT symmetry-corrected as of this
     commit.** #365 adds element-blocked assignment but is still open; until it merges this
     metric pairs atoms positionally and systematically OVERSTATES error on symmetric ligands
     (measured on real 1gpk poses: up to 1.66 Å, enough to move a pose across the 2.0 Å bar).
     Do not describe a gate number as symmetry-corrected before #365 lands.
  2. **Offline reference scorer** — `benchmarks/astex_repro/score_reference.py`: spyrmsd
     graph-isomorphism, heavy atoms, pose selection on PDB serial ≥ 90000, element-blocked
     Hungarian fallback only when spyrmsd raises (logged). Treat this as the strongest
     instrument.
  3. **Offline permissive scorer** — `benchmarks/astex_repro/score_offline.py`: element-blocked
     Hungarian, no graph isomorphism. **The repo's own audit records it as over-permissive:
     1HP0 reads success under `score_offline.py` and failure under `score_reference.py`**
     (`docs/audit/26h-swarm/9971dff7e.md`). Do not quote it as a docking-power number.
  - **Neither offline scorer is wired into anything.** Searched: `.github/workflows/`,
    `benchmarks/` (including `run.py`), and `python/`. The only invocation recorded anywhere is
    a manual one — `benchmarks/astex_repro/MONITORING.md:8` says "SCORING: python3
    score_reference.py". So the strongest instrument is the one nothing runs, and the gate runs
    the one no document described until this section.
  - `score_offline.py`'s partial `poster_metric_results.csv` is still in-tree beside
    `score_reference.py`'s. The audit's standing recommendation is to deprecate the permissive
    one; until that happens, **check which CSV you are reading.**
  4. **Engine Hungarian** — `LIB/calc_rmsd.cpp::calc_Hungarian_RMSD`, its own assignment
     implementation. Called via `calc_rmsd(..., Hungarian)` from `BindingMode.cpp:779,785,922,928`.
     **This is what writes the `REMARK RMSD` line inside every pose PDB.**
  5. **DatasetRunner Hungarian** — `LIB/DatasetRunner.cpp:436 dataset::hungarian_rmsd`, a
     SEPARATE implementation with its own solver (`munkres_solve`, `:397`). Feeds
     `compute_pose_ligand_rmsd` / `pose_pose_rmsd`. **This is what `result.csv` carries.**
  - 🔴 **4 and 5 are two independent implementations of the same algorithm, each with its own
    test file (`tests/test_hungarian_rmsd_bounds.cpp` and `tests/test_dataset_runner.cpp`), and
    NO test compares them to each other.**
    They also partition atoms differently — 4 blocks by SYBYL **type**, 5 by **element** — so
    they are not required to agree, and neither dominates the other: type refines element for
    typed atoms but COARSENS it for atoms sharing one type across elements (see the DUMMY and
    aliasing notes below). Two "correct" numbers can differ in either direction. A pose PDB can therefore carry a `REMARK RMSD` from
    one while the `result.csv` row beside it carries a number from the other. Until a
    cross-check exists, **never mix a REMARK RMSD and a CSV RMSD in the same table**, and say
    which file a quoted number came from.

- **Known element-typing limitations (affect scoring, not only RMSD).** `atom.type` is used
  directly as a row/column index into the VCT energy matrix, so a mistyped atom is *scored*
  against the wrong element, and it also changes which atoms the engine's Hungarian may swap:
  - **Unknown elements → DUMMY (39).** `SdfReader.cpp:88` / `read_lig.cpp:218`. DUMMY is one
    bucket spanning every element, so the engine may pair N with O and report a *flattered*
    RMSD. `check_types_assigned` (`TargetValidation.cpp:145`) would detect this but **nothing
    calls it**, and `read_lig.cpp:227`'s stderr warning is discarded because the run succeeds.
  - **Deliberate aliasing.** `I → 25 (Br)` and `Se → 18 (S.3)`, identically in all THREE
    typing sites — `SdfReader.cpp:78,83`, `Mol2Reader.cpp:78,83` and the receptor path
    `read_coor.cpp:63,64`. The canonical table reserves `26=I` and `27=Se`; **no site emits
    either**, so the substitution is total and symmetric across ligand and receptor. It is a
    modelling choice, not a bug — the I and Se matrix rows are near-empty. But an iodine is then indistinguishable from a bromine downstream:
    no diagnostic can catch it, and I/Br in one ligand share a type bucket.
  - **Scope, checked:** 187 `*.sdf`/`*.mol2` under `benchmarks/` contain **no** I, Se, B, Si,
    As, Pt or Ru, and every element present maps to a real VCT row (`normalize_element` fixes
    the `CL`→`Cl` case before the lookup). **These are reachable trapdoors with no current
    input that opens them** — they become live the day a halogen series or a metal-organic
    ligand is benchmarked, and they will be silent when they do.

- **Why a reported number needs a discrepancy, not an agreement.** Every defect found in the
  RMSD work of 2026-08-01 came from someone declining to round off a mismatch: a suite count
  that differed by nine, a pose count that differed by one, an inequality that held in one
  direction and not the other, a `strcmp` against an input nobody had traced. **None came from
  a review that agreed.** When two people's numbers differ, the difference is the finding —
  reconcile it in the open before either number is published.

---

## 0.1 The CI gate and the campaign path are NOT the same experiment

Two harnesses invoke the same engine with different configuration. Neither is wrong; they are
simply different runs, and a number from one does not transfer to the other.

| | **CI / tier1 gate** | **Campaign** |
|---|---|---|
| entry point | `python -m benchmarks.run` (`.github/workflows/benchmark-tier1.yml`) | `benchmark_datasets` → `DatasetRunner` |
| config file | **none** — the engine takes its compiled-in defaults | writes `dock_config.json` |
| `permeability` | 1.0 (`top.cpp:601`) | 0.9 |
| `normalize_area` | 0 (`top.cpp:585`) | true |
| `intermolecular_clash_ratio` | 0.0 | 0.75 |
| `coarse_init.enabled` | **OFF** | ON (hardcoded, `DatasetRunner.cpp:6036`) |
| `mif_enabled` | **OFF** | ON (hardcoded, `DatasetRunner.cpp:6012`) |
| retained poses | 10 (`ga_constants.h:16` `GA_DEFAULT_NUM_PRINT`, applied `gaboom.cpp:315`) | 50 per restart (`DatasetRunner.cpp:86` `kBenchmarkPoseLimit`, emitted as `max_results` at `:5959`) |

**The consequence that matters:** with `mif_enabled = 0` and `grid_prio_percent = 100.0`, the
guard at `top.cpp:1836` makes `initialize_direct_mif` return immediately. **The gate docks
without ever building its pocket field.** Every green tier1 tick to date passed under that
configuration.

Rules that follow:

- **Never compare an RMSD from one path against an RMSD from the other** without stating the
  divergences above. On 1mq6 the two paths differ by 7.5 Å and the cause is not yet isolated —
  four of the six divergences can move a centroid.
- A cross-path result is not evidence until the divergence responsible has been ablated.
- Changing either harness's defaults is a methodology change: it lands here first.

## 0.2 Provenance — what a receipt must capture

Three tiers, by how a parameter can be recovered after the fact:

| tier | example | recoverable from |
|---|---|---|
| **FIXED** | anything written to `dock_config.json` | the artifact |
| **VARIABLE** | the COM floor, via `CF.com` in the pose | the poses themselves |
| **LOST** | `permeability`, `pb_pocket_weight`, `pb_clash_weight` | **nothing** |

**The rule, in one line: anything that goes through `getenv` and nowhere else is gone the moment
the shell exits.** Those three are read from the environment and never written to any config or
receipt, so a completed run cannot be told apart from one at different weights.

- **Every run MUST emit its `getenv`-only scoring environment beside its results** (sidecar or
  `RUN_RECEIPT`). This blocked two separate investigations that had the artifacts in hand.
- The two harnesses currently capture **disjoint** provenance fields, so cross-path comparison
  has no common provenance even where both wrote something.
- 🔴 **`_sweep_config.json` records one field that is not a setting.** The tier-1 workflow
  offers a `tier1_subset_size` dispatch input (default `"5"`) and writes it into the provenance
  sidecar as `tier1_subset_size_input` (`benchmark-tier1.yml:225`). **Nothing consumes it** —
  `benchmarks.run` has no subset flag, and the size actually used comes from the dataset YAML
  (`astex_diverse.yaml:35` → `2`, read at `runner.py:252`). So dispatching a run with `5`
  docks 2 targets and writes a sidecar claiming 5. The sidecar's own comment says it records
  "exactly which switches this run used"; this field records a switch that was not used, which
  is worse than not recording it. Until it is wired or removed, **do not read
  `tier1_subset_size_input` as provenance** — count the entries in the run log instead.
- `ops/reference_config.env` pins `FLEXAIDDS_PARALLEL_RESTARTS=0` but **not**
  `FLEXAIDDS_RESTARTS` — a campaign run against it silently takes the compiled-in default of 5,
  not the published Astex 10. Pin it explicitly in any run you intend to publish.

---

## 1. Reproducibility / parity gate (run before ANY merge)

Purpose: prove a change is bit-identical to the baseline under default flags.

1. Build the candidate engine and the baseline engine (main) separately; record both md5s.
2. Dock 1G9V, `FLEXAID_SEED=12345`, `OMP_NUM_THREADS=1`, config `/tmp/parity.json`
   (2000 gen / pop 1000, no crystal-pose seed: `pose_seed_enabled=false, seed_fraction=0`).
3. Assert: elected CF equal AND all 10 elected poses byte-identical between candidate and main.
4. PASS = default-flag behavior unchanged. Any intended behavior change must be opt-in behind an
   env flag that defaults OFF, and parity must hold with the flag OFF.

Reference baseline engine md5 (main @ 7f1f10a0…): `7f1f10a0f10b682b33a76622a40f1a60`.

---

## 2. Determinism check (multi-thread)

For changes touching parallel regions (GA eval, cleft detection):

- **Cleft grid:** dock a fixed seed at `OMP_NUM_THREADS=1` and `=4`, twice each; assert the emitted
  `.rrg` grid-cache file is byte-identical across thread counts AND run-to-run.
- **GA population:** with the parallel-reproduce flag ON, assert the gen-0 order-independent CF
  checksum (Σ cf.com, Σ cf.wal over the population) is identical run-to-run at 1 AND 4 threads, and
  that all 10 elected poses are byte-identical across two 4-thread runs.

---

## 3. Astex-85 accuracy A/B (the science gate)

Purpose: prove a determinism/perf change does not regress docking accuracy.

- **Protocol:** autonomous (blind) mode — `--mode autonomous` — which exercises SURFNET cleft
  detection (top.cpp "Always run SURFNET"). Seed 12345, `FLEXAIDDS_NO_SEC=1`, 2000 gen / pop 1000.
- **Engine selection:** `FLEXAIDDS_BINARY=<engine>`; `benchmark_datasets --benchmark astex
  --only-codes <list>` subsets targets. (Direct-CLI per-target is acceptable when the harness's
  oracle-site guard blocks a blind run — document which was used.)
- **Metric:** top-1 rank-0 RMSD, 2.0 Å, in-place. Name the instrument (§0: in-repo metric vs
  `score_reference.py`) — an unlabelled RMSD is not reportable.
- **Acceptance:** candidate within noise of baseline; **no target flips success→fail** attributable
  to the change. A full landing decision uses the full 85 × 10-restart protocol; a fast pre-check
  may use a documented subset.
- **Published anchors (for context, not gates):** Gaudreault & Najmanovich 2015 JCIM **Table 2**
  Astex native FLRP: **top-1 = 45.2%**, **top-10 = 66.7%** (do not swap labels). Morency 2017
  3Dsig/poster FlexAID ~66% / entropy ~69% are **S_top10-family** medians — not the JCIM top-1
  figure. Report the protocol next to any number; never mix top-1, top-10, and S_top10 without labels.

### 3.1 Comparative three-arm goal (JCIM FlexAID vs first entropy vs FlexAIDdS)

For **fair comparison** of (A) JCIM-era CF FlexAID, (B) first entropy FlexAID, and (C) current
FlexAIDdS under frozen fairness axes, do **not** invent a parallel protocol here. Execute:

- **Goal design / phases / G1–G9 fulfillment:** `docs/implementation/COMPARATIVE_GOAL_METHODOLOGY.md`
- **Arm specs + confound checklist:** `docs/implementation/COMPARATIVE_BENCHMARK_METHODOLOGY.md`
- **Source commits / paths:** `docs/implementation/arm_pins.json`
- **Deck knobs / FO MinPts:** `docs/implementation/3dsig_red_pair_protocol.md`

Primary comparative statistic for that goal is **S_top10** (10k bootstrap median), not the §3
autonomous S1 gate used for merge-time accuracy A/B. Serial arms only; local-first + thin iCloud
mirror (`docs/ICLOUD_BENCHMARK_STORAGE.md`).

---

## 4. ctest

`cd build && /opt/homebrew/bin/ctest --output-on-failure` — expect **11/11**. A stale binary can
cause a false PoseBustTests failure; rebuild before trusting a red.

---

## 5. Branch / merge discipline (OPS owns the merge)

- SWE agents implement on a **feature branch** and push; they do **not** merge to `main`.
- OPS (Benchmarker/Monitor/CI-CD) runs §1–§4 and the §6 audit, then merges.
- Every commit is reviewed (see §6). No commit lands on `main` unreviewed.

---

## 6. Commit review & audit (OPS/CI)

Every commit by any agent is reviewed against: (a) §1 parity, (b) §2 determinism if it touches
parallel code, (c) §3 accuracy if it touches scoring/cleft/GA, (d) §4 ctest, (e) code-conduct in
`AGENTS.md`, (f) no leaked instrumentation / debug prints / commented-out dead code introduced.
The review is run through the strongest reviewer model available in the runtime and its verdict is
recorded. NOTE: model availability is runtime-dependent — if a specifically requested reviewer model
(e.g. Fable 5) is not present in `host.list_models()`, OPS uses the strongest available and says so;
it never silently relabels a different model as the requested one.

---

## 7. Model routing for SWE tasks (guidance)

Scientific subtlety → stronger model; mechanical wiring → lighter model.
- Deep concurrency / numerics / scoring physics → Opus 4.8 (or Fable 5 where available).
- Test wiring, mechanical refactors, doc plumbing → Sonnet 5.
State the model used in the commit/PR description so the run is reproducible.
