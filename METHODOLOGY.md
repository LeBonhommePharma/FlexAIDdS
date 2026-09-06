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
  a nominal population-generation product of 2,000,000 per restart). Do not change for accuracy runs.
  This product is not an observed evaluation count; throughput requires a witnessed evaluator-boundary counter.
- **Restarts:** `FLEXAIDDS_RESTARTS=<n>`. Published Astex protocol = 10 restarts; a fast A/B may
  use 1–3 but MUST state it.
- **Security channel off for benchmarks:** `FLEXAIDDS_NO_SEC=1`.
- **Thread rule:** workers × omp-threads ≤ physical P-cores. On this host use `--threads 1
  --omp-threads 4` (or 6). Never oversubscribe.
- **RMSD reporting:** report **rank-0 elected pose RMSD** (the elected `_0.pdb`), heavy-atom,
  2.0 Å cutoff, **in-place in the receptor frame — never superposed**. Success ⇔ rank-0 in-place
  RMSD **`<= 2.0 Å`**. NEVER report seed-elitism
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
     Do not describe a *gate* number as symmetry-corrected: this metric is still
     positional. #365 landed the symmetry correction on the CLAIM path only
     (`scripts/rmsd_symmcorr.py` → S1/S2 in `scripts/aggregate_claim_metrics.py`),
     not inside this in-repo metric and not inside the engine's `success_rmsd`.
  2. **Offline reference scorer** — `benchmarks/astex_repro/score_reference.py`: spyrmsd
     graph-isomorphism, heavy atoms, pose selection on PDB serial ≥ 90000, element-blocked
     Hungarian fallback only when spyrmsd raises (logged). Treat this as the strongest
     instrument.
  3. **Offline permissive scorer** — `benchmarks/astex_repro/score_offline.py`: element-blocked
     Hungarian, no graph isomorphism. **The repo's own audit records it as over-permissive:
     1HP0 reads success under `score_offline.py` and failure under `score_reference.py`**
     (`docs/audit/26h-swarm/9971dff7e.md`). Do not quote it as a docking-power number.
  - **The strongest instrument is now wired into the claim path (#365).** It used to be
    true that "the strongest instrument is the one nothing runs": neither offline scorer was
    invoked by `.github/workflows/`, `benchmarks/` (including `run.py`), or `python/`, and the
    only recorded invocation was manual (`benchmarks/astex_repro/MONITORING.md:8`).
    `scripts/rmsd_symmcorr.py` closes that gap. It is **not a sixth implementation** — it is
    method 2's invocation contract made callable (method 2 is a top-level script that executes
    on import, so it cannot be imported as a library), calling the same
    `spyrmsd.rmsd.symmrmsd` with the same arguments: crystal SDF bond block, heavy atoms,
    `center=False, minimize=False`, ligand selected on PDB serial ≥ 90000.
    `tests/test_rmsd_symmcorr.py` pins the contract. **The CI gate still runs method 1**, which
    remains positional — so a gate number and a claim number are still different quantities.
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
    test file (`tests/test_hungarian_rmsd_bounds.cpp` and `tests/test_dataset_runner.cpp`).**
    They group atoms differently (SYBYL type vs element), so a REMARK RMSD and a `result.csv`
    RMSD can legitimately differ. The five-way harness (`scripts/rmsd_five_way_crosscheck.py`
    plus `RmsdCrossCheck` / `RmsdClaimCutoff` in `tests/test_dataset_runner.cpp`) pins that
    **claim success is rank-0 in-place RMSD `<= 2.0 Å`** and compares 4 vs 5 on shared
    assignments. Still never mix an unlabelled REMARK RMSD and a CSV RMSD in the same table,
    and say which file a quoted number came from.

### 0.0 Which RMSD is the claim? (#365)

The three quantities are ordered by construction, not by accident:

```
hungarian  ≤  symmcorr  ≤  serial
```

* **serial** (`rmsd_to_crystal`, engine) pairs atoms by position — the identity mapping.
* **symmcorr** (`rmsd_symmcorr`, `scripts/rmsd_symmcorr.py`) minimises over the ligand's
  **graph automorphisms**. The identity mapping is one of them, so symmcorr can never exceed
  serial. Equality holds exactly when the identity mapping is already optimal.
* **hungarian** (`rmsd_hungarian`, engine) minimises over all **same-element bijections** — a
  *superset* of the automorphisms, including chemically invalid ones. So it can score *below*
  the true value: it is over-permissive, not merely different. Measured: it inflated the pool
  ceiling from 48.8% to 57.8%.

**The claim metric is symmcorr.** Two consequences follow directly from the ordering:

1. Correcting serial → symmcorr moves targets **fail → PASS only**. No banked number can
   shrink. Verified on the 84-target parent campaign: 76 elected poses scored, **0 direction
   violations**, 1 flip (1TZ8, 6.8360 → 1.0871 Å).
2. Where symmcorr is unavailable, falling back to serial is **safe but conservative** — it
   under-counts successes and never over-counts. Falling back to *hungarian* would not be safe,
   which is why no code path does.

**Still inheriting the serial definition:** the engine's `success_rmsd` gates on
`rmsd_to_crystal`, and `success_pb := success_rmsd ∧ pb_pass`, and `claim_ready` requires
`success_pb`. So the **STRICT headline is still a serial number** and is a conservative lower
bound. Repointing the engine gate is the remaining half of #365 and requires a rebuild.
`aggregate_claim_metrics.py` states this in its JSON report under
`strict_metric_inheritance`.

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
| **what a timed-out dock records** | **nothing.** `runner.py:1415` logs and `continue`s without touching the crash count or the exit-code map, so a timeout is indistinguishable **in the artifact** from "ran and found nothing". The only witness is one `logger.error("Docking timed out")` line in the CI job log. | **three places.** `wait_with_timeout` returns `-1` (`DatasetRunner.cpp:335`) → `docking_completed=false` (`:6856`) → every scoring stage skipped and `result.pb_failed_keys = "docking_incomplete"` written (`:7553-7555`), plus `[TIMEOUT]` in the captured `stderr.log`. |

🔴 **The last row is a divergence of a different kind.** The other rows are about what the engine
*computes*; that one is about what the harness *admits went wrong*. It matters because
`#326`'s liveness gate exists precisely to distinguish "the engine did not run" from "the engine
ran and produced nothing" — and on the gate path a timeout is the former reported as the latter.
The `OSError` handler three lines below (`runner.py:1418`) does record both, with a comment
explaining that omitting it would make the run "look like 'executed, 0 poses' (productivity)
rather than 'engine did not run' (liveness)." That reasoning applies verbatim to
`TimeoutExpired` and was not applied to it. (Found by Honey; campaign side traced by both of us.)

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
- `ops/reference_config.env` pins `FLEXAIDDS_PARALLEL_RESTARTS=0` but **not**
  `FLEXAIDDS_RESTARTS` — a campaign run against it silently takes the compiled-in default of 5,
  not the published Astex 10. Pin it explicitly in any run you intend to publish.

---

## 0.3 Admission identity, missingness, and repair evidence

The CSV aggregator validates receipt fields; it does not independently witness engine or
validator execution. Its report must identify that evidence level and must not describe
string equality as verification of raw artifacts. STRICT requires complete, consistent
protocol, matrix, finite serial RMSD, PoseBusters, score-pose consistency, tENCoM/Eigen,
and syntactically valid pose-linked hash fields. A stored `claim_ready` or success flag
cannot override contradictory or missing measurements. This hardening does not change
the RMSD instrument or pose election described in §0.0.

The frozen manifest is mandatory for a primary rate: validate its schema, unique codes,
declared count, and sorted-code digest. Never silently substitute admitted row count.
Each target contributes at most once. Existing one-row-per-target producer output is
explicitly a **single-observation** analysis. Repeated seeds require an explicit expected
seed list; a target passes only on a strict majority of that list, with absent/failed
seeds counting as failures. Reject duplicate observation identities and implicit mixing
of arms or endpoints; any explicit arm selection must be recorded. Do not infer a union
or majority across unrelated experiments. A diagnostic mode may expose legacy data but
must not emit a primary STRICT rate without this contract.

Reject duplicate CSV headers, inconsistent widths, ambiguous source layouts, and sidecar
joins lacking matching pose identity. Report S1/S2 diagnostics over their declared eligible
population independently of filtering on STRICT success. Empty measurement populations
have unavailable statistics and an explicit valid count, never an invented zero RMSD.
Process completion and scientific success are separate: a wrapper must retain runtime
failure state and signal incomplete docking without calling a completed inaccurate pose
a process failure.

Memory-ownership repairs are checked against allocation and workspace invariants. Correct
paths must retain their outputs; paths that use invalid pointers are not a scientifically
valid parity baseline. Document those intentional invalid-path corrections and test the
production ownership transitions, including repeated and parallel workspace use. Do not
hide an invalid-pointer repair behind an option that leaves undefined behavior as default.
Changes to valid scoring, ranking, clustering, or thermodynamic models still require the
feature flags and science gates below. No claim of historical crash attribution follows
from a source fix alone.

## 1. Reproducibility / parity gate (run before ANY merge)

Purpose: compare the baseline and candidate's emitted scientific results under default
flags, with exact input and executable provenance. Publication acceptance is separate.

1. Build main and candidate in separate pinned checkouts using the same compiler/options.
   Preserve source commits and any diagnostic-only patch, compiler commands, binary MD5
   and SHA-256, runtime-data hashes, input hashes and exact argv/environment.
2. Dock 1G9V, `FLEXAID_SEED=12345`, `OMP_NUM_THREADS=1`, population 1000 and generations
   2000, without crystal/native pose seeding (`pose_seed_enabled=false, seed_fraction=0`).
   Use fresh output directories keyed by run ID and source identity, never engine basename.
3. Require child exit zero, completed fresh artifacts, all ten emitted ranks and a valid
   fresh `.rrg` grid. Missing or stale output is a failed gate. A printed FAIL must return
   nonzero. Never compare a file overwritten by the second run with itself.
4. Compare all ten scientific PDB payloads and elected CF fields at their emitted precision.
   Preserve raw files/hashes. Across builds, normalize **only** the commit and dirty values
   in the exact `REMARK FLEXAID.commit=... FLEXAID.dirty=... FLEXAID.seed=...` provenance line.
   Keep seed, every other REMARK, coordinates, atom identities, scores, ranking and grid
   order unchanged. Report this as provenance-normalized byte equality, not raw hash identity
   or proof of equality below the instrument's serialized precision.
5. Exact initial-population gene/score receipts (§2) complement the emitted-output comparison.
   PASS means the declared observations agree. Invalid-pointer baseline execution is handled
   by §0.3; it cannot establish a valid scientific reference. Changes to valid scoring,
   ranking, clustering or thermodynamic models still require their own gated validation.

`ops/gate_parity.sh` and `ops/engine_repro_gate.py` implement the fail-closed invocation
and comparison interface. The obsolete hard-coded baseline MD5 is not a current engine pin.

---

## 2. Determinism check (multi-thread)

For changes touching parallel regions (GA eval, cleft detection), repeat candidate runs
at `OMP_NUM_THREADS=1` and `=4` twice each with `FLEXAIDDS_PARALLEL_REPRODUCE=1`.
Default-flag parity runs are separate and cannot silently substitute for flag-ON runs.

- Use an OpenMP-enabled build; record actual compiler flags and actual worker/team
  participation. Disable dynamic team sizing, reject thread-limit overrides and leave
  `FLEXAID_DETERMINISTIC` unset so the tested path is not silently serialized.
- **Cleft grid:** generate a fresh grid for each run; compare valid `.rrg` files in their
  original record order across runs/thread counts. A shared cached grid is not evidence
  of independently reproduced grid construction.
- **Initial population:** enable the opt-in `FLEXAIDDS_GEN0_RECEIPT` observer at a unique
  path. It snapshots stored values immediately after the initial population is evaluated
  and sorted, before reproduction, without rescoring, RNG calls or population mutation.
  Require successful engine completion as well as a complete receipt. Record actual seed,
  population count, exact gene/ring identity and exact stored score bits. Check complete
  order-independent gene/score record multisets; derive CF.com/CF.wal checksum summaries
  from those records. Equality of two sums alone cannot exclude compensating differences.
  Later generation traces or re-scored terminal populations are not this observation.
- Compare the two four-thread elected outputs as in §1. Keep the one-thread repeats and
  cross-thread comparisons explicit; never hide a failure by sorting grid or pose ranks.
- A matching diagnostic-only observer may be applied to baseline and candidate with its
  exact patch recorded. The observer is off by default and is tested for nonmutation and
  I/O failure. Failed writes, duplicate output paths, missing rows or wrong population
  counts fail the gate. The enabled observer does not make a run a publication campaign.

---

## 3. Astex-85 accuracy A/B (the science gate)

Purpose: prove a change to valid scoring, search, cleft selection, ranking or thermodynamic
models does not regress docking accuracy.

**Applicability to correctness repairs:** ownership, receipt, reporting and diagnostic-observer
repairs that preserve valid scientific operations use §§0.3, 1, 2, 4 and 6 for merge acceptance.
They do not require a new 85×10 publication campaign merely because the repaired code is in a
GA source file. Any unexplained valid-path scoring/output difference must be investigated;
if resolution changes the valid scientific model or search behavior, this full accuracy A/B
applies. This exception does not turn merge-validation runs into docking-success evidence.

- **Protocol:** autonomous (blind) mode — `--mode autonomous` — which exercises SURFNET cleft
  detection (top.cpp "Always run SURFNET"). Seed 12345, `FLEXAIDDS_NO_SEC=1`, 2000 gen / pop 1000.
- **Engine selection:** `FLEXAIDDS_BINARY=<engine>`; `benchmark_datasets --benchmark astex
  --only-codes <list>` subsets targets. (Direct-CLI per-target is acceptable when the harness's
  oracle-site guard blocks a blind run — document which was used.)
- **Metric:** top-1 rank-0 RMSD, 2.0 Å, in-place. Success ⇔ rank-0 in-place RMSD **`<= 2.0 Å`**.
  Name the instrument (§0: in-repo metric vs
  `score_reference.py`) — an unlabelled RMSD is not reportable.
- **Blind republish protocol (no % without a receipt):**
  `scripts/blind_astex85_receipt_protocol.py`. Fixed 85, `native_pose_seeded=0`,
  `seed_echo=0`, matrix MD5 `72d7c7396702331d96ff12d18f831796`. Default
  `SEED_ELITISM=0` / `NATIVE_SEED_FRAC=0`. `claim` refuses to print a success %
  without `RUN_RECEIPT.json`. Do not treat `--oracle-ceiling` as docking power.
  `scripts/reproduce_astex85.sh --dry-run` writes a receipt and does not dock.
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

`cd build && /opt/homebrew/bin/ctest --output-on-failure`. All configured tests must pass;
record the actual test count from this build instead of relying on a historical count. A
stale binary can cause a false PoseBustTests failure; rebuild before trusting a red.

---

## 5. Branch / merge discipline (OPS owns the merge)

- SWE agents implement on a **feature branch** and push; they do **not** merge to `main`.
- OPS (Benchmarker/Monitor/CI-CD) runs §1–§4 and the §6 audit, then merges.
- Every commit is reviewed (see §6). No commit lands on `main` unreviewed.

---

## 6. Commit review & audit (OPS/CI)

Every commit by any agent is reviewed against: (a) §1 parity, (b) §2 determinism if it touches
parallel code, (c) §3 accuracy when its scientific applicability above is met, (d) §4 ctest, (e) code-conduct in
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
