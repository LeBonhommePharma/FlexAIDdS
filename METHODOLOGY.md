# METHODOLOGY.md — FlexAIDdS canonical benchmarking, determinism & CI methodology

**Single source of truth for HOW work is measured, validated, and gated — for every agent
(Claude, Claude Code, Codex, Grok, GPT) and every human.** `AGENTS.md` governs coding conduct;
this file governs *methodology*. When a skill, script, or agent describes a benchmark, a
determinism check, or a merge gate, it must **defer to this file**, not restate it. If a
procedure here changes, it changes here first and everything else re-reads it. The goal is
**agent-independence**: the same task run by any agent produces the same result because all of
them execute the identical procedure below.

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
- **RMSD reporting:** report **rank-0 elected pose RMSD** (the elected `_0.pdb`), symmetry-corrected,
  heavy-atom, 2.0 Å cutoff. NEVER report seed-elitism / `_INI.pdb` RMSD as the result.
- **RMSD engine:** spyrmsd 0.9.0 in the `python` conda env
  (`/Users/lp.more/.claude-science/conda/envs/python/bin/python`). Fallback: element-blocked
  Hungarian only if spyrmsd raises (must be logged).

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
- **Metric:** top-1 rank-0 RMSD, spyrmsd, 2.0 Å. Report per-target and aggregate success.
- **Acceptance:** candidate within noise of baseline; **no target flips success→fail** attributable
  to the change. A full landing decision uses the full 85 × 10-restart protocol; a fast pre-check
  may use a documented subset.
- **Published anchors (for context, not gates):** Gaudreault & Najmanovich 2015 JCIM **Table 2**
  Astex native FLRP: **top-1 = 45.2%**, **top-10 = 66.7%** (do not swap labels). Morency 2017
  3Dsig/poster FlexAID ~66% / entropy ~69% are **S_top10-family** medians — not the JCIM top-1
  figure. Report the protocol next to any number; never mix top-1, top-10, and S_top10 without labels.

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
