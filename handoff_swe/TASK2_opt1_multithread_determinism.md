# TASK 2 — Make Opt1 parallel-reproduce bit-deterministic at >1 thread
Role: software engineer (deep concurrency + numerics). Model: **Opus 4.8** (or Fable 5) —
scientifically subtle, high-value; NOT a mechanical task.
Repo: /Users/lp.more/Projects/FlexAIDdS   Base: branch opt1/stale-status-fix (commit 9b514ca82)

## Established by OPS (do not re-derive; instrumentation already reverted)
Env: FLEXAID_SEED=12345 FLEXAIDDS_NO_SEC=1, 1G9V, gen-0 population CF checksum:
  - flag OFF, 4 threads:  -302083.242271  reproducible
  - flag ON,  1 thread:   -302096.627183  reproducible (stale-status fix works)
  - flag ON,  4 threads:  -302130.80 / -302302.87  DIVERGES  <-- the bug
Per-chromosome decomposition (flag ON, 4-thread, run1 vs run2):
  998/1000 chromosomes identical CF; **2/1000 genuinely differ** (one with large `wal`
  swing ~1016 -> a rare clash/penalty branch). Offspring genes are created SERIALLY in
  reproduce() (gaboom.cpp ~1760-1970), so population contents are deterministic; only
  calculate_fitness is parallel. ThreadSanitizer: GA eval path (gaboom/vcfunction/
  Vcontacts/ic2cf) is RACE-CLEAN. cleftgrid is read-only. tl_fa/tl_res/tl_atoms are
  deep-ish copies; resid.bond and other FA pointer members alias shared memory but are
  READ-only in eval. DEE list is !omp_in_parallel()-guarded (inactive under Opt1).

## The remaining bug (find it)
A rare, branch/boundary-sensitive numeric nondeterminism reached only under multi-thread
batch eval, in the per-thread tl_fa/tl_vc/tl_atoms COPY/INIT path — TSan-invisible
(ordering / rare uninitialised-read, not a coarse race). Candidate hunt order:
  1. Diff the 2 divergent chromosomes' full CF decomposition (com/wal/sas) 1-thread vs
     4-thread; identify the branch (the wal=1016 one is a clash path).
  2. Audit every FA_Global pointer member NOT re-pointed after `tl_fa(n_thr,*FA)`
     (gaboom.cpp ~2499): confirm each is read-only OR give it a per-thread buffer.
  3. Check tl_vc / tl_calc scratch for a field consumed before it is written on a rare branch.

## Deliverables
- A fix that makes flag-ON 4-thread gen-0 checksum reproducible AND full-dock 10/10 poses
  byte-identical run-to-run at 4 threads (1G9V, seed 12345).
- Keep flag OFF by default. Push branch `fix/opt1-mt-determinism`. Do NOT merge to main.

## OPS acceptance gate (I verify)
- gen-0 checksum reproducible at 1 and 4 threads, flag ON.
- 10/10 elected poses byte-identical across two 4-thread runs.
- Default-flag parity vs main preserved (CF + 10/10 poses byte-identical, flag OFF).
- ctest 11/11 green.
