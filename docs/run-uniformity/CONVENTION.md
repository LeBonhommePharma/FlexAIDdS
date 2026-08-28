# Run-directory convention

What a run directory must contain, by tier, and why each entry is there.

Scope: this describes the target. Adoption is Phase 3 and forward-only — no
existing driver is retrofitted, and nothing in Phase 0 is sourced or invoked by
anything currently running.

---

## The problem this is fixing

60 result-bearing trees under `~/flexaidds_results`, in **53 distinct structural
shapes**. Only 8 follow `<campaign>_<TS>/arm<N>_<name>/run/<TARGET>/`.

Coverage does not split by how important a file is. It splits by **who writes
it**:

| Artifact | Written by | Present |
|---|---|---|
| `dock_config.json` | engine | 88% |
| `result.csv` | engine | 86% |
| `RUN_RECEIPT.json` | engine | 83% |
| `DONE` | shell driver | 40% |
| `claim.log` | shell driver | 36% |
| `provenance.txt` | shell driver | 36% |
| `STATUS.md` | shell driver | 21% |
| `inputs_at_launch` | shell driver | 18% |
| `inputs_at_end` | shell driver | 11% |

Engine-written ~85%. Driver-written 11–40%.

The gap is not carelessness. The drivers that exist are careful —
`run_t13_twotarget.sh` verifies the staged binary hash, freezes its
pre-registration next to the data, and runs a per-target void test. The gap is
structural: **11% is what "the driver author remembers" achieves**, and it does
not improve with effort, because each new driver re-decides the same five things
from scratch. `run_arm15_dofscale.sh`, `run_t13_twotarget.sh` and
`run_cachefix_2targets.sh` name their run directories three different ways
(`arm15_dofscale`, `run_t13_s374`, `13_intragenes`), use two different disk
floors, and only one of the three verifies that the binary it copied is the
binary it hashed.

The fix is not more discipline. It is to move the decision out of the driver:
`scripts/driver_layout.sh` makes the layout the default and the omission the
work.

Neither `inputs_at_launch` nor `inputs_at_end` is written by any of those three
drivers. At 18% and 11% they are the two least-present artifacts in the corpus,
which is what an artifact looks like when nobody's template includes it.

---

## Tier T — test

A small, throwaway-but-attributable run. Something a seat launches to answer one
question.

| File | Writer | Purpose |
|---|---|---|
| `KIND` | `new_run_dir` | What this directory is. See `KIND_SCHEMA.md`. |
| `RUN_RECEIPT.json` | engine | Configuration + engine/matrix hashes. See `RUN_RECEIPT_CONTRACT.md`. |
| `dock_config.json` | engine | Resolved per-target parameters as the engine actually saw them. |
| `provenance.txt` | driver | Human-readable intent: what varies, against what baseline, what is pinned. |
| `DONE` | driver | Terminal state. See below. |

## Tier A — arm

A run that will be **compared against another run**. Everything in tier T, plus
what a comparison needs.

| File | Writer | Purpose |
|---|---|---|
| *(everything in tier T)* | | |
| `result.csv` | engine | The scored output. The thing being compared. |
| `claim.log` | driver | Full stdout+stderr. The only place a `Fatal error:` survives. |
| `STATUS.md` | driver | Rolling progress, written atomically (`.tmp` + `mv`) while running. |
| `inputs_at_launch` | driver | Hashes of receptor / ligand / site / matrix as staged. |
| `inputs_at_end` | driver | The same hashes re-taken at completion. |

**Why an arm needs `inputs_at_launch` *and* `inputs_at_end`.** An arm's claim is
that one variable differed. That claim is only as good as the assertion that
nothing else did. Inputs on this machine are not immutable: caches get rebuilt,
site directories get regenerated, a `TMPDIR` sweep can remove a directory
mid-run — that last one cost arm14 eight of 85 targets and is the reason
`driver_preamble.sh` exists. Hashing once proves what was staged. Hashing twice
proves it stayed staged. A single hash at launch cannot detect a mid-run
substitution, and a single hash at the end cannot tell you the run started from
what you thought.

They are the two rarest artifacts in the corpus (18%, 11%) and the two that a
between-arm comparison most depends on. That is the coverage-by-author problem
in one line.

**Why `claim.log` is not optional at tier A.** Return codes compress a run to one
integer. `run_arm15_dofscale.sh:72` greps `claim.log` for `Fatal error`
specifically because a run can emit fatal per-target errors and still exit 0 —
the engine writes a sentinel row and continues. Without the log the arm looks
clean.

---

## Why `DONE` earns its place

`DONE` is the single most valuable driver-written artifact and it is present 40%
of the time.

**The argument, from first principles.** Consider using directory mtime instead —
the obvious "free" alternative, since the filesystem maintains it for you. mtime
is one scalar: the time of the last write. The terminal state of a run is a
categorical variable with at least four values:

1. finished normally
2. killed — SIGKILL, OOM, laptop slept, someone hit the power
3. never started — disk floor, precondition refused
4. still running, just slow

All four produce a directory with a last-write time. mtime is therefore not an
injective function of terminal state: **the four states map onto the same
observable.** No amount of care reading mtime recovers the distinction, because
those bits were never written to begin with. You cannot recover information that
was never recorded — you can only arrange for it to be recorded.

That is what `DONE` is: a token that can only be created by code that ran *after*
the thing it attests to. Its existence is not a description of the run, it is
**evidence of reaching a specific line**. A killed run leaves everything a
finished run leaves, minus `DONE`, and that asymmetry is the entire signal.

**This is live right now, not hypothetical.** In
`gan2vsq5_20260828_162000/`, `S1_1N2V` is void — the runner died after r1 — and
`S1_1N2V_b` is the relaunch, currently mid-flight. Both have recent mtimes.
Neither has `DONE`. By mtime alone, a dead arm and a running arm are the same
observation. That is exactly the discrimination `DONE` exists to make, and it is
why `scripts/backfill_kind.sh` refuses to walk a batch that is being written.

**Corollaries, which are conditions on `DONE` being worth anything:**

- **Written last, and only on the success path.** A `DONE` written before the
  work is a `DONE` that means nothing.
- **`SKIPPED` is a separate token, not `DONE` with a flag.** All three drivers
  write `"never started, SKIPPED not killed"`. "Never started" and "died" have
  different scientific consequences: one is a missing observation, the other is a
  possibly-biased one. Collapsing them loses the distinction that matters most.
- **Existence is necessary, not sufficient — content carries the void test.**
  `run_t13_twotarget.sh:243-244` records that arm10 once wrote
  `rc=0 85/85` while 83 of 85 targets held zero poses and `-1.0000` sentinels. A
  `result.csv` count is not a void test. `DONE` must record per-target poses and
  a real RMSD, which is why the t13 line reads
  `rc=$rc targets_with_poses=${nreal}/$NT poses=... sentinel=...` rather than
  just `rc=0`.
- **`DONE` and `KIND.status` are not redundant.** `DONE` is raw evidence written
  by the driver. `status=` is the interpretation, in a fixed vocabulary a census
  can read without parsing free text. `DONE` is the primary record; `status` is
  the index. When they disagree, `DONE` wins and the index is wrong.

**The cost asymmetry, which is the whole reason this is an easy call.** Writing
`DONE` costs one `echo` at the end of a run that already took hours of CPU.
Not writing it costs the ability to ever distinguish a completed arm from a
killed one — retroactively, permanently, for every run in that shape. The
`--redock` receipt gap on 2026-08-28 is the same trade already resolved the wrong
way once: the missing write was free, and its absence was unrecoverable.
