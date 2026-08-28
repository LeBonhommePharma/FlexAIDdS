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

---

# `DONE` — normative specification

**This section is normative.** Everything above it is rationale; this is the
contract an implementer conforms to. Frozen 2026-08-28 ahead of its first real
implementer.

The keywords MUST, MUST NOT and SHOULD carry their usual force.

## Location and name

The file MUST be named exactly `DONE`, with no extension, at the **run root** —
the same directory that holds `bin/` and `run/`, i.e. the directory `new_run_dir`
returns.

## Contents

Plain text, one `key=value` per line, LF-terminated, UTF-8. Same lexical rules as
`KIND` (see `KIND_SCHEMA.md`): no quoting, no comments, no blank lines, order not
significant, readers MUST ignore unknown keys.

```
status=ok
rc=0
targets_done=2
targets_total=2
finished_utc=2026-08-28T21:14:07Z
engine_sha=5ecbb89eebede8cba9271cbdd386496583bfbe2178cd9ded3db0f5512f11b511
run=gan2vsq5_20260828_162000/S1_1SQ5
```

| Key | Required | Value | Meaning |
|---|---|---|---|
| `status` | yes | `ok` \| `partial` \| `failed` | The runner's verdict on its own work |
| `rc` | yes | integer | Runner process exit code, verbatim |
| `targets_done` | yes | integer | Targets that produced **real output** — see below |
| `targets_total` | yes | integer | Targets the run was asked to do |
| `finished_utc` | yes | ISO-8601 UTC, `%Y-%m-%dT%H:%M:%SZ` | When the runner reached its end |
| `engine_sha` | yes | 64 lowercase hex | SHA-256 of the binary **actually invoked** |
| `run` | yes | string | Batch/arm identifier, e.g. `<batch>/<arm>` |

### `targets_done` counts poses, not rows

`targets_done` MUST count targets that produced at least one real pose with a
non-sentinel RMSD. It MUST NOT be a count of `result.csv` files or of rows in
them.

This is not pedantry. arm10 once recorded `rc=0 85/85` while 83 of the 85 targets
held zero poses and `-1.0000` sentinels — a `result.csv` existed for every
target and every one of them was empty. `run_t13_twotarget.sh:243-252` implements
the correct test per target; a conforming writer does the same.

### `status` values

| Value | Condition |
|---|---|
| `ok` | `rc == 0` **and** `targets_done == targets_total` |
| `partial` | The runner finished, but `targets_done < targets_total` |
| `failed` | `rc != 0`, or the run produced no usable output at all |

`status=ok` with `targets_done < targets_total` is a contract violation. A reader
that sees it MUST treat the file as untrustworthy rather than believing the
`status` field, because the two claims cannot both be true and the count is the
one derived from evidence.

## The three rules

These are the specification. The field list above is just how it is spelled.

**1 · `DONE` MUST be written last**, after every other output the run produces.

If it can be written before the work completes, it degenerates into
`RUN_RECEIPT.json` — which this codebase writes at `DatasetRunner.cpp:5395`,
ten lines before docking starts, making it a statement of intent that proves
nothing about completion (see `RUN_RECEIPT_CONTRACT.md` §5). That failure already
exists here. Do not reproduce it.

**2 · `DONE` MUST be written atomically** — write a temp file in the same
directory, then `rename()` it into place.

`rename(2)` within a filesystem is atomic; a plain `>` redirect is not. A crash
midway through a non-atomic write leaves a truncated `DONE` that still parses as
far as a naive reader gets, which reads as completion. A partial completion
marker is worse than none, because none is honest.

**3 · `DONE` MUST be written whenever the runner reaches its end, including on
failure.** It is not a success marker.

Therefore, and this is the only thing a reader needs to remember:

> **Absence means "died or still running." Presence means "the runner finished —
> see `status` for the outcome."**

Those are different claims, and conflating them is what made this document
necessary. A writer that emits `DONE` only on success collapses them again: its
absence would then mean "died, still running, *or* finished badly", which is not
a signal.

## Why nothing else can substitute

`DONE` is the **only durable killed-versus-finished signal**. Every cheaper
candidate fails the same way — it is a fact about the *artifacts*, and completion
is a fact about the *process*:

| Candidate | Why it fails |
|---|---|
| directory `mtime` | One scalar. Finished, killed, never-started and still-running all produce a last-write time. Not injective — the four states share the observable. |
| `result.csv` present | Written per target as the run proceeds. Present long before the run ends, and present in full for arm10's 83 empty targets. |
| receipt `binary_sha256` | Identifies *which engine*, written before docking starts. Says nothing about whether it ran to completion. |
| pose count | Cannot distinguish "all poses written" from "killed after the last one it happened to write". |
| process absent | True of both a finished run and a killed one. Also unobservable retroactively. |

Grok Bot currently uses `result.csv` + receipt sha as a stand-in for exactly this
reason: they are the best available proxies, and they are proxies. They answer
"what artifacts exist", never "did the runner reach its end". Only a token
written *by code that ran after the work* answers that, because only such a token
is evidence of reaching a line.

## Relationship to `KIND` and `SKIPPED`

Three files, three different claims. They are not redundant:

- **`DONE`** — primary evidence, written by the runner. What the process did.
- **`KIND.status`** — the index, in a fixed vocabulary a census reads without
  parsing free text. What the directory is worth.
- **`SKIPPED`** — never started. A missing observation, not a biased one.

Mapping, for a driver sealing `KIND` after writing `DONE`:

| `DONE` | → `KIND.status` |
|---|---|
| `status=ok` | `ok` |
| `status=partial` | `partial` |
| `status=failed`, `targets_done == 0` | `void` |
| `status=failed`, `targets_done > 0` | `partial` |
| *(no `DONE`, `SKIPPED` present)* | `unknown` |
| *(no `DONE`, no `SKIPPED`)* | `unknown` — died or still running |

The two vocabularies differ on purpose. `DONE` has no `void` or `unknown`:
writing it is itself the evidence, so `unknown` is unreachable, and `void` is a
judgement about scientific worth rather than about what the runner did. `KIND`
has no `failed` for the same reason in reverse — a census asks whether results
are usable, not how the process exited. **When they disagree, `DONE` wins and the
index is stale.**

## Compatibility

`ga1jd0_20260828_005342/DONE` and `ga1jd0_20260828_012729/DONE` predate this
specification and use an earlier informal single-line form
(`DONE <ISO8601>`). They are grandfathered, not retrofitted. A reader
encountering a `DONE` with no `=` on its first line MUST treat it as the legacy
form: presence still means the runner finished, and no other field is available.
