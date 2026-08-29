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
| `status` | yes | `ok` \| `partial` \| `failed` \| `unverified` | The verdict on the work. `unverified` is observer-only — see *Attestation* |
| `rc` | yes | integer, or `unknown` | Runner process exit code, verbatim. `unknown` **only** in a reconstructed `DONE` — see *Attestation* |
| `targets_done` | yes | integer | Targets that produced **real output** — see below |
| `targets_total` | yes | integer | Targets the run was asked to do |
| `finished_utc` | yes | ISO-8601 UTC, `%Y-%m-%dT%H:%M:%SZ` | When the runner reached its end |
| `engine_sha` | yes | 64 lowercase hex | SHA-256 of the binary **actually invoked** |
| `run` | yes | string | Batch/arm identifier, e.g. `<batch>/<arm>` |
| `source` | no | `runner` \| `reconstructed` | Who wrote this file. Absence means `runner` — see *Attestation* |

### `targets_done` counts poses, not rows

`targets_done` MUST count targets that produced at least one real pose with a
non-sentinel RMSD. It MUST NOT be a count of `result.csv` files or of rows in
them.

This is not pedantry. arm10 once recorded `rc=0 85/85` while 83 of the 85 targets
held zero poses and `-1.0000` sentinels — a `result.csv` existed for every
target and every one of them was empty. `run_t13_twotarget.sh:243-252` implements
the correct test per target; a conforming writer does the same.

### `status` values

| Value | Condition | Writer |
|---|---|---|
| `ok` | `rc == 0` **and** `targets_done == targets_total` | runner only |
| `partial` | The runner finished, but `targets_done < targets_total` | either |
| `failed` | `rc != 0`, or the run produced no usable output at all | either |
| `unverified` | `targets_done == targets_total`, and `rc` is unrecoverable | **observer only** |

`status=ok` with `targets_done < targets_total` is a contract violation. A reader
that sees it MUST treat the file as untrustworthy rather than believing the
`status` field, because the two claims cannot both be true and the count is the
one derived from evidence.

`ok` is unreachable for a reconstructed `DONE`: `ok` requires `rc == 0`, and a
reconstructed `DONE` MUST carry `rc=unknown`, which is not `0`. `unverified` is
the value that fills the cell `ok` cannot occupy — all the work is present, and
nobody observed how the process exited. See *Attestation*.

`unverified` is **not** a weaker `ok` for the same evidence. It is the correct
value for *different* evidence: `ok` rests on a reported exit code, `unverified`
rests on a recount of artifacts. The two are not on a quality ladder, they are
claims of different kinds, which is why `rc` and `source` must be read alongside.

### Why the value is named `unverified` and not `complete`

`complete` was the obvious candidate and it is the wrong word here, on three
independent counts, all of them checkable in this repository:

1. **This project already sorts `complete` into the claim-asserting bucket.**
   `scripts/validate_thermo_claims.py:29-45` maintains two sets:
   `CLAIM_STATUS_WORDS = {complete, completed, published, reproducible,
   validated, verified}` and `NONCLAIM_STATUS_WORDS = {candidate, draft, example,
   pending, planned, unverified}`. A status word that means "nobody checked the
   exit" belongs in the second set by construction, and `unverified` is already
   in it. Naming the value `complete` would put a claim word on a non-claim, in
   a codebase that has a validator whose entire job is to catch exactly that.
2. **`complete` is live in a neighbouring `status` namespace.**
   `benchmarks/m3pro/dashboard/fleet_monitor.py:557,587,930` compares
   `status == "complete"`, and `benchmarks/re-dock/orchestrator.py:313` sets
   `status = "completed"`. Both mean *finished and fine*. Reusing the token for
   *finished, disposition unknown* is a collision inside one repository.
3. **`COMPLETE` is already a leading token in a `DONE` on disk.**
   `vctent_20260828_034330/DONE` opens `COMPLETE <ts> — NOT truncated.` — used
   there to make the *strong* claim. The word is spoken for, and it is spoken for
   in the opposite direction.

There is also a plain-English hazard independent of this repository: `complete`
reads as a stronger claim than `ok`, when the value is weaker than `ok`. An
enum whose ordering inverts on a casual read is a bad enum. `unverified` cannot
be misread that way — it announces its own gap.

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

## Attestation — testimony versus reconstruction

*Amendment ratified 2026-08-28, after the freeze above. Normative.*

Everything in this document so far assumes the runner wrote its own `DONE`. That
assumption is what makes the token worth anything: a `DONE` is **testimony** — a
claim by the process that did the work, and evidence that a specific line of code
was reached. Rule 1 is not really about ordering, it is about that. A file
written by an observer after the fact, reading the artifacts on disk and
inferring what must have happened, is **inference**. It looks identical on disk.

Both are legitimate. A reconstruction is often the only record obtainable for a
run whose driver died, and refusing to write one loses information. What is not
legitimate is silently formatting inference as testimony, because the reader has
no way to tell them apart and the whole discrimination this file exists to make
collapses.

The `rc` rule follows directly, and is the reason `rc` is singled out below.
Every other field has an observable counterpart: `targets_done` can be recounted
from poses, `engine_sha` re-hashed from the binary, `run` read off the path.
`rc` has none. An exit code is a property of a process that no longer exists and
leaves no trace in the artifacts. It is therefore the one field where
reconstruction is not weak inference but pure invention, and a reconstructed
integer would read to every downstream consumer as a reported one.

**The rule:**

A `DONE` written after the runner has already exited — reconstructed by an
observer from the artifacts on disk — MUST carry `source=reconstructed`. A `DONE`
written by the runner itself MAY carry `source=runner` or omit the field; absence
of the field means runner-written.

A reconstructed `DONE`:

- MUST set `rc=unknown`. The exit code is not recoverable from artifacts and a
  reconstructed value would read as reported.
- MUST NOT assert `status=ok` on inference alone. It may assert only what is
  independently checkable from the artifacts.
- MUST use `status=unverified` when, and only when, its recount finds
  `targets_done == targets_total` and `rc` is unrecoverable. This is the value
  `ok` cannot occupy; rounding up to `ok` or down to `partial` are both
  misreports, in opposite directions.
- MAY set `targets_done` from a recount of real poses with non-sentinel RMSD,
  since that is directly observable.

A reader MUST treat `source=reconstructed` as weaker evidence than a
runner-written `DONE`.

### `unverified` — normative

*Ratified by LP, 2026-08-29. Normative.*

**A runner MUST NOT write `status=unverified`.** A process knows its own exit
code; `rc` is available to it by construction, at the only moment it is ever
available to anyone. A runner emitting `unverified` is either lying about what it
can see or is not the process that did the work. There is no third case, so there
is no legitimate one.

This is the field's defining property and the reason it is stated rather than
left to inference: **`unverified` is entangled with `source` in a way the other
three values are not.** `ok`, `partial` and `failed` are meaningful from either
author. `unverified` is meaningful from exactly one.

The consequent rules:

- A `DONE` carrying `status=unverified` MUST also carry `source=reconstructed`
  and `rc=unknown`. All three or none of them; the combination is one claim
  written in three fields, not three independent facts.
- `status=unverified` with `rc` set to an integer is a contract violation. The
  writer either had the exit code, in which case the status is wrong, or invented
  it, in which case the `rc` is.
- `status=unverified` with `targets_done < targets_total` is a contract
  violation. The short count is observable and `partial` reports it; `unverified`
  claims completeness it does not have.
- `status=unverified` with `source=runner`, or with `source` absent (which means
  runner-written), is a contract violation.

A reader encountering any of those four combinations MUST treat the file as
untrustworthy in whole, not repair it in part. The fields disagree, and which one
is wrong is not recoverable from the file.

**Worked example — all three fields, not one.** A conforming reconstruction of a
run whose recount finds every target present:

```
status=unverified
rc=unknown
targets_done=1
targets_total=1
finished_utc=2026-08-28T21:03:39Z
engine_sha=5ecbb89eebede8cba9271cbdd386496583bfbe2178cd9ded3db0f5512f11b511
run=<batch>/<arm>
source=reconstructed
```

`status`, `rc` and `source` move together. A retrofit that writes
`status=unverified` while leaving `rc` or `source` as they were has produced one
of the four violations above, not a partial improvement.

**What a writer records instead, when it cannot use `unverified`.** A
reconstruction that finds a short count uses `partial`. One that finds no usable
output uses `failed` — whose second disjunct ("produced no usable output at all")
is observable from artifacts and therefore reachable by an observer without
inventing an `rc`. `unverified` is for the one cell those two do not cover.

**`aborted` is not a value.** `gan2vsq5_20260828_162000/S1_1N2V/DONE` carries
`status=aborted`, which has never been in this vocabulary. The evidence recorded
in that same file — `targets_done=0`, no `result.csv`, no elected pose — is
exactly the `failed` condition, and `failed` was available. Inventing a value at
write time is how a vocabulary stops being one. That file is grandfathered as a
record of what happened; it is not a precedent.

### Consequences a writer should expect

`source=reconstructed` does **not** relax rule 2 (atomic write) or the field
list. A reconstruction is still a `DONE` and is still written with `rename()`.

Rule 3 does not apply to a reconstruction, because there is no "runner reaching
its end" for the observer to hook. A reconstruction is written when someone
chooses to write one.

The justification given under *Relationship to `KIND` and `SKIPPED`* for `DONE`
having no `unknown` status — "writing it is itself the evidence, so `unknown` is
unreachable" — holds only for runner-written files. For a reconstruction, writing
it is *not* evidence of reaching a line, so that argument lapses.

The amendment as first ratified left that consequence unresolved: a
reconstruction observing `targets_done == targets_total` with an unrecoverable
`rc` had no status value that fits. `unverified` is that value, added 2026-08-29.
The diagnosis behind it is that `status` was carrying two independent axes —
**completeness**, observable from artifacts, and **exit disposition**, knowable
only at termination. For a runner-written file the two always co-occur, so one
enum could carry both without anyone noticing. Admitting observer-written files
decouples them and exposes the missing cell.

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
| `status=unverified` | `partial` — see below |
| `status=failed`, `targets_done == 0` | `void` |
| `status=failed`, `targets_done > 0` | `partial` |
| *(no `DONE`, `SKIPPED` present)* | `unknown` |
| *(no `DONE`, no `SKIPPED`)* | `unknown` — died or still running |

**`unverified` → `partial` is deliberately conservative, and it under-reports.**
A full-count reconstruction is not "stopped short", so `partial` is not literally
right. It is chosen because KIND's `ok` is defined as *"Ran to completion. `DONE`
records `rc=0`"* (`KIND_SCHEMA.md`), and a reconstruction records `rc=unknown`.
Promoting `unverified` to KIND `ok` would require redefining KIND `ok` away from
`rc`, which is a separate ruling and a separate document. Until then the index
under-claims, which is the direction an index should err. Note that `partial` is
also what `scripts/backfill_kind.sh` already infers for such a file without any
change, so the conservative mapping costs nothing to adopt.

The two vocabularies differ on purpose. `DONE` has no `void` or `unknown`:
writing it is itself the evidence, so `unknown` is unreachable, and `void` is a
judgement about scientific worth rather than about what the runner did. `KIND`
has no `failed` for the same reason in reverse — a census asks whether results
are usable, not how the process exited. **When they disagree, `DONE` wins and the
index is stale.**

That precedence assumes testimony. A driver sealing `KIND` from a `DONE` carrying
`source=reconstructed` is reading inference, not evidence, and MUST NOT treat it
as outranking its own observations: where a reconstructed `DONE` disagrees with
what the directory shows, neither wins automatically and the disagreement is the
finding.

## Compatibility

`ga1jd0_20260828_005342/DONE` and `ga1jd0_20260828_012729/DONE` predate this
specification and use an earlier informal single-line form
(`DONE <ISO8601>`). They are grandfathered, not retrofitted. A reader
encountering a `DONE` with no `=` on its first line MUST treat it as the legacy
form: presence still means the runner finished, and no other field is available.

`source` has the same problem one layer up, and it is worth stating rather than
leaving to be discovered. Every `DONE` written before the *Attestation* amendment
lacks the field, including any that were in fact reconstructed, so the
"absence means runner-written" default reads them all as testimony. The default
is right going forward and wrong retroactively. A reader MUST NOT infer
runner-written from absence alone for a `DONE` whose `finished_utc` precedes this
amendment; for those files, provenance is established out of band or not at all.
Existing files are grandfathered, not retrofitted — retrofitting `source` onto
one is a deliberate act by whoever knows how it was written, not a migration.

### What `unverified` touches — measured 2026-08-29

Census over `~/flexaidds_results`, read-only, no live runner:

```
DONE files                                             61
  in this document's key=value form                     5
  pre-spec free text (8 distinct informal dialects)     56
status= present                                         5
  status=ok                                             4
  status=aborted                                        1   <- never a valid value
source=reconstructed                                    3
rc=unknown                                              3
KIND files                                              0   <- sidecar not yet adopted
```

**Two files are the migration set**, and both are in violation of the
*Attestation* amendment as it stands today:

| File | Now | Under this amendment |
|---|---|---|
| `gan2vsq5_20260828_162000/S0_1N2V/DONE` | `status=ok rc=unknown source=reconstructed` | `status=unverified` |
| `gan2vsq5_20260828_162000/S1_1N2V_b/DONE` | `status=ok rc=unknown source=reconstructed` | `status=unverified` |

Both assert `status=ok` on inference alone, which the amendment forbids, and both
carry `rc=unknown`, which makes `ok` unreachable by the value table regardless.
They are not edge cases discovered by this work — they are the files that raised
the question. `unverified` is the value they should have had.

A third file, `gan2vsq5_20260828_162000/S1_1N2V/DONE`, carries `status=aborted`
and should be `failed`; see *Attestation*. Rewriting any of the three is a
deliberate act by their author, not a migration this document performs.

The remaining 56 `DONE` files carry no `status=` field at all and are unaffected:
they are the grandfathered informal forms, read by presence and by free text.

### Readers, and the one that needs fixing first

Nothing in the repository or the corpus reads `status=` today. That is why
`status=aborted` has sat on disk unremarked. The inventory:

| Reader | Reads | Behaviour on `unverified` |
|---|---|---|
| `scripts/check_run_receipt.py` | receipts only — never opens `DONE` | unaffected |
| `scripts/benchmark_ops_monitor.py` | pose counts, process liveness | unaffected — never opens `DONE` |
| `scripts/backfill_kind.sh:150-151` | `DONE`, via `grep -q 'rc=0'` | infers `partial` — correct, but see below |
| corpus chain gates (`run_arm1[2-4]*.sh`, `run_combo.sh`, `run_priority.sh`, `run_boom.sh`, `run_temp2016.sh`) | `[ -f DONE ]` | unaffected — presence only, `status` never read |
| corpus status tables (`monitor_priority.sh`, `run_t13_twotarget.sh`, `run_followon.sh`, …) | `cat DONE` for display | unaffected |
| `workorders/run_pshare_v2.sh:82` | `grep -q "^VALID"` | unaffected — anchored, fails closed |

**`backfill_kind.sh:151` is the one to fix.** `grep -q 'rc=0'` is unanchored and
matched against the whole file, not against a `^rc=` line. It gives the right
answer for every file on disk today, and the wrong one for a plausible next one:
a reconstruction whose free-text `evidence=` field quotes a per-target log line
containing `rc=0` would be indexed `status=ok` — an inferred run promoted to
"records `rc=0`" by a substring. The same defect already fires on multi-target
`DONE` files, where `g43_mutgran_20260827_090753/run_ctrl_s12345/DONE` is
inferred `ok` from a *per-target* `rc=0` on line 2, with no run-level `rc` in the
file at all. `grep -qE '^rc=0$'` fixes both. Not applied here — this document
does not change code.
