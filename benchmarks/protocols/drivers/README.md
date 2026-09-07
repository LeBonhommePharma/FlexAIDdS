# Campaign drivers — and the failure each gate was built from

These are the shell drivers that produced the 2026-08/09 Astex-84 campaign arms and the
determinism investigation. They are committed **verbatim, exactly as they ran** — not
re-pathed, not tidied. That is deliberate: they are the record of what produced the
numbers, so a committed script that never ran would be worse than a slightly awkward one
that did.

**Portability caveat, stated up front.** Each carries 1–2 hardcoded references to a
results root (`/Users/lp.more/flexaidds_results`). They are session drivers, not portable
tools. To reuse one, change the root at the top; the gates below are the transferable part.

## Why this directory exists

Every gate in these scripts is here because something failed once and cost hours. The
gates are the actual deliverable — more durable than any single arm's numbers — so each
is documented with the failure that motivated it.

### The resume guard (`prescrub`) — `run_water_ablation_w4.sh`

A driver was killed mid-cell, leaving `run/W0/1igj` with `RUN_RECEIPT.json` and
`provenance.json` written but **no `DONE`**. The resume logic skipped only cells carrying
`DONE`, so it re-ran into that dirty directory — and `DatasetRunner` honoured the stale
receipt instead of `--only-codes`, docking the **entire dataset inside one cell**: 72
target directories, 9,034 poses, 12.4 hours of wall, while the driver's other two slots
sat idle.

Two gates came out of it:

1. **Pre-scrub.** Any cell directory with partial output and no `DONE` is renamed aside
   with a timestamp (`<cell>.dirty_<TS>`) before anything runs. Renamed, never deleted —
   see the deletion convention below.
2. **Post-check.** Any cell whose output holds more than one target directory fails with
   `rc=90`. This is the check that would have caught the above in minute one rather than
   hour twelve.

A related trap, same root cause: `posehash()` globs every `*.pdb` in the cell directory,
so orphan poses from a killed run get hashed together with fresh ones. The pre-scrub fixes
both.

### The overlap gate — `run_contention_probe_flex.sh`, `run_contention_probe_confirst.sh`

A contention experiment whose cells did not actually run concurrently measures nothing.
These record per-cell start/end times and compute real pairwise overlap, refusing to
report a verdict if the arms did not overlap. Written after realising the isolated and
concurrent arms could serialise and still "pass".

### The flexibility gate — `run_contention_probe_flex.sh`

Earlier probes were reported as engine reproducibility when they had all run **rigid**,
while the campaign arm ran with flexible sidechains. The gate refuses to continue unless
the emitted `dock_config.json` confirms flexibility actually engaged — because a
"flexible" probe that silently runs rigid repeats the exact error it exists to correct.

### The stability gate — `run_overnight_chain_v2.sh`

The first chain required `load1 < 4.0` before launching. Measured minimum on this machine
across 18 polls was **7.26**, with nothing writing to disk — so the gate was below the
box's floor and would have burned its whole deadline launching nothing. Replaced with a
gate on **stability** (coefficient of variation over a rolling window) rather than
silence, which is what actually matters for a timing measurement: a steady baseline still
leaves a real manipulation, and every receipt records `load1` so the contrast is
verifiable rather than assumed.

Lesson generalised: *a gate that cannot pass is as broken as one that cannot fail.*

### Prediction-before-launch — `run_flexdet_probe.sh`

The expected outcome table is written into the driver **before** it runs, so the result
cannot be reinterpreted after the fact. The arm that landed matched the row the table
named.

### Non-vacuity — all probes

Several gates in this project's history reported PASS on empty output — in one case the
compared hash was `e3b0c442…`, the SHA-256 of the empty string, in all three arms. Every
gate here checks that it had something to compare (pose counts, comparison counts) and
announces that count **before** any verdict, so "no differences" can never be confused
with "nothing was examined".

### Wake-on-exit — `monitor_v3.sh`

A background monitor only notifies the session when it **exits**. An earlier version had a
wake branch that logged and continued looping, so it could never wake anything. Every wake
condition here exits, and the exit-branch count is verified before arming. It also lives
in a file rather than an inline heredoc, so the loop is not tied to the kernel that
launched it.

## Conventions these drivers follow

- **Never delete a directory carrying state or provenance.** Rename it aside with a
  timestamp. Verified mechanically: these drivers contain zero `rm` verbs.
- **Never edit a script a running process is still reading.** Write a new file with a new
  name — which is why `w3` and `w4` both exist here, and why the chain has a `v2`.
- **A stop lever, not a kill.** Each long-running driver polls a sentinel file
  (`CHAIN2_STOP` and similar) so it can be halted cleanly without signalling a detached
  process, and without writing skip-markers that would fake completion and cascade-release
  successors.
- **Failure receipts are honest.** A killed or failed cell records `rc != 0`; it is never
  marked complete to unblock a queue.
- **bash 3.2 only.** macOS ships 3.2, so `wait -n`, `declare -A`, `mapfile` and `${x,,}`
  are unavailable. All 11 drivers parse clean under `/bin/bash -n`.
- **A syntax check is not a behavioural test.** One patch to `w4` inserted a function's
  own call inside its loop, creating infinite recursion; `bash -n` passed it and only a
  timeout-bounded behavioural test caught it by hanging.

## Inventory

| driver | what it produced |
|---|---|
| `run_optresfix_requeue_v2.sh` | the 345-cell OptRes-fix requeue campaign |
| `run_water_ablation_w3.sh` | the water-retention arms (W0/W2), 3-wide |
| `run_water_ablation_w4.sh` | w3 + the resume guard and one-target gate |
| `run_igj_redo.sh` | single-cell recompute of the contaminated `1IGJ` cell, off the critical path |
| `run_determinism_probe.sh` | thread-count × keyed-jitter reproducibility matrix |
| `run_confound_probe.sh` | campaign-binary × GA-budget confound isolation |
| `run_contention_probe_confirst.sh` | process-concurrency effect, concurrent arm first |
| `run_contention_probe_flex.sh` | the same on the flexible arm, with the flexibility gate |
| `run_flexdet_probe.sh` | the decisive flexible-path determinism arms |
| `run_overnight_chain_v2.sh` | stability-gated chain: contention probe → water ablation |
| `monitor_v3.sh` | campaign monitor, exit-on-wake |
