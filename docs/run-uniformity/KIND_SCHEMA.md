# `KIND` — the run-directory sidecar

One file, per run directory, that says what the directory *is*. Small enough
that writing it is never the expensive step, and greppable enough that a census
never needs a JSON parser.

```
kind=arm
tier=A
status=unknown
name=15_dofscale
created=2026-08-28T20:41:07Z
by=claude-code
engine_sha=5ecbb89eebede8cba9271cbdd386496583bfbe2178cd9ded3db0f5512f11b511
```

---

## Format

- Plain text, one `key=value` per line, LF-terminated, UTF-8.
- No quoting, no escaping, no comments, no blank lines, no sections.
- Keys are lowercase ASCII. Values contain no newline and no leading or trailing
  whitespace.
- **Order is not significant.** Every field is addressable as
  `grep -m1 '^status=' KIND | cut -d= -f2-`. The order below is for humans.
- Unknown keys must be ignored by readers, so the schema can grow.
- `cut -d= -f2-`, not `-f2` — values may contain `=` in future fields.

Why not JSON: the file has to be readable by a `find`/`grep` one-liner across
tens of thousands of directories, from any seat, with no dependencies. The
project already has a JSON provenance format — `RUN_RECEIPT.json` — and its
existence is not what is missing. Coverage is what is missing.

---

## Fields

| Key | Required | Values | Meaning |
|---|---|---|---|
| `kind` | yes | `test` \| `arm` | What the directory is for |
| `tier` | yes | `T` \| `A` | Which artifact set is expected (see `CONVENTION.md`) |
| `status` | yes | `ok` \| `void` \| `partial` \| `unknown` | What the contents are worth |
| `name` | yes | slug `[a-z0-9][a-z0-9_-]*` | Short identifier, no prefix |
| `created` | yes | ISO-8601 UTC, `%Y-%m-%dT%H:%M:%SZ` | When the directory was created |
| `by` | yes | `grok` \| `science` \| `claude-code` \| `dispatch` \| `lp` \| `unknown` | Which seat created it |
| `engine_sha` | yes | 64 lowercase hex, or `unknown` | SHA-256 of the binary that will run in it |

### `engine_sha` — the field that justifies the file

On 2026-08-28 a binary was relinked mid-sequence. `--redock` writes no receipt,
the inodes were gone, and the project builds with Unix Makefiles, so there was no
build ledger either. The engine identity for those runs is not recoverable —
not "hard to recover", **gone**.

`scripts/driver_layout.sh` therefore treats this as a precondition, not a field:
`new_run_dir` refuses to create a directory at all unless
`FLEXAIDDS_ENGINE_BIN` points at a hashable file. A run directory that cannot be
attributed to an engine is not cheaper than no directory — it is a directory
whose results cannot be used later, produced at full CPU cost.

The literal `unknown` is permitted **only** for backfilled records, where it is
the honest recording of exactly this gap.

### `status` — enumerable is not the same claim as valid

`find -name KIND` is meant to enumerate runs. If stamping a sidecar also implied
the run was good, the first void directory to be backfilled would poison every
census built on it — and there is a void directory on disk right now
(`gan2vsq5_20260828_162000/S1_1N2V`, runner died after r1).

So the two claims are kept separate. `KIND` makes a tree *enumerable*; `status`
says what it is *worth*.

| Value | Means |
|---|---|
| `ok` | Ran to completion. `DONE` records `rc=0`. |
| `partial` | Ran and stopped short: non-zero `rc`, or results with no `DONE`. |
| `void` | Ran and produced nothing usable — zero poses, `-1.0000` sentinels, or an explicit void marking. Distinct from `partial`: a void arm is unscoreable, not merely incomplete. |
| `unknown` | No positive evidence either way. **The default.** Also covers `SKIPPED` — never started is not void. |

`new_run_dir` writes `status=unknown` at creation and never anything else: a
directory that has produced nothing yet has not earned `ok`. Sealing it is the
driver's job at end of run, alongside `DONE`. Backfill defaults to `unknown` and
promotes only on positive evidence.

`DONE` uses a different vocabulary (`ok|partial|failed`) on purpose — it reports
what the *runner* did, where `status` reports what the *directory is worth*. The
normative mapping between them is in `CONVENTION.md`, under
"`DONE` — normative specification". Where the two disagree, `DONE` is primary
evidence and `status` is a stale index.

### `by` — carried, pending LP's ruling

LP has not ruled on whether seat attribution belongs in the schema. It is
included on the assumption that it is wanted; it is one line, one field, and one
`case` arm in `driver_layout.sh` to remove if not.

The argument for keeping it: the coverage split that motivates this whole
exercise is a split by *author*, not by importance. Engine-written artifacts land
~85% of the time and driver-written ones 11–40%. Recording which seat produced a
tree makes that measurable per seat instead of in aggregate.

### `kind` and `tier` are currently 1:1

`kind=test` implies `tier=T`; `kind=arm` implies `tier=A`. As written, one of the
two is redundant.

Both are kept because they answer different questions — `kind` is *what this is*,
`tier` is *what must be in it* — and the second is the one a checker reads. If a
T-tier arm or an A-tier test ever becomes real, they stop being redundant.
`driver_layout.sh` refuses a mismatch unless
`FLEXAIDDS_KIND_ALLOW_TIER_MISMATCH=1` is set. **If that override is never used,
delete one of the fields** — that is the deciding experiment, and the override
exists to run it.

---

## Worked examples

New arm, created by `new_run_dir arm A 15_dofscale`:

```
kind=arm
tier=A
status=unknown
name=15_dofscale
created=2026-08-28T20:41:07Z
by=claude-code
engine_sha=5ecbb89eebede8cba9271cbdd386496583bfbe2178cd9ded3db0f5512f11b511
```

Same directory after the driver seals it at end of run:

```
status=ok
```

Backfilled onto an old tree whose engine identity is unrecoverable:

```
kind=arm
tier=A
status=partial
name=arm13_intragenes
created=2026-08-21T14:03:00Z
by=unknown
engine_sha=unknown
```

A backfilled record is identifiable by its `unknown` values. That is a feature:
the census can report "how much of the corpus has recoverable provenance" as a
single `grep -c`.

---

## Census one-liners

```sh
# every enumerable run
find "$FLEXAIDDS_RESULTS" -name KIND

# only the ones that finished
find "$FLEXAIDDS_RESULTS" -name KIND -exec grep -l '^status=ok$' {} +

# runs with unrecoverable engine identity
find "$FLEXAIDDS_RESULTS" -name KIND -exec grep -l '^engine_sha=unknown$' {} +

# which seat produced what
find "$FLEXAIDDS_RESULTS" -name KIND -exec grep -h '^by=' {} + | sort | uniq -c

# every run of one engine build
find "$FLEXAIDDS_RESULTS" -name KIND -exec grep -l "^engine_sha=$SHA$" {} +
```

The last one is the query that was impossible on 2026-08-28.

---

## Writers

| Writer | When | `status` |
|---|---|---|
| `scripts/driver_layout.sh` → `new_run_dir` | at directory creation | always `unknown` |
| the driver itself | at end of run | seals to `ok`/`partial`/`void` |
| `scripts/backfill_kind.sh` | retroactively, create-only | inferred, default `unknown` |

Backfill never overwrites an existing `KIND` — the create is `O_EXCL`, so a
native record always wins over an inferred one.
