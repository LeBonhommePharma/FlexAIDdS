# Launcher defect: `launch_seed2.sh` splices a comment onto its `env` line

**Status: UNFIXED as of 2026-09-23.** The script has not been repaired. Nothing
in this note has landed in the campaign directory. If you are reading this to
decide whether the `c1_off` arm can be launched as-is: it cannot.

Script: `~/flexaidds_results/campaigns/wall_paired_85_seed2/launch_seed2.sh`
(mtime 2026-09-20 17:03 local). That tree is research data. Read it; change
nothing in it. Owner of launches: Claude Science.

Everything below was re-verified against that file and against `main` on
2026-09-23. Where an earlier chat-transcript claim did not survive
verification, the section says so.

---

## 1. The defect

Lines 113–128 of the script are one intended command written as a
backslash-continued block, with a comment block inserted in the middle:

```bash
env FLEXAIDDS_NO_SEC=1 FLEXAIDDS_ALLOW_LOCAL_OUT=1 \
    # TMPDIR MUST NOT POINT INTO THE AGENT SESSION WORKSPACE. ...
    # ... (seven more comment lines, the last one NOT ending in a backslash)
    TMPDIR=/Users/lp.more/flexaidds_results/tmp_engine TMP=... TEMP=... \
    OMP_NUM_THREADS=2 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    FLEXAIDDS_RESTARTS="$RESTARTS" FLEXAIDDS_MAX_CONCURRENT_RESTARTS=1 \
    FLEXAIDDS_SEED_BASE="$SEED_BASE" $WALENV \
    "$BIN" --benchmark astex_diverse --only-codes "$T" \
    --cache "$CACHE" --output "$d" --threads 1 --omp-threads 2 \
    > "$d/run.log" 2>&1
```

Bash removes the backslash-newline on line 113, so the first comment line is
joined onto the `env` line and becomes a trailing comment. The comment block's
last line has no continuation, so the shell sees **two commands**:

1. `env FLEXAIDDS_NO_SEC=1 FLEXAIDDS_ALLOW_LOCAL_OUT=1` with no program to
   run. `env` with only assignments prints the **entire process environment
   to stdout** and exits 0. That dump lands in the shard log
   (`w_*.log`, `log_rerun_*.log`, `log_fill_*.log` all contain it).
2. `TMPDIR=... OMP_NUM_THREADS=2 ... FLEXAIDDS_SEED_BASE=... $WALENV "$BIN" ...`
   as a plain assignment-prefixed command, with the `> run.log 2>&1` redirect.

So `TMPDIR`, the thread pins, `FLEXAIDDS_RESTARTS`,
`FLEXAIDDS_MAX_CONCURRENT_RESTARTS` and `FLEXAIDDS_SEED_BASE` **are** applied.
Only `FLEXAIDDS_NO_SEC` and `FLEXAIDDS_ALLOW_LOCAL_OUT` fall off the command.

Reproduced 2026-09-23 with a stub binary in place of `$BIN`, same block shape:

```
--- c1_on, WALENV=''
FAKEBIN ran NO_SEC=unset LOCAL_OUT=unset WAL=unset TMPDIR=/x
rc=0
--- c1_off, WALENV='FLEXAIDDS_WAL_C1=0'
line 6: FLEXAIDDS_WAL_C1=0: command not found
```

## 2. Why the two arms fail differently

`$WALENV` is expanded **after** the shell has decided which leading words are
assignments. Its expansion is therefore never treated as an assignment; it is
the first non-assignment word, i.e. the command name.

| Arm | `WALENV` | First non-assignment word | Result |
|---|---|---|---|
| `c1_on` | `""` (expands to nothing) | `"$BIN"` | runs, `rc=0` |
| `c1_off` | `FLEXAIDDS_WAL_C1=0` | `FLEXAIDDS_WAL_C1=0` | `command not found`, `rc=127` |

This is why `c1_on` has a result set on disk and the bug was not noticed.
Verified on disk 2026-09-23: 77 target directories under `c1_on/`, 72 with a
top-level `result.csv`, every `CELL_PROVENANCE.json` reading
`"wall_env":"default_c1_on"`. There is no `c1_off/` directory: that arm has
never been launched.

## 3. What a `c1_off` launch would produce

With the script as it stands, each `c1_off` cell would be a directory with
`run.log` containing one line (`command not found`), `rc=127`, zero poses, and
a `CELL_PROVENANCE.json` asserting `"wall_env":"FLEXAIDDS_WAL_C1=0"` as though
the condition had been applied. The receipt would attest to a treatment that
never reached the engine. Treat any `c1_off` cell with `rc=127` as this
defect, not as an engine failure.

## 4. The trap: the obvious repair changes the design

The obvious fix is to delete the comment block (or move it above line 113) so
`env` reattaches to the command. **Before doing that, read this section.**

The earlier transcript claimed that `c1_on` "never had" `FLEXAIDDS_NO_SEC=1`
and `FLEXAIDDS_ALLOW_LOCAL_OUT=1`, so reattaching `env` would add them to
`c1_off` only and put a second difference between the arms.

**That claim did not survive verification.** The `c1_on` arm *did* run with
both variables set. They came from the launching session's environment, not
from the script. Two independent pieces of evidence:

- The environment dump that the orphaned `env` writes into every shard log
  (section 1) records the parent environment at launch time. In
  `w_2hr7_20260922T002518Z.log`, `log_rerun_s0.log` and
  `w_20260921T164406Z_s1.log` it contains `FLEXAIDDS_NO_SEC=1` and
  `FLEXAIDDS_ALLOW_LOCAL_OUT=1`, alongside `OMP_NUM_THREADS=8` and a `TMPDIR`
  inside a `~/.claude-science/.../workspaces/.../.tmp` directory. The launches
  were made from a Claude Science session whose shell already exported both
  flags.
- The engine confirms it took effect. `LIB/gaboom.cpp` prints
  `[SEC] All entropy-convergence early exits DISABLED (FLEXAIDDS_NO_SEC=1)`
  to stderr when the flag is present. That banner is in the restart
  `stdout.log`s of 75 of the 76 `c1_on` cells that have provenance (the one
  without is `2HR7`, an `rc=2`, zero-pose cell).

So the principle behind the trap stands, but the direction is inverted:

- `FLEXAIDDS_NO_SEC` is a **presence** gate (`env_present_any` in
  `LIB/ProtocolConfig.cpp`). Any value, including empty, disables every
  entropy-convergence early exit and runs the GA to `max_generations`. It
  changes the search budget. An arm with it and an arm without it are not
  paired.
- `c1_on` was run **with** it. Therefore `c1_off` must also run **with** it,
  and with `FLEXAIDDS_ALLOW_LOCAL_OUT=1` (a storage-policy gate read by
  `scripts/require_icloud_out.sh`; harmless to the engine but part of the
  recorded environment).
- Whether `c1_off` gets them today depends on **what the launching shell
  happens to export**, not on the script. A launch from a clean shell would
  silently produce a `c1_off` arm with early exits enabled, paired against a
  `c1_on` arm that had them disabled. Neither the provenance JSON nor the
  shard header records either flag, so the mismatch would be invisible.

The transcript's proposed "correct form" (drop `env`, put only the thread,
temp-dir, restart, seed and wall assignments on the prefix) would therefore
**also** be wrong: it would omit `FLEXAIDDS_NO_SEC` from `c1_off` on purpose,
and the arms would differ in search budget by construction.

### What the block should look like

Pin both flags explicitly, for both arms, on the same assignment prefix as
everything else, and keep no comment lines inside the continuation:

```bash
# (comment block moved ABOVE this command; nothing between the lines below)
FLEXAIDDS_NO_SEC=1 FLEXAIDDS_ALLOW_LOCAL_OUT=1 \
TMPDIR=/Users/lp.more/flexaidds_results/tmp_engine TMP=/Users/lp.more/flexaidds_results/tmp_engine TEMP=/Users/lp.more/flexaidds_results/tmp_engine \
OMP_NUM_THREADS=2 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
FLEXAIDDS_RESTARTS="$RESTARTS" FLEXAIDDS_MAX_CONCURRENT_RESTARTS=1 \
FLEXAIDDS_SEED_BASE="$SEED_BASE" FLEXAIDDS_WAL_C1="$WALC1" \
"$BIN" --benchmark astex_diverse --only-codes "$T" \
--cache "$CACHE" --output "$d" --threads 1 --omp-threads 2 \
> "$d/run.log" 2>&1
```

with the `case` changed from a whole-word `WALENV` to a value:

```bash
case "$ARM" in
  c1_on)  WALC1=1 ;;   # the shipped default, made explicit
  c1_off) WALC1=0 ;;
esac
```

`env` is unnecessary here; an assignment prefix on a simple command exports
those names to that command only. Do not put `$WALENV` back as a bare word.
Record `FLEXAIDDS_NO_SEC` and `FLEXAIDDS_WAL_C1` in `CELL_PROVENANCE.json`
from the value actually passed, and after the first `c1_off` cell finishes,
confirm the `[SEC] ... DISABLED` banner is in its `stdout.log`s before
launching the rest.

Whether to re-run `c1_on` under the fixed script, or accept the existing
`c1_on` cells (which had the correct environment by accident), is a design
decision for the campaign owner. The existing cells are usable only if the
fixed `c1_off` reproduces the same flag set, which the banner check verifies.

## 5. Invocation contract

```
launch_seed2.sh ARM SHARD [NSHARD]
```

- `ARM` is `c1_on` or `c1_off`. Required.
- `SHARD` is `0` or `1`. Required (`${2:?}`); the script aborts without it.
- `NSHARD` defaults to `2`.

Target `i` (1-based, sorted cache order) runs only if `(i-1) % NSHARD == SHARD`.
**One invocation covers half the 85 targets.** A full arm is two concurrent
invocations, `SHARD=0` and `SHARD=1`, which is what `relaunch.sh` does. A single
invocation completes, writes `SHARD_<n>_DONE`, and leaves 42 or 43 targets
untouched with no warning.

The resume rule is `[ -f "$d/result.csv" ]` per target directory, so a re-run
skips completed cells and fills the rest.

## 6. Semantics that make `c1_off` meaningful

Verified on `main` at 30ad3ac3 and on the campaign build reference (`afe8f58a`
per provenance):

- `LIB/vcfunction.cpp:310`
  `static const bool wal_c1 = flexaids::env_bool("FLEXAIDDS_WAL_C1", true);`
  The transcript cited line 257; that line is now inside the flip-record
  comment above the definition. Unset is equivalent to `=1`. The comment at
  lines 304–308 records that the old byte-0 read treated `FLEXAIDDS_WAL_C1=""`
  as ENABLED and that `=0` is the documented way to reach the legacy wall.
- `LIB/EnvFlags.h:15–17` and `env_bool_str` (line 31): the default is
  returned when the variable is **unset, empty, or unparseable**, after
  whitespace trimming. `c1_off` must therefore pass literally `0` (or
  `false`/`no`/`off`). An empty value, a stray space-only value, or a typo
  silently yields `c1_on`.

## 7. Other things wrong in the script (beyond the transcript)

- **Malformed provenance on zero-pose cells.** Line 131:
  `np=$(find "$d" -name '*.pdb' | grep -vc _INI || echo 0)`. When there are no
  poses, `grep -vc` prints `0` and exits 1, so `|| echo 0` prints a second `0`.
  `$np` becomes `0<newline>0` and the JSON is invalid. On disk:
  `c1_on/2HR7/CELL_PROVENANCE.json` reads `"poses":0` followed by `0,` on the
  next line. Any parser walking the campaign will choke on that cell.
- **The environment dump is a leak.** Every iteration of the loop writes the
  launcher's full environment into the shard log. Those logs are in the
  campaign tree; anything exported in the launching shell is now recorded
  there.
- **The provenance does not record the flags that matter.** `wall_env`
  records the string in `$WALENV`, not what the engine saw, and nothing
  records `FLEXAIDDS_NO_SEC` at all. Section 4 explains why that hides an
  arm mismatch.
- **Sibling script `launch_seed2_redock.sh`** (mtime 2026-09-23 02:09) has
  `env` correctly attached and no comment splice, but still passes `$WALENV`
  as a bare word, so its `c1_off` arm would fail with the same
  `command not found`.
- The 15 `rc=2` cells in `c1_on` are the `temp_directory_path` sweep failure
  described in the script's own comment block, not this defect. They are
  separate and already understood.

## 8. Checklist before any launch of this campaign

1. Confirm the `env` line and the comment block have been separated
   (`grep -n '^ *env ' launch_seed2.sh` should show nothing, or a line whose
   command is on the same logical line).
2. Confirm `FLEXAIDDS_NO_SEC=1` is on the command prefix for **both** arms.
3. Confirm `FLEXAIDDS_WAL_C1` is passed as `NAME=value`, never as `$WALENV`.
4. Launch both shards.
5. After the first `c1_off` cell: `run.log` is not one line, `rc` is not 127,
   and the `[SEC] ... DISABLED` banner is present in its `stdout.log`s.
6. Check `CELL_PROVENANCE.json` parses as JSON for a zero-pose cell.
