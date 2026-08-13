# FlexAIDdS swarm pack — 2026-08-13

**Canonical copy: `docs/swarm/2026-08-13/` in the FlexAIDdS git repository.**
Lane briefs, the frozen pose CSV, and `score_canonical.py` are sourced there.

`$FLEXAIDDS_LOCAL_ROOT/workorders/` (default `~/flexaidds_results/workorders/`)
is a **mirror**, not the source. If the two disagree, trust the in-repo tree
and refresh the mirror. Do not edit the mirror and copy back.

| If you are… | Start here |
|-------------|------------|
| **A new Cursor Desktop local session with multitask** | **[`CURSOR_LOCAL_SESSION.md`](CURSOR_LOCAL_SESSION.md)** — paste the marked block, spawn five lanes, do not merge |
| Claude Science / any other seat | [`HANDOFF_README.md`](HANDOFF_README.md) |

Do not edit the `SWARM_*.md` files to “improve” them. `score_canonical.py` was
ported off a hardcoded `/Users/<name>/...` cache default: it now reads
`FLEXAIDDS_CACHE_V2` or `$FLEXAIDDS_RESULTS/cache_v2/astex_diverse`, and `--run`
fails closed if neither is set. `--frozen` is unchanged. Re-run
`sha256sum -c SHA256SUMS.txt` after pulling.
