# FlexAIDdS swarm pack — 2026-08-13

Canonical Claude Science handoff. Lane briefs match the tarball (`SHA256SUMS.txt`).
`score_canonical.py` is the same referee with a portable cache default (see below).

| If you are… | Start here |
|-------------|------------|
| **A new Cursor Desktop local session with multitask** | **[`CURSOR_LOCAL_SESSION.md`](CURSOR_LOCAL_SESSION.md)** — paste the marked block, spawn five lanes, do not merge |
| Claude Science / any other seat | [`HANDOFF_README.md`](HANDOFF_README.md) |

Do not edit the `SWARM_*.md` files to “improve” them. `score_canonical.py` was
ported off a hardcoded `/Users/<name>/...` cache default: it now reads
`FLEXAIDDS_CACHE_V2` or `$FLEXAIDDS_RESULTS/cache_v2/astex_diverse`, and `--run`
fails closed if neither is set. `--frozen` is unchanged. Re-run
`sha256sum -c SHA256SUMS.txt` after pulling.
