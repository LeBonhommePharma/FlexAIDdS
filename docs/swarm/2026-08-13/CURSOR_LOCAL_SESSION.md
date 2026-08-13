# Cursor local multitask — paste this into the new session

This file is the adapter between the Claude Science swarm pack (this directory,
byte-identical to the 2026-08-13 tarball; SHA256SUMS.txt verifies) and a **Cursor
Desktop local session with multitask / parallel agents**.

The cloud agent that prepared this pack could not be that session: it is a single
background Cloud Agent on one worktree, the Task-tool enum here has no
`ci-investigator`, and five lanes cannot share one git index. **You** are that
session. Spawn the five lanes. Do not re-plan the swarm.

---

## PASTE INTO THE NEW CURSOR LOCAL CHAT (from the `---` below to the end of §SPAWN)

---

You are the **parent agent** of a Cursor local **multitask** session on FlexAIDdS.

**Do not implement any lane yourself. Do not merge. Do not rebase one lane onto
another. Do not launch a docking campaign.**

Git is the orchestrator. `score_canonical.py` is the referee. Claude Science owns
the box and every campaign. Claude Code (a different seat) owns the merge queue
later. You only **spawn five parallel lane agents** and then stop.

### 0. Checkout

Work from **`main`**, not from PR #419 (`cursor/audit-mt-reproducibility-science-1c2f`).
That PR is a docs-only multithreading/reproducibility audit. It is a different
job. Ignore its `claude-review` red check: `cursor[bot]` is not on
`allowed_bots` in `.github/workflows/claude-code-review.yml`. Unrelated, not a flake.

Pack root (repo-relative):

```
docs/swarm/2026-08-13/
```

Read, in this order, then spawn:

1. `docs/swarm/2026-08-13/SWARM_COMMON_PREAMBLE.md`
2. `docs/swarm/2026-08-13/SWARM_ASSIGNMENT.md`
3. `docs/swarm/2026-08-13/SWARM_ORCHESTRATION.md`
4. the five lane files listed in §SPAWN
5. this file (path resolution + spawn recipe)

### 1. Referee (must print these numbers before anyone codes)

From the repo root:

```bash
python3 docs/swarm/2026-08-13/score_canonical.py \
  --frozen docs/swarm/2026-08-13/ASTEX84_FROZEN_POSE_BENCHMARK.csv
```

It MUST print:

```
min-CF election   26/84 =  31.0%
pool ceiling      41/84 =  48.8%
selection gap     15 targets
```

If it prints anything else, **stop**. The checkout or the frozen file is wrong.
`--frozen` does not need spyrmsd (RMSDs are already in the CSV). `--run` does,
and refuses to fall back.

Science-box live data (campaigns, cache, sites) still lives under
`$FLEXAIDDS_RESULTS` (on LP's Mac that is the iCloud/local-first tree documented
in `AGENTS.md` / `docs/ICLOUD_BENCHMARK_STORAGE.md`). Do **not** bake
`/Users/<name>/...` into new committed files. For `--run`:

```bash
python3 docs/swarm/2026-08-13/score_canonical.py \
  --run "$FLEXAIDDS_RESULTS/<campaign_dir>" \
  --cache "${FLEXAIDDS_CACHE_V2:-$FLEXAIDDS_RESULTS/cache_v2/astex_diverse}"
```

cmake on that Mac is `/opt/homebrew/bin/cmake` (often not on PATH).

### 2. Hard rules (from the preamble — non-negotiable)

- **R1** No dataset campaign. Single-target probes for a gate are fine.
- **R2** One metric: spyRMSD graph automorphism. Never `rmsd_to_crystal` or `rmsd_hungarian`.
- **R3** Denominator N=84 (2HR7 excluded).
- **R4** `result.csv` "success" means the docking *ran*.
- **R5** Every feature **env-gated OFF**. A PR that moves a default is rejected.
- **R6** Never sum gains across fixes. One combined re-election number.
- **R7** Measure; do not infer from source.
- **R8** Stay in the lane's OWNED files. Forbidden paths stay forbidden.
- **R9** One branch per lane, names below. Do not merge.
- **R10** If a gate cannot be met, say so. "Could not measure" is useful.
- **R11** New `FLEXAIDDS_*` boolean gates use `flexaids::env_bool` (`LIB/EnvFlags.h`) only.
- **R12** Before editing, read/write `$FLEXAIDDS_LOCAL_ROOT/workorders/CLAIMED_<lane>.txt`. If another claim holds the files, stop. Do not commit those files.

Also: Apache-2.0 only. No GPL. Parity per `METHODOLOGY.md` §1 for anything that
could move a score. Separate CF scoring proxy from thermodynamic ΔG language.

### 3. Merge order (you do not merge; this is so agents do not fight)

`D → A → B → E (rebase onto D) → C`

A before B is required: `LIB/vcfunction.cpp:727` metal block reads `FA->use_elec`
and the charge fields lane A populates. E rebases onto D: both edit
`LIB/DatasetRunner.cpp`, 146 lines apart at closest approach.

Three gates before *anyone* merges later: (a) lane acceptance gates green,
(b) number from this scorer, (c) `--frozen` still 31.0% / 48.8% with the new
flag unset.

### 4. §SPAWN — emit all five Task/subagent calls in **one** assistant message

Each subagent prompt is **exactly three items**, in this order, then the
verbatim sentence:

1. Full text of `docs/swarm/2026-08-13/SWARM_COMMON_PREAMBLE.md`
2. Full text of that seat's lane file
3. Absolute or repo-relative path to `docs/swarm/2026-08-13/score_canonical.py`
   and to `docs/swarm/2026-08-13/ASTEX84_FROZEN_POSE_BENCHMARK.csv`

Then this sentence, **verbatim**:

> Work only this lane. Do not launch a docking campaign — Claude Science owns the box. Every number you report must come from score_canonical.py. Ship the feature env-gated OFF so defaults do not move the baseline. Open a PR on your own branch; do not merge.

| lane | seat (original assignment) | lane file | git branch |
|------|----------------------------|-----------|------------|
| D waste + selection | Cursor Pro | `SWARM_D_selection.md` | `lane/d-selection-waste` |
| A partial charges | Codex | `SWARM_A_charges.md` | `lane/a-partial-charges` |
| B metal coordination | Grok | `SWARM_B_metal.md` | `lane/b-metal-coord` |
| E autoflex | Claude Code Max | `SWARM_E_autoflex.md` | `lane/e-autoflex` |
| C thermo path | Claude Code 2 | `SWARM_C_thermo.md` | `lane/c-thermo-gate` |

Branch each lane **from `main`**. Do not branch from PR #419. Do not share a
worktree: if Cursor multitask does not isolate git, use separate worktrees
(`git worktree add`).

Tell every agent:

- OWNED / FORBIDDEN paths are in its lane file. Honour them.
- Env gates (all default OFF):
  - D: `FLEXAIDDS_DROP_SENTINEL_POSES`, `FLEXAIDDS_REP_MINCF`, `FLEXAIDDS_VOID_ONLY_IF_ALL_FAIL`, `FLEXAIDDS_SPREAD_GUARD_ALWAYS`
  - A: `FLEXAIDDS_PARTIAL_CHARGES` **and** existing `scoring.electrostatics_enabled` (both required)
  - B: `FLEXAIDDS_METAL_FIX` — diagnose upstream of `vcfunction.cpp`; do not patch the 688–742 window until A has landed
  - E: `FLEXAIDDS_FLEX_SIDECHAINS` must dominate; do not flip `autoflex_enabled`'s default
  - C: `FLEXAIDDS_THERMO_FIX`; banner may always print; **do not wire TdS_vib into election**
- Drop results as `<LANE>_RESULT.md` (D/A/B/E/C) in the local workorders dir if
  `$FLEXAIDDS_RESULTS/workorders` exists, otherwise at
  `docs/swarm/2026-08-13/results/<LANE>_RESULT.md` (create that directory).
- Open a PR. Do not merge.

After the five spawn calls: parent confirms the frozen scorer still prints
31.0% / 48.8% on this checkout, then waits. Do not "orchestrate" by rewriting
lane reports.

---

## Why the previous Cloud Agent session could not be this

- One Cloud Agent = one branch (`cursor/audit-mt-reproducibility-science-1c2f`) and one worktree.
- Cursor Desktop **multitask / multi-agent** is a local-session feature; it is not
  exposed as a mode switch inside a running cloud background agent.
- Nested Task `ci-investigator` is not in this environment's subagent enum
  (`generalPurpose`, `explore`, …). That only affected the Checks-tab
  investigation of `claude-review`, which is **unrelated** bot policy.

## What is already shipped (do not redo)

| Item | Where |
|------|--------|
| Multithreading / reproducibility / science audit | PR https://github.com/LeBonhommePharma/FlexAIDdS/pull/419 |
| Test-coverage quality (green `ctest` ≠ F1/F5/F13 coverage) | same PR, `docs/audit/2026-08-13_test_coverage_quality.md` |
| CI STRICT orphans (`EnvFlags.h`, `ShannonBinning.h`) | same PR, `build_sources.ignore` |
| This swarm pack in-repo | `docs/swarm/2026-08-13/` on branch `cursor/swarm-local-handoff-1c2f` |

Do **not** start implementing audit findings F1–F24 in the swarm lanes. The
swarm is Astex-84 landscape / selection / charges / metal / thermo-observability.
F1 (`lazy_thread_rng`) etc. stay on the audit track unless LP reassigns.

## Ladder (for context; numbers only from the scorer)

```
as-run T=300 election         15/84 = 17.9%
min-CF election               26/84 = 31.0%   <- lane D, bounded +15 targets
reweighting cap, proven               32.5%
published bar (FlexAID-2015)          45.2%
pool ceiling, 10 restarts     41/84 = 48.8%   <- cap of any selection fix
ceiling at 30 restarts        50/84 = 59.5%   (4/19 measured, extrapolated)
ceiling if A+E convert half   62/84 = 73.8%   (projection)
```

Ranking alone cannot reach the published bar. That is why A and E exist.
