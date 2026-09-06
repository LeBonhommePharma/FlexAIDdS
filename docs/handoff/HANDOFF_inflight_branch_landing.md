# Handoff — land the five remaining in-flight branches on `main`

**For:** Claude Dispatch / Claude Science
**From:** Claude Code session 2026-09-03
**Repo:** `/Users/lp.more/Projects/FlexAIDdS`
**Authoritative rules:** `AGENTS.md` (conduct), `METHODOLOGY.md` (§1 parity, §4 ctest, §5 feature-branch rule, §6 commit review)
**Original plan:** `local://inflight-branch-landing-plan.md` — steps 3–7 remain; step 1 (baseline) and step 2 are superseded by this note.

---

## 1. State you are inheriting

`main` = `7975fc33`, clean, level with `origin/main`, pushed.

```
7975fc33 Merge 2a615c3d: ParEvalWSIsolation TLS barrier
2a615c3d Tests: keep ParEvalWSIsolation TLS objects alive with a barrier
4b6c0d1a Merge origin/main: dependabot fast-uri bump
0e9e7466 Merge fix/baseline-test-gates-3: restore red baseline test gates
5a59773f Tests: --top assertion matches the design-system header
6959aaac Tests: add validator fixture with required docking_mode
73eee781 Fix: stub rot_gene_index so test_cf_aggregator links
```

**Why this exists:** the planned step-1 baseline on `c372a94d` was **red** (3 failures), so no branch was merged.
A separate repair branch `fix/baseline-test-gates-3` fixed the baseline, was gated, and landed. The five
originally-planned merges have **not** happened.

What was red on `c372a94d`, and the fix:

| failure | cause | fixed by |
|---|---|---|
| `test_cf_aggregator` link error, `rot_gene_index` undefined | `94d2d0ba` added the clamp; `tests/cf_aggregator_stubs.cpp` was never updated | `73eee781` (range-safe stub) |
| `test_cli.py::test_top_limits_table_rows` | asserted obsolete `"Binding modes: 3"`; CLI emits `"binding modes=3"` (design-system header, `26f8ee17`) | `5a59773f` |
| `test_validate_benchmark_results.py::…uses_manifest_comparator` | fixture YAML lacked required `docking_mode`; `DatasetConfig.from_yaml` fail-closes | `6959aaac` (also force-added the file — it was hidden by `.git/info/exclude`) |
| `ParallelDockTests` (found later, flaky) | `t0` could recycle its `thread_local` address before `t1` published; `EXPECT_NE(p0,p1)` racy | `2a615c3d` (`std::barrier`, both threads `arrive_and_wait`) |

### Gate evidence on `7975fc33` (do not re-derive; re-run only if you distrust it)

- **§1 parity — PASS.** `origin/main` `53715546` vs `7975fc33`, clean worktrees, both built
  `-DFLEXAIDS_GIT_COMMIT_OVERRIDE=parity-s1`. 1G9V, `/tmp/parity.json` (pop 1000 / 2000 gen,
  `pose_seed_enabled=false`, `seed_fraction=0`), unique output dirs, sanitized env.
  `d_0..d_9.pdb` **10/10 raw-byte identical** (`filecmp shallow=False`); elected CF `-237663.32736` equal.
  Record: `/tmp/gate_parity2/section1_result.txt`, env `/tmp/gate_parity2/parity.env`.
- **§4 — ctest 98/98; pytest 1275 passed, 72 skipped.** `ParallelDockTests --repeat until-fail:20` → 20/20.
- **§6 — reviewed** by `completion(model=slow)` (Fable 5 absent in this runtime, not relabeled).
  `73eee781`/`6959aaac`/`5a59773f` PASS_WITH_CAVEATS; `2a615c3d` PASS. Record: `/tmp/gate_parity/section6_review.txt`.

`/tmp/*` is volatile — copy those two records into the repo if you need durable provenance.

**Use `7975fc33` as your baseline.** It is green. Do not re-establish a baseline on `c372a94d`.

### Untouched, do not disturb

- `feat/ff14sb-lumped-charges` `235f43c4` — commit titled "WIP", 3 `LIB/` files, no CMake wiring, no tests. Not ready.
- `wip/stash0-astex85-receipt-20260828` `e56ff5cb` — parked, **do not drop** (Claude Science / receipt protocol).
- 13 stashes. Never `git stash drop`/`clear`.
- Local-only, unpushed, superseded: `fix/baseline-test-gates` and `203dce74`/`a0ef949c` (a first attempt that
  bundled two logical changes; rebuilt as three clean commits rather than rewritten). Ignore; do not push.

---

## 2. Work remaining — five branches

| branch | tip | delta vs `7975fc33` | verdict entering this handoff |
|---|---|---|---|
| `fix/stale-cli-summary-assertion` | `972939bb` | `python/tests/test_cli.py` +36/−1 | **superseded, do not merge** — see 2.1 |
| `fix/loud-spyrmsd-gate` | `3e188260` | `tests/test_rmsd_symmcorr.py`, `.github/workflows/ci.yml` +61/−4 | ready; gate was green in worktree |
| `chore/coverage-uplift` | `901d1c4c` | 3 new `tests/*.cpp` + `CMakeLists.txt`, +1146 | ready; needs **fresh configure** |
| `docs/benchmark-methods-v2` | `c0e78851` | `benchmarks/protocols/METHODS_v2.md` +435 | **blocked** — see 2.4 |
| `chore/installability-matrix` | `d5d88216` | 13 files, +1144/−27 (6 commits) | ready; **new conflict risk** — see 2.5 |

### 2.1 `fix/stale-cli-summary-assertion` — superseded, do not merge (verified)

`main` already carries this exact fix via `5a59773f`. Verified this session:

```bash
git diff main:python/tests/test_cli.py fix/stale-cli-summary-assertion:python/tests/test_cli.py
```

The only differences are **comment prose and a docstring** — every assertion is byte-identical:
`assert "binding modes=3" in out`, `assert len(data_rows) == 1`, and the `test_top_none_shows_all_rows`
counterpart all already exist on `main`. There is no behavioural delta.

**Disposition: do not merge. Record as superseded by `5a59773f` and leave the branch in place.**
The branch's longer comments (they cite `26f8ee17` and explain why the literal moved) are arguably better
documentation; if you want them, cherry-pick the comment hunk onto a feature branch as a docs-only change.
Do not merge the branch to obtain them — it would be a no-op merge commit.

### 2.2 `fix/loud-spyrmsd-gate`

Fail-closed symmcorr gate. Verified green in a throwaway worktree this session: `tests/test_rmsd_symmcorr.py`
**17 passed**, `ci.yml` parses. `spyrmsd 0.9.0` is installed here.

```bash
python3 -m pytest tests/test_rmsd_symmcorr.py -q
python3 -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml')); print('ci.yml ok')"
```

A failure *with* spyrmsd present blocks. A loud failure *without* it is the intended design — never soften the test.

### 2.3 `chore/coverage-uplift`

Adds `test_cleft_determinism`, `test_binary_snapshot`, `test_instream_clustering` (+ `CMakeLists.txt`
registration). Only step touching the build system → **fresh configure, not incremental**.

```bash
rm -rf build && cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 11
ctest --test-dir build -N | grep -E "CleftDeterminism|BinarySnapshot|InStreamClustering"
ctest --test-dir build --output-on-failure
```

Expect 98 → 101 tests. Duplicate CMake target name = configure error → abort, report the colliding name,
never drop a target to silence it. `main` already pins FetchContent GoogleTest — do not add a second source.
Note: its worktree had no `build/` this session, so the three targets are **unbuilt and unproven** — this is
the one branch whose tests have never actually run.

### 2.4 `docs/benchmark-methods-v2` — blocked as-is

`benchmarks/protocols/METHODS_v2.md` restates and **forks `METHODOLOGY.md` numbering** ("§7.1", its own
defect table, its own gate definitions). `METHODOLOGY.md` is the single source of truth and forbids exactly
this. Hygiene passes; that is not the problem.

Two acceptable outcomes — **not** a plain merge:
1. Change `METHODOLOGY.md` first, then reduce `METHODS_v2.md` to a pointer at the changed sections; or
2. Rewrite it to cite `METHODOLOGY.md §N` instead of restating, then merge.

Escalate to the user for the call. Do not merge the fork.

### 2.5 `chore/installability-matrix` — new conflict since the plan was written

It touches `tests/test_parallel_dock.cpp` (+9), and `main` just changed that file (`2a615c3d`, the TLS barrier).
**Expect a conflict there — and it is a benign one: both sides fix the same race.** Verified this session,
`git diff main...chore/installability-matrix -- tests/test_parallel_dock.cpp` shows the branch fixes
`ParEvalWSIsolation` independently with `std::latch both_published(2)` + `arrive_and_wait()` in each thread,
plus a comment noting it was "observed on macOS under `ctest -j4` (2 of 4 runs); never in isolation".
`main` uses `std::barrier sync(2)`. **Semantically equivalent** for a two-thread single-phase rendezvous
(`latch` is the tighter fit — single use).

Resolution: **keep exactly one rendezvous, not both.** Preferred — take the branch's `<latch>` version and its
explanatory comment, since it documents the observed failure mode; then delete `main`'s `<barrier>` include and
`sync` object so no duplicate synchronisation remains. Either choice is defensible; a tree with both is not.
After resolving, re-run `ctest --test-dir build -R ParallelDockTests --repeat until-fail:20` — 20/20 required.
Never resolve by dropping both and restoring the racy original.

```bash
rm -rf build && cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 11 && ctest --test-dir build --output-on-failure
pip install -e ./python && pytest python/tests -q
./build/FlexAIDdS --version
python3 -c "import json,glob; p=glob.glob('build/**/build-provenance.json',recursive=True)[0]; d=json.load(open(p)); print(d); assert d"
```

Known, expected, **not** a regression: `python3 -c "import yaml; yaml.safe_load(open('conda/meta.yaml'))"`
**fails** — `conda/meta.yaml` is a Jinja template (`{% set … %}`), invalid plain YAML. Verify it with
`conda-build`/Jinja rendering, or just confirm the file parses as a template. Container and conda install
paths stay unexercised locally; say so in the final report.

---

## 3. Procedure — non-negotiable

1. **Feature branch per repair (METHODOLOGY.md §5).** Fixes go on a branch, get gated, then merge. Never commit
   a fix directly to `main`.
2. **Test in a throwaway worktree before touching `main`:**
   `git worktree add /tmp/land-<slug> <branch>` … `git worktree remove /tmp/land-<slug>`.
3. **Merge** `git merge --no-ff <branch> -m "Merge <branch>: <reason>"`. One logical change per commit.
4. **Gate the merge tip, then push.** Zero failures before any push. A red suite blocks even if the failure predates the merge.
5. **§1 parity is required for anything touching `LIB/`, `CMakeLists.txt`, or scoring.** Two clean worktrees, same
   `-DFLEXAIDS_GIT_COMMIT_OVERRIDE`, **unique output dirs**, sanitized env (strip inherited `FLEXAIDDS_BINARY`,
   `_BUILD`, `_RESULTS`, `_RUNNER`, `_ENGINE_SHA256`, `_RUNNER_SHA256`, `_POSEBUSTERS_BIN`), compare raw bytes of
   all 10 poses + elected CF. Copy `/tmp/gate_parity2/section1_result.txt` as your template. ~35 min/engine dock.
   Test-only branches (2.1, 2.2) are exempt — record "test-only, §1 N/A".
6. **§6 review every commit** with the strongest model actually present in your runtime. If Fable 5 is absent,
   say which model reviewed — never relabel.
7. **Never** `git push --force`, `reset --hard`, `rebase`, delete branches, or rewrite history. Bad unpushed merge →
   `git revert -m 1 <sha>` after telling the user.
8. **Never edit a test or CI file to make a gate pass.** `git config core.fsmonitor` hangs git → `kill $(pgrep -f git)`.
9. `python3 scripts/check_repo_hygiene.py` after any doc/skill/agent-instruction change.
10. GitHub identity must be `LeBonhommePharma` — `gh api user --jq .login` before push work.

**Suggested order:** 2.1 decision → 2.2 → 2.3 → 2.5 → escalate 2.4. Small→large so the `CMakeLists.txt` and
`test_parallel_dock.cpp` conflicts land last, on a tree already proven green.

## 4. Done means

```bash
git status --porcelain=v1 -b     # clean, main level with origin/main
for b in fix/loud-spyrmsd-gate chore/coverage-uplift chore/installability-matrix; do
  echo "$b: $(git log --oneline main..$b | wc -l) unmerged"; done   # expect 0
ctest --test-dir build --output-on-failure && pytest python/tests -q
python3 scripts/check_repo_hygiene.py
```

Plus: `feat/ff14sb-lumped-charges` and `wip/stash0-astex85-receipt-20260828` still at `235f43c4` / `e56ff5cb`,
13 stashes intact, and an explicit written disposition for `fix/stale-cli-summary-assertion` (2.1) and
`docs/benchmark-methods-v2` (2.4) — including "not merged, and why" if that is the answer.

## 5. Traps this session actually hit

- A blank-line insertion silently deleted `def test_top_default_is_none`; an include edit dropped `<vector>`.
  **Re-read the file after every structural edit** — do not trust the edit echo alone.
- `ctest --repeat until-pass:20` stops at the **first pass** (1 run). For flake-hunting use `--repeat until-fail:20`.
- `python/tests/test_validate_benchmark_results.py` was listed in `.git/info/exclude`; `git add` silently no-ops on
  such paths. It is force-added and tracked now — check `git ls-files` before believing a file is committed.
- Engine md5s differ between two identical trees purely from the configure-time `built_utc` stamp. §1 asserts
  pose/CF equality, **not** engine-hash equality.
