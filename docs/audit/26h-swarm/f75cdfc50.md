# Audit: f75cdfc50 — Fix: Production CloudDocs anti-hang I/O for ops and agents

**Commit:** `f75cdfc50a2322a76e5f9acc3cda575ca39ecefa`  
**Short:** `f75cdfc50`  
**AuthorDate:** 2026-07-15 06:30:14 -0400  
**Audit date:** 2026-07-15  
**Scope:** git show + static review + live repro of isolation/timeout/reaper semantics. **Report only — no source edits.**

## Summary

This commit introduces a production “anti-hang” toolkit for iCloud Drive (FileProvider) abuse: `scripts/icloud_safe_io.py` (timeout-bounded read/md5/materialize), a local-first rewrite of `scripts/benchmark_ops_monitor.py`, a wall-clock wrapper in `scripts/run_benchmark_ops_monitor.sh`, and `scripts/reap_hung_icloud_walkers.sh` to kill stuck CloudDocs walkers while claiming never to kill docking. The **direction is correct** and aligns with `AGENTS.md` local-first / thin-iCloud policy.

However, the **core isolation primitive is broken in two independent ways**: (1) `ProcessPoolExecutor` + `fut.result(timeout=…)` does **not** kill a stuck worker — `__exit__` → `shutdown(wait=True)` re-blocks the parent, so the advertised hard timeout is false under real FileProvider stalls; (2) CloudDocs paths for `safe_read_bytes` / `materialize` / `safe_glob_result_csvs` submit a **lambda** or **nested function** that cannot be pickled under spawn, so those APIs fail immediately (swallowed as `None`/`[]`) rather than performing isolated I/O. The reaper’s “never kill docking” claim holds for normal claim binaries (`bin/C/FlexAIDdS`, `benchmark_datasets`, `claim_icloud_sync_loop`) but has pattern gaps and can reap legitimate `python …CloudDocs…` tools including `icloud_safe_io.py` itself. Residual unguarded `is_dir`/`open` on iCloud fallbacks in the ops monitor can still stick the process; the outer 90s wall-clock is the real safety net, and the macOS non-GNU-timeout path is weaker (TERM only, no KILL).

## Severity: HIGH

| Dimension | Rating | Notes |
|-----------|--------|-------|
| Overall | **HIGH** | Timeout contract is false under hang; CloudDocs materialize/read/glob paths non-functional via pickle |
| Docking / ranking impact | **None** | No `LIB/` / scoring / GA / StatMech changes |
| Ops / agent hang risk | **HIGH** | Broken isolation + residual direct CloudDocs I/O |
| Accidental kill of live compute | **LOW–MEDIUM** | Claim binaries protected; broad `*python*CloudDocs*` can reap safe tools / agent work |
| Test coverage in this commit | **Missing** | No tests in `f75cdfc50`; later `tests/test_icloud_safe_io.py` only covers local/path-string cases |

## Scope reviewed

| Path | Change | Role |
|------|--------|------|
| `scripts/icloud_safe_io.py` | **Added** (342 lines) | Timeout isolation, pin-cache materialize, CloudDocs detect |
| `scripts/reap_hung_icloud_walkers.sh` | **Added** (103 lines) | Kill hung find/md5/rglob on CloudDocs |
| `scripts/run_benchmark_ops_monitor.sh` | Modified | Local-first env, `--reap-walkers`, 90s wall-clock |
| `scripts/benchmark_ops_monitor.py` | Modified | Local campaign scan, no `**/` rglob, shallow dock_config |
| `docs/ICLOUD_BENCHMARK_STORAGE.md` | Modified | Documents toolkit + agent rules |

## What the commit fixed well

1. **Policy alignment** with local APFS live work + thin iCloud mirror; removes hard dependency on `use_icloud_benchmark_storage.sh` for the monitor wrapper.
2. **Ops monitor local-first campaign roots** (`~/flexaidds_results/campaigns/<name>` before iCloud `results/…`).
3. **Shallow globs only** for campaign CSV (`*/result.csv`) and dock_config (`dock_config.json`, `r*/dock_config.json`) — eliminates `**/` rglob from the monitor.
4. **Avoids `Path.resolve()`** on CloudDocs roots in several path helpers (resolve itself can hang).
5. **Reaper default posture** protects `bin/C/FlexAIDdS`, `benchmark_datasets`, `run_C0_claim*`, `claim_icloud_sync_loop`, and build tools; age gate `MIN_AGE_SEC=90`.
6. **Outer monitor wall-clock** (`FLEXAIDDS_MONITOR_TIMEOUT_SEC`, default 90) is the right last line of defense for cron.
7. **Docs** state hard agent rules: no `find`/`rglob` under `Mobile Documents/`, hash via safe CLI, prefer `$FLEXAIDDS_LOCAL_ROOT`.
8. **Scratch/logs** moved off machine-specific `/var/folders/.../grok-goal-...` default to `$FLEXAIDDS_LOCAL_ROOT/logs/ops_monitor`.

## Findings

### F1. CRITICAL — `_run_isolated` timeout does not kill stuck workers (timeout is a lie under hang)

**File:** `scripts/icloud_safe_io.py` (`_run_isolated`)

```python
def _run_isolated(fn, arg: str, timeout_s: float, default=None):
    try:
        with ProcessPoolExecutor(max_workers=1) as ex:
            fut = ex.submit(fn, arg)
            return fut.result(timeout=timeout_s)
    except FuturesTimeout:
        return default
    except Exception:
        return default
```

**Mechanism:** When `fut.result(timeout=…)` raises `TimeoutError`, the `with` block still runs `ProcessPoolExecutor.__exit__` → `shutdown(wait=True)` **before** the outer `except FuturesTimeout` can return. A worker blocked in an uninterruptible/long FileProvider `open`/`read` keeps the parent blocked indefinitely. Docstring claims “kill on timeout”; implementation does not call `Process.kill` / `terminate` / `shutdown(wait=False, cancel_futures=True)` + explicit kill.

**Repro (session):** Module-level hang worker under `ProcessPoolExecutor`; `fut.result(timeout=0.5)` then context exit — parent still running after 2.0–2.5s wall; outer `subprocess` TimeoutExpired. Comment in repro: *ProcessPoolExecutor context does not return after FuturesTimeout while worker hangs*.

**Impact:** Any CloudDocs path that actually stalls I/O will still hang the calling agent/ops process for the full stall duration (plus orphaned pool workers if shutdown is eventually forced from outside).

**Fix direction (not applied):** Prefer `multiprocessing.Process` + `join(timeout)` + `terminate`/`kill`, or `subprocess` with `timeout=` running a tiny worker module; never rely on `ProcessPoolExecutor` context exit for kill semantics. Always reap the child PID on timeout.

---

### F2. CRITICAL — CloudDocs `safe_read_bytes` / `materialize` / `safe_glob_result_csvs` use non-picklable callables

**File:** `scripts/icloud_safe_io.py`

| API | Submitted callable | Spawn pickle |
|-----|-------------------|--------------|
| `safe_read_bytes` (CloudDocs) | `lambda s: _worker_read(s, max_bytes)` | **Fails** (`PicklingError`) |
| `safe_glob_result_csvs` (CloudDocs) | nested `_glob` | **Fails** (`PicklingError`) |
| `safe_md5` / `safe_sha256` / `safe_exists` | module-level `_worker_*` | OK if module importable |

Bare `except Exception: return default` **swallows** pickle failures, so callers see `None` / `[]` in ~50–150 ms and believe a “timeout,” while **no I/O was attempted**.

**CLI chain breakage:** `python3 scripts/icloud_safe_io.py md5|materialize <CloudDocs path>` calls `materialize` → `safe_read_bytes` (lambda) → always `TIMEOUT_OR_ERROR` for CloudDocs even for missing paths (session: `materialize` → `TIMEOUT_OR_ERROR`, real ≈ 0.15s). The documented agent hash path is non-functional for real CloudDocs files.

**Top-level workers alone are not enough:** even after pickling is fixed, F1 still applies on true hangs.

---

### F3. HIGH — Ops monitor still performs unguarded CloudDocs syscalls on iCloud fallback

**File:** `scripts/benchmark_ops_monitor.py`

When local campaign dir is absent, root becomes iCloud:

- `scan_campaign`: `root.is_dir()` — **not** isolated; can block forever.
- After glob (if it returned paths): `parse_result_csv` → `path.open()` — **not** isolated.
- `pid_file(q / "logs/…")` when queue is still under CloudDocs: `path.is_file()` + `read_text()` — unguarded.
- `queue_root()` may still point at `$FLEXAIDDS_ICLOUD/queues/...` when local queue staging is missing.

Local-first reduces frequency, but cold machines / missing local layout re-expose hang. Wrapper 90s timeout (F6) is the only backstop.

Also: `safe_glob_result_csvs` on CloudDocs is broken (F2), so iCloud campaign scan returns **empty results** without distinguishing “no data” vs “isolation failed.”

---

### F4. MEDIUM — Reaper protection matrix: docking mostly safe; false positives and gaps

**File:** `scripts/reap_hung_icloud_walkers.sh`

**Simulated matrix (patterns from commit):**

| Process class | Protected? | Classified walker? | Action |
|---------------|------------|--------------------|--------|
| `…/bin/C/FlexAIDdS …` | YES | no | skip |
| `…/benchmark_datasets …` | YES | no | skip |
| `caffeinate … benchmark_datasets` | YES | no | skip |
| `claim_icloud_sync_loop.sh` | YES | no | skip |
| `find …/Mobile Documents/…` | no | YES | **REAP** |
| `md5 …/CloudDocs/…` | no | YES | **REAP** |
| `python3 …/icloud_safe_io.py md5 '…CloudDocs…'` | no | YES (`*python*CloudDocs*`) | **REAP** after age |
| `FlexAIDdS` without `bin/C/` path | no | no (unless CloudDocs+md5-like) | ignore |
| `sync_claim_local_to_icloud.sh` / `rsync …CloudDocs…` | no | no | ignore (not walker) |
| `find ~ -path '*CloudDocs*'` without marker substrings | no | **no** | **missed hang** |

**Docking safety (AGENTS.md):** Normal claim layout uses `$FLEXAIDDS_LOCAL_QUEUE/bin/C/FlexAIDdS` and `benchmark_datasets` — **protected**. Reaper does **not** match bare docking without CloudDocs walker patterns. **Does not kill claim caffeinate runners** when cmdline includes FlexAIDdS/benchmark_datasets.

**Gaps / risks:**

1. Protection is path-shape based (`*bin/C/FlexAIDdS*`), not binary-name based — nonstandard install paths are unprotected *if* they somehow also look like walkers (unlikely for pure dock argv).
2. `*python*CloudDocs*` is **very broad**: any Python with CloudDocs in argv older than 90s is killed — including the safe I/O tool, agent notebooks, and future Python sync helpers that are **not** named `claim_icloud_sync_loop`.
3. Does not protect `sync_claim_local_to_icloud.sh` by name (usually OK — rsync not classified as walker).
4. Misses hung `find` that never puts `Mobile Documents` / `com~apple~CloudDocs` in the argv string.
5. Always `kill -TERM` then 0.5s later `kill -KILL` with no check that PID still matches the same command (classic TOCTOU; low practical risk for 0.5s window).
6. `ps` parsing via `awk`/`sed` is fragile for exotic `etime` but standard macOS `[[dd-]hh:]mm:ss` works; bare `"59"` would mis-parse as 59 minutes (not produced by `ps etime`).

**Verdict on “never kills docking”:** **Mostly true for production claim binaries;** not a formal guarantee for all FlexAIDdS invocations or all protected-class relatives.

---

### F5. MEDIUM — macOS wall-clock path in `run_benchmark_ops_monitor.sh` is incomplete

**File:** `scripts/run_benchmark_ops_monitor.sh`

```bash
TIMEOUT_SEC="${FLEXAIDDS_MONITOR_TIMEOUT_SEC:-90}"
if command -v timeout >/dev/null 2>&1; then
  exec timeout "$TIMEOUT_SEC" python3 ... ${@//--reap-walkers/}
else
  python3 ... &
  mon_pid=$!
  ( sleep "$TIMEOUT_SEC"; kill -TERM "$mon_pid" 2>/dev/null || true ) &
  wait "$mon_pid" 2>/dev/null || true
fi
```

Issues:

1. Stock macOS often **lacks** GNU `timeout` → fallback path is default on this platform.
2. Fallback sends **TERM only**, never **KILL** — insufficient if Python is stuck in non-alertable I/O wait.
3. Killer `sleep` subshell is not tracked/cleaned; races if monitor exits early.
4. `wait … \|\| true` forces exit 0 even on timeout/failure — cron cannot detect stuck-killed runs.
5. `${@//--reap-walkers/}` is unquoted → word-splitting; also substring-strips inside other args.
6. Reaper runs **before** monitor without its own wall-clock; usually fast (`ps` scan only).

GNU `timeout` path is fine where coreutils is installed.

---

### F6. LOW–MEDIUM — `safe_exists` CloudDocs worker only checks files

`_worker_exists` uses `os.path.isfile` only, while local branch uses `is_file() or is_dir()`. CloudDocs **directories** report non-existent via `safe_exists` even when present. Misleading for callers that gate on dir existence.

---

### F7. LOW — Pin-cache unbounded growth; silent provenance

`materialize` writes `$FLEXAIDDS_LOCAL_ROOT/pins/materialize/{sha16}_{name}` with no TTL/eviction. Sidecar `.src.txt` write failures are ignored. Collisions on `name` are avoided by content-key prefix; same path rewrite without `force=True` returns stale pin if CloudDocs content changed (no mtime/size check).

---

### F8. LOW — `is_clouddocs` is substring-only

Markers: `"Mobile Documents/com~apple~CloudDocs"` only. No `resolve()`, good. False negatives: other FileProvider roots, unicode/normalization variants, paths already materialized elsewhere. False positives: unlikely path that embeds the marker string. Acceptable heuristic; should not be the only gate for “safe to block.”

---

### F9. LOW — Swallowed exceptions hide diagnostics

`_run_isolated` maps all failures to `default` with no logging. Operators cannot tell pickle bug vs timeout vs ENOENT vs permission. CLI only prints `TIMEOUT_OR_ERROR`.

---

### F10. LOW — No tests in this commit for the anti-hang contract

`f75cdfc50` adds no tests. Later tree has `tests/test_icloud_safe_io.py` (local paths, marker detection, one-level glob) — **does not** exercise:

- hung worker + wall-clock return,
- CloudDocs materialize/read success,
- pickle of workers,
- reaper dry-run protection matrix,
- monitor 90s timeout.

## Ranking / Repro / Tests

### Ranking / science

- **No change** to pose ranking, CF scoring, clustering, StatMech, or DatasetRunner election.
- Safe to treat as pure ops/infra from a thermodynamics perspective.

### Repro performed this session

| Check | Result |
|-------|--------|
| `ProcessPoolExecutor` + hung worker + `result(timeout=0.5)` | Parent still blocked >2s on context exit — **F1 confirmed** |
| Lambda / nested submit pickling | `PicklingError` — **F2 confirmed** |
| `safe_read_bytes` / `materialize` / `safe_glob` on CloudDocs-marker path | Immediate `None`/`[]`/`TIMEOUT_OR_ERROR` (~0.05–0.15s) |
| CLI `is-cloud` | `yes` / `no` correct for marker vs `/tmp` |
| Reaper pattern matrix (shell simulation) | Claim dock protected; `icloud_safe_io.py md5 CloudDocs` would REAP; bare FlexAIDdS not protected but not walker |
| etime parser | `01:02`, `1:02:03`, `2-03:04:05` OK |

### Tests

- **In commit:** none.
- **Suggested gates (for a follow-up fix commit):**
  1. Unit: picklable top-level workers for read/glob with max_bytes arg as tuple or `functools.partial` of module function.
  2. Unit: fake hang worker; assert parent returns within `timeout_s + ε` **and** child PID is dead.
  3. Integration: materialize from a local path rewritten to contain CloudDocs marker via symlink **or** inject fake `open` in worker (no real iCloud required).
  4. Shell: `reap_hung_icloud_walkers.sh --dry-run` fixture lines must never mark `bin/C/FlexAIDdS` / `benchmark_datasets` / `claim_icloud_sync_loop` as REAP.
  5. Wrapper: macOS fallback must escalate to KILL and non-zero exit on timeout.

## Residual hang paths (checklist)

| Path | Guarded? |
|------|----------|
| Monitor local `*/result.csv` | Yes (local APFS) |
| Monitor iCloud `root.is_dir` | **No** |
| Monitor iCloud `parse_result_csv` open | **No** |
| Monitor iCloud pid_file | **No** |
| `icloud_safe_io` CloudDocs read under true stall | **No** (F1); currently often fails earlier (F2) |
| Wrapper `timeout` (GNU) | Yes |
| Wrapper macOS fallback | Partial (TERM only) |
| Reaper itself | `ps` only — OK |

## Verdict

| Question | Answer |
|----------|--------|
| Intent correct? | **Yes** — local-first ops + kill only CloudDocs walkers |
| Implementation meets “hard timeout / never hang”? | **No** — F1 + F2 invalidate the isolation primitive |
| Safe for live docking when reaper enabled? | **Mostly yes** for standard claim binaries; not a formal proof |
| Ranking/thermo risk? | **None** |
| Production-ready as sole anti-hang control? | **No** — keep outer process managers; **do not trust** `icloud_safe_io` alone on CloudDocs until F1/F2 fixed |
| Ship recommendation | **Fix-forward required** before advertising agent/cron CloudDocs I/O as safe; local-first monitor changes are still net-positive as hang *reduction* |

**One-line verdict:** Valuable policy and local-first ops refactor, but the CloudDocs isolation core is **incorrect under hang and non-functional for materialize/read/glob**, so the commit overclaims production anti-hang safety.
