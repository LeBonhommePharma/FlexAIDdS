# Canary handoff — post PR #296 (CMA-ES merge)

**Audience:** Claude Science (relaunch owner) and any agent resuming canary work  
**Written:** 2026-07-21 17:09 local (America/Toronto)  
**Repo HEAD:** `7c901a687` — *Merge pull request #296 from LeBonhommePharma/merge/cmaes-into-main*  
**Branch:** `main` (matches `origin/main`)

---

## Status snapshot

| Item | State |
|------|--------|
| PR #296 (CMA-ES into main) | **MERGED** |
| Stale pre-merge canaries | **Killed** earlier (Claude session `5ebf01f0-…`; binary lacked gate strings) |
| Post-merge `build/FlexAIDdS` | **Rebuilt + pinned** (gates present) |
| Machine pin | `~/.flexaidds/active_build.json` + `~/.flexaidds_env` updated |
| Relaunch | **Claude Science owns relaunch** — do not double-launch from this handoff |

At handoff write time, live dock processes were already using the pinned path:

```text
…/build/FlexAIDdS … -o …/canary_pbclash_w1_20260721_post_grok/1G9V/…
```

(OUT root under `Documents/PhD/Programs/FlexAIDdS/results/canary_pbclash_w1_20260721_post_grok/`.)

---

## Pinned binary (authoritative)

| Field | Value |
|-------|--------|
| Engine path | `build/FlexAIDdS` (repo-relative) |
| Engine SHA256 | `8b85ba77a936d24e6b896598fc906a677d82cc018e90ae9a5559d7d0ffee535f` |
| Runner path | `build/benchmark_datasets` |
| Runner SHA256 | `d21af08953df778c01a7d3b482e7f6a418b4451d8077a75f88536a2dc2661d07` |
| mtime (local) | 2026-07-21 17:06:26 |
| git_head at pin | `7c901a687bf8e1d4f36edf1b08f743c4744ba757` |
| fresh vs sources | **true** |

### Required gate / feature strings (verified via `strings`)

| Marker | Count in pinned engine |
|--------|------------------------|
| `FLEXAIDDS_PB_CLASH_GRID_HOIST` | 1 |
| `FLEXAIDDS_WAL_COERCIVE` | 1 |
| `cmaes_run_dock` | ≥1 (present) |

**Do not relaunch on:**

| Binary | Why |
|--------|-----|
| Pre-17:06 `build/FlexAIDdS` (old 16:41 artifact) | STALE — hoist/WAL markers were 0; was running killed canary |
| `build_v135/FlexAIDdS` | STALE (Jul 20) |
| `.swarm/cmaes/.../build_fast/FlexAIDdS` | STALE / pre-merge swarm |
| Old pin `build_lto` SHA `d4937d1c…` | Tree gone; pin superseded |

**Smoke-only reference (not the active pin):**

| Path | SHA256 | Notes |
|------|--------|--------|
| `build_merge_smoke/FlexAIDdS` | `1a79f6ae7b58085e6e92b3524ada9acdc9abf86eecb634988689fed5b19d56f4` | Post-merge smoke tree (17:02); gates OK; different SHA than `build/` (expected) |

---

## How the pin was applied (local machine)

```bash
# From repo root after post-merge rebuild of build/FlexAIDdS + benchmark_datasets
ENGINE_SHA=$(shasum -a 256 build/FlexAIDdS | awk '{print $1}')
unset FLEXAIDDS_ENGINE_SHA256 FLEXAIDDS_BUILD FLEXAIDDS_BINARY FLEXAIDDS_BUILD_DIR
export FLEXAIDDS_BUILD="$(pwd)/build"
python3 .grok/skills/flexaidds/scripts/resolve_build.py \
  --pin-sha "$ENGINE_SHA" --write-pin --sync-env --json
python3 .grok/skills/flexaidds/scripts/resolve_build.py --check
```

Expect:

```text
OK: build resolved
  build_dir:     …/FlexAIDdS/build
  engine_sha256: 8b85ba77a936d24e6b896598fc906a677d82cc018e90ae9a5559d7d0ffee535f
  source:        pin  fresh=True
```

Pin files (machine-local, **not** committed):

- `~/.flexaidds/active_build.json`
- `~/.flexaidds_env` → `FLEXAIDDS_BINARY`, `FLEXAIDDS_ENGINE_SHA256`, etc.

---

## Pre-flight for any further canary work

1. Confirm HEAD is still `7c901a687` (or a descendant that rebuilt the same pin intentionally).
2. Re-check pin before launch:

   ```bash
   source ~/.flexaidds_env
   python3 .grok/skills/flexaidds/scripts/resolve_build.py --check
   test "$(shasum -a 256 "$FLEXAIDDS_BINARY" | awk '{print $1}')" = "$FLEXAIDDS_ENGINE_SHA256"
   strings "$FLEXAIDDS_BINARY" | grep -E 'FLEXAIDDS_PB_CLASH_GRID_HOIST|FLEXAIDDS_WAL_COERCIVE'
   ```

3. If the binary was rebuilt, **re-pin** with `--write-pin --sync-env` before claiming runs.
4. Prefer local OUT (`$FLEXAIDDS_LOCAL_ROOT` / claim staging) per `AGENTS.md` benchmark storage rules; thin-mirror to iCloud only after local success.

---

## What Claude Science owns

- Canary / sweep **relaunch** (already under way as of this handoff: `canary_pbclash_w1_20260721_post_grok`).
- Gate monitoring (PB clash hoist behavior, WAL coercive, RMSD/PoseBusters success criteria).
- Any re-pin after a further rebuild on `main`.

## What this handoff deliberately does **not** do

- Start or kill canary processes.
- Change scoring, ranking, or campaign configs.
- Force-push or rewrite history.

---

## Related evidence

| Artifact | Path |
|----------|------|
| Merge smoke summary | `validation_evidence/merge_smoke/SUMMARY.txt` |
| Smoke engine SHA file | `validation_evidence/merge_smoke/FlexAIDdS.sha256` |
| CMA-ES orchestrator notes | `validation_evidence/orchestrator/SUMMARY.md` |
| A–G validation pack | `validation_evidence/build_ab/` |
| Pin script | `.grok/skills/flexaidds/scripts/resolve_build.py` |

### Known non-blockers

- `test_cmaes_search` under g++-16 + Homebrew GTest: **ABI link FAIL** (libstdc++ vs libc++).
- Clang standalone CMA-ES mock tests: **PASS 5/5** (see merge_smoke `mock_tests.log`).

---

## One-line pin for chat paste

```text
PIN: build/FlexAIDdS @ 8b85ba77a936d24e6b896598fc906a677d82cc018e90ae9a5559d7d0ffee535f
     HEAD 7c901a687 (PR #296 merged) | gates: PB_CLASH_GRID_HOIST + WAL_COERCIVE
     Relaunch: Claude Science | Grok: pin+handoff only
```
