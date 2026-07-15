# Audit: 9dbbd9fa9 — Fix: Admission metrics fail-closed seeds + three-engine rmsd_top1

## Summary (2–4 sentences)

Commit `9dbbd9fa9edc36543e6b3966f259870559df5b17` correctly hardens **claim admission** and **S1/S2/S3 evaluation** in `scripts/aggregate_claim_metrics.py`: missing/blank `seed_echo` / `native_pose_seeded` now fail (fail-closed), elected RMSD prefers `rmsd_hungarian` → `rmsd_top1` → `rmsd_to_crystal`, and **finite numeric RMSD always overrides** stale `success_s1` / `success_s3` flags. Multi-source CLI is rejected and matrix pins must be 32 hex chars; **12 unit tests pass** and are wired into the skill CI job. Residual risks are **medium doc drift** (contract field map still omits `rmsd_top1` and still documents flag-first S1) and **negative RMSD sentinels** (`_f` treats `-1` as real, unlike later DPFO aggregator). **Ranking risk is NONE** (post-hoc claim aggregator only; no engine election change).

## Severity: **LOW–MEDIUM** (claim-table integrity fix; residual sentinel/doc gaps)

## Verdict (claim science first)

| Question | Verdict |
|----------|---------|
| Missing/blank seed columns fail claim admission? | **YES** — fail-closed (parent was fail-open) |
| Explicit `seed_echo=0` / `0.0` still admit? | **YES** |
| Elected RMSD: hungarian → `rmsd_top1` → crystal? | **YES** |
| Finite RMSD overrides `success_s1`? | **YES** (unit-tested) |
| Finite BCR overrides `success_s3`? | **YES** in code; **no dedicated unit test** |
| S2 ⊆ S1; S3 diagnostic-only; headline s3 gated? | **YES** (unchanged, still correct) |
| Multi-source CLI rejected? | **YES** (probed exit 2) |
| Matrix pin 32-hex validated? | **YES** (probed) |
| Contract doc fully synced with code? | **NO** — field map + S1 order lag |
| Negative RMSD sentinels (`-1`) treated as missing? | **NO** — `_f` accepts finite negatives (M2) |
| Engine ranking / GA / CF changed? | **NO** |
| Unit tests green? | **YES** — 12/12 |

**Overall:** **ACCEPT / ship.** Direction is correct for claim science. Follow up with contract sync + sentinel handling aligned to DPFO (`v >= 0`).

---

## Scope audited

| Path | Role in commit |
|------|----------------|
| `scripts/aggregate_claim_metrics.py` | `_flag0` fail-closed; `elected_rmsd` + `rmsd_top1`; `is_s1`/`is_s3` finite-first; pin normalize; multi-source CLI |
| `tests/test_aggregate_claim_metrics.py` | +4 tests (12 total) |
| `benchmarks/protocols/admission_metrics_contract.md` | Fail-closed seed note only (partial) |
| `.github/workflows/ci.yml` | Wire `tests/test_aggregate_claim_metrics.py` into skill job |
| Cross-read (not edited) | Parent `668cc3095`; three-engine schema §5; `parse_flexaid_arm_results.py`; DatasetRunner `-1` RMSD sentinels; later `c93b866a4` DPFO `_f` |

**Parent:** `668cc3095734d6f074c50814fd117c3f70c8eb0c` — *Add: Enforceable S1/S2/S3 admission + metrics claim contract*

**Author / date:** LP \<lp@thebonhomme.com\> · 2026-07-15 00:57:35 -0400

**Diffstat:** 4 files, +141 / −23

---

## Verification (this session)

```text
python3 -m pytest tests/test_aggregate_claim_metrics.py -v
# 12 passed in 0.19s
```

Behavioral probes (import module; not committed):

- `_flag0` missing/blank → False; `"0"`/`"0.0"`/`"false"` → True; `"0.00"`/`"FALSE"` → False
- `elected_rmsd`: priority hungarian > top1 > crystal; `h=NA,top1=1.2` → 1.2; **`h=-1,top1=1.2` → -1.0** (shadow)
- `is_s1(rmsd=5, success_s1=1)` → False; `is_s1(rmsd=1, success_s1=0)` → True; boundary 2.0 inclusive
- `is_s3(BCR=5, success_s3=1)` → False; `is_s3(rmsd_bcr=1.5)` → True
- CLI multi-source → exit 2; bad pin `zzzz` → exit 2; good 32-hex pin → exit 0

---

## Findings

### F1. Fail-closed seed admission — PASS (primary fix)

**Parent (`668cc3095`):**
```python
# missing key or blank → True  (claim-pass; treat as 0)
```

**This commit:**
```python
def _flag0(row, key) -> bool:
    """Fail-closed: missing or blank keys fail admission (return False)."""
    if key not in row or row.get(key) is None:
        return False
    s = str(row.get(key, "")).strip()
    if s == "":
        return False
    return s in ("0", "0.0", "False", "false", "NO", "no")
```

`is_claim_eligible` still gates on `_flag0(seed_echo)` and `_flag0(native_pose_seeded)` before matrix pin and optional `protocol_claim_eligible`.

**Why this matters:** Claim tables for three-engine / C0 must never treat “no seed columns” as “no seed”. Parent’s fail-open path could inflate N_claim on incomplete CSVs.

**Intentional break (I1):** Legacy flat CSVs without seed columns → N_claim=0, process exit 1. Correct for claim science; ops must re-export rows with explicit zeros.

**Test:** `test_missing_seed_columns_fail_closed` — PASS.

**Minor (L1):** `"0.00"` and `"FALSE"` rejected. Engine writers use `0`/`1` integers as strings — production risk low. Drop reason always `seed_echo!=0` even when **missing** (L2 — ops noise).

---

### F2. `rmsd_top1` elected RMSD mapping — PASS with sentinel caveat

```python
def elected_rmsd(row):
    """Preferred elected RMSD: Hungarian → three-engine rmsd_top1 → crystal."""
    return _f(row, "rmsd_hungarian", "rmsd_top1", "rmsd_to_crystal")
```

| Producer | Elected field | BCR | Covered? |
|----------|---------------|-----|----------|
| DatasetRunner | `rmsd_hungarian` / `rmsd_to_crystal` | `best_cluster_rmsd` | Yes |
| Three-engine / `parse_flexaid_arm_results.py` | `rmsd_top1` (empty if unknown) | `rmsd_bcr` | Yes |
| Protocol `three_engine_entropy_comparison.md` §5 | lists `rmsd_top1`, `rmsd_bcr` | — | Code yes; contract table **no** |

**Test:** `test_rmsd_top1_three_engine_schema` — PASS (S1=1, S3=1 on top1=1.2, bcr=0.8).

**Priority is scientifically correct:** symmetry-corrected Hungarian > three-engine top-1 > serial crystal.

**Report string lag (L3):** JSON still defines S1 as `"elected RMSD <= 2.0 A (Hungarian preferred)"` without mentioning `rmsd_top1`.

---

### F3. S1 finite-RMSD override — PASS (integrity fix)

**Parent bug:** `success_s1` honored **before** numeric recompute → a stale `success_s1=1` with `rmsd_hungarian=5` would count as S1 success.

**This commit:**
```text
if seed_echo truthy → False
rh = elected_rmsd
if finite(rh) → (0 ≤ rh ≤ 2.0)   # ALWAYS wins
else success_s1 → success_rmsd → success → False
```

| Case | S1 | Correct? |
|------|----|----------|
| rmsd=5, success_s1=1 | False | Yes |
| rmsd=1, success_s1=0 | True | Yes |
| no RMSD, success_s1=1 | True | Yes (fallback) |
| rmsd=2.0 | True | Yes (≤) |
| rmsd=2.0001 | False | Yes |
| seed_echo=1, rmsd=0.5 | False | Yes |

**Test:** `test_success_s1_flag_cannot_override_high_rmsd` — PASS.

---

### F4. S2 correctness — PASS (unchanged this commit)

```text
require S1
success_pb if present else pb_pass else False
```

- S2 never exceeds S1.
- Missing PB columns → S2 false (fail-closed secondary). Three-engine rows often leave `pb_pass` empty → S2=0 until PB wired. Expected.
- Parent election-gap test still asserts S2 ⊆ S1.

---

### F5. S3 finite-BCR override + diagnostic role — PASS with sentinel/test gap

```text
bc = _f(best_cluster_rmsd, rmsd_bcr)
if finite(bc) → (0 ≤ bc ≤ 2.0)
elif success_s3 → truthy
else False
```

| Case | S3 | Correct? |
|------|----|----------|
| BCR=5, success_s3=1 | False | Yes (finite wins) |
| BCR=1, success_s3=0 | True | Yes |
| only success_s3=1 | True | Yes |
| rmsd_bcr=1.5 | True | Yes |
| empty | False | Yes |
| BCR=-1 sentinel, success_s3=1 | **False** | Flag unreachable (M2); parent checked flag **first** |

Headline guard unchanged: `--headline s3` without `--diagnostic-only` → exit 2. Election gap = S3 ∧ ¬S1 on claim rows only.

**No unit test** for S3 finite override (symmetric to S1 test) — L4.

---

### F6. Negative RMSD sentinels — MEDIUM residual (M2)

`_f` returns first **finite** float; negatives are finite:

```python
# DatasetRunner: uncomputed RMSD / BCR often -1.0f
# elected_rmsd({rmsd_hungarian: "-1", rmsd_top1: "1.2"}) → -1.0
# is_s1 → False (range), never falls through to top1 or flags
```

Three-engine writer emits **empty** (not `-1`) for missing RMSD → pure three-engine rows OK.

DatasetRunner-only with only `-1` → S1 False (fail-closed under-count; safer than over-claim).

**Sibling policy (later `c93b866a4` `aggregate_dpfo_pilot.py`):**
```python
if math.isfinite(v) and v >= 0.0:
    return v
```

**Fix recommendation:** Align claim aggregator `_f` with `v >= 0.0`; then re-test mixed schemas and flag fallback when only sentinels present.

---

### F7. Matrix pin + multi-source CLI — PASS

```python
def _normalize_matrix_pin(md):
    s = str(md).strip().lower()
    if len(s) != 32 or any(c not in "0123456789abcdef" for c in s):
        raise ValueError(...)
    return s
```

- Applied to CLI, RUN_RECEIPT, provenance, default constant.
- `n_sources > 1` → exit 2; zero sources → help + exit 2.
- Probed: multi-source and bad pin both exit 2.

**No unit tests** for these paths (L4) — logic is simple and probed manually.

---

### F8. Contract documentation drift — MEDIUM (M1)

| Contract statement | Code after commit | Match? |
|--------------------|-------------------|--------|
| Fail-closed missing/blank seeds | `_flag0` | **Yes** |
| Elected RMSD: hungarian → crystal | hungarian → **top1** → crystal | **No** |
| S1: success_s1 first, else recompute | recompute first if finite | **No** |
| S3 never primary | `apply_headline` | **Yes** |
| Default pin `72d7c739…` | `DEFAULT_MATRIX_MD5` | **Yes** |

Only a short “Fail-closed seed flags” appendix was added. Normative §1 field mapping is stale relative to the enforcement script.

---

### F9. Tests + CI — PASS with structure smell

| Behavior | Tested? |
|----------|---------|
| Seeded rows dropped | Yes (parent) |
| Matrix mismatch | Yes (parent) |
| S1/S3 election gap + rates | Yes (parent) |
| Headline S3 contract | Yes (parent) |
| CLI happy path / `--csv` | Yes (parent) |
| Missing seed columns fail-closed | **Yes (new)** |
| `rmsd_top1` three-engine schema | **Yes (new)** |
| Finite RMSD overrides success_s1 | **Yes (new)** |
| `"0.0"` seed spelling | **Yes (new)** |
| Multi-source CLI / bad pin | **No** |
| S3 finite override | **No** |
| Negative sentinel shadowing | **No** |

New tests sit **after** `if __name__ == "__main__"`. Pytest file-path collection still finds all 12 (verified). Style smell only.

CI: skill job pytest list includes `tests/test_aggregate_claim_metrics.py`. Appropriate (pure Python, no C++ deps).

---

### F10. Scientific / ranking impact — NONE on engine; POSITIVE on claims

| Guardrail | Status |
|-----------|--------|
| No-seed claim admission | **Strengthened** |
| S1 headline KPI | Preserved |
| S3 not abstract success | Preserved |
| Three-engine CSV aggregation | **Enabled** |
| Pose ranking / GA / CF / clustering | **Untouched** |
| Matrix pin integrity | **Strengthened** |

Pure post-hoc aggregator. Safe under AGENTS.md “preserve current ranking” (no ranking change).

---

## Findings table (actionable)

| ID | Sev | Finding | Fix |
|----|-----|---------|-----|
| M1 | Medium | Contract §1 field map omits `rmsd_top1`; S1 order wrong | Sync `admission_metrics_contract.md` with code |
| M2 | Medium | `_f` accepts RMSD `< 0` (sentinel) → shadows keys / blocks flags | `v >= 0` like DPFO; unit tests |
| L1 | Low | `_flag0` misses `0.00` / `FALSE` | Normalize carefully or document accepted spellings |
| L2 | Low | Drop reason `seed_echo!=0` for missing columns | Distinct reasons (`missing` vs `nonzero`) |
| L3 | Low | JSON S1 definition omits top1 | Update definition string |
| L4 | Low | Tests after `__main__`; missing multi-source/pin/S3/sentinel tests | Reorder + expand |
| I1 | Info | Breaking change for seed-column-less legacy CSVs | CLI help / ops note |

---

## Diff intent vs delivery

| Commit message claim | Delivered? | Evidence |
|----------------------|------------|----------|
| Missing/blank seed_echo and native_pose_seeded fail claim | **Yes** | `_flag0`; fail-closed test |
| Elected RMSD hungarian → rmsd_top1 → crystal | **Yes** | `elected_rmsd`; three-engine test |
| Finite RMSD overrides success_s1 / success_s3 | **Yes** (S1 tested; S3 logic yes) | `is_s1` / `is_s3` |
| Reject multi-source CLI | **Yes** | `n_sources > 1`; probed |
| Validate 32-hex matrix pin | **Yes** | `_normalize_matrix_pin`; probed |
| Wire tests into skill CI; 12 unit tests | **Yes** | `ci.yml`; 12 passed |

---

## Final recommendation

**Keep / ship.** Parent fail-open seed admission was the more dangerous claim bug; this commit’s direction is correct.

**Follow-up (small, test-gated PR):**

1. Sync contract §1 (M1).
2. Align `_f` with `v >= 0` sentinel policy (M2).
3. Add tests: multi-source CLI, bad pin, S3 finite override, `h=-1` + `top1` fallback (L4).

No revert. No engine rebuild required for this commit alone.
)
