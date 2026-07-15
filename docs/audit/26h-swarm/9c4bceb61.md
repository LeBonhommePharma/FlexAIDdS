# Deep code audit — `9c4bceb61` (PR #274)

| Field | Value |
|-------|--------|
| **Short** | `9c4bceb61` |
| **Full SHA** | `9c4bceb6132606bf7c59a7f2e4a24b2a966ed277` |
| **Subject** | Merge pull request #274 science-admission-metrics-contract |
| **PR** | [#274](https://github.com/LeBonhommePharma/FlexAIDdS/pull/274) — *Add: Enforceable S1/S2/S3 admission + metrics claim contract* |
| **Branch** | `fix/science-admission-metrics-contract` → `main` |
| **Merged** | 2026-07-15T04:58:08Z |
| **Parents** | `439a2be8c` (main / PR #273) + `9dbbd9fa9` (feature tip) |
| **Merge type** | GitHub merge commit (no squash) |
| **Net diff (vs first parent)** | 5 files, **+1055 / −2** |
| **Audit date** | 2026-07-15 |
| **Audit scope** | S1/S2/S3 admission + metrics claim contract only. **No source edits** in this audit pass. |
| **Verdict** | **ACCEPT with residual doc/test gaps** (P2/P3). Contract is scientifically aligned, fail-closed where it matters, CI-gated, and does not touch ranking/election. |

---

## 1. Executive summary

PR #274 lands an **enforceable claim-table contract** for three-engine / TIER-1 science:

| Layer | Artifact |
|-------|----------|
| Normative prose | `benchmarks/protocols/admission_metrics_contract.md` |
| Protocol cross-link | `benchmarks/protocols/three_engine_entropy_comparison.md` §1.4–§5 |
| Enforcement | `scripts/aggregate_claim_metrics.py` (556 lines) |
| Unit tests | `tests/test_aggregate_claim_metrics.py` (12 tests; all pass) |
| CI wire | skill job in `.github/workflows/ci.yml` adds the new test file |

**Science intent (correct):**

- **S1** = elected (top-1) pose RMSD ≤ 2.0 Å → **primary / headline KPI**
- **S2** = S1 ∧ PoseBusters pass → secondary
- **S3** = BCR / any-pose RMSD ≤ 2.0 Å → **diagnostic only** (never abstract success)
- Claim admission requires `seed_echo==0`, `native_pose_seeded==0`, and `matrix_md5` pin match
- Election gap = S3 ∧ ¬S1 reported separately (sampling vs election)

**Critical non-goals (honoured):** no changes to pose ranking, clustering, GA, CF scoring, or docking C++. Oracle / native-seeded ceilings remain on a separate track (`scripts/aggregate_oracle_ceiling.py`).

The feature branch is two commits:

1. `668cc3095` — initial contract + aggregator + 8 tests  
2. `9dbbd9fa9` — **hardening follow-up**: fail-closed seeds, `rmsd_top1` path, finite-RMSD overrides flags, multi-source CLI reject, 32-hex pin validation, CI wiring, +4 tests → **12 total**

That follow-up is the difference between a soft aggregator and a claim gate. Merge quality is clean (disjoint paths from base; no conflict markers; no accidental reverts of #273 content into the merge tree vs first parent).

---

## 2. Commit / merge topology

```
439a2be8c  Merge PR #273 (astex-apo-strip-validation)     ← first parent (main)
    \
     9c4bceb61  Merge PR #274
    /
9dbbd9fa9  Fix: fail-closed seeds + rmsd_top1             ← second parent
    |
668cc3095  Add: S1/S2/S3 admission + metrics claim contract
```

| Check | Result |
|-------|--------|
| Second parent is ancestor of merge | Yes (`9dbbd9fa9` ⊂ `9c4bceb61`) |
| Diff first-parent…merge equals branch tip contents | Yes (same 5-file, +1055/−2 stat) |
| C++ / `LIB/` / DatasetRunner election touched | **No** |
| Secrets / hardcoded user paths introduced | **No** (`check_repo_hygiene.py` clean on current tree) |
| Executable bit on aggregator | `100755` (consistent with other scripts) |

PR body “test proof” quoted **8 passed**; tip of branch (and merge tree) has **12**. Residual documentation lag only — not a functional defect at merge.

---

## 3. File-by-file review

### 3.1 `benchmarks/protocols/admission_metrics_contract.md` (+130)

**Role:** Normative contract for claim aggregation and abstract/headline rates.

**Strengths**

- Explicit metric table (S1 primary, S2 secondary, S3 diagnostic).
- Hard rules: never headline S3; always report S1/S2/S3 separately; election gap as count, not folded into S1.
- Claim admission triple gate: seeds + matrix pin.
- Default matrix pin `72d7c7396702331d96ff12d18f831796` (`MC_st0r5.2_6.dat`) with override order: CLI → `RUN_RECEIPT.json` → `provenance.json` → default.
- Explicit separation of seeded / oracle-ceiling campaigns.
- “What this contract does **not** change” protects ranking/GA/CF.
- Fail-closed seed paragraph (post-`9dbbd9fa9`).

**Gaps (doc ↔ code drift)**

| Doc says | Code does | Severity |
|----------|-----------|----------|
| Elected RMSD: `rmsd_hungarian` → `rmsd_to_crystal` | `rmsd_hungarian` → **`rmsd_top1`** → `rmsd_to_crystal` | **P2** — three-engine schema needs the middle key in the normative field map |
| S1 flag order: `success_s1` first, then recompute | Finite elected RMSD **always wins**; flags only if no finite RMSD | **P2** — doc still implies flag-first |
| S3 field: `best_cluster_rmsd` | Also accepts `rmsd_bcr` | **P3** — missing alternate name |
| Drop reason language | Missing seed reported as `seed_echo!=0` | **P3** — true failure mode but message is imprecise |

None of these reverse the science contract; they can mislead a human re-implementer or an agent that reads only the markdown.

### 3.2 `benchmarks/protocols/three_engine_entropy_comparison.md` (+15/−1)

- Related-docs line points at the new contract + aggregator.
- §5 admission block gains enforceable CLI examples and default pin.
- Aligns with existing §1.4 metric definitions and McNemar / election-gap stats.
- Three-engine CSV schema already listed `rmsd_top1`, `rmsd_bcr`, `success_s1/2/3` — aggregator follow-up correctly tracked that schema; the **standalone contract field map did not fully catch up**.

### 3.3 `scripts/aggregate_claim_metrics.py` (+556)

#### Architecture

| Component | Behaviour |
|-----------|-----------|
| `_f` | First finite float among keys; skips empty/`NA` |
| `_truth` | Explicit true spellings only (`1`, `True`, …) |
| `_flag0` | **Fail-closed** explicit zero (`0`, `0.0`, `false`, …); missing/blank → fail |
| `load_matrix_pin` | CLI / receipt / provenance / default; **32-hex normalize** |
| `load_campaign_rows` | Prefer `*/result.csv` (first row per target), else flat summary names |
| `elected_rmsd` | Hungarian → top1 → crystal |
| `is_s1` / `is_s2` / `is_s3` | Finite RMSD/BCR override flags; S2 gated on S1 |
| `is_claim_eligible` | Seeds + matrix + optional `protocol_claim_eligible` |
| `aggregate_rows` | N_raw / N_claim / dropped reasons / S1–S3 ids / election_gap / headline default S1 |
| `apply_headline` | `--headline s3` without `--diagnostic-only` → exit **2** + `CONTRACT VIOLATION` |
| `main` | Exactly one of dir / `--csv` / `--c0-full85`; exit 0 if N_claim>0 else 1 |

#### Admission logic (core)

```text
claim ⇔ _flag0(seed_echo)
      ∧ _flag0(native_pose_seeded)
      ∧ row_matrix_ok (empty matrix_md5 → campaign pin; else exact pin)
      ∧ (protocol_claim_eligible missing OR truthy)
```

**Correct scientific choices**

1. **Fail-closed seeds** after `9dbbd9fa9` — older “missing = 0” behaviour would have admitted legacy CSVs that never recorded seed provenance. That was unsafe for claim tables.
2. **Finite RMSD overrides `success_s1` / `success_rmsd`** — prevents stale engine flags from minting S1 on high-RMSD poses (tested).
3. **S3 finite BCR overrides `success_s3`** — same integrity rule on the diagnostic axis.
4. **S1 requires `0.0 ≤ RMSD ≤ 2.0`** — rejects negative sentinels and 999-class placeholders.
5. **`--headline s3` hard-fail** — process-level enforcement, not mere prose.
6. **Matrix pin hex validation** — rejects garbage pins before aggregation.
7. **Multi-source CLI reject** — avoids ambiguous campaign vs CSV vs C0 path mixes.
8. **S2 never exceeds S1** under the aggregator’s S1 definition (even if `success_pb` is present).
9. **Oracle path not mixed** — docstring and contract point to `aggregate_oracle_ceiling.py`.

#### Edge cases verified at audit time (live interpreter)

| Case | Result |
|------|--------|
| RMSD 999 / −1 sentinels | S1=false, S3=false |
| Blank `seed_echo` | Not claim-eligible (`seed_echo!=0` reason) |
| `protocol_claim_eligible=0` | Dropped even if seeds clean |
| Finite RMSD 1.5 with `success_s1=0` | S1=true (finite wins) |
| Only `rmsd_top1=1.1` | S1=true |
| Both hungarian=1.0 and top1=5.0 | Elected = 1.0 (Hungarian first) |
| `seed_echo="0.0"` | Accepted by `_flag0` |

#### Residual risks in the implementation

| ID | Finding | Severity | Notes |
|----|---------|----------|-------|
| R1 | **Row `matrix_md5` empty always passes** (campaign pin only). | P3 / intentional | Documented. Risk: heterogeneous matrices inside one campaign without per-row MD5 go unnoticed. Prefer engine always emit pin. |
| R2 | **`protocol_claim_eligible` missing is pass** (unlike seeds). | P3 / intentional | Contract says missing → recompute from seed flags. Asymmetric fail-closed vs seeds is OK but worth one line in contract “why”. |
| R3 | **`load_campaign_rows` keeps only first CSV row** per `*/result.csv`. | P3 | Fine for DatasetRunner single-row summaries; multi-row per target would silently drop restarts. |
| R4 | **No `pdb_id` de-duplication** across rows. | P3 | Duplicate dirs or flat CSVs with repeated targets inflate N. |
| R5 | **Drop reason `seed_echo!=0` for missing/blank**. | P3 | Operational noise only. Prefer `seed_echo_missing_or_nonzero`. |
| R6 | **`resolve_c0_full85_dir` default walks into `~/Library/Mobile Documents/...`**. | P2 ops | Only `is_dir` / path join (no `rglob`), but CloudDocs existence checks can still FileProvider-stall. Prefer requiring `FLEXAIDDS_RESULTS` / fail fast without home iCloud probe when env unset. Aligns with local-first / thin-iCloud policy. |
| R7 | **`--headline s2` allowed without extra flag**. | OK | Secondary KPI as headline is intentional; only S3 is restricted. |
| R8 | **S1 still consults `seed_echo` via `_truth` after admission**. | OK | Defence in depth; admission already fail-closed. |
| R9 | **No test for multi-source CLI / invalid matrix pin / `protocol_claim_eligible=0`**. | P3 | Behaviour present; untested corners. |
| R10 | **Tests appended after `if __name__ == "__main__"`**. | P3 style | Pytest still collects 12 functions; reorder for readability. |
| R11 | **Does not replace older claim aggregators** (`v48_selector_official`, `v26_v27_analysis`, ops monitor). | P2 process | Agents can still quote non-contract scripts for “success %”. Contract is normative; **ops/docs should point claim work only at this CLI**. |

### 3.4 `tests/test_aggregate_claim_metrics.py` (+354)

| Test | Intent | Status |
|------|--------|--------|
| `test_claim_filter_drops_seeded_rows` | seed_echo / native_seeded drops | Pass |
| `test_matrix_md5_pin_filters_wrong_matrix` | pin mismatch | Pass |
| `test_s1_vs_s3_diverge_election_gap` | election gap + S2 subset | Pass |
| `test_headline_s3_rejected_without_diagnostic_flag` | apply_headline API | Pass |
| `test_cli_headline_s3_exits_nonzero` | process exit 2 + message | Pass |
| `test_cli_happy_path_json` | JSON write + quiet | Pass |
| `test_flat_summary_csv` | flat summary path | Pass |
| `test_cli_csv_flag` | `--csv` | Pass |
| `test_missing_seed_columns_fail_closed` | fail-closed seeds | Pass |
| `test_rmsd_top1_three_engine_schema` | three-engine column names | Pass |
| `test_success_s1_flag_cannot_override_high_rmsd` | flag vs finite RMSD | Pass |
| `test_seed_echo_0_0_accepted` | `0.0` spelling | Pass |

**Audit verification (this session):**

```text
$ python3 -m pytest tests/test_aggregate_claim_metrics.py -q --tb=line
............                                                             [100%]
12 passed in 0.25s
```

Coverage quality is strong on the science-critical paths (admission, S1≠S3, S3 headline ban, three-engine schema). Gaps listed under R9.

### 3.5 `.github/workflows/ci.yml` (+1/−1)

Skill unit-test step now includes `tests/test_aggregate_claim_metrics.py` alongside skill/resolve_build tests. Placement is appropriate (pure Python, no C++ build). No other jobs modified.

---

## 4. Scientific guardrail alignment (`AGENTS.md`)

| Guardrail | Assessment |
|-----------|------------|
| Verify with execution | Unit suite green at audit; CLI contract paths exercised in tests |
| Separate CF proxy from thermodynamics | Aggregator speaks only RMSD / PB / BCR / seeds — **no ΔG / free-energy claims** |
| Do not alter pose ranking | **No engine ranking/election changes** |
| No silent seed success | Fail-closed `seed_echo` + `native_pose_seeded` |
| Matrix identity | Pin required at campaign level; per-row mismatch drops |
| S3 not abstract success | Hard exit 2 + report labels `diagnostic_only` + warning strings |
| Apache-2.0 / no GPL | Pure Python stdlib + pytest; no new deps |
| Local-first / thin-iCloud | Soft concern: C0 default path probes iCloud home (R6); claim computation itself is local CSV read |

---

## 5. Security / hygiene / ops

- No network calls, no shell interpolation of CSV fields, no `eval`.
- CSV read via `csv.DictReader` only.
- CLI paths are `Path` existence checks; no `subprocess` of untrusted binaries.
- Hygiene script: no tracked secrets / no hardcoded agent user paths in this merge.
- Exit code contract documented and tested for the S3 violation path.

---

## 6. Merge-risk assessment

| Risk class | Level | Rationale |
|------------|-------|-----------|
| Functional regression in docking | **None** | No `LIB/` / election changes |
| Claim inflation (false high S1) | **Low** | Finite RMSD wins; seeds fail-closed; sentinels excluded |
| Claim deflation (false low N_claim) | **Low–med** | Fail-closed seeds will drop incomplete historical CSVs — **correct for claims**, may surprise ops on pre-seed-column dumps |
| Headline S3 abuse | **Low** | Process-level ban |
| Matrix confusion | **Low** | Pin validation + drop on mismatch; empty per-row still soft (R1) |
| Agent misuse of alternate scripts | **Med process** | Older aggregators still exist (R11) |
| Doc drift | **Med doc** | Field map missing `rmsd_top1` / flag precedence (P2) |

**Overall merge risk: LOW** for code correctness; **MEDIUM residual process/doc** for agent/operator compliance outside this CLI.

---

## 7. Findings catalogue

### P0 — none

### P1 — none

### P2 — follow-ups recommended (non-blocking for this merge)

1. **Sync normative field map** in `admission_metrics_contract.md` with code:
   - Elected RMSD: `rmsd_hungarian` → `rmsd_top1` → `rmsd_to_crystal`
   - S3: `best_cluster_rmsd` → `rmsd_bcr`
   - S1: finite RMSD first; flags only if no finite elected RMSD
2. **C0 path resolution:** when neither `FLEXAIDDS_RESULTS` nor `FLEXAIDDS_ICLOUD` is set, fail with exit 2 rather than probing `Mobile Documents` (R6).
3. **Ops pointer:** make claim/abstract pipelines and handoff docs default to `aggregate_claim_metrics.py` only (R11).

### P3 — nice-to-have

1. Clearer drop reasons for missing seed columns.
2. Tests: multi-source CLI, invalid pin, `protocol_claim_eligible=0`, blank seed reason text.
3. Reorder tests above `if __name__ == "__main__"`.
4. Optional: warn if N_claim rows lack per-row `matrix_md5` when campaign has heterogeneous sources.
5. Optional: de-dupe by `pdb_id` with explicit conflict report.

---

## 8. What this merge gets right (keep)

1. **Normative + enforceable** pair (markdown contract + failing CLI).
2. **S1 primary / S3 diagnostic** enforced at exit codes, report roles, and human text (`[NOT abstract success]`).
3. **Election gap** as first-class output — directly supports three-engine H2/H3 analysis without polluting S1.
4. **Fail-closed seed provenance** after the follow-up commit — the correct default for claim science.
5. **Three-engine schema awareness** (`rmsd_top1`, `rmsd_bcr`) so FlexAID A/B0/B arms and DatasetRunner CSVs share one aggregator.
6. **CI attachment** so regressions on admission cannot silently land.
7. **Zero ranking surface area** — pure post-hoc metrics discipline.

---

## 9. Reproduction / verification commands

```bash
# Unit suite (merge tip / main after PR)
python3 -m pytest tests/test_aggregate_claim_metrics.py -q --tb=line

# Contract violation (expect exit 2)
python3 scripts/aggregate_claim_metrics.py <campaign_dir> --headline s3

# Happy path
python3 scripts/aggregate_claim_metrics.py <campaign_dir> --json /tmp/claim.json

# C0 full85 (requires env)
source ~/.flexaidds_env
python3 scripts/aggregate_claim_metrics.py --c0-full85

# Merge contents
git show 9c4bceb6132606bf7c59a7f2e4a24b2a966ed277 --stat
git log --oneline 439a2be8c..9dbbd9fa9   # 668cc3095, 9dbbd9fa9
```

---

## 10. Verdict

**ACCEPT** merge `9c4bceb61` / PR #274 as a **correct, low-risk science-governance merge**.

It successfully operationalizes the S1/S2/S3 admission contract without contaminating the docking engine or ranking path. The follow-up commit (`9dbbd9fa9`) is essential: fail-closed seeds, finite-RMSD flag override, and three-engine `rmsd_top1` make the aggregator trustworthy for claim tables.

Ship residual work as **docs/ops P2** (field-map sync, C0 env-hard-fail, exclusive claim CLI pointers) and **test P3** — not as blockers to this merge.

| Dimension | Score |
|-----------|-------|
| Scientific correctness | Strong |
| Enforcement fidelity | Strong (after follow-up) |
| Test coverage of critical paths | Strong (12/12 green) |
| Doc completeness | Good with P2 drift |
| Engine safety (no ranking touch) | Excellent |
| Merge cleanliness | Excellent |

---

*Audit method: full read of both feature commits and merge tree; line-level review of aggregator + tests + contract; live pytest + interpreter edge probes; cross-check against `three_engine_entropy_comparison.md` §1.4–§5, `AGENTS.md` scientific guardrails, and sibling `aggregate_oracle_ceiling.py`. No production source modified for this report.*
