# Scoring provenance schema (Wave 0.3)

**Purpose:** Fail-closed audits of docking campaigns. Empty `git_commit` or missing
scoring env (as on the UNCITABLE com-cap CAP=-130 run) must not pass claim review.

**Checker:** `python3 scripts/check_run_receipt.py <campaign_dir>`

---

## Required documents (any one may be primary)

| File | Role |
|------|------|
| `RUN_RECEIPT.json` | Campaign-level claim receipt (preferred) |
| `provenance.json` | DatasetRunner per-run provenance |
| `campaign.log` / `provenance.txt` | Fallback text lines (weaker) |

---

## Required keys (JSON receipts)

| Key | Meaning |
|-----|---------|
| `matrix_md5` | Must be **`9dc93717dfed0698006d88dd6a9627bc`** for claim-style 9dc9 runs |
| `binary_sha256` or `binary_path` | Engine identity |
| `git_commit` | **Non-empty** commit of the **built** binary when available |
| `seed_elitism` | Must be 0 / false for claim-style fair rates |

## Strongly recommended scoring env (JSON object `scoring_env` or top-level)

Record actual process environment used at launch (not ambient shell after the fact):

| Env | Notes |
|-----|--------|
| `FLEXAIDDS_ACF_STRICT` | E1b — default unset/0 |
| `FLEXAIDDS_COM_BURIAL_CAP` | Prefer unset; never cite CAP=-130 without full probe set |
| `FLEXAIDDS_COM_FLOOR` | Soft com clamp if used |
| `FLEXAIDDS_VCT_NORM` | Intensive com if used |
| `FLEXAIDDS_SOFTBETA_ELECTION` | DatasetRunner Softβ S1 (default OFF) |
| `FLEXAIDDS_ELECTION_ENTROPY` | Legacy alias if used |
| `FLEXAIDDS_WAL_STIFF` / wall knobs | When E2 wall experiments run |
| `OMP_NUM_THREADS` | Per worker |
| workers / `--threads` | Campaign concurrency |

## Optional but useful

`pop`, `gen`, `restarts`, `mode` (autonomous vs defined-cleft), `temperature_K` / TEMPER,
`matrix_path`, `runner_sha256`.

---

## Checker exit codes

| Code | Meaning |
|------|---------|
| 0 | Required keys present and non-empty (claim-eligible shape) |
| 1 | Missing required keys or empty `git_commit` / matrix pin mismatch when `--require-matrix-9dc9` |
| 2 | Usage / path not found |

Ops launchers should write `scoring_env` into RUN_RECEIPT or a sibling `SCORING_PROVENANCE.json`
before any rate is cited.
