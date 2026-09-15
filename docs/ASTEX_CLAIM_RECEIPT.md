# Astex Diverse-85 claim receipt (scaffold)

This file is a **receipt template**, not a measured success rate.

Do **not** fill the rate fields with invented numbers. A live claim is allowed
only from a completed on-disk `result.csv` plus `RUN_RECEIPT` produced by a
real engine run of `scripts/reproduce_astex85.sh` **without** `--dry-run`.

`--dry-run` writes a protocol receipt and exits. It does **not** dock 85
targets and must not be quoted as docking power.

## Denominator (normative)

- **N = 85** (Hartshorn 2007 Astex Diverse).
- **2HR7 is included.** A miss (RMSD > 2.0 Å or PoseBusters fail, or absent
  from the campaign) is a **failure**, not an exclusion.
- There is **no** `astex85_codes_84.txt` roster.
- Frozen codes: `benchmarks/protocols/astex85_target_manifest.json`.

## What 0.70 is (and is not)

`expected_baselines.docking_power_top1: 0.70` in the Astex YAML is
**`expected_baselines_role: ci_gate_only`**.

It is **not**:

- a published FlexAIDdS docking-power rate
- a receipted campaign result
- interchangeable with Gaudreault & Najmanovich 2015 JCIM Table 2 (0.452,
  different protocol) or the 3Dsig 2017 deck top-10 bootstrap medians

Do not quote “70%” as a FlexAIDdS result.

## Placeholders (UNSET until a real run lands)

| Field | Value |
|-------|-------|
| n_targets | 85 |
| git_commit | UNSET |
| binary_path | UNSET |
| binary_sha256 | UNSET |
| matrix_md5 | UNSET |
| claim_ready_count | UNSET — no invented rate |
| claim_ready_rate | UNSET — no invented rate |
| success_rmsd_count | UNSET — no invented rate |
| success_rmsd_rate | UNSET — no invented rate |
| success_pb_count | UNSET — no invented rate |
| success_pb_rate | UNSET — no invented rate |
| 2HR7_outcome | UNSET (counts in the 85) |

## How to populate for real

```bash
# Protocol-only (CI / this scaffold): no docking, no rates.
bash scripts/reproduce_astex85.sh --dry-run

# Live 85-target run (local-first OUT; see AGENTS.md benchmark storage).
# Only that path may replace UNSET above.
bash scripts/reproduce_astex85.sh
```

Success language still requires RMSD ≤ 2.0 Å **and** PoseBusters (STRICT also
requires tENCoM/Eigen). See `.grok/skills/flexaidds/SKILL.md` deception-proof
claim contract and `METHODOLOGY.md`.
