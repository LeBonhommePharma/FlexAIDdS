# P0 — Claim Contract Repair (fixed denominator + PB 0.6.5 fixtures)

**Commit:** 410e63cb7 (origin/main) · **Date:** 2026-07-19 · **Role:** OPS/CI (Algorithmic Optimizer)
**Scope:** instrument fix only — NO change to the docking engine or scoring function.

## Problem (Codex P0)
The claim aggregator computed success rate as `k / len(claim_eligible_rows)`. Any target
dropped in admission (seed contamination, missing receipt, error) or absent from the CSV
left the denominator entirely — so **dropping hard targets mechanically inflated the rate**.
Separately, PoseBusters was **not installed in any environment**, so the PB backend could only
schema-reject (the "0/85 from schema rejection" failure mode).

## Fix
1. **Frozen denominator.** `benchmarks/protocols/astex85_target_manifest.json` pins the 85
   Astex codes (cross-verified: dataset YAML ∩ strict-audit CSV, identical sets; sha256
   `da89650afd79…`). Rates are now `k / 85`. Missing/dropped targets count as **failures**
   and never leave the denominator. `max(manifest_N, observed)` so extra rows can't shrink it.
2. **Real PoseBusters 0.6.5** installed in `benchmarker` env; `bust` CLI present (backend
   `bust_cli`, required for claims).
3. **Transparency fields** in every report: `N_denominator`, `N_denominator_source`,
   `N_missing_from_manifest`, `missing_targets`.

## Kill/promote gate (both pass)
- `tests/p0_claim_contract/test_fixed_denominator.py` — **10/10**: denominator=85 not
  len(eligible); 10 strict successes over 85 report 11.76% not 100%; hash-mismatch row fails
  closed; missing seed columns fail closed; off-manifest extra can't shrink denom.
- `tests/p0_claim_contract/test_posebusters_fixtures.py` — real PB 0.6.5 discriminates:
  clean pose **12/12**, broken pose **8/12** failing exactly `bond_lengths`, `bond_angles`,
  `internal_steric_clash` (the same checks dominating the audit failure tallies). Proves the
  harness is not schema-rejecting everything and not passing everything.

## Effect (same file, before vs after)
| File | Old summary | Under fixed contract |
|---|---|---|
| `oracle_ceiling_restore_v43proto_r3/astex_diverse_results.csv` | 83/85 = **97.65%** | **0/85 = 0.00%** STRICT |

All 85 rows are `native_pose_seeded=1` / `protocol_claim_eligible=0` → correctly dropped, with
the denominator held at 85 and per-target drop reasons emitted. The contract now refuses to let
a 90%-native-seeded campaign report as a claim.

## What P0 does NOT do
It fixes the *ruler*, not the *science*. It does not change any docking result, and it does not
tell you which layer (objective / search / election) is broken. Next: native-CF oracle on current
defaults (minutes) to discriminate objective-broken vs search-broken before spending CPU-hours.
