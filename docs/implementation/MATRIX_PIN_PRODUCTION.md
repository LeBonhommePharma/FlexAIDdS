# Production matrix pin (MC_st0r5.2_6.dat)

**Decision (LP, 2026-07-23): `9dc9` is production; `72d7` stands down.**

| MD5 | Role |
|-----|------|
| **`9dc93717dfed0698006d88dd6a9627bc`** | **Production** — baseline-validated JCIM matrix; what clean `main` ships and what claim launchers must assert |
| `72d7c7396702331d96ff12d18f831796` | **Not production** — packing-sweetened fork (7 cells differ). Branch-only on `fix/matrix-repin-72d7`; do **not** merge those matrix commits to `main` |

## What differs (exactly 7 / 820 cells)

1. **Packing sweetened** (more attractive burial): `C.2×C.ar` (2-4), `N.ar×O.2` (10-13)
2. **O-desolvation cut** (cheaper to bury O): `O.2/O.3/O.co2 × SOLVENT` (13/14/15-40)
3. **Zn coord zeroed**: `N.ar×Zn` (10-35), `N.pl3×Zn` (12-35)

No measured claim benefit (e.g. 1G9V ~neutral); direction pushes toward confirmed over-burial pathology.

## Enforced by (hard fail on mismatch)

| Path | Constant |
|------|----------|
| `scripts/generate_flexaid_inp.py` | `MATRIX_MD5_PIN` |
| `scripts/run_flexaid_arm_pilot8.sh` | `MATRIX_PIN` |
| `scripts/run_3dsig_red_pair_{serial,full85}.sh` | `MATRIX_PIN` |
| `scripts/aggregate_claim_metrics.py` | `DEFAULT_MATRIX_MD5` |
| Protocol docs under `benchmarks/protocols/` | fallback / default pin |

## Fleet / launch rules

1. **Do not** merge matrix commits `41a385ed` / `5ee49ba1` to `main`.
2. Cherry-pick only non-matrix code fixes (water, coarse_init, hard-clash, probe_cf, …) onto 9dc9-`main` when landing.
3. Any full-85 / claim run must use a data dir whose `MC_st0r5.2_6.dat` MD5 is **9dc9** (e.g. clean `main` worktree or isolated data dir with that blob). **`build/` on a 72d7 branch may still hold the wrong bytes** — preflight MD5, never assume.
4. Historical one-shot launchers (`scripts/launch_v40.py` / `v41` / `v42`) keep their frozen 72d7 expectations for re-verification of those campaign directories only; they are **not** production claim paths.
5. Historical audit notes under `docs/audit/` may still mention 72d7; treat them as period snapshots, not the live pin.

## Verify

```bash
md5 -q MC_st0r5.2_6.dat
# expect: 9dc93717dfed0698006d88dd6a9627bc

# claim scripts (must match pin above)
rg -n 'MATRIX_(MD5_)?PIN|DEFAULT_MATRIX_MD5' \
  scripts/generate_flexaid_inp.py \
  scripts/run_flexaid_arm_pilot8.sh \
  scripts/run_3dsig_red_pair_serial.sh \
  scripts/run_3dsig_red_pair_full85.sh \
  scripts/aggregate_claim_metrics.py
```
