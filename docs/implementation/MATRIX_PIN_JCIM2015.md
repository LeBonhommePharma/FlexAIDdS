# Interaction matrix pin — FlexAID JCIM 2015 lineage

**Rule (operator + agents):** The pairwise contact matrix used for claim / three-engine / 3Dsig-style reproduction **must not change** relative to the FlexAID JCIM 2015 production pin. Never “improve” matrix numbers mid-campaign.

## Canonical pin (admission + v137)

| Field | Value |
|--------|--------|
| Filename | `MC_st0r5.2_6.dat` |
| **MD5** | **`72d7c7396702331d96ff12d18f831796`** |
| SHA256 | `2265b1ba08e763887ab49b0bfbcfacf0505fd1c4b612fcbd3ec1fc851accd71c` |
| Size | **18043** bytes |
| Role | Three-engine / claim admission pin (`scripts/aggregate_claim_metrics.py` `DEFAULT_MATRIX_MD5`) |

Preflight (fail if mismatch):

```bash
test "$(md5 -q "$MATRIX_PATH")" = "72d7c7396702331d96ff12d18f831796"
```

## Observed drift (2026-07-15 audit)

A **second** file with the same name exists in the FlexAIDdS repo and many build trees:

| Field | Value |
|--------|--------|
| Path (example) | `WRK/MC_st0r5.2_6.dat`, `build_v137/MC_st0r5.2_6.dat`, repo root `MC_st0r5.2_6.dat` |
| **MD5** | **`9dc93717dfed0698006d88dd6a9627bc`** |
| Size | **18040** bytes |
| vs pin | **≠** claim pin; `cmp` reports thousands of differing bytes |

**Do not** use `9dc93717…` for claim or JCIM-comparative benchmarks.

Historical note: matrix edit experiments are documented in `docs/VCT_MATRIX_AUDIT.md` / enthalpy audits. Those edits are **research only** and must never silently replace the JCIM pin in production OUT.

## Staging for live claim (anti-hang local)

Immutable copy + live data dir (operator machine):

```text
~/flexaidds_results/pins/MC_st0r5.2_6.dat.JCIM2015_claim_pin   # md5 72d7c739…
~/flexaidds_results/three_engine_entropy_q1/data/MC_st0r5.2_6.dat
~/flexaidds_results/three_engine_entropy_q1/bin/C/MC_st0r5.2_6.dat
```

Source of the pin copy on this machine:

```text
$FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1/data/MC_st0r5.2_6.dat
```

`FLEXAIDDS_DATA_DIR` for claim must point at a directory containing **only** the `72d7c739…` matrix (or the file at that path must hash to that MD5).

## Reproducibility receipt keys

Every campaign `RUN_RECEIPT.json` / `provenance.json` must include:

```json
"matrix_filename": "MC_st0r5.2_6.dat",
"matrix_md5": "72d7c7396702331d96ff12d18f831796",
"matrix_sha256": "2265b1ba08e763887ab49b0bfbcfacf0505fd1c4b612fcbd3ec1fc851accd71c"
```

`aggregate_claim_metrics.py` drops rows whose `matrix_md5` ≠ pin.

## What to do if repo WRK is wrong

1. Do **not** auto-overwrite git history of `WRK/MC_st0r5.2_6.dat` without an explicit operator decision.  
2. For all live docks: stage `72d7c739…` into `FLEXAIDDS_DATA_DIR` as above.  
3. Optionally open a separate PR to restore repo `WRK/` / root matrix to the pin **after** confirming against FlexAID 2015 archival media / lab golden file.

## Related

- `benchmarks/protocols/three_engine_entropy_comparison.md` §1.1  
- `benchmarks/protocols/admission_metrics_contract.md`  
- `docs/implementation/v137_3dsig_clean_run.md` (campaign contract)  
