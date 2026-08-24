# Blind Astex-85 receipt protocol

Wave 4 protocol for a **receipted blind** Astex Diverse 85 republish. This
document does not report a success rate.

## Pins

| Field | Value |
|---|---|
| N | 85 |
| `native_pose_seeded` | 0 |
| `seed_echo` | 0 |
| Matrix | `MC_st0r5.2_6.dat` MD5 `72d7c7396702331d96ff12d18f831796` |
| Default seed | `FLEXAIDDS_SEED_ELITISM=0`, `FLEXAIDDS_NATIVE_SEED_FRAC=0` |
| Claim metric | rank-0 in-place RMSD **≤ 2.0 Å** (METHODOLOGY.md §0) |

`SEED_ELITISM=1` / `NATIVE_SEED_FRAC=0.90` is the oracle ceiling and is **not**
the default. Cite METHODOLOGY.md §3.

## Commands (no 85-target dock)

```bash
python3 scripts/blind_astex85_receipt_protocol.py --help
python3 scripts/blind_astex85_receipt_protocol.py validate-defaults
python3 scripts/blind_astex85_receipt_protocol.py emit --out /tmp/astex85_receipt --dry-run \
    --git-commit "$(git rev-parse HEAD)" --binary-path unspecified
python3 scripts/blind_astex85_receipt_protocol.py claim --dir /tmp/astex85_receipt
# REFUSE: receipt OK, but no result.csv. Not printing a success %.

bash scripts/reproduce_astex85.sh --dry-run
```

`claim` prints an S1 percentage only when `RUN_RECEIPT.json` passes the blind
gates **and** a 85-row results CSV is present. Missing receipt → no %.

A live 85-target dock is **not** part of this PR and must not be launched from
the dry-run / unit-test path.
