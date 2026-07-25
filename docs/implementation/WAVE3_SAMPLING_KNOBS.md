# Wave 3 sampling knobs (default-OFF)

Env-gated diversity / niche controls for BCR raisers. **No default path change.**

| Env | Default | Effect |
|-----|---------|--------|
| `FLEXAIDDS_BOOM_INTERVAL` | unset → JSON (100) | Generations between BOOM random injection; `0` disables |
| `FLEXAIDDS_BOOM_FRAC` | unset → JSON (1.0) | Fraction of worst half replaced each injection |
| `FLEXAIDDS_SIGMA_SCALE` | unset → 1.0 | Multiplies `GB->scale` (niche radius / sharing) |
| `FLEXAIDDS_COARSE_ORIENTATIONS` | unset → JSON (64) | Coarse-init orientation count (1–4096) |
| `FLEXAIDDS_MEMETIC` | off | **Ignored** unless `FLEXAIDDS_WALL_PILOT_PASS=1` (W2 wall oracle PASS) |
| `FLEXAIDDS_WALL_PILOT_PASS` | unset | Set only after score-only wall oracle PASS |

## Pilot recipe (after live baseline finishes)

```bash
# One variable at a time; workers=2; 1J3J + 1K3U; restarts=1 diversity bar
export FLEXAIDDS_BOOM_INTERVAL=50
# vs control unset
export FLEXAIDDS_SIGMA_SCALE=0.5   # or 2.0 — factorial later
export FLEXAIDDS_COARSE_ORIENTATIONS=256
```

## Memetic (W3.4)

**SKIP** unless wall oracle records `wall_pilot_pass: true`. Enabling memetic before wall uncap risks CF-minimizing into over-burial (Opus E5).

## Cartesian niche (W3.2)

Not fully landed as gene0-decoupled distance in this wave — tracked as follow-up.
`SIGMA_SCALE` is the interim niche-radius lever.
