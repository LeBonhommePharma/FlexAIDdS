# Soft-wall uncap (k_wal · o²)

**Status:** production default (soft_wall_cutoff = 0.40 Å).  
**CF/scoring-proxy change:** deep-clash ranking (overlap o ≳ 1 Å) differs from pre-uncap soft-core.

## Physics

Soft-core path (`soft_wall_cutoff > 0`):

```text
o = cr − d
E = k_wal · max(o, 0)²     (uncapped, C¹ at o=0)
```

- `cr = permeability · (r_i + r_j)`
- Default `k_wal = 50` (overridable)
- Pure quadratic is C∞ for o > 0 and matches deep-wall value **and** slope (no Hermite kink)

Legacy path (`soft_wall_cutoff == 0`):

```text
E = min(KWALL · (d⁻¹² − cr⁻¹²), WAL_CONTACT_CAP=50)
```

## Why uncap

The previous soft-core form applied `min(E, WAL_CONTACT_CAP)` on the soft-core branch.
With k_wal=50 that bound activates at o ≥ 1 Å, **flattening** the wall and zeroing the GA
fitness gradient into buried poses (attractive CF.com could win unbounded).

## Knobs

| Source | Key | Default |
|--------|-----|---------|
| JSON | `flexibility.soft_wall_cutoff` | 0.40 |
| JSON | `flexibility.k_wal` | 50 |
| CONFIG.inp | `SOFTWA` / `KWAL` | 0.40 / 50 |
| Env | `FLEXAIDDS_SOFT_WALL`, `FLEXAIDDS_K_WAL`, `FLEXAID_KWAL` | — |
| ProtocolConfig / RUN_RECEIPT | `soft_wall_cutoff`, `k_wal`, `soft_wall_uncapped` | 0.40 / 50 / true |

## Startup log

```text
[soft_wall] soft_wall_cutoff=0.400 A  k_wal=50.000  soft_wall_uncapped=true
```

## Claim / receipt requirements

Any post-merge claim campaign must record in RUN_RECEIPT / protocol_config:

- `soft_wall_uncapped: true`
- `soft_wall_cutoff`, `k_wal`
- engine binary SHA

Language: this changes the **CF/contact-function scoring proxy**, not thermodynamic ΔG.

## Backend parity

CPU `vcfunction` / `Vcontacts` and batch `cpu_eval` / CUDA / Metal WAL terms use the same soft-core formula when configured.
