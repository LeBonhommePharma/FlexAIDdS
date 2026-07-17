# v135 — Crystal-blind basin recovery election (BCR-proxy)

**Status:** implemented behind flags (default OFF). Does **not** change claim ranking unless enabled.  
**Code:** `LIB/ProtocolConfig.*`, `LIB/DatasetRunner.cpp` (`select_pose_freq_gated_pooled`)  
**Launcher:** `scripts/launch_v135.py`

## Problem

Campaign diagnostics often show **election gaps**:

- **BCR** (best cluster-head RMSD ≤ 2 Å) finds a near-native pose  
- **S1** (elected head) is a deep **CF scoring-proxy** false minimum  

On 1G9V (claim): BCR CF ≈ +25 a.u., Frequency=1 vs elected CF ≈ −68 a.u., Frequency=69.  
CF is **arbitrary units**, not kcal/mol (`UNITS_CORRECTION.md`).

BCR itself is **oracle** (uses crystal). v135 is a **crystal-blind** approximation of the *idea* “consider the full head pool, not only min-CF among freq>1”.

## What v135 changes (when `FLEXAIDDS_ELECTION_V135=1`)

| Knob | Env | Default under v135 | Effect |
|------|-----|--------------------|--------|
| Master | `FLEXAIDDS_ELECTION_V135` | off | Enables package |
| Score τ | `FLEXAIDDS_ELECTION_SCORE_TAU` | **25** (CF a.u.) | Replaces legacy τ=0.592 mixed scale |
| Singletons | `FLEXAIDDS_ELECTION_INCLUDE_SINGLETONS` | **on** | Frequency=1 heads stay in the pool |

Composite (same form as before, correct unit language):

```text
score = Z(τ) · exp(−α H) · log1p(N)
τ in CF arbitrary units  — not physical kT
```

Also enumerates FastOPTICS dual-suffix poses: `<prefix>_<minPts>_<rank>.pdb`.

## What v135 does **not** claim

- Not “true ΔG” or thermodynamic Boltzmann weights  
- Not guaranteed to elect BCR when |ΔCF| ≫ τ (e.g. +93 a.u. on 1G9V still loses)  
- Not a substitute for fixing the CF false-minimum physics  

When ΔCF is huge, only **scoring** (or multi-restart geometry consensus) can recover; v135 mainly fixes **selection bookkeeping** (freq gate + unit-mixed τ).

## Claim path

Default claim (`FLEXAIDDS_ELECTION_V135` unset): **unchanged** ranking (AGENTS.md).  
Ablation / research runs: use `scripts/launch_v135.py` or export the three env vars on a **new** OUT namespace.

## Offline validation (no re-GA)

On finished pose pools:

```bash
# compare elected under legacy vs v135 knobs after a flagged rebuild
python3 scripts/aggregate_claim_metrics.py ...
```

Log tag when active: `[V135-NNBR] crystal-blind basin recovery election: ...`
