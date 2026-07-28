# Benchmark self-eval contract (a priori + a posteriori)

**Status:** Canonical checklist for Phase-4+ FlexAIDdS experiments.  
**Defers to:** [`METHODOLOGY.md`](../METHODOLOGY.md) (numbers, RMSD, matrix, claim language) and  
[`PHASE4_GATES_ACTUALIZED.md`](PHASE4_GATES_ACTUALIZED.md) (L1–L4, magnitude floors, class-matched panels).  
Do **not** fork budgets or success floors here — re-read those files.

---

## Why this exists (a posteriori learning)

| Failure mode observed | Rule encoded |
|----------------------|--------------|
| Wrong panel class (pb_clash on SEARCH-MISS) | Class-matched lever only |
| No magnitude floor (noise as PASS) | Sampling: mean ΔBCR ≤ −0.5 Å **or** ≥1 BCR&lt;2 |
| Unwired knob (BOOM_INTERVAL, WAL) | L2: env must beat claim JSON; L4 log fire |
| Early-stop truncating every lever | `FLEXAIDDS_NO_SEC=1` on Phase-4 docks |
| Averaging 1J3J SCORING_PULL with 1N1M near-miss | Two-tier SEARCH-MISS reporting |
| Emission heads as full-pop ceiling | DUMP_POP for any-pose claims |
| Dual-dock / lock thrash | Sol #9 preflight / release |

---

## Status enum (a posteriori only)

| Status | Meaning |
|--------|---------|
| **PASS** | L1–L4 live + magnitude ACCEPT floors met |
| **FAIL** | L1–L4 live + floors missed (honest negative) |
| **PASS_LIVENESS** | Knob fires; no magnitude claim |
| **VOID** | Wrong instrument / structurally inert |
| **INVALID** | Unwired or multi-var confounded |
| **MISSING_OUT** | Numbers only from workorder; OUT gone |
| **IN_PROGRESS** | Dock running; do not invent BCR |
| **NOT_RUN** | Never launched |

Never re-label VOID/INVALID as FAIL to “rank” levers.

---

## A priori checklist (fill **before** launch)

Copy into OUT as `APRIORI.json` (or fill via `scripts/benchmark_self_eval.py preflight`).

| Field | Required value / rule |
|-------|------------------------|
| `one_variable` | Exactly one intentional delta vs matched control |
| `panel_class` | `SEARCH_MISS` \| `SCORING_LOCKED` \| `NEAR_MISS` \| `GROSS_MISS` |
| `codes` | Explicit list; near-miss default `1N1M,1L7F` |
| `l1_l4_plan` | How knob is read, not JSON-stuck, can act, log marker |
| `magnitude_floor` | Sampling: PHASE4 floors; scoring: \|dCF\|≥1 + sign flip |
| `matched_control` | Same binary SHA, matrix, R, pop/gen, only var differs |
| `matrix_pin` | **9dc9** (`md5` of `MC_st0r5.2_6.dat`) |
| `no_sec` | `true` for Phase-4 sampling docks |
| `sol9` | `benchmark_coord.py preflight` token held |
| `workers` | ≤4; prefer 2 on this host |
| `restarts` | State R; claim Astex uses 10 (METHODOLOGY); pilots may use 1–5 |
| `report_tiers_separately` | Never mean gross-miss into near-miss lever read |
| `forbid` | dual-dock, silent re-run of matched FAIL, burial re-panel |

---

## A posteriori checklist (fill **after** arm completes)

| Field | Rule |
|-------|------|
| `status` | Enum above from on-disk OUT only |
| `l4_evidence` | Paths + counts of log markers (`[BOOM]`, `[NICHE-CART]`, …) on **stderr + r\*** (not stdout-only) |
| `per_target_metrics` | BCR, S3, elect RMSD, gens_reached if available |
| `delta_vs_control` | mean ΔBCR (treatment − control) on **same codes** |
| `wipeout` | false unless gen~300 + CF≈0 signature |
| `verdict_reason` | One line: scientific FAIL vs instrument VOID |
| `flip_order_applied` | Path to `campaign_flip_order.py` output if used |

### S2 closed-gate pin pack (required for publication audit)

Every **closed** OUT root must pass:

```bash
python3 scripts/benchmark_self_eval.py validate-pins --out OUT
```

| Pin | Rule |
|-----|------|
| `evidence/accept.txt` | Non-empty; machine-readable accept/status lines (`ACCEPT_*=True/False`, `status=…`) |
| Per-arm **binary_sha256** | `evidence/arm_pins.json` with `arms.<name>.binary_sha256`, **or** hashable `arm_<name>/bin/FlexAIDdS.stamped` |
| Matrix | Prefer `arm_pins.json` `matrix_pin` = **9dc9** / full md5 |
| Shared binary | If all arms share one stamp, set `"shared_binary": true` in `arm_pins.json` and still list each arm’s SHA |

`arm_pins.json` schema (minimal):

```json
{
  "matrix_pin": "9dc93717dfed0698006d88dd6a9627bc",
  "shared_binary": true,
  "arms": {
    "control": {"binary_sha256": "<64 hex>", "git_tip": optional},
    "mut_gran": {"binary_sha256": "<64 hex>", "git_tip": optional}
  }
}
```

Missing accept.txt or any arm SHA → **PINS_FAIL** (exit 2). Does not re-dock or rewrite scientific status.

---

## Publication residual path (unblock sequence)

1. Phase-4 **sampling ACCEPT** on **near-miss** class (magnitude floor).  
2. Optional: non-burial scoring arm for SCORING-LOCKED (separate).  
3. Full-85 per `METHODOLOGY.md` § claim protocol (R=10, spyrmsd, matrix 9dc9, autonomous).  
4. Claim language: CF/contact-function proxy vs full thermo only if STRICT (PB + tENCoM) validated (`AGENTS.md`).

---

## Agent enforcement

```bash
python3 scripts/benchmark_self_eval.py preflight --apriori APRIORI.json --write-out OUT
# after dock:
python3 scripts/benchmark_self_eval.py posteriori --control CTRL_OUT --treatment NAME PATH ...
python3 scripts/benchmark_self_eval.py validate-pins --out OUT
python3 scripts/campaign_flip_order.py g4_1 ...   # or g4_3
```
