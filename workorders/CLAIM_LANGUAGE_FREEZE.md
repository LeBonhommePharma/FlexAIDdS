# Claim language freeze (CF proxy vs thermodynamics)

**Status:** FREEZE for Methods / SI / agent claims  
**Date:** 2026-07-28  
**Aligned with:** `AGENTS.md` Scientific Guardrails, `METHODOLOGY.md`,  
`BENCHMARK_SELF_EVAL_CONTRACT.md` publication residual path.

---

## Allowed (default for Phase-4 / null-stack results)

| Phrase | When |
|--------|------|
| **CF/contact-function scoring proxy** | Always for VoronoiCF / Vcontacts ranks |
| **Elected pose / top CF pose** | Rank-0 under production LOCCLF |
| **best_cluster_rmsd (BCR)** | Pool/sampling ceiling metric |
| **Ensemble-derived free energy estimate** | Only if StatMech layer ran and is labeled as estimate |
| **RMSD to crystal ≤ 2.0 Å** | Geometry only — not “binding success” alone |
| **PASS_LIVENESS** | Knob fired; no magnitude claim |
| **FAIL (null magnitude)** | Honest negative vs a priori floor |
| **SEARCH-MISS / SCORING-LOCKED** | Per inversion map class only |

## Forbidden without STRICT validation

| Phrase | Why |
|--------|-----|
| **True binding free energy ΔG** | Requires full partition function + validated thermo stack |
| **Computed ΔG / experimental-grade affinity** | Overclaim |
| **Docking success** from RMSD alone | Need RMSD≤2 **and** PoseBusters (and STRICT: tENCoM/Eigen) |
| **pb_clash / WAL fixed scoring** | Oracles VOID/FAIL; do not claim unlock |
| **Phase-4 sampling ACCEPT** | Near-miss stack is **null** (freeze table) |
| **Full-85 claim rates** | Blocked until sampling ACCEPT |

## STRICT path (only when all true)

Allowed only if documented on-disk: production dock + `result.csv` / RUN_RECEIPT +  
RMSD≤2.0 Å **and** PoseBusters pass **and** tENCoM/Eigen status OK. Then:

- “STRICT-validated pose” / “thermodynamic ledger (F, H, −TS, Cv) with tENCoM correction”  
- Still separate CF proxy ranking from thermo ledger unless product integrates them.

## Agent / paper checklist

1. Cite CF proxy when discussing ranks from G4.1 / election / G4.3.  
2. Cite freeze table for sampling nulls — do not invent ACCEPT.  
3. SCORING-LOCKED text points to SI package, not more BOOM.  
4. Never set memetic unlock flags from burial/pb_clash fails.

## Non-claims for this freeze doc

- Not a new scientific result.  
- Not a dock.  
- Not permission to weaken magnitude floors.
