# BOOM small-fraction A/B — PASS

**Written:** 2026-07-25T23:14:46.035267+00:00  
**One variable (Arm B):** `FLEXAIDDS_BOOM_FRAC=0.1` (interval left at claim JSON **100**)  
**Target:** 1N1M · R=2 · gen=2000 · autonomous · matrix 9dc9  
**OUT:** `/Users/lp.more/flexaidds_results/boom_frac_ab_20260725_184135`

**Verdict:** **PASS** — B live (JSON 0 + [BOOM] inject frac≈0.1); A silent; no CF≈0 wipeout signature

## Prior liveness (required)

Smoke `boom_liveness_smoke_20260725_183333`: JSON frac=0, `[BOOM]` gen 50/100 with n_inject=50/1000 → **env wins**.

Claim path hardcodes `boom_inject_fraction: 0.0` **deliberately** (frac=1.0 wiped blind GA). Small frac only.

## Results

| Check | Arm A (control) | Arm B (frac=0.1) |
|-------|-----------------|------------------|
| JSON boom_inject_fraction | 0.0 | 0.0 |
| JSON boom_inject_interval | 100.0 | 100.0 |
| [BOOM] line count | **0** | **11** |
| n_inject / pop (first) | None/None | 50/1000 |
| implied frac | n/a | 0.100 |
| max SMFREE gen (any restart) | 950 | 750 |
| restart max gens | {'1N1M': 550, 'r1': 950} | {'1N1M': 450, 'r1': 750} |
| elected CF | -99.314 | -99.314 |
| rmsd_hungarian | 5.66 | 5.66 |
| best_cluster_rmsd | 4.55 | 4.55 |
| seed_echo | 0.0 | 0.0 |

### Arm B [BOOM] samples
```
[BOOM] injection #1 at gen 100: re-randomized worst 50/1000 chromosomes (fresh random, better half preserved)
[BOOM] injection #2 at gen 200: re-randomized worst 50/1000 chromosomes (fresh random, better half preserved)
[BOOM] injection #3 at gen 300: re-randomized worst 50/1000 chromosomes (fresh random, better half preserved)
[BOOM] injection #4 at gen 400: re-randomized worst 50/1000 chromosomes (fresh random, better half preserved)
[BOOM] injection #1 at gen 100: re-randomized worst 50/1000 chromosomes (fresh random, better half preserved)
[BOOM] injection #2 at gen 200: re-randomized worst 50/1000 chromosomes (fresh random, better half preserved)
```

### Termination notes
**A:**  
- 1N1M: TIMING SUMMARY: 599 gens timed, avg 199.06 ms/gen, ~99.53 us/eval (2x-pop est), est 398.1 s for 2000-gen x 1000-chrom run
- r1: TIMING SUMMARY: 999 gens timed, avg 204.77 ms/gen, ~102.38 us/eval (2x-pop est), est 409.5 s for 2000-gen x 1000-chrom run
- 1N1M/stdout: GA terminated: CF stagnant for 400 gens (best_CF=-98.5618) with gene-space collapsed
- 1N1M/stdout: GA terminated early by fitness stagnation
- r1/stdout: GA terminated: CF stagnant for 300 gens (best_CF=-99.3141) with gene-space collapsed
- r1/stdout: GA terminated early by fitness stagnation

**B:**  
- 1N1M: TIMING SUMMARY: 489 gens timed, avg 199.45 ms/gen, ~99.73 us/eval (2x-pop est), est 398.9 s for 2000-gen x 1000-chrom run
- r1: TIMING SUMMARY: 799 gens timed, avg 199.01 ms/gen, ~99.50 us/eval (2x-pop est), est 398.0 s for 2000-gen x 1000-chrom run
- 1N1M/stdout: Entropy convergence at generation 490 (H=2.0340 nats, stable for 5 checks)
- 1N1M/stdout: GA terminated early by entropy convergence
- r1/stdout: GA terminated: CF stagnant for 500 gens (best_CF=-97.2527) with gene-space collapsed
- r1/stdout: GA terminated early by fitness stagnation

## Interpretation

- **Lever validity (primary):** BOOM_FRAC=0.1 is **instrumented and live** on claim JSON-0 path. Control stays silent.
- **Docking success (secondary):** not a campaign success-rate gate; 1N1M N=1×R=2 is noise for RMSD.
- **Early stop:** both arms may terminate under 2000 gens (SEC/stagnation) — same class as audit L2. Watch that B does not add CF≈0 wipe at ~300.
- **Do not** set frac=1.0. **Do not** treat STEP 3 INTERVAL-only pilot as BOOM evidence.

## Cadence

- Phase: BOOM lever re-test (post B1 + OPS caveat)  
- One variable: `BOOM_FRAC=0.1`  
- Genuine metrics secondary on this smoke  
- **PASS**

## Next

1. Optional: panel of 5–8 with frac=0.1 only if lever PASS and no wipe signature.  
2. Still open: **pb_clash burial oracle** (wall STEP 2 replacement).  
3. No full-85; no memetic; no WALL_PILOT_PASS.
