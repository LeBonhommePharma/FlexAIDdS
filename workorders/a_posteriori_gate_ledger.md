# A posteriori gate ledger (live)

## G4.1 BOOM near-miss — CLOSED
- OUT: g4_1_boom_near_miss_20260726_200953
- L4: PASS (control 0, treatments 236 [BOOM] on stderr after scanner fix 70ed4f51)
- Magnitude: FAIL null (best mean_dBCR=−0.0192 at frac010; floor −0.5)
- accept_g4_1: False
- flip: election_fix_P0 (merged with 1N1M DUMP_POP offline)
- SCRATCH: flip_order_decision.txt, g4_1_posteriori_read.txt, g4_1_fixed_l4/

## ELECTION_V135 — IN FLIGHT
- OUT: election_v135_near_miss_20260726_225823
- One var: FLEXAIDDS_ELECTION_V135=1 vs control
- Codes: 1N1M,1L7F · R=5 · matrix 9dc9 · NO_SEC
- Floor: 1N1M elect ≤2.5 OR elect gap shrink ≥1.0 Å; no wipeout
- Evaluator: OUT/evaluate_on_complete.sh
- Sol#9 lock: daa3e200…

## Still blocked for full-85
Phase-4 ACCEPT not yet (BOOM null; election pending).
