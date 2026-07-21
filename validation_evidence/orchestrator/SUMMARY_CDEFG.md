# Orchestrator summary — live docks C–F + G local (2026-07-21)

## CLOSED: A B C D E F (6/7) · OPEN: G (.sif)

| Item | Key number(s) | Path |
|------|----------------|------|
| C | CMA smoke best_cf=457.54, DOF=10, n_snap=32, evals=3000 | `validation_evidence/build_ab/C_cmaes_smoke/` |
| D short | GA CF=-37.12 RMSD=10.11; CMA CF=643.05 RMSD=12.13 (5000 evals) | `validation_evidence/build_ab/D_1G9V_dual/` |
| E 2e6 | GA CF=-68.56 RMSD=5.56 (355s); CMA CF=-20.19 RMSD=6.52 evals=2e6 (526s) | `validation_evidence/build_ab/E_ab_2e6/E_ab_summary.txt` |
| F | fingerprint sha256=6b051f0d…0cac; 2000 gens; best_cf_end=-20.19 | `validation_evidence/build_ab/F_trace/` |
| G | manifest OK; apptainer MISSING | `validation_evidence/build_ab/G_harness/G_local.txt` |

Binary sha256: `404b3ccddc22c12bf3cfaced9b0eaf996faa16d5caa501d7ff691282d66f9eb1`

No fabricated poses. All RMSDs = ordered heavy-atom ligand (RQ3) vs crystal SDF.
