# Expected Outputs — Astex Non-Native Tier-1

This directory will contain reference results after the first validated
cross-docking benchmark run. Expected files per pair:

- `{target_pdb}_x_{ligand_pdb}/best_pose.pdb` — Top-ranked cross-docked pose
- `{target_pdb}_x_{ligand_pdb}/binding_modes.json` — Binding modes + thermodynamics

Aggregated metrics (computed by DatasetRunner) must match the baselines
defined in `../manifest.yaml` within `baseline_tolerance`.

## Success Criterion

RMSD < 2.0 Å vs the co-crystal pose in the ligand source PDB (`ligand_pdb`).
This is harder than native docking: the receptor was crystallized with a
*different* ligand, so the binding-site conformation may not match.
