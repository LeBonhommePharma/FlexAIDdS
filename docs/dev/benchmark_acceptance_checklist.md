# Thermodynamic Feature Benchmark Acceptance Checklist

Use this checklist before promoting any thermodynamic feature from experimental → validated.

For each feature, answer **Yes / No / N/A** and provide evidence (PR, test, doc link).

## General Questions (apply to all features)

- [ ] Does the feature change any ranking, scoring, or pose selection behaviour?
- [ ] Are all units explicit in field names and documentation?
- [ ] Are all sign conventions documented (especially for energies and entropy)?
- [ ] Is the reference state (standard state concentration, temperature, etc.) explicit?
- [ ] Is calibration required before the output can be treated as physical?
- [ ] Is the output fully backward compatible with existing JSON/CSV/REMARK consumers?
- [ ] Do C++ and pure-Python paths produce matching results (within tolerance) on toy systems?
- [ ] Are experimental / diagnostic / model-derived fields clearly labelled in code, JSON, and docs?
- [ ] Are all scientific claims backed by unit tests or a reproducibility bundle?
- [ ] Has the feature been exercised in at least one end-to-end test or benchmark workflow?

## Feature-Specific

### ThermodynamicBreakdown + Component Averages
- [ ] G_config, H_eff, S_config, Cv, σ_E identities verified on single-state, two-equal, two-unequal, and multiplicity cases?
- [ ] Component means reduce correctly to per-pose values on 1-pose ensembles?
- [ ] `components_complete` flag behaves correctly when decomposition is partial?

### Diagnostic Metrics (compensation_score etc.)
- [ ] Zero-denom and near-zero safety verified?
- [ ] Metrics never appear in any ranking or selection path?
- [ ] Documentation explicitly forbids their use for affinity claims?

### Joint Receptor–Ligand Ensemble
- [ ] Fallback (no receptor conformer IDs) correctly sets S_receptor=0 and I=0?
- [ ] Probabilities sum to 1 in all tested cases?
- [ ] Mutual information is zero for independent distributions?

### Affinity Calibration
- [ ] Round-trip ΔG° ↔ Kd validated?
- [ ] Invalid T and invalid Kd are rejected with clear errors?
- [ ] `calibrated=false` prevents emission of `predicted_Kd_M` as real affinity?

### Temperature Scan + ΔCp
- [ ] T ≤ 0 rejected?
- [ ] Scan at T_ref matches normal `compute()` within tolerance?
- [ ] Fit refuses <4 points?
- [ ] Output always carries `model_derived=true` and `experimental=true`?

### Cleft Annotation + Flexible Residue Selector
- [ ] Fixed residues are never included?
- [ ] Forced-flexible residues are always included (when valid)?
- [ ] No duplicates and deterministic ordering guaranteed?
- [ ] Selector is pure preprocessing (no scoring side effects)?

## Sign-off

Feature: _______________________________

Reviewed by: __________________________

Date: ________________________________

Evidence links (PRs, tests, docs): _______________________________

Decision: [ ] Promote to validated   [ ] Keep experimental   [ ] Needs more work

Notes: ______________________________________________________________