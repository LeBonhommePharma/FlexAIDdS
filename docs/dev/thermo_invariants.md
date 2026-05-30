# Thermodynamic Invariants

This file is the sign and unit ledger for FlexAIDdS thermodynamic output.
Contact-function energy is an effective docking scoring-energy proxy, not a
calorimetric enthalpy unless a calibration layer explicitly says otherwise.

## Canonical Ensemble

- `beta = 1 / (kB * T)`
- `logZ = logsumexp_i(log multiplicity_i - beta * E_i)`
- `G_config = -kB * T * logZ`
- `H_eff = <E> = sum_i p_i * E_i`
- `p_i = exp(log multiplicity_i - beta * E_i - logZ)`
- `S_config = (H_eff - G_config) / T`
- `-T*S_config = G_config - H_eff`
- `Cv = Var(E) / (kB * T^2)`
- `sigma_E = sqrt(Var(E))`

Units:

- Energies and free energies: `kcal/mol`
- Entropies: `kcal/mol/K`
- Heat capacity: `kcal/mol/K`
- Temperature: `K`
- `logZ`: dimensionless natural log

## Ledger Fields

- `G_config_kcal_mol`: canonical configurational free energy from the sampled scoring-energy ensemble.
- `H_eff_kcal_mol`: Boltzmann-weighted effective energy average.
- `S_config_kcal_mol_K`: configurational entropy.
- `minus_T_S_config_kcal_mol`: `G_config_kcal_mol - H_eff_kcal_mol`.
- `G_vib_kcal_mol`: optional vibrational correction, currently heuristic/model-scale unless calibrated metadata says otherwise.
- `G_natural_kcal_mol`: optional NATURaL co-translational correction.
- `G_other_kcal_mol`: reserved explicit correction bucket; must not hide unknown physics.
- `G_total_kcal_mol`: sum of configurational and explicitly present corrections.

Legacy `Thermodynamics.free_energy` is preserved for compatibility. In
`BindingMode::get_thermodynamics()` it may include vibrational and NATURaL
corrections. New code should use `ThermodynamicBreakdown` when it needs to
distinguish configurational and correction terms.

## Grand Canonical Ensemble

- `Xi = 1 + sum_i z_i * Z_i`
- `z_i = c_i / c0`
- `p_i(bound) = z_i * Z_i / Xi`
- `p_empty = 1 / Xi`

Intrinsic selectivity must ignore concentration. Apparent selectivity includes
concentration.

## Affinity Boundary

- `DeltaG_standard = R*T*ln(Kd / c0)`
- `Kd = c0 * exp(DeltaG_standard / (R*T))`

Do not expose uncalibrated docking scores as real `Kd`, `Ki`, or affinity.
`IC50` is not `Kd` without an explicit conversion model.

## Forbidden Ranking Metric

Do not use `(DeltaH + T*DeltaS) / DeltaG` as a ranking score.

Allowed compensation diagnostics must remain diagnostic-only and must not feed
sorting, scoring, docking, pose selection, or optimization.

