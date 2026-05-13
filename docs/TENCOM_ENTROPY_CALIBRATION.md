# tENCoM Eigenvalue-to-Frequency Calibration

## Executive Rule

tENCoM eigenvalues are stiffness-like model quantities. They are useful for
ranking flexible modes, sampling torsional perturbations, B-factor-like
fluctuation profiles, and relative rigidity changes. They are not automatically
physical angular frequencies.

An absolute vibrational entropy or `-T*S_vib` free-energy claim requires one of
the following:

1. A generalized normal-mode problem with an explicit mass/inertia metric:
   `K q = omega^2 M q`.
2. An empirical calibration bundle that maps model eigenvalues to angular
   frequency:
   `omega_i = eigenvalue_to_omega * sqrt(lambda_i)`.

Without that bundle, FlexAIDdS reports tENCoM vibrational entropy as
`model_scale_heuristic`. Those numbers can still be useful for relative
comparisons under the same model settings, but they are not absolute
thermodynamic observables.

## Why This Exists

The torsional Hessian has the form:

```text
V(theta) = 1/2 * delta_theta^T K delta_theta
```

For tENCoM, `K` is assembled from contact springs and torsional Jacobians. Its
eigenvalues measure stiffness in the torsional coordinates. A real vibrational
frequency also needs inertia. In first-principles language, stiffness alone is
like knowing how hard a spring is without knowing the mass attached to it.

The classical harmonic-oscillator entropy formula used here is:

```text
S_mode = kB * [1 + ln(kB*T / (hbar*omega_i))]
```

That equation is dimensionally valid only when `omega_i` is in rad/s. If
`omega_i = sqrt(lambda_i)` is just a model-scale proxy, the formula produces a
relative entropy score, not an absolute entropy.

## Upstream ENCoM/NRGTEN Context

The original ENCoM literature and NRGTEN use vibrational entropy differences as
predictive model features. NRGTEN's current public documentation explicitly
states that ENCoM is pseudo-physical and that its vibrational entropy value has
no definite units. NRGTEN also preserves the older ENCoM rigid-rotor entropy path
for reproducibility and computes eigenfrequencies from `sqrt(eigenvalue)`.

FlexAIDdS therefore keeps the useful relative model but now records calibration
status instead of silently promoting model eigenvalues into SI-looking absolute
frequencies.

Primary references:

- Frappier V, Najmanovich RJ. 2014. ENCoM model. PLoS Computational Biology.
  https://doi.org/10.1371/journal.pcbi.1003569
- Frappier V, Najmanovich RJ. 2015. Vibrational entropy differences and protein
  engineering. Protein Science. https://doi.org/10.1002/pro.2592
- Mailhot O, Najmanovich R. 2021. NRGTEN. Bioinformatics.
  https://doi.org/10.1093/bioinformatics/btab189
- NRGTEN documentation, vibrational entropy notes:
  https://nrgten.readthedocs.io/en/latest/enm.html

## Calibration Bundle Contract

Calibration metadata should be stored as JSON matching
`docs/tencom_calibration/schema.json`. The minimum fields are:

```json
{
  "schema_version": 1,
  "model": "tENCoM",
  "status": "calibrated",
  "eigenvalue_to_omega_rad_s_per_sqrt_unit": 1000000000000.0,
  "label": "example-qm-fit",
  "provenance": "Fit against QM or MD normal-mode frequencies for the stated training set",
  "absolute_entropy_claims_allowed": true
}
```

Use `status: "uncalibrated"` and `absolute_entropy_claims_allowed: false` when
the scale is only the default model scale.

## CLI Use

Default behavior:

```bash
tENCoM ref.pdb target.pdb -f all
```

Output REMARK/JSON/CSV metadata will say:

```text
S_VIB_STATUS=model_scale_heuristic
EIGENVALUE_TO_OMEGA=1
CALIBRATION_LABEL=model-scale
```

Calibrated behavior:

```bash
tENCoM ref.pdb target.pdb \
  --omega-scale 1.0e12 \
  --calibration-label qm-frequency-fit \
  --calibration-provenance "fit against stated QM normal-mode training set" \
  -f all
```

Only calibrated runs should be described as absolute vibrational entropy or
absolute `-T*S_vib` free-energy corrections.

## Practical Use Of Stiffness

Stiffness is valuable in computational chemistry flexibility simulation:

- high stiffness modes identify constrained torsions and rigid cores;
- low stiffness modes identify collective flexibility directions;
- Boltzmann sampling can use `sigma^2 = kB*T/lambda` in the model coordinates;
- B-factor-like profiles can use `1/lambda` weighting for relative fluctuations;
- differential stiffness can flag ligand-induced rigidification or loosening.

The mistake is not using stiffness. The mistake is calling stiffness-derived
model frequencies physical frequencies without the missing inertia/calibration
layer.
