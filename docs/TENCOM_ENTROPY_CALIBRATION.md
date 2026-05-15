# tENCoM Entropy Calibration Status

tENCoM and ENCoM-derived vibrational entropy are currently supported as a
relative, protocol-fixed heuristic in FlexAIDdS.

The important boundary is simple: ENCoM eigenvalues are model-scale quantities.
They are not physical SI frequencies unless an eigenvalue-to-frequency
calibration is supplied and documented for the benchmark system. Therefore:

- differential comparisons such as `Delta S_vib = S_vib(holo) - S_vib(apo)` are
  acceptable when apo/holo are processed with the same protocol
- absolute `S_vib` and `-T*S_vib` magnitudes must be labeled heuristic unless a
  calibration artifact is present
- benchmark tables must not present tENCoM vibrational entropy as absolute
  thermodynamic truth without calibration provenance

## Required calibration artifact

Before promoting absolute vibrational entropy claims, commit a bundle with:

- reference structures and preprocessing commands
- ENCoM/tENCoM parameters
- eigenvalue-to-frequency scale or normal-mode reference
- expected `S_vib`, `Delta S_vib`, and `-T*S_vib` values
- metric scripts and tolerances
- git SHA, compiler, backend, thread count, and temperature

Until then, use labels such as `S_vib_heuristic` or `relative Delta S_vib`.
