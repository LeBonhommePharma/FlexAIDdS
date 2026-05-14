# tENCoM Calibration Bundles

This directory stores machine-readable calibration metadata for converting
tENCoM model eigenvalues into angular frequencies.

- `schema.json` defines the required fields.
- `uncalibrated_model_scale.json` documents the default model-scale behavior.

Until a calibrated bundle is added, outputs tagged
`model_scale_heuristic` must be treated as relative flexibility scores, not
absolute vibrational entropy measurements.
