# Experimental Capabilities

This file lists visible repository surfaces that should be treated as experimental until they are explicitly promoted into `docs/VALIDATED_CAPABILITIES.md`.

Experimental means one or more of the following:

- incomplete CI coverage
- unstable API or UX contract
- incomplete installation documentation
- incomplete reproducibility artifacts
- performance path exists but is not release-gated

## Experimental surfaces

At the current stage, the following should be treated as experimental:

- Swift packages and Apple-platform integration layers
- TypeScript, PWA, dashboard, and browser-facing viewers
- Bonhomme Fleet and iCloud-driven distributed execution
- NATURaL and related co-translational or co-transcriptional workflows.
  Default `FlexAIDdS` docking is **fail-closed**: DualAssembly growth does not
  run unless `--natural` / `FLEXAIDDS_NATURAL=1` / `advanced.enable_natural=true`.
  Without that opt-in, `FA->natural_deltaG` stays 0.
- PoseHelixThermoRewrite / pose→helix ΔH/ΔS rewrite (2014 NATURAL seminar glue; default OFF; not Astex / not `G_natural`). DualAssemblyRunner / DualAssemblyEngine feed GA atom coords only when `enable_pose_helix_rewrite` or `FLEXAIDDS_POSE_HELIX_THERMO_REWRITE=1` is set; missing coords fail closed. See `docs/POSE_HELIX_THERMO_REWRITE.md`.
- PoseLocalThermoRewrite (`LIB/NATURaL/PoseLocalThermoRewrite.{h,cpp}`): experimental
  docking-pose rewrite of local secondary-structure ΔH/ΔS. RNA, DNA, protein
  α-helix, and protein β-sheet share the same algebra and use **class-specific**
  increment tables (Xia 1998 RNA NN ≠ SantaLucia 1998 PNAS Table 2 DNA ≠
  Scholtz 1991 / Zavrtanik 2026 helix ≠ Meier–Seelig 2008 sheet midpoint). Default OFF
  (`DualAssemblyConfig::enable_pose_local_thermo_rewrite`,
  `FLEXAIDDS_POSE_LOCAL_THERMO_REWRITE`). DualAssemblyRunner prefers GA atom
  coords (`GAResult::receptor_nts` / `ligand_poses`); labelled `pose_rewrite_poses`
  remain a fallback. Missing both sources fail closed. Motifs without a calorimetric
  consensus stay unset (`TODO(burgundy)` cite slots) — no invented placeholders.
  Does not feed `G_natural` / Astex / FlexADS claim contracts. See
  `docs/POSE_LOCAL_THERMO_REWRITE.md`.
- backend-specific acceleration paths not required by the Core 1.0 support matrix
- benchmark claims not yet backed by a repository reproducibility bundle

## Promotion rule

A capability should remain experimental until it has:

1. documentation sufficient for external use
2. automated validation or release validation coverage
3. an explicit place in the support matrix

## Experimental Thermodynamic Features

The following are intentionally kept experimental (see `docs/thermodynamics.md` and individual task PRs):

- Joint receptor–ligand ensemble analysis (`JointEnsembleResult`, mutual information)
- Standard-state affinity calibration and Kd conversion (safe utilities only; `calibrated=false` guard)
- Temperature scan + model-derived ΔCp fitting (explicitly labelled `model_derived` / `experimental`)
- Cleft annotation and flexible residue selection (preprocessing only)
- Advanced compensation / enthalpy-entropy diagnostic metrics (for analysis only)
- PoseHelixThermoRewrite (`LIB/NATURaL/PoseHelixThermoRewrite`) — docking-pose mixture that
  rewrites RNA decision-helix ΔH/ΔS (loop ΔS, stem stacking ΔH). Experimental; default OFF;
  DualAssemblyEngine::run() does not apply it. DualAssemblyRunner / DualAssemblyEngine
  sidecars accept GA atom coords only behind `enable_pose_helix_rewrite` /
  `FLEXAIDDS_POSE_HELIX_THERMO_REWRITE`. See `docs/POSE_HELIX_THERMO_REWRITE.md`.
- PoseLocalThermoRewrite (`LIB/NATURaL/PoseLocalThermoRewrite`) — same algebra generalized
  to RNA/DNA/helix/sheet with class-specific tables. Experimental; default OFF.
  DualAssemblyRunner opt-in: `enable_pose_local_thermo_rewrite` /
  `FLEXAIDDS_POSE_LOCAL_THERMO_REWRITE`. See `docs/POSE_LOCAL_THERMO_REWRITE.md`.

These features are fully implemented with tests and JSON exposure but require additional benchmarking or calibration data before promotion to validated status.
4. an unambiguous ownership in the product boundary defined by `PRODUCT.md`
