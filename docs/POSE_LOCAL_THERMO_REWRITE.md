# Pose-local secondary-structure thermodynamic rewrite

**Status: experimental. Not claim-ready. Default OFF.** Draft ΔH/ΔS increments
are tunable placeholders. This path must not be cited as Astex success, FlexADS
claim-ready ranking, or a validated binding free energy.

This document specifies the physics, class-specific tables, and API for
`LIB/NATURaL/PoseLocalThermoRewrite.{h,cpp}`.

## NATURAL lineage

**NATURAL** = Native Assembly of Transcriptionally / Translationally Unified
Receptor and Ligand.

- **Transcriptional** → RNA/DNA + ligand. Cotranscriptional RNA secondary-structure
  kinetics: Zhao, Zhang, Chen, *J. Chem. Phys.* **135**, 245101 (2011),
  doi:10.1063/1.3671644 ("Cotranscriptional folding kinetics of ribonucleic acid
  secondary structures").
- **Translational** → protein + ligand. Ribosome elongation master equation:
  Zhao et al., *J. Phys. Chem. B* **115**, 3987 (2011). That JPCB paper is
  **ribosome/protein only** — do not retarget it as an RNA/DNA folding citation.

2014 seminar (LP Morency R2, 2014-07-25) Perl NATURAL: clustering, population
inheritance, pause at nt repeats, docking poses rewrite helix ΔH/ΔS. 2026
DualAssembly (`LIB/NATURaL/`) is the broader C++ home. This module generalizes
the 2014 ligand coupling from protein helix rewrite to RNA, DNA, α-helix, and
β-sheet.

`G_natural` on `ThermodynamicBreakdown` is an optional StatMech flag. This PR
does **not** feed pose patches into `G_natural` / `has_natural`. Optional later
work may do so only on a labelled path.

This rewrite is **not** `k_fold` in `RibosomeElongation`. Arrhenius coupling of
rewritten ΔG into elongation/folding rates is a later PR.

## Physics contract

The **same calculation** applies to every polymer:

1. Annotate decision elements (RNA/DNA stem-loop indices, or protein residue
   range plus optional sheet partner).
2. Docking poses contribute geometry-weighted contacts on those elements.
3. Each contact looks up **its own class table** and adds `ΔH += coeff_H · w`,
   `ΔS += coeff_S · w`.
4. `ΔG = ΔH − TΔS` with Turner/SantaLucia units (below).
5. Mix the top-k poses with Boltzmann weights from the rewritten local ΔG.

α-helix and β-sheet use the **same method** as an RNA/DNA stem-loop and
**different proper values**.

RNA hairpin ≠ DNA duplex/hairpin ≠ protein helix ≠ protein sheet.

NucleationDetector remains sequence-only (Turner / Chou-Fasman seeds). It does
not perform pose→local ΔH/ΔS rewrite.

## Units

Nearest-neighbour convention (Turner RNA, SantaLucia DNA):

| Symbol | Unit |
|--------|------|
| ΔH | kcal mol⁻¹ |
| ΔS | cal mol⁻¹ K⁻¹ (e.u.) |
| T | K |
| ΔG | kcal mol⁻¹ = ΔH − T·(ΔS / 1000) |

`kB = 0.001987206 kcal mol⁻¹ K⁻¹` (matches `statmech::kB_kcal`) is used only
for the Boltzmann mixture, not for converting ΔS.

## Class-specific increment tables (draft, experimental)

`SSClass` members:

- `RNA_STEM_LOOP` — Turner-shaped defaults
- `DNA_STEM_LOOP` — SantaLucia-shaped defaults — **never** an alias of the RNA table
- `PROTEIN_HELIX`
- `PROTEIN_SHEET`
- `PROTEIN_OTHER` — no rewrite unless annotated **and** that table is filled

Draft coefficients, geometry weight `w`:

| Class | Contact | ΔH (kcal mol⁻¹) | ΔS (cal mol⁻¹ K⁻¹) |
|-------|---------|-----------------|---------------------|
| RNA | stem stack | −1.2 w | 0 |
| RNA | loop H-bond | 0 | −3.0 w |
| DNA | stem stack | −1.0 w | 0 |
| DNA | loop H-bond | 0 | −2.5 w |
| DNA | intercalation-like | −1.5 w | −1.0 w |
| PROTEIN_HELIX | backbone H-bond | −1.5 w | −1.0 w |
| PROTEIN_HELIX | packing | −0.8 w | −2.0 w |
| PROTEIN_HELIX | clash | +2.0 w | 0 |
| PROTEIN_SHEET | bridge H-bond | −1.8 w | −1.5 w |
| PROTEIN_SHEET | hydrophobic pack | −1.0 w | −2.5 w |
| PROTEIN_SHEET | shear | +2.5 w | +1.0 w |

An RNA `ContactKind` never reads the DNA or protein tables (and conversely).
Shared algebra, separate numbers.

## API

Pure functions in `namespace natural` — no GA / `FA_Global` / DualAssemblyEngine
required.

```cpp
DecisionElement { SSClass ss; label; RNA/DNA na_start/na_end
                  OR protein res_start/res_end + optional sheet partner }
PoseView { score_kcal; contacts[] }   // score used only to pick top-k
PoseThermoRewriteConfig { ThermoIncrementTable per class; T; top_k;
                          require_decision_element_only }
LocalThermoMixture rewrite_local_thermo_from_poses(poses, elements, cfg);
```

- `require_decision_element_only` (default true): contacts outside every
  annotated element of the matching class are dropped. If nothing remains, the
  mixture is empty (finite zeros, `n_poses_used == 0`).
- Boltzmann mixture: sort by `score_kcal` ascending, keep `top_k`,
  `p_i ∝ exp(−ΔG_i / kT)` via log-sum-exp. Ranking of docking poses is not
  changed; this is a diagnostic local-SS mixture.
- No NaN in the returned mixture.

## DualAssembly hook (default OFF)

`DualAssemblyConfig::enable_pose_local_thermo_rewrite` defaults to `false`.
Environment `FLEXAIDDS_POSE_LOCAL_THERMO_REWRITE` (1/true/yes/on) may opt in.
CLI: `dual_assembly --pose-local-thermo-rewrite`.

When the flag is on **and** labelled `pose_rewrite_elements` +
`pose_rewrite_poses` are supplied, `DualAssemblyRunner` stores the mixture in
`CheckpointOutcome` diagnostic fields (`pose_local_*`). It does **not** overwrite
`dG_A_kcal` / `dG_B_kcal` or set `has_natural`. The 21-column trajectory CSV is
unchanged.

See also `docs/DUAL_ASSEMBLY_COTRANSLATIONAL.md` §10 and
`docs/EXPERIMENTAL_CAPABILITIES.md`.
