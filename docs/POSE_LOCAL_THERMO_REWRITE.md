# Pose-local secondary-structure thermodynamic rewrite

**Status: experimental. Not claim-ready. Default OFF.** Filled ΔH/ΔS increments
are peer-reviewed NN / helix–coil means (cited below). Motifs without a
calorimetric consensus stay `{0,0}` experimental gaps — they are **not**
invented placeholders. This path must not be cited as Astex success, FlexADS
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
4. `ΔG = ΔH − TΔS` with nearest-neighbour units (below).
5. Mix the top-k poses with Boltzmann weights from the rewritten local ΔG.

α-helix and β-sheet use the **same method** as an RNA/DNA stem-loop and
**different proper values**.

RNA hairpin ≠ DNA duplex/hairpin ≠ protein helix ≠ protein sheet.

NucleationDetector remains sequence-only (Turner / Chou-Fasman seeds). It does
not perform pose→local ΔH/ΔS rewrite.

## Units

Nearest-neighbour convention (Xia/Turner RNA, SantaLucia DNA):

| Symbol | Unit |
|--------|------|
| ΔH | kcal mol⁻¹ |
| ΔS | cal mol⁻¹ K⁻¹ (e.u.) |
| T | K |
| ΔG | kcal mol⁻¹ = ΔH − T·(ΔS / 1000) |

`kB = 0.001987206 kcal mol⁻¹ K⁻¹` (matches `statmech::kB_kcal`) is used only
for the Boltzmann mixture, not for converting ΔS.

## Class-specific increment tables (published vs experimental gaps)

`SSClass` members:

- `RNA_STEM_LOOP` — Xia 1998 INN-HB WC-stack **mean** (sequence-specific NN is future work)
- `DNA_STEM_LOOP` — SantaLucia 2004 NN-stack **mean** — **never** an alias of the RNA table
- `PROTEIN_HELIX` — Scholtz 1991 calorimetric helix-formation ΔH; packing/clash unset
- `PROTEIN_SHEET` — all contact kinds unset until a per-bridge consensus exists
- `PROTEIN_OTHER` — no rewrite unless annotated **and** that table is filled

Geometry weight `w` multiplies the increment. An RNA `ContactKind` never reads the
DNA or protein tables (and conversely). Shared algebra, separate numbers.

### Filled (peer-reviewed)

| Class | Contact | ΔH (kcal mol⁻¹) | ΔS (cal mol⁻¹ K⁻¹) | Source |
|-------|---------|-----------------|---------------------|--------|
| RNA | stem stack | −10.756 w | −27.9026… w | Xia et al. *Biochemistry* **37**:14719 (1998), doi:10.1021/bi9809425. Unweighted mean of the 10 unique WC **propagation** stacks (GC/CG, CC/GG, GA/CU, CG/GC, AC/UG, CA/GU, AG/UC, UA/AU, AU/UA, AA/UU) as restated in Zuber et al. *Nucleic Acids Res.* **50**:5251 (2022) Table 1A “1998 Model”, doi:10.1093/nar/gkac261. Initiation, AU-end, and symmetry terms are excluded. ΔS is recovered from each stack’s published ΔH and ΔG°37 at 310.15 K: ΔS = (ΔH − ΔG)/T. |
| DNA | stem stack | −8.33 w | −22.28 w | SantaLucia & Hicks, *Annu. Rev. Biophys. Biomol. Struct.* **33**:415 (2004), doi:10.1146/annurev.biophys.32.110601.141800, Table 1 (1 M NaCl). Unweighted mean of the 10 WC **propagation** stacks (AA/TT … GG/CC). Initiation, terminal AT, and symmetry are excluded. |
| PROTEIN_HELIX | backbone H-bond | −1.3 w | *gap* (0) | Scholtz, Marqusee, Baldwin et al., *PNAS* **88**:2854 (1991), doi:10.1073/pnas.88.7.2854. Calorimetric helix **unfolding** ΔH ≈ +1.3 kcal mol⁻¹ residue⁻¹ (best estimate; lower limit 0.9 if ΔCp = 0). The formation increment is the negative. Companion helix–coil CD fit: Scholtz, Qian, Baldwin, *Biopolymers* **31**:1463 (1991), doi:10.1002/bip.360311304 (ΔH° ≈ −0.955 kcal mol⁻¹ residue⁻¹). Per-residue ΔS is **not** filled — there is no calorimetric consensus independent of ΔCp assumptions; do not invent ΔS from *s*. |

RNA mean ≠ DNA mean (required). Sequence-specific stack lookup is future work;
the default is the unweighted propagation-stack mean, not a single invented
“typical stack”.

### Experimental gaps (`TODO(burgundy)` cite slots)

These slots are `{0,0}` in `default_pose_thermo_rewrite_config()`. They are
**labelled experimental**, not draft numbers. Burgundy is harvesting literature;
align the harvested set here and in `LIB/NATURaL/PoseLocalThermoRewrite.{h,cpp}`.

| Class | Contact | Why unset | Canonical reference (wrong quantity or wrong scale — do not paste as the increment) |
|-------|---------|-----------|-------------------------------------------------------------------------------------|
| RNA | loop H-bond | No calorimetric consensus for a **ligand–loop** H-bond ΔH/ΔS | Mathews, Sabina, Zuker, Turner, *J. Mol. Biol.* **288**:911 (1999), doi:10.1006/jmbi.1999.2700 — hairpin/internal-loop **initiation** ΔG, not this contact |
| DNA | loop H-bond | Same mismatch of physical quantity | SantaLucia 2004 hairpin initiation tables (doi:10.1146/annurev.biophys.32.110601.141800) |
| DNA | intercalation | Ligand-specific; no generic calorimetric mean | Leave cite slot for a harvested intercalator review (e.g. Chaires) once a consensus increment exists |
| PROTEIN_HELIX | packing | No per-contact calorimetric consensus distinct from the backbone H-bond term | `TODO(burgundy)` |
| PROTEIN_HELIX | clash | No calorimetric consensus | `TODO(burgundy)` |
| PROTEIN_SHEET | bridge H-bond | No per-bridge consensus | Maynard, Sharman, Searle, *J. Am. Chem. Soc.* **120**:1996 (1998), doi:10.1021/ja9726769 is **whole 16-residue β-hairpin** folding (endothermic / entropy-driven in water, ΔH ≈ +7 kJ mol⁻¹ for the peptide). Wrong scale and often opposite sign versus a per-bridge increment — **do not reuse** |
| PROTEIN_SHEET | hydrophobic pack | No per-contact consensus | `TODO(burgundy)` |
| PROTEIN_SHEET | shear | No calorimetric consensus | `TODO(burgundy)` |

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

## Literature (DOIs)

- Xia, SantaLucia, Turner et al. (1998) *Biochemistry* **37**:14719, doi:10.1021/bi9809425
- Zuber, Schroeder, Kennedy, Turner (2022) *Nucleic Acids Res.* **50**:5251, doi:10.1093/nar/gkac261
- SantaLucia & Hicks (2004) *Annu. Rev. Biophys. Biomol. Struct.* **33**:415, doi:10.1146/annurev.biophys.32.110601.141800
- Scholtz, Marqusee, Baldwin et al. (1991) *PNAS* **88**:2854, doi:10.1073/pnas.88.7.2854
- Scholtz, Qian, Baldwin (1991) *Biopolymers* **31**:1463, doi:10.1002/bip.360311304
- Mathews, Sabina, Zuker, Turner (1999) *J. Mol. Biol.* **288**:911, doi:10.1006/jmbi.1999.2700 — loop **initiation** (not a ligand–loop increment)
- Maynard, Sharman, Searle (1998) *J. Am. Chem. Soc.* **120**:1996, doi:10.1021/ja9726769 — whole-hairpin folding (not a per-bridge increment)
