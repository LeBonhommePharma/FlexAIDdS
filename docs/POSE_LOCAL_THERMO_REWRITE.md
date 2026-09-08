# Pose-local secondary-structure thermodynamic rewrite

**Status: experimental. Not claim-ready. Default OFF.** Filled ΔH/ΔS increments
are peer-reviewed NN / helix–coil values (cited below). Motifs without a
calorimetric consensus stay `{0,0}` experimental gaps — they are **not**
invented placeholders. This path must not be cited as Astex success, FlexADS
claim-ready ranking, or a validated binding free energy.

This document specifies the physics, class-specific tables, and API for
`LIB/NATURaL/PoseLocalThermoRewrite.{h,cpp}` and `LIB/NATURaL/NnMotifTables.h`.
It **generalizes** PR #471 `PoseHelixThermoRewrite` (RNA decision-helix, atom
coords) to RNA/DNA/helix/sheet class tables. `#471` files stay; this module is
the cited-table API.

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

- `RNA_STEM_LOOP` — Xia 1998 INN-HB WC stacks (mean default; `PoseContact::motif` selects a dimer)
- `DNA_STEM_LOOP` — SantaLucia 1998 PNAS Table 2 — **never** an alias of the RNA table
- `PROTEIN_HELIX` — Scholtz 1991 calorimetric helix-formation ΔH + Zavrtanik 2026 backbone ΔS
- `PROTEIN_SHEET` — Meier–Seelig 2008 labelled experimental midpoint on the bridge H-bond only
- `PROTEIN_OTHER` — no rewrite unless annotated **and** that table is filled

Geometry weight `w` multiplies the increment. An RNA `ContactKind` never reads the
DNA or protein tables (and conversely). Shared algebra, separate numbers.

### Filled (peer-reviewed)

| Class | Contact | ΔH (kcal mol⁻¹) | ΔS (cal mol⁻¹ K⁻¹) | Source |
|-------|---------|-----------------|---------------------|--------|
| RNA | stem stack | Xia dimer, or −10.756 w mean | recovered ΔS, or −27.90… w mean | Xia et al. *Biochemistry* **37**:14719 (1998), doi:10.1021/bi9809425. Ten unique WC **propagation** stacks (GC/CG … AA/UU) as restated in Zuber et al. *Nucleic Acids Res.* **50**:5251 (2022) Table 1A “1998 Model”, doi:10.1093/nar/gkac261. Empty `motif` → unweighted mean. Named `motif` (e.g. `GC/CG`) → that row. TSV copy: `LIB/NATURaL/data/xia_1998_rna_wc_stacks.tsv`. |
| RNA | loop (hairpin initiation) | Lu ΔH°(*n*) · w | *gap* (0) | Lu, Turner, Mathews, *Nucleic Acids Res.* **34**:4912 (2006), doi:10.1093/nar/gkl472, Table 1. Default *n*=4 → +4.8. `loop_unpaired` selects *n* (3:1.3; 4:4.8; 5:3.6; 6:−2.9; 7:1.3; 8:−2.9; ≥9:5.0). This is hairpin **initiation enthalpy**, not a ligand–loop H-bond. ΔS is not tabulated (extra length cost is treated as entropic). TSV: `LIB/NATURaL/data/lu_2006_rna_hairpin_initiation.tsv`. |
| DNA | stem stack | Table 2 dimer, or −8.36 w mean | Table 2 dimer, or −22.37 w mean | SantaLucia Jr, *PNAS* **95**:1460 (1998), doi:10.1073/pnas.95.4.1460 (PMC19045), **Table 2** in 1 M NaCl. Ten WC propagation stacks (AA/TT … GG/CC); Init G·C / Init A·T / Symmetry are stored but are not the stem-stack default. Empty `motif` → WC-stack mean. Named `motif` (e.g. `AA/TT` = −7.9 / −22.2) → that row. 2004 review is a reprint, not the table source. TSV: `LIB/NATURaL/data/santalucia_1998_pnas_table2.tsv`. |
| PROTEIN_HELIX | backbone H-bond | −1.3 w | −5.2 w | ΔH: Scholtz, Marqusee, Baldwin et al., *PNAS* **88**:2854 (1991), doi:10.1073/pnas.88.7.2854 (calorimetric helix **unfolding** +1.3 → formation −1.3). ΔS: Zavrtanik, Lah, Hadži, *Biophys. J.* **125**:305 (2026), doi:10.1016/j.bpj.2025.11.2689 (PubMed 41318999). Helix→coil ΔS_BB = 5.2 ± 0.3 cal mol⁻¹ K⁻¹ per peptide unit → formation −5.2. |
| PROTEIN_SHEET | bridge H-bond | −0.4 w | −1.0 w | **Labelled experimental midpoint.** Meier & Seelig, *J. Am. Chem. Soc.* **130**:1017 (2008), doi:10.1021/ja077231r — membrane coil⇄β, ΔH_fold ≈ −0.2 to −0.6 kcal mol⁻¹ residue⁻¹, TΔS_fold ≈ −0.1 to −0.5; mid ΔH = −0.4, ΔS ≈ −1.0 e.u. at 298 K. Caveat: Deechongkit et al., *Nature* **430**:101 (2004), doi:10.1038/nature02611 — β H-bond energetics are context-dependent. Not a universal NN sheet table. |

RNA mean ≠ DNA mean (required). `w` multiplies the matched published term. An RNA `ContactKind` never reads the DNA or protein tables (and conversely). Shared algebra, separate numbers.

### Experimental gaps (`TODO(burgundy)` cite slots)

These slots are `{0,0}` in `default_pose_thermo_rewrite_config()`. They are
**labelled experimental**, not draft numbers. Burgundy is harvesting literature;
align the harvested set here and in `LIB/NATURaL/PoseLocalThermoRewrite.{h,cpp}`.

| Class | Contact | Why unset | Canonical reference (wrong quantity or wrong scale — do not paste as the increment) |
|-------|---------|-----------|-------------------------------------------------------------------------------------|
| DNA | loop H-bond | No calorimetric consensus for a **ligand–loop** H-bond ΔH/ΔS | SantaLucia 1998/2004 hairpin initiation tables are loop-initiation ΔG, not this contact |
| DNA | intercalation | Ligand-specific; no generic calorimetric mean | Leave cite slot for a harvested intercalator review (e.g. Chaires) once a consensus increment exists |
| PROTEIN_HELIX | packing | No per-contact calorimetric consensus distinct from the backbone H-bond term | `TODO(burgundy)` |
| PROTEIN_HELIX | clash | No calorimetric consensus | `TODO(burgundy)` |
| PROTEIN_SHEET | hydrophobic pack | No per-contact consensus | `TODO(burgundy)` |
| PROTEIN_SHEET | shear | No calorimetric consensus | `TODO(burgundy)` |

Maynard, Sharman, Searle, *J. Am. Chem. Soc.* **120**:1996 (1998), doi:10.1021/ja9726769 is **whole 16-residue β-hairpin** folding — wrong scale versus a per-bridge increment; **do not reuse**.

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
- NNDB: Turner & Mathews (2010) *Nucleic Acids Res.* **38**:D280, doi:10.1093/nar/gkp892 — https://rna.urmc.rochester.edu/NNDB/
- Zuber, Schroeder, Kennedy, Turner (2022) *Nucleic Acids Res.* **50**:5251, doi:10.1093/nar/gkac261
- Lu, Turner, Mathews (2006) *Nucleic Acids Res.* **34**:4912, doi:10.1093/nar/gkl472
- SantaLucia Jr (1998) *PNAS* **95**:1460, doi:10.1073/pnas.95.4.1460 (PMC19045) — DNA Table 2
- SantaLucia & Hicks (2004) *Annu. Rev. Biophys. Biomol. Struct.* **33**:415, doi:10.1146/annurev.biophys.32.110601.141800 — review reprint, not the DNA table source
- Scholtz, Marqusee, Baldwin et al. (1991) *PNAS* **88**:2854, doi:10.1073/pnas.88.7.2854
- Zavrtanik, Lah, Hadži (2026) *Biophys. J.* **125**:305, doi:10.1016/j.bpj.2025.11.2689 (PubMed 41318999)
- Meier & Seelig (2008) *J. Am. Chem. Soc.* **130**:1017, doi:10.1021/ja077231r — labelled experimental sheet midpoint
- Deechongkit et al. (2004) *Nature* **430**:101, doi:10.1038/nature02611 — β H-bond context caveat
- Maynard, Sharman, Searle (1998) *J. Am. Chem. Soc.* **120**:1996, doi:10.1021/ja9726769 — whole-hairpin folding (not a per-bridge increment)
