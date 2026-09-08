# Pose→helix ΔH/ΔS rewrite (experimental)

**Status:** experimental. Not claim-ready. Does not elect poses, does not change
Astex CF ranking, and does not write `FA->natural_deltaG` or DatasetRunner claim
metrics.

## What this is

`LIB/NATURaL/PoseHelixThermoRewrite.{h,cpp}` is the 2026 C++ glue for piece (4) of
the **2014 NATURAL** seminar:

> **NATURAL** = Native Assembly of Transcriptionally-Unified RNA And Ligand
> (LP Morency, master’s seminar R2, 2014-07-25)

That Perl implementation followed Zhao, Zhang, Chen, *Cotranscriptional folding
kinetics of ribonucleic acid secondary structures*, J. Chem. Phys. **135**, 245101
(2011), doi:10.1063/1.3671644, and added:

1. clustering
2. population inheritance across elongation
3. slower transcription at nucleotide repeats
4. **docking poses that rewrite helix ΔH/ΔS** (this module)

The 2026 DualAssembly acronym is broader: Native Assembly of co-Transcriptionally /
co-Translationally Unified Receptor–Ligand (`NATURaLDualAssembly.h`).

## Distinct Zhao 2011 papers

| Paper | Role |
|-------|------|
| Zhao, Zhang, Chen, J. Chem. Phys. **135**, 245101 (2011) doi:10.1063/1.3671644 | RNA cotranscriptional folding kinetics. 2014 NATURAL lineage. |
| Zhao et al., J. Phys. Chem. B **115**, 3987 (2011) | Ribosome master equation / **protein cotranslation ONLY**. `RibosomeElongation.h`. |

Do not cite JPCB as the RNA folding-kinetics source, and do not cite JCP as the
ribosome ODE.

## Physics (seminar)

Decision-helix patches only (`HelixSegment` list):

- Ligand H-bond / base-pair-like contacts to **loop** nucleotides → ΔS < 0
  (loop entropy reduction). Default −3.0 cal K⁻¹ per unique loop nt.
- Aromatic stacking on the **stem** → more negative ΔH (stabilizing stacks).
  Default −1.2 kcal per unique stem nt. Ring-normal angle ≤ 40° (parallel or
  antiparallel); a zero normal skips the angle test.
- Mixture over top-*k* poses with Boltzmann weights from docking **CF scores**
  (scoring proxy, not ΔG_bind).
- ΔG [kcal] = ΔH [kcal] − T · (ΔS [cal/K] / 1000).

`NucleationDetector` hairpin ΔG is sequence-only (Turner). DualAssemblyEngine
CF + Shannon instantaneous ΔG is **not** this ligand coupling. Optional StatMech
`G_natural` is not fed by these patches.

## API

Pure function (no GA):

```text
rewrite_helices_from_poses(helices, rna_nts, poses, cfg)
```

Coordinate types (`ReceptorNtCoord`, `LigandPose`) are lightweight and independent
of `atom_struct`, so unit tests do not need a full engine.

## DualAssembly hook (default OFF)

`NATURaLConfig.enable_pose_helix_rewrite` defaults to **false**.

- `DualAssemblyEngine::run()` never applies patches to CF or `final_deltaG_`.
- `DualAssemblyEngine::compute_pose_helix_rewrite(rna, poses)` is a sidecar: it
  returns empty unless the flag is on **and** the receptor is nucleic acid.
- `DualAssemblyRunner` does not call the rewrite (injected GA results currently
  carry pose RMSDs, not atom coords). Follow-up: wire when a real-GA backend
  emits snapshots.

## Not claim-ready

Do not use this module for Astex success, PoseBusters gates, or METHODS claim
contracts. Treat outputs as a helix-parameter diagnostic for RNA DualAssembly
experiments only.

## Class-table generalization

`LIB/NATURaL/PoseLocalThermoRewrite.{h,cpp}` (see `docs/POSE_LOCAL_THERMO_REWRITE.md`)
extends this 2014 ligand coupling to **RNA, DNA, protein helix, and protein sheet**
with class-specific published tables (Xia 1998, Lu/Mathews 2006 / NNDB,
SantaLucia 1998 PNAS Table 2, Scholtz 1991, Meier–Seelig 2008). Same algebra
(`ΔG = ΔH − TΔS`, geometry weight `w`, Boltzmann mixture); **different tables**.
PoseLocal does **not** use this module’s seminar defaults (−1.2 kcal/stack,
−3.0 e.u./loop H-bond). Both DualAssembly flags default **OFF**. Neither path
feeds `G_natural` or Astex claim metrics.
