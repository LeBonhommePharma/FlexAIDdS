# Audit: dea70ea88 — Add: no-seed cognate docking stack + DirectLigandIC geometry

## Summary (2–4 sentences)
Commit `dea70ea885490acd9b8bf7d5acdc48d8b7d11edc` rewires **cognate / known-site redocking** into a genuine **no-seed protocol**: DatasetRunner empties `reference_ligand.file`, forces `seed_fraction=0` / `pose_seed_enabled=false`, blinds ligand orientation (and recenters on site centroid for defined-cleft / oracle when a site file exists), enables coarse LDS pocket scan, and marks ORACLE_CEILING as non-native-seeded for claim eligibility. It also introduces **DirectLigandIC** (topology-chosen GPA + bridge rotatable bonds for Mol2/SDF), hardens emission score↔geometry consistency (exact re-score before cluster/thermo + REMARK CF from emitted pose), switches default PoseBusters claim backend to **upstream `bust` CLI**, and expands retained pose scan to 50. Scientifically this is a large, mostly-correct fairness upgrade with **real ranking-path changes** (cluster ACF when T>0, exact re-score, cofactor radii, coarse-init ranking) and a **stale mode-contract** that still documents ORACLE_CEILING as seed_elitism ON / blinding OFF.

## Severity: HIGH

## Findings

### F1. ORACLE_CEILING semantics flipped; headers still claim crystal-seed ceiling (HIGH)
- Evidence: Runtime JSON emission now always writes:
  - `"file": ""`, `"seed_fraction": 0.0`, `"pose_seed_enabled": false`, `hetatm_fallback: false`
  - `coarse_init.enabled` true for AUTONOMOUS | DEFINED_CLEFT_REDOCK | **ORACLE_CEILING**
  - boom injection disabled for the same trio
  - pose blinding active for all publishable modes; `oracle_direct_active` only for legacy **UNSET**
  - `result.native_pose_seeded` false for ORACLE_CEILING; `protocol_claim_eligible` gated on that
- But `LIB/DatasetRunner.h` still documents:
  ```text
  ORACLE_CEILING → seed_elitism ON, blinding OFF; ceiling measurement with crystal IC.
  ```
  and `LIB/benchmark_datasets.cpp` help still prints `oracle-ceiling  seed_elitism=ON, blinding=OFF (ceiling)`.
- Why it matters for science/repro: Historical “oracle-ceiling” runs with **~90% native IC inheritance** are **not comparable** to post-commit “oracle-ceiling” (known-site, no pose seed). Operators / meta aggregators that trust the enum name will silently mix ceiling and no-seed series. Violates AGENTS.md “inspect first, claim never” for protocol labels.
- Fix recommendation: Rename the new protocol (e.g. `known-site-redock` / keep `DEFINED_CLEFT_REDOCK` as the claim path) **or** update all contracts (header, CLI help, receipt fields, campaign docs) to “ORACLE_CEILING = known-site no-seed.” Pin a `protocol_version` / receipt field distinguishing pre- vs post-`dea70ea88` oracle. Never compare RMSD success rates across that boundary without labeling.

### F2. “Pose-blind” is orientation/centroid-blind, not conformation-blind (HIGH–MEDIUM)
- Evidence: `write_blinded_ligand` applies a deterministic rigid rotation about the heavy-atom centroid and optionally **translates** the centroid to a site/cleft centroid (`target_center`). Bond lengths, angles, and **crystal torsions are preserved** in the input SDF/MOL2.
- DirectLigandIC then builds IC from that geometry (`buildic` after GPA/tree), so gen-0 torsional IC still encodes the **bound conformer** until the GA mutates dihedrals.
- Why it matters: Standard for **cognate self-docking** (many public benchmarks start from bound conformation and scramble pose), but the commit message / comments overclaim “pose-blind” / “crystal coordinates must NOT enter the GA.” Crystal **internal geometry** still enters as the IC reference. For high-DoF ligands this is a strong prior vs fully randomized torsion starts; for rigid ligands impact is mostly orientation.
- Fix recommendation: Document explicitly as **cognate conformation-preserving, pose-scrambled redock**. If a fully no-structure claim is required, add optional torsion scrambling (with seed) and report `conformation_seeded=true/false` in RUN_RECEIPT. Keep current behavior as default for Astex-style cognate redock.

### F3. Cluster ranking formula changes when T>0 (HIGH ranking impact)
- Evidence: `LIB/cluster.cpp` removes global partition function + Shannon-style per-member term:
  - **Before:** \(A_j = \sum_i\bigl(P_i E_i + T P_i\log P_i\bigr)\) with \(P_i\) from **global** \(Z=\sum_k e^{-\beta(E_k-E_0)}\)
  - **After:** cluster-local log-sum-exp free energy  
    \(A_j = E_{\min}^{(j)} - \frac{1}{\beta}\ln Z_{\mathrm{local}}^{(j)}\)  
    (or plain `app_evalue` when T=0 / β≤0)
- QuickSort still orders emission by ACF when T>0 (“Classic FlexAID contract”). Clus_GAPOP correctly **stops** being swapped with cluster ranks (GAPOP is chromosome→rep map keyed by TOP chromosome indices — swapping it was a latent membership bug).
- Why it matters: Soft-β / TEMPER>0 campaigns (e.g. entropy arms) **change rank-0 and emission order** even with identical GA snapshots. Not a pure “no-seed” change — it redefines basin scoring. Language “after considering cluster’s entropy” in the sort comment is now misleading (local Helmholtz-like F, not the old global Shannon mix).
- Fix recommendation: Version the ACF policy in receipts (`cluster_score=local_logsumexp_v1`). Add a unit/sim test with synthetic multi-basin chromosomes proving order vs old formula. Keep T=0 path byte-stable for pure CF campaigns where possible.

### F4. Exact-pose re-score before clustering + emitted REMARK CF (HIGH ranking / audit)
- Evidence: `LIB/top.cpp` re-evaluates every retained chromosome with `eval_chromosome` / `get_cf_evalue`, **always overwrites** `evalue`/`app_evalue` (including clash penalties — no HBOND clash guard retention), then `QuickSort`s. Post-GA thermo samples use the new scores. `cluster.cpp` sets `REMARK CF=` to **emitted** `get_cf_evalue`, and adds `CF.search`, `CF.pose_score_delta`, `CF.pose_score_consistent`.
- DatasetRunner re-reads those REMARKs, recomputes RMSD against the elected artifact SHA, and requires `score_pose_consistent` for `claim_ready`.
- Why it matters: Correct scientific practice (score the coordinates you emit). Previously REMARK CF could report search evalue while geometry re-scored differently; clashy poses could keep soft search scores. Ranking of snapshot chromosomes and cluster representatives can change when search vs exact diverge (OMP races, GPA rebuild order, hbond rank path).
- Fix recommendation: Keep. Log aggregate `Exact-pose score audit: inconsistent=…` into RUN_RECEIPT. Do **not** claim ΔG from the still unitless/soft-β post-GA “Free energy F” print (pre-existing; still present).

### F5. DirectLigandIC geometry: sound core + residual DoF gaps (MEDIUM)
- Evidence: New `LIB/DirectLigandIC.h`:
  - Rotatable = single bond + heavy deg≥2 both ends + **bridge** (acyclic) — rings excluded correctly.
  - GPA triad chosen by rigid-edge preference + frame quality (cross product / edge lengths), not file atom order.
  - BFS reconstruction tree sets `rec[0..2]`; `configure_rotatable_bonds` wires `ligand.bond[]` / `fdih` and cyclic shift chain for multi-control atoms.
  - Mol2/SDF readers replace first-three-atom GPA + ad-hoc BFS with this shared path.
  - Test `TopologyDerivedFrameAndTorsionPreserveLocalGeometry` asserts bonded GPA, rebuild preserves bonds/angles after `buildlist`/`buildcc`, and torsion degrees of freedom.
  - `buildlist.cpp` drops the rigid-body span bypass; relies on acyclic IC tree (fallback warns on cycles).
- Residual risks:
  1. **Amide / partial-double bonds** often appear as order 1 in SDF → counted rotatable (industry-common overcount; expands DoF and GA budget via dihedral scale).
  2. Torsions whose control child is **GPA1/GPA2** are skipped (`child == tree.gpa[1|2] continue`) — may drop edge torsions.
  3. Bridge BFS is O(E) per edge → O(E²) perception (OK for ligands; not a science bug).
  4. Mol2 heavy test uses `sybyl[0] != 'H'` (fine for `H`/`H.spc`; odd types need monitoring).
- Why it matters: Wrong GPA/rec[] historically broke rigid-body rebuild and produced non-physical geometry. This is the right fix direction; residual DoF miscounts change search budget and can change success rates without being “scoring bugs.”
- Fix recommendation: Keep geometry tests as release gates. Add amide / peptide-bond exclusion and a golden IC round-trip RMSD test on real Astex ligands (buildic→mutate dihedral→buildcc→RMSD). Log `[DIRECT-IC] fdih=N` into claim receipts for DoF audit.

### F6. Coarse-init LDS + apparent-CF ranking (MEDIUM ranking / init impact)
- Evidence: `coarse_init.cpp`:
  - Halton-like radical inverse per gene dimension + per-restart random shift; **cos-θ** sampling for polar angle (gene 1, typ==1) for isotropic orientation.
  - Ranks candidates by `get_apparent_cf_evalue` not bare `cf.com` (avoids deep-penetration attractors).
  - Restores atom/residue/optres/VC state after each trial and after the full scan (reproducible gen-0).
  - Without reflig nearest grids, covers full cleft via **voxel representatives** (not first 50 storage-order points).
  - Default `n_orientations` 16 → **64**; DatasetRunner enables coarse_init for oracle-ceiling too.
- Why it matters: Replaces crystal IC seeds with a **CF-scored pocket prior**. Fair for known-site protocols if the cleft is the only prior; changes gen-0 distribution vs old random/cleft-Gaussian. Removing gaboom “cleft-biased GPA0” Box-Muller fallback means non-MIF, non-seed paths no longer densify near grid index 0 — DatasetRunner sets `mif_enabled: true`, so production cognate path uses **MIF Boltzmann gene[0]** + coarse seeds.
- Fix recommendation: Treat coarse_init + MIF as part of the protocol fingerprint. Receipt should record `coarse_init.n_orientations`, `mif_enabled`, and binary SHA. Optional ablation: coarse_init off to measure prior contribution (not for claim parity).

### F7. Residual native-channel / centroid leakage edges (MEDIUM)
- Evidence:
  1. **ORACLE_CEILING without readable site/cleft centroid** still blinds orientation but keeps the **native heavy-atom centroid** (`target_center` null → rotate in place). DEFINED_CLEFT fails closed; oracle only logs `random blind (no site centroid for recenter)`.
  2. `FLEXAIDDS_SCORE_NATIVE=1` is still prefixed for ORACLE_CEILING (and UNSET) when an RMSD reference exists. Implementation is diagnostic-only (`score_native_pose` restores coordinates; GA continues) — **not a ranking leak**, but with blinded input the “native” CF is scoring the **processed/blinded** frame, so the label is misleading.
  3. `top.cpp` gates reflig nearest-grid seeding on `reflig_file[0] || reflig_hetatm_fallback`; DatasetRunner sets both empty/false — good fail-closed against HETATM crystal inheritance.
- Why it matters: Cognate fairness requires no native **pose** inheritance. Native **site** is allowed by design. Oracle-without-site is a weaker protocol than the commit narrative.
- Fix recommendation: Fail closed for ORACLE_CEILING without site centroid (same as defined-cleft), or rename that path. Only inject SCORE_NATIVE under explicit diag env for all claim modes. Rename diagnostic tag if input is blinded.

### F8. Receptor cofactor typing / radii fix (MEDIUM ranking when HETATMs retained)
- Evidence: `read_coor.cpp` maps element → **canonical VCT type indices** (C=3, N=11, … Fe=37) and stores `atom.element`. `assign_radii_types` skips only `FA->resligand`, not all `type==1` residues — cofactors/ions no longer keep radius 0.
- Test: `RetainedCofactorsHaveCanonicalTypesAndNonzeroRadii`.
- Why it matters: Previously invisible cofactors under-scored steric/contact terms near metals/HEM — silent scoring error on real Astex/PDB complexes. This **does** change CF landscapes vs parent when such atoms are present (correctly).
- Fix recommendation: Keep. Note in campaign diffs when comparing pre/post radii.

### F9. Intermolecular hard clash (PoseBusters radii) default-off (LOW–MEDIUM)
- Evidence: `intermolecular_clash_ratio` default **0.0** (disabled). When >0, protein–direct-ligand pairs (`atom.number >= 90000`) that violate relative RDKit/PoseBusters vdW cutoff set `clash_value = CLASH_THRESHOLD`. Soft-wall path otherwise unchanged. Unit test for boundary math.
- Why it matters: Optional alignment of search exclusion with PoseBusters intermolecular distance — good for claim parity **if** enabled and documented. Default 0 → no ranking change for current DatasetRunner JSON (unless config sets it).
- Fix recommendation: If claim campaigns enable it, pin the ratio in ARM_SPEC / receipt; never silently flip default mid-campaign.

### F10. PoseBusters default → upstream `bust` CLI (MEDIUM claim-gate)
- Evidence: `PoseBust/Engine.cpp` `resolve_backend_from_env()` default `Backend::BustCli` (was Native). DatasetRunner treats NativePoseQC as parity diagnostic; `claim_ready` requires `pb_backend == "bust_cli"`, hash equality of elected pose across RMSD/PB/tENCoM, and non-empty `posebusters_input_sha256`.
- Why it matters: Aligns claims with official PoseBusters — correct per AGENTS.md (RMSD∧PB). Environments without `bust` installed will fail PB closed (good) but can look like docking regressions if misread as engine failures.
- Fix recommendation: Document required `FLEXAIDDS_POSEBUSTERS_BIN` / PATH in claim runbooks. Keep NativePoseQC dual log for debug.

### F11. Vcontacts thread_local indexing + OpenMP race fix (MEDIUM correctness)
- Evidence: `indexed` / `prev_box` become `thread_local`; box allocation path simplified when `!FA->vindex`.
- Why it matters: Shared mutable contact boxes under parallel eval are a classic non-reproducible scoring race. Fix improves determinism; may change multi-thread results vs racy parent (reproducibility win).
- Fix recommendation: Keep. Prefer serial exact re-score path (already added) as ground truth for emission.

### F12. Stale GAPOP-in-sort + thermo language (LOW–INFO)
- Removing GAPOP from `QuickSort_Clusters` is a **correctness fix** for membership maps, not a ranking sort-key change.
- Post-GA print still says “Free energy F … Entropy S” on CF samples at `FA->temperature` — soft-β / proxy; not true ΔG (AGENTS.md). Commit does not fix naming.
- `ligand_tencom_pose` tool + eigenvalue rigid-mode cutoff (`max_λ·1e-8`) improve H_vib hygiene for validators — good; still Shannon spectrum entropy, not full tENCoM+solvent thermo.

### F13. Tests adequacy (MEDIUM gap)
- Present: SDF topology/geometry round-trip; cofactor radii; soft_wall relative vdW; extractor V5 stamp; intermolecular_clash default in JSON config.
- Missing for this commit’s claim surface:
  - DatasetRunner no-seed contract test (emitted JSON: seed_fraction 0, empty reflig, coarse_init on, blinding fail-closed for defined-cleft).
  - ORACLE_CEILING mode-doc vs runtime assertion.
  - Cluster local-logsumexp ordering golden test.
  - End-to-end seed_echo / native_pose_seeded flags under blinded cognate.
  - Coarse-init restore_scan_state isolation test.
- Fix recommendation: Add pure unit tests for config emission + cluster ACF; one short integration smoke with tiny pop/gen.

### F14. Security / hygiene (INFO)
- No secrets. Paths resolved from entry / env. Sibling binary auto-set of `FLEXAIDDS_BINARY` is build-tree local. Apache-2.0 only. Absolute user paths not introduced in shared scripts by this commit.

## Ranking/scoring impact: YES

| Subsystem | Default claim path impact |
|-----------|---------------------------|
| No crystal pose seed + blinding + site recenter | Search start distribution (not CF formula) |
| Coarse-init LDS + apparent CF + n_orient=64 | Gen-0 seeds |
| MIF gene[0] (still on) + removal of cleft-Gaussian fallback | Gen-0 translation prior |
| Exact re-score + QuickSort snapshot | Chromosome order into clustering |
| Cluster ACF local log-sum-exp (T>0) | **Emission rank-0 / order** |
| REMARK CF = emitted pose CF | Downstream election CF columns |
| Cofactor VCT types/radii | CF when HETATMs retained |
| `intermolecular_clash_ratio` default 0 | None unless enabled |
| PoseBust default bust_cli | Claim gate, not CF ranking |
| max_results / pose limit 50 | Election pool size |

Per AGENTS.md: ranking path changed without an isolated feature flag for the cluster ACF rewrite — treat as a **new protocol revision**, not a drop-in engine patch.

## Reproducibility impact: YES

**Positive:** empty reflig, fail-closed blinding (defined-cleft), elected-pose SHA chain for RMSD/PB/tENCoM, score_pose_consistent gate, sibling binary pin, V5 ligand cache stamp, thread_local Vcontacts, coarse-init state restore, validator_provenance.json.

**Negative:** ORACLE_CEILING name/docs disagree with runtime; conformation prior not labeled; T>0 ranking formula change without version pin; SCORE_NATIVE still auto-injected for oracle; bust CLI dependency for claims.

## Tests adequate: PARTIAL

Geometry and cofactor unit tests are strong for DirectLigandIC/read_coor. Protocol/ranking changes (no-seed JSON, cluster ACF, blinding recenter) lack dedicated automated gates in this commit.

## No-seed fairness verdict

| Check | Status |
|-------|--------|
| Crystal pose IC flood (`seed_fraction` / `pose_seed_enabled`) | **Closed** (0 / false) |
| `reference_ligand.file` empty + hetatm_fallback false | **Closed** |
| Seed elitism override for publishable modes | **Off** (0) for oracle + defined-cleft + autonomous |
| Orientation scramble | **Yes** (deterministic PDB-id seed) |
| Native centroid leakage | **Removed** when site/cleft centroid available; residual for oracle-without-site |
| Crystal torsion / bound conformer prior | **Still present** (cognate-standard; must label) |
| Known-site prior (cleft / MIF / coarse CF) | **Intentional** for cognate stack |
| Claim eligibility vs native_pose_seeded | **Wired** |
| seed_echo still fails success_rmsd | **Yes** |

## DirectLigandIC geometry verdict

| Check | Status |
|-------|--------|
| GPA topology-derived, not record order | **Yes** (tested) |
| IC tree acyclic BFS | **Yes** |
| Rotatable = bridge single bonds | **Yes** (amide overcount residual) |
| Local geometry preserved under rebuild | **Yes** (bond/angle test) |
| fdih wired for flexible SDF (was rigid risk) | **Yes** (shared path) |
| buildlist cycle infinite-loop | **Mitigated** (warn + deterministic append) |

## Verdict: MERGE_WITH_FIX

Ship as the foundation of the **no-seed cognate stack**, but do not treat historical oracle-ceiling numbers as comparable and do not publish T>0 rankings without acknowledging the new local free-energy cluster score.

**Required before claim CSVs labeled “no-seed cognate / oracle-ceiling”:**
1. Fix F1 mode contracts (docs/CLI/header) to match runtime — or rename mode.
2. Label F2 conformation-preserving redock in protocol/receipts.
3. Version F3 cluster score in RUN_RECEIPT; regression-test T>0 ordering.
4. Fail-closed or document F7 oracle-without-centroid.
5. Add F13 protocol emission tests.

**Accept as-is for engineering merge** if the above are tracked and campaigns use a single post-`dea70ea88` binary SHA with explicit protocol notes.
