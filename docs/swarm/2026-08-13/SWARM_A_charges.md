# LANE A — PARTIAL CHARGES (electrostatics has no input data)
BRANCH: `lane/a-partial-charges`      SEAT: Codex / ChatGPT Plus
DEPENDS ON: nothing.   BLOCKS: lane B (see §OWNERSHIP).   MACHINE TIME: zero.

## THE DEFECT, MEASURED
Electrostatics is fully implemented and completely inert.
  LIB/vcfunction.cpp:692   if (FA->use_elec && qA != 0.0 && qB != 0.0) { ... }
  LIB/vcfunction.cpp:697   E_elec = KCOULOMB * qA * qB / (FA->dielectric * dist * dist)
  LIB/vcfunction.cpp:726-731  reads atoms[i].resp_charge when has_resp, else atoms[i].charge
It never fires because THERE ARE NO PARTIAL CHARGES:
  LIB/SdfReader.cpp:193    a.charge = 0.0f;                  <- zeroed
  LIB/SdfReader.cpp:195    only reads "M  CHG" = FORMAL charge (0 for a neutral ligand)
  LIB/CifReader.cpp:338    receptor gets pdbx_formal_charge only
MEASURED: CF.elec is exactly 0.0 on all 35,873 frozen poses. Zero. Not small — zero.
CONSEQUENCE: setting scoring.electrostatics_enabled=true is a NO-OP. Never ship that alone.

## OWNERSHIP — read this before editing
YOU OWN:
    LIB/SdfReader.cpp
    LIB/Mol2Reader.cpp        (already propagates MOL2 partial charges at :319 — reference)
    LIB/CifReader.cpp
    a NEW file for the receptor charge template (e.g. LIB/charge_template.h/.cpp)
    LIB/config_parser.cpp  **lines 116-125 ONLY** (the use_elec block)
YOU CO-OWN, WITH CARE:
    LIB/vcfunction.cpp lines 688-742 — this window contains BOTH the elec blocks (692-698,
    726-731) AND lane B's metal block (721-739), interleaved in ONE loop. Metal sits
    BETWEEN the two elec blocks.
FORBIDDEN:
    LIB/config_parser.cpp lines 108-115 (lane B's metal knobs — five lines from yours)
    LIB/metal_coordination.h, LIB/DatasetRunner.cpp, LIB/gaboom.cpp, LIB/read_input.cpp
CRITICAL COUPLING: at vcfunction.cpp:727 the METAL block reads `FA->use_elec` and the same
charge fields you are populating. So **lane B's behaviour changes when your charges land.**
You land FIRST; lane B rebases onto you. Coordinate by PR, not by editing each other's lines.

## OUT OF SCOPE — do not touch GIST
GIST desolvation is HARD-DISABLED on purpose (config_parser.cpp:132-145, audit 2026-07-17):
`FA->use_gist = 0; FA->gist_weight = 0.0; FA->gist_evaluator = nullptr;` with the comment
"re-enable only behind a new validated gate + tests; never via gist_enabled alone until that
repair lands." That is a separate, pre-existing repair. **Do not enable it. Do not touch it.**
Earlier swarm drafts said "charges + desolvation" — desolvation is withdrawn from this lane.

## THE WORK
1. Ligand partial charges at load. Gasteiger is adequate and dependency-light; use AM1-BCC
   only if a toolchain already exists in-tree (check before adding a dependency — this
   project values dependency minimalism). Populate atoms[].charge for the ligand.
2. Receptor partial charges: a residue+atom-name template table for the 20 standard residues
   (AMBER or CHARMM per-atom-type charges) is sufficient and needs no new dependency.
3. Resolve a real inconsistency at vcfunction.cpp:695-697 before shipping: the comment says
   `E_elec = 332.0637*qA*qB/(eps*r)` with "distance-dependent dielectric: eps = dielectric*r",
   but the CODE divides by `dielectric * dist * dist`. Decide which is intended, state your
   reasoning in the PR, and do not silently pick one.
4. Gate it: `FLEXAIDDS_PARTIAL_CHARGES` (default OFF) plus the existing
   scoring.electrostatics_enabled. Both must be ON for anything to change.

## ACCEPTANCE GATES (all offline, no campaign)
  G1. Load one Astex ligand: assert sum(|q|) > 0 AND net charge == the formal charge.
  G2. probe_cf on ONE campaign pose **with --config** and both flags ON: CF.elec != 0.
      Also run with flags OFF and confirm CF.elec == 0 and CF_total is UNCHANGED (R5).
  G3. Rescore ~200 frozen poses across 5 targets; report Spearman(CF.elec, rmsd_spyrmsd).
      **A term with the WRONG SIGN is worse than a dead one.** Check the sign. If elec
      correlates positively with RMSD (worse poses score better), say so and stop.
  G4. score_canonical.py --frozen with the flag OFF must still print 31.0% / 48.8%.

## WHAT TO AVOID
  * Do NOT claim a success-rate gain. Electrostatics changes the SEARCH landscape; only a
    re-dock measures that, and only Claude Science launches one.
  * Do NOT add a heavy cheminformatics dependency to get charges. State the cost first.
  * Do NOT "fix" the ~200x probe_cf magnitude discrepancy by patching probe_cf — it is the
    missing --config, already diagnosed.
  * Do NOT touch the metal block even though it is 20 lines away.
