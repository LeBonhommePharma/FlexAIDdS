# LANE B — METAL COORDINATION: enabled in config, firing on 2.2% of poses
BRANCH: `lane/b-metal-coord`      SEAT: Grok / Super Grok
DEPENDS ON: lane A for the vcfunction window (see §OWNERSHIP).   MACHINE TIME: zero.

## THE DEFECT, MEASURED
The campaign's dock_config.json carries `"metal_coord_enabled": true`.
CF.metal is nonzero on only **783 of 35,873 frozen poses (2.2%)** — on a benchmark full of
zinc metalloenzymes (carbonic anhydrases, MMPs). 1JD0 is a carbonic anhydrase whose
sulfonamide coordinates the catalytic zinc directly; its elected pose carries CF.metal = 0.0.
Enabled-and-not-firing is a DEFECT, not a setting. The compute site EXISTS and is reached:
  LIB/vcfunction.cpp:721   if (FA->use_metal_coord) {
  LIB/vcfunction.cpp:736   E_mc = metal_coord::compute_metal_coord_energy(..., mc_weight, sigma)
  LIB/vcfunction.cpp:739   cfs->metal_coord += E_mc;
So the cause is almost certainly UPSTREAM of the compute site, not in it. Start upstream.

## OWNERSHIP
YOU OWN:
    LIB/metal_coordination.h
    LIB/config_parser.cpp **lines 108-115 ONLY** (metal_coord knobs)
    the receptor-construction path in LIB/DatasetRunner.cpp around the cofactor whitelist
      (cofactor_blacklist() at :2008, call sites :3249 and :3312, write_receptor_without_ligand's
       keep_catalytic) — coordinate with lane D, which owns DatasetRunner election/emission
       regions (~1005-1024, ~1266, ~6433-6466). Yours are 2000-3400; textually disjoint.
    atom typing where metals are classified (find it; likely LIB/ProcessLigand/SybylTyper.cpp
      or read_input's typing path)
FORBIDDEN UNTIL LANE A LANDS:
    LIB/vcfunction.cpp lines 688-742. This window is shared: elec at 692-698 and 726-731,
    your metal block at 721-739, interleaved in ONE loop.
CRITICAL COUPLING: vcfunction.cpp:727 — YOUR metal block reads `FA->use_elec` and the charge
fields lane A is populating. **Your term's behaviour will change when lane A lands.** Do your
diagnosis first (it is upstream anyway), then rebase onto lane A before patching vcfunction.
FORBIDDEN OUTRIGHT:
    LIB/SdfReader.cpp, LIB/CifReader.cpp, config_parser.cpp:116-145, LIB/gaboom.cpp

## THE WORK — diagnosis first, in this order
1. Which targets even HAVE a coordinating metal? Compute, from the crystal structures: how
   many of the 84 have a metal ion within 3.0 A of the crystal ligand. **That set is the
   denominator for this defect** — if it is 12 targets, a perfect fix is worth at most 12.
2. For 3 of those, check whether the metal atom SURVIVES into the receptor the engine loads.
   Receptor construction uses a POSITIVE whitelist (`keep_catalytic`) plus a
   `cofactor_blacklist()`. A stripped ZN cannot be coordinated. This is the top hypothesis.
3. If the metal is present, check whether it is TYPED as a metal by the atom typer. A metal
   typed as a generic heavy atom will never enter the metal branch.
4. Only then look at geometry/cutoffs (`metal_coord_sigma` default 0.45,
   `metal_coord_cn_weight` default 0.5) and the compute function itself.
5. Gate any change: `FLEXAIDDS_METAL_FIX` (default OFF).

## ACCEPTANCE GATES (all offline)
  G1. Print the count and list of targets with a metal within 3.0 A of the crystal ligand.
      State that count as the maximum possible value of this lane.
  G2. For 3 of them, show whether the metal exists in the receptor PDB the engine loads —
      with the grep/parse command and its output, not a source-reading argument (R7).
  G3. Root cause named as file:line PLUS a runtime observation confirming it.
  G4. After the fix, probe_cf **with --config** on a known coordinating pose gives
      CF.metal != 0; with the gate OFF it gives CF.metal == 0 and an unchanged CF_total.
  G5. score_canonical.py --frozen with the gate OFF still prints 31.0% / 48.8%.

## WHAT TO AVOID
  * Do NOT start in vcfunction.cpp. The compute site is reached; the input is missing.
  * Do NOT touch lane A's lines even though they interleave with yours.
  * Do NOT strip or add cofactors globally to "help" — 1TW6 shipped its own ligand as ATOM
    records and 1T9B was missing its FAD; receptor content is audited and fragile.
  * Do NOT claim a success-rate gain without a re-dock, which you may not launch.
  * Do NOT assume metal_coord_enabled=true means the term ran. It was true and it did not.
