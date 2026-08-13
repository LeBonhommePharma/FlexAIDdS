# LANE E — AUTOFLEX IS DEFAULT-ON AND SILENTLY INERT (highest upside)
BRANCH: `lane/e-autoflex`      SEAT: Claude Code (Max), primary seat
DEPENDS ON: lane D (rebase onto it).   MERGES LAST.   MACHINE TIME: zero to diagnose.

## THE DEFECT, VERIFIED THREE WAYS
  LIB/top.cpp:458   FA->autoflex_enabled = 1;   // auto-flex key binding residues BY DEFAULT
  LIB/top.cpp:459   FA->autoflex_max     = 5;
  LIB/read_input.cpp:635   if (FA->autoflex_enabled && FA->mif_energies && FA->mif_count > 0)
Evidence it never ran in the completed campaign:
  1. "AUTOFLEX" appears 0 times in claim.log and 0 times in ANY per-restart log
  2. "MIF" appears 0 times in claim.log
  3. every pose carries exactly 1 "REMARK optimizable residue" (the ligand alone)
=> zero flexible side chains on all 84 targets. A DEFAULT-ON feature has been silently inert,
and every campaign this project has ever run was rigid-receptor.

TRAP: dock_config carries `"receptor_rotamer_prep": true` and it is NOT runtime flexibility.
It is STATIC file prep (DatasetRunner.cpp:6292, `prep_receptor_rotamers()` writing
<PDB>_prepped.pdb). Do not mistake it for side-chain sampling. The runtime paths are:
  * FLEXSC directive — format `FLEXSC <resnum> <chain> <resname>`, parsed at
    LIB/read_flexscfile.cpp:43, collected at read_input.cpp:184
  * autoflex — gated on mif_count > 0 as above
DatasetRunner emits NEITHER. That is the hole.

## WHY THIS IS THE HIGHEST-UPSIDE LANE
The pool ceiling is 48.8% and the bar is 45.2%, so selection cannot get there. Receptor
flexibility is the mechanism most likely to RAISE the ceiling:
  * 62% of sub-2 A poses carry wall penalty wal > 20 (median 29.4, max 367)
  * 1Z95: true 1.75 A poses score CF -9 while a 7.5 A decoy wins at CF -187
Near-native geometry is being punished for clashing with un-relaxed side chains.
Corroborating (measured 2026-08-13): 30 restarts converted 4/19 ceiling-miss targets to
sub-2 A, so the ceiling IS movable — but median gain was only +0.02 A, i.e. more sampling of
the SAME rigid landscape is not the answer. Changing the landscape is.

## OWNERSHIP
YOU OWN:
    LIB/read_input.cpp        (the autoflex gate ~635, FLEXSC collection ~184)
    LIB/read_flexscfile.cpp
    LIB/top.cpp               (the autoflex defaults ~458-459)
    LIB/DatasetRunner.cpp **~5949-6320 ONLY** (.inp emission + receptor_rotamer_prep)
    the MIF computation path, wherever it is invoked from
FORBIDDEN, in DatasetRunner.cpp:
    ~1005-1024, ~1266, ~6433-6466  election / guard / voiding      -> LANE D (merges before you)
    ~2008, ~3249, ~3312  cofactor whitelist              -> LANE B
FORBIDDEN OUTRIGHT:
    LIB/vcfunction.cpp (A/B), LIB/gaboom.cpp (C), LIB/SdfReader.cpp (A), config_parser.cpp (A/B)
You merge LAST and rebase onto lane D — you both touch DatasetRunner.cpp, 146 lines apart at
closest approach, so a clean rebase is expected but must be verified.

## THE WORK
1. Determine why mif_count is 0. Is MIF computation gated off, failing silently, or never
   invoked in the dataset-runner path? Name it file:line WITH runtime evidence (R7).
2. Give DatasetRunner a way to request flexible pocket side chains — EITHER fix the autoflex
   prerequisite, OR emit FLEXSC lines for residues within ~5 A of the binding-site centroid,
   capped at 5 to match autoflex_max. State which you chose and why.
3. Confirm the rotamer library resolves at runtime (rotlib / rotobs.lst path). A missing
   library will make this fail silently too — that is this defect's whole family.
4. Gate it: `FLEXAIDDS_FLEX_SIDECHAINS` (default OFF). NOTE: autoflex_enabled is currently
   default 1 but inert. If your fix makes it live, flipping the default ON would change every
   future campaign silently — so your gate must dominate, and with the gate OFF behaviour must
   be bit-identical to today (R5).

## ACCEPTANCE GATES
  G1. Single-target run on 1Z95 emits AUTOFLEX (or FLEXSC) lines AND more than 1
      "REMARK optimizable residue" in the pose file. Show the grep output.
  G2. On that target, report the wall-penalty distribution of sub-2 A poses before vs after.
      The prediction under test: flexible side chains REDUCE wal on near-native poses.
      If wal does not drop, say so — that falsifies the mechanism and is a valuable result.
  G3. Report the per-target POOL CEILING before vs after on 1Z95 (spyRMSD, canonical scorer).
      This is the number that matters: does flexibility let the search FIND what it could not?
  G4. With the gate OFF, score_canonical.py --frozen still prints 31.0% / 48.8%, and a
      single-target rerun is bit-identical to the parent campaign's pose for that target.
  G5. Report the RUNTIME COST multiplier per target. Flexible side chains add degrees of
      freedom; if it is 5x, a full campaign is 45 h and that changes the whole plan.

## WHAT TO AVOID
  * Do NOT run a full campaign. Hand the binary to Claude Science for a ~20-target pilot.
  * Do NOT flip autoflex_enabled's default or remove the mif gate without a dominating env gate.
  * Do NOT touch lane D's DatasetRunner regions; rebase onto D instead.
  * Do NOT confuse receptor_rotamer_prep (static) with runtime flexibility.
  * Do NOT skip G5. An unbounded runtime cost makes a working fix unusable on one laptop.
