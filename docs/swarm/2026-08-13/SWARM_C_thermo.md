# LANE C — THE THERMODYNAMICS PATH NEVER EXECUTED
BRANCH: `lane/c-thermo-gate`      SEAT: Claude Code (2nd seat)
DEPENDS ON: nothing.   BLOCKS: nothing.   MACHINE TIME: zero to diagnose.

## THE DEFECT, MEASURED
In campaign astex84_dG_20260809_141245:
  * ZERO "[THERMO]" lines in claim.log
  * result.csv TdS_vib, G_bind, H_vct, TdS_shannon, D_vib_thermo all 0/84 nonzero
  * thermo_n_heavy = 0/84
  * BUT eigen_n_modes and elected_H_vib ARE nonzero for 76/84 -> the eigen/vibrational
    machinery RAN. Only the TdS_vib TERM never entered the ledger or the election.
The THERMO block is around LIB/gaboom.cpp:1605-1665.

## WHY THIS MATTERS MORE THAN THE RATE
LP presented 66% (69% WITH ENTROPY) at 3Dsig. This campaign's "entropy" was the soft-beta CF
Shannon term over duplicate-inflated cluster members — NOT vibrational entropy. Therefore
**the 69%-with-entropy result has never been reproduced by this codebase**, and no recent
campaign is comparable to it. This lane is what makes that comparison possible at all. It is
about claim integrity, not about moving top-1.

## OWNERSHIP
YOU OWN:
    LIB/gaboom.cpp  (the THERMO block region ~1605-1665 and its gate)
    the thermo columns' emission path
FORBIDDEN:
    LIB/vcfunction.cpp (lanes A and B), LIB/DatasetRunner.cpp (lanes D and E),
    LIB/SdfReader.cpp, LIB/config_parser.cpp, LIB/read_input.cpp
NOTE: gaboom.cpp is 4,962 lines and no other lane owns it — you have it to yourself. Do not
touch the four GA hot-loop pragmas (3221, 3338, 3449, 4022): they are already order-
independent and static scheduling there was measured to be a NON-fix. That request is withdrawn.

## THE WORK
1. Find the gate that skipped the THERMO block on every target. Name it with a file:line AND
   show the runtime evidence that it is closed (R7). Candidates to TEST, not assume: a
   protocol flag, a build gate, an eigen prerequisite that silently failed, a weight of zero
   (dock_config carries tencom_weight = 0.0 — check whether that alone disables the block).
2. Establish whether TdS_vib was ever intended to enter the ELECTION objective, or only the
   reported ledger. These are different claims with different blast radii. Say which, with
   evidence.
3. Make it observable: emit a one-line banner stating thermo ON/OFF **and why**, so a future
   campaign cannot silently run without it. This is the highest-value part of the lane —
   the defect class here is "a default-on feature silently inert", which has now bitten this
   project three times (thermo, autoflex, metal).
4. Gate any behaviour change: `FLEXAIDDS_THERMO_FIX` (default OFF). The banner may always print.

## ACCEPTANCE GATES (all offline)
  G1. Name the gate and show the runtime evidence it is closed.
  G2. Single-target probe with thermo enabled emits [THERMO] lines and a nonzero TdS_vib.
  G3. State explicitly whether enabling it changes the ELECTION or only the REPORT.
  G4. score_canonical.py --frozen with the gate OFF still prints 31.0% / 48.8%.
  G5. The banner prints on a target where thermo is OFF, and names the reason.

## WHAT TO AVOID
  * Do NOT claim any success-rate change. Hand back for a pilot; you may not launch one.
  * Do NOT wire TdS_vib into the election in this PR even if you find you can. Enabling the
    term and changing the objective are two decisions; the second needs LP's ruling.
  * Do NOT touch the CF Shannon / soft-beta election path — that is lane D.
  * Do NOT re-derive the soft-beta T=300 finding; it is settled (|T*S|/|H| median 17x,
    elects the largest cluster in 73/76 targets).
