# SWARM ASSIGNMENT — five lanes, one referee, one merge queue   (v2, 2026-08-13)

## Files to hand each seat (three, in this order)
    1. SWARM_COMMON_PREAMBLE.md   (identical for all five)
    2. <the seat's lane file>
    3. score_canonical.py         (or its absolute path — it must be the ONLY scorer used)

| lane | file | seat | why that seat |
|------|------|------|----------------|
| E autoflex | SWARM_E_autoflex.md | Claude Code (Max) | deep single-file tracing; longest context carry. START HERE — highest upside |
| A partial charges | SWARM_A_charges.md | Codex / ChatGPT Plus | long mechanical C++ across 3 readers, tight compile loop |
| D waste + selection | SWARM_D_selection.md | Cursor Pro | many small surgical edits, fastest edit-verify loop |
| B metal coordination | SWARM_B_metal.md | Grok / Super Grok | whole-codebase gate hunt on an unfamiliar surface |
| C thermo path | SWARM_C_thermo.md | Claude Code (2nd seat) | archaeology on a path that never ran |

## FILE COLLISION MAP — measured, not assumed
TWO pairs collide, not one:

  vcfunction.cpp lines 688-742  — lanes A and B, INTERLEAVED IN ONE LOOP
      elec   692-698 and 726-731   (lane A)
      metal  721-739               (lane B)  <- sits BETWEEN the two elec blocks
      COUPLING: line 727 the METAL block READS FA->use_elec and the charge fields lane A
      populates. So B is FUNCTIONALLY dependent on A, not merely adjacent.

  config_parser.cpp — lanes A and B, five lines apart
      108-115 metal knobs (lane B)   116-125 use_elec (lane A)

  DatasetRunner.cpp — lanes D, E and B, textually disjoint
      1005-1024 / 1266 / 6433-6466   election, spread guard, voiding   (lane D)
      2008 / 3249 / 3312   cofactor whitelist                (lane B)
      5949-6320            .inp emission, receptor prep      (lane E)
      closest approach D<->E is 146 lines

  gaboom.cpp — lane C ALONE (4,962 lines, no other lane touches it)

## MERGE ORDER — D → A → B → E (rebased on D) → C
  1. D first  — smallest, fully offline-verifiable, no rebuild needed to evaluate
  2. A next   — MUST precede B (B reads use_elec + the charge fields at vcfunction.cpp:727)
  3. B        — rebase onto A before touching the vcfunction window
  4. E        — largest blast radius; rebase onto D
  5. C        — independent; land any time
An earlier version of this file said "D → B → A → E → C". That was WRONG: it merged B before
its prerequisite. Use the order above.

## THREE GATES before any merge
  (a) the lane's own acceptance gates green
  (b) its number came from score_canonical.py
  (c) after rebuild, `score_canonical.py --frozen` still prints 26/84 = 31.0% and
      41/84 = 48.8% with the new feature OFF
A PR that moves a DEFAULT is an automatic reject. A knob documented as OFF being ON is exactly
what cost 13 points on the last campaign.

## Ranked by expected value
  E autoflex  — the only lane that can raise the ceiling by changing the landscape.
                34 of 42 ceiling-misses already place their best pose in the RIGHT pocket
                (centroid <=3 A) at the wrong geometry — several within 0.04 A of threshold.
                That is landscape-limited, which is what flexibility addresses.
  A charges   — restores a whole missing physics term; also a landscape change.
  D selection — smallest but CERTAIN, bounded at +15 targets, zero compute.
  B metal     — bounded by however many targets actually have a coordinating metal (lane B's
                first gate is to measure that count).
  C thermo    — unblocks comparability to the 69%-with-entropy claim, not the rate itself.

## Who owns what
  Claude Science : owns the box, launches every docking campaign, holds the single ledger,
                   re-verifies every lane's number against the frozen benchmark.
  Claude Code    : owns the merge queue (repo access, runs git directly).
  LP             : ~10 min/day — says "merge D" when gates are green, forwards each lane's
                   number to Claude Science, refuses any number not from the scorer.
