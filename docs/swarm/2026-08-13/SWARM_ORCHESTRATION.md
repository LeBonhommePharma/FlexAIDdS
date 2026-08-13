# SWARM ORCHESTRATION — git is the orchestrator, the scorer is the referee   (v3, 2026-08-13)

## There is no AI orchestrator, and that is deliberate
Do NOT appoint an agent to relay between lanes. An AI in the middle adds a lossy
re-description hop, and that is precisely how this project lost hours: three different RMSD
metrics in play at once, gains double-counted twice in one night, and a "the fix landed"
report that pointed at an unmerged branch. A merge queue plus a shared scorer cannot make
those mistakes.

Spawn-time needs no coordination at all — the five lanes are independent, so the briefs go
out in any order. Orchestration exists only at MERGE time.

## The referee — score_canonical.py
    python3 /Users/lp.more/flexaidds_results/workorders/score_canonical.py \
        --frozen /Users/lp.more/flexaidds_results/workorders/ASTEX84_FROZEN_POSE_BENCHMARK.csv

Validated: prints 26/84 = 31.0% (min-CF) and 41/84 = 48.8% (pool ceiling), 15-target
selection gap. It refuses to run without spyrmsd rather than silently degrading to a weaker
metric, drops sentinel poses (CF > 1e3), and prints the anti-double-count warning.
RULE: a number that did not come out of this script does not enter any report, PR
description, or message.

## Collision map — MEASURED (two pairs collide, not one)
    vcfunction.cpp 688-742   lanes A and B INTERLEAVED IN ONE LOOP
        elec  692-698, 726-731 (A)   metal 721-739 (B)
        line 727: the METAL block READS FA->use_elec and the charge fields lane A populates
        => B is FUNCTIONALLY dependent on A. A merges first.
    config_parser.cpp        A at 116-125, B at 108-115 — five lines apart
    DatasetRunner.cpp        D at 1005-1024 / 1266 / 6433-6466
                             B at 2008 / 3249 / 3312
                             E at 5949-6320
                             closest approach D<->E is 146 lines
    gaboom.cpp               lane C alone (4,962 lines, no other lane touches it)

## Merge order — D → A → B → E (rebased on D) → C
  1. D — smallest, fully offline-verifiable, no rebuild needed to evaluate
  2. A — MUST precede B (the vcfunction.cpp:727 coupling above)
  3. B — rebase onto A before touching the vcfunction window
  4. E — largest blast radius; rebase onto D
  5. C — independent; land any time
A prior version of this file said D → B → A → E → C. That merged B before its prerequisite
and was wrong.

## Three gates before any merge
  (a) the lane's acceptance gates green
  (b) its number came from score_canonical.py
  (c) after rebuild, --frozen still prints 31.0% / 48.8% with the new feature OFF
Reject any PR that moves a DEFAULT. Reject any PR whose number came from elsewhere.

## Who owns what
  Claude Science — owns the box, launches every docking campaign, holds the single ledger,
                   re-verifies every lane's number against the frozen benchmark.
  Claude Code    — owns the merge queue (repo access, runs git directly).
                   Standing instruction:
                   "You own the merge queue for FlexAIDdS. Merge order is D, A, B, E (rebased
                    on D), C. Before merging any lane verify: (a) its acceptance gates are
                    green, (b) its number came from score_canonical.py, (c) after rebuild
                    score_canonical.py --frozen still reproduces 31.0% / 48.8% with the new
                    feature OFF. Reject any PR that changes a default. One branch per lane.
                    Report to LP after each merge."
  LP             — ~10 min/day: say "merge D" when gates are green, forward each lane's
                   number to Claude Science, refuse any number not from the scorer.

## Note on Shannon
LP's hub (~/Projects/Shannon) already encodes this swarm: a role map naming claude_code /
codex / grok_build / opencode, a "sole heavy-arm owner … refuses dual heavy owners" rule, and
"emits result / benchmark updates with cf/rmsd (never invent)". Use it for the four code seats
if you want lifecycle messaging. Claude Science CANNOT reach the gate — all socket families
are blocked in its sandbox, outbound included (AF_INET connect to a closed local port returns
EPERM, not ECONNREFUSED) — so it participates by file drop into this workorders/ directory.
Before spawning through the gate, land the staged gate_socket_up fix: on origin that function
is stat-only, so a crashed gate reports healthy forever and greenlights a spawn that fails.
