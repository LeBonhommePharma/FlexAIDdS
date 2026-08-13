# FlexAIDdS SWARM HANDOFF — start here   (2026-08-13, validated)

Everything below was checked this session: every path resolves, every source line reference
was re-grepped against the current checkout, and the scorer was re-run to confirm its numbers.

## 1. What to paste, per seat
Three items, in this order:
    1. SWARM_COMMON_PREAMBLE.md      (identical for all five seats)
    2. the seat's lane file          (see the table)
    3. score_canonical.py            (or its absolute path)
Then this sentence, verbatim:

    "Work only this lane. Do not launch a docking campaign — Claude Science owns the box.
     Every number you report must come from score_canonical.py. Ship the feature env-gated
     OFF so defaults do not move the baseline. Open a PR on your own branch; do not merge."

| seat | lane file | branch |
|------|-----------|--------|
| Claude Code (Max)      | SWARM_E_autoflex.md   | lane/e-autoflex        |
| Codex / ChatGPT Plus   | SWARM_A_charges.md    | lane/a-partial-charges |
| Cursor Pro             | SWARM_D_selection.md  | lane/d-selection-waste |
| Grok / Super Grok      | SWARM_B_metal.md      | lane/b-metal-coord     |
| Claude Code (2nd seat) | SWARM_C_thermo.md     | lane/c-thermo-gate     |

Start with E — it is the only lane that can raise the ceiling by changing the landscape.

## 2. The one command every seat must be able to run
    python3 /Users/lp.more/flexaidds_results/workorders/score_canonical.py \
        --frozen /Users/lp.more/flexaidds_results/workorders/ASTEX84_FROZEN_POSE_BENCHMARK.csv

It MUST print:
    min-CF election   26/84 =  31.0%
    pool ceiling      41/84 =  48.8%
    selection gap     15 targets
If it prints anything else, the checkout or the frozen file is wrong — stop and say so.
Requires spyrmsd; the script refuses to fall back to a weaker metric by design.

## 3. Merge queue — Claude Code owns it
Order: D → A → B → E (rebased on D) → C.
A before B is not cosmetic: vcfunction.cpp:727 has the METAL block reading FA->use_elec and
the charge fields lane A populates. Merging B first merges it before its prerequisite.
Three gates per merge: acceptance gates green; number from the scorer; after rebuild the
scorer still prints 31.0% / 48.8% with the new feature OFF.

## 4. Where the numbers come from
    campaign      : /Users/lp.more/flexaidds_results/astex84_dG_20260809_141245/
                    84/84 complete, INPUT INTEGRITY OK, engine dfc065ac…, repo aa15464e
    deep-restart  : /Users/lp.more/flexaidds_results/deeprestart_20260813_000719/
                    20 targets x 30 restarts; 4/19 ceiling-misses converted to sub-2 A
    frozen bench  : ASTEX84_FROZEN_POSE_BENCHMARK.csv — 35,873 poses, 83 targets,
                    artifact 3cc422aa-6ba6-4691-a060-3f705fd63c14
                    version  559a2075-6a6f-4429-b871-f1ac86ec0192
    inputs        : /Users/lp.more/flexaidds_results/cache_v2/astex_diverse/<PDB>/
    sites         : /Users/lp.more/flexaidds_results/astex85_sites_clean/<PDB>/
    cmake         : /opt/homebrew/bin/cmake   (NOT on PATH)

## 5. The ladder, and where each lane sits
    as-run T=300 election         15/84 = 17.9%
    min-CF election               26/84 = 31.0%   <- lane D territory, bounded at +15 targets
    reweighting cap, proven               32.5%   (leave-one-target-out CV — a hard cap)
    published bar                         45.2%
    pool ceiling, 10 restarts     41/84 = 48.8%   <- cap of ANY selection fix
    ceiling at 30 restarts        50/84 = 59.5%   (extrapolated from 4/19 measured)
    ceiling if A+E convert half   62/84 = 73.8%   (projection, not a measurement)
    hard cap                      79/84 = 94.0%   (4 never-explored targets + 1UNL are out)

Why A and E matter most: of the 42 ceiling-miss targets, 34 already place their best pose in
the RIGHT pocket (centroid <= 3 A of the crystal ligand) at the wrong geometry — 1IGJ at
2.04 A, 1N2J 2.04, 1W1P 2.05, 1JJE 2.07. That is landscape-limited, not search-limited.
Only 4 targets were never explored (centroid > 6 A) and those are a pocket-finding problem.

Framing number: our pool ceiling is the best of ~430 poses per target and reaches 48.8%; the
published TOP-10 is the best of just 10 poses and reaches 66.7%. The 2015 engine found
near-native geometry on more targets with ten poses than current FlexAIDdS does with its whole
pool. This is a regression to recover, not a frontier to advance.

## 6. Files in this pack
    HANDOFF_README.md             this file
    SWARM_ASSIGNMENT.md           seat map, collision map, merge order, ranking
    SWARM_ORCHESTRATION.md        why git orchestrates, the referee, ownership, Shannon note
    SWARM_COMMON_PREAMBLE.md      the 10 hard rules + refuted list (paste to every seat)
    SWARM_A_charges.md            lane A
    SWARM_B_metal.md              lane B
    SWARM_C_thermo.md             lane C
    SWARM_D_selection.md          lane D
    SWARM_E_autoflex.md           lane E
    score_canonical.py            the referee
    ASTEX84_FROZEN_POSE_BENCHMARK.csv   the unit of experiment
    deeprestart_result.json       per-target 30x vs 10x ceilings
