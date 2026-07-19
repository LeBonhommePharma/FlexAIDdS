# FlexAIDdS SWE handoff — OPS ↔ external coding agent

**Constraint (stated by OPS agent):** OPS *cannot* delegate to Claude Code / Grok Build /
Codex from its interface — no connector or terminal bridge exists. All SWE tasks below are
forwarded by the human, who boosts the chosen agent. OPS retains Benchmarker / Monitor /
OPS / CI-CD and owns every acceptance gate.

## Split
- **SWE (external agent):** implements TASK1, TASK2 on feature branches. Never merges to main.
- **OPS (me):** parity gate, ctest, Astex-85 A/B, merge decision, benchmarks, monitoring.

## Model routing (rule: scientific subtlety -> stronger model; mechanical -> Sonnet)
| Task | Nature | Model |
|---|---|---|
| TASK1 cleft deterministic order | mechanical + test wiring | **Sonnet 5** |
| TASK2 Opt1 >1-thread determinism | deep concurrency + numerics | **Opus 4.8** (fallback **Fable 5**) |

## Paste-in command lines (human runs one, in repo root)
Grok 4.5:
  grok --model grok-4.5 --repo /Users/lp.more/Projects/FlexAIDdS \
    --task-file handoff_swe/TASK1_cleft_determinism.md
Codex 5.6 "Sol":
  codex exec --model sol-5.6 --cd /Users/lp.more/Projects/FlexAIDdS \
    "$(cat handoff_swe/TASK1_cleft_determinism.md)"
Claude Code:
  claude -p "$(cat handoff_swe/TASK2_opt1_multithread_determinism.md)" \
    --model opus-4.8 --add-dir /Users/lp.more/Projects/FlexAIDdS
(Adjust flags to each CLI's actual syntax — these are templates; the task FILE is the contract.)

## Environment invariants every agent MUST honor
- cmake: /opt/homebrew/bin/cmake (NOT on PATH — export PATH or call full path)
- Build: cd build && /opt/homebrew/bin/cmake --build . --target FlexAIDdS -j4
- Determinism env: FLEXAID_SEED=12345  (NOT FLEXAIDDS_SEED_BASE)
- Benchmark env: FLEXAIDDS_NO_SEC=1 FLEXAIDDS_RESTARTS=<n> FLEXAIDDS_DATA_DIR=$PWD/build
- GA benchmarking norm: 2000 generations, pop 1000 (do not change for accuracy runs)
- ALWAYS use the energy matrix (MC_st0r5.2_6.dat); never let provenance fall back to empty.
- Report rank-0 pose RMSD (elected _0.pdb), never seed-elitism/_INI RMSD.
- Feature branch only; push, do not merge. OPS gates the merge.

## OPS gate scripts (I run these; agents can pre-run to self-check)
- Parity (default flags == main): dock 1G9V seed 12345 both engines, assert CF + 10/10 poses byte-identical.
- ctest: cd build && /opt/homebrew/bin/ctest --output-on-failure  (expect 11/11).
- Astex-85 A/B accuracy: autonomous 2000-gen, spyrmsd top-1 @2Å (python env: spyrmsd 0.9.0).
