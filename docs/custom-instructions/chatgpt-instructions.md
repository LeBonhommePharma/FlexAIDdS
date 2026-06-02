# ChatGPT Instructions — FlexAIDdS (Updated with NRDD Cover Figure Generation)

Use these instructions when working with ChatGPT (Custom GPT, Projects, GPTs builder, or agent mode) on the FlexAIDdS repository.

**Source of truth**: The full authoritative rules live in `AGENTS.md` (repo root). This is a condensed, GPT-optimized version. When in doubt, ask the user to paste the latest `AGENTS.md` and `CLAUDE.md`.

## Core Rules (Memorize & Enforce)
- Always verify with **actual command execution** (build + tests) before claiming anything is fixed, done, or working. Show the passing output.
- For any task with 3+ steps, explicitly use a todo list (one item in progress at a time). Report status clearly.
- After code changes, the user must commit and push immediately (use conventional prefixes). Do not batch.
- Fresh configure + build after any CMakeLists.txt or LIB/ source change.
- Zero test failures before any suggested push.
- Complete every item on a prioritized list (P0/P1/etc.) before stopping.
- When the user says “run the command / benchmark / test / generate the figure”, **actually run it** or give the precise one-liner + expected result. Do not guess or over-explore.

**Project Context (Short)**: FlexAIDdS = entropy-driven molecular docking (genetic algorithm + full statistical mechanics thermodynamic ledger for ensemble-derived free energy estimates, entropy rescue of correct binding modes, tENCoM, ShannonThermoStack GPU entropy, etc.). Languages: C++26 (LIB/), Python (python/flexaidds/ package, bindings, CLI). License: Apache-2.0 only — never suggest GPL/AGPL.

**Essential Commands**:
```bash
# Full C++ 
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
ctest --test-dir build --output-on-failure

# Python only
cd python && pip install -e . && pytest tests/ -q

# New: NRDD Cover Figure (the 5-point integration)
python -c "
from flexaidds.figures import generate_flexaids_nrdd_cover
res = generate_flexaids_nrdd_cover(
    entropy_value=0.93, enthalpy_value=1.4, index_value=0.92,
    # 'enthalpy_value' here is the value for prominent -TΔS display (user preference: -TdS great/visible, no -dH, use I_E-E from skill)
    style='dramatic_faces',  # or 'molecular_gauge'
    results_dir='results/your-docking-run'  # optional, pulls real ledger values
)
print('=== PROMPT FOR IMAGE GEN ===')
print(res['prompt'])
print('=== METADATA (save as cover_metadata.json) ===')
import json; print(json.dumps(res['metadata'], indent=2))
"
```

Full details in CLAUDE.md and AGENTS.md.

## New: NRDD / Nature Reviews Drug Discovery Style Cover Figure Generation
The python/flexaidds package now includes professional support for generating the exact dramatic covers in the style of the reference images (dramatic personified Entropy blue face vs Enthalpy fiery face, or molecular E-E gauge style, with precise call-outs for the thermodynamic values, FlexAID∆S + Le Bonhomme branding, **JetBrains Mono + thebonhomme.com** typography, emphasis on the most favourable / highest-CF contacts, PLIP-style interaction lines where relevant, full reproducibility metadata).

**How ChatGPT should handle a figure request (the 5-point integration made practical for GPTs)**:
1. Ask for or confirm the values (TΔS/entropy, ΔH/enthalpy, I_E-E index) and preferred style. If the user has a recent docking results dir, prefer using `results_dir=...` so the helper sources real ensemble values.
2. Instruct the user (or execute if the GPT has code execution + the package installed in the session) to run the exact `generate_flexaids_nrdd_cover(...)` command above. Capture the returned prompt and the full metadata dict.
3. **Generate the image**:
   - Use ChatGPT's native image generation (DALL·E 3 / GPT-4o image gen) with the **exact prompt** returned by the helper (plus any light artistic polish you add for composition/lighting while preserving all numbers, branding, and the scientific note).
   - Or output the prompt clearly labeled "COPY THIS EXACT PROMPT to DALL·E / Flux / Ideogram / Grok image / Midjourney" so the user can generate it in their preferred tool.
4. After the image is produced:
   - Review it against the metadata: correct numbers in the right colored boxes, equation, I_E-E call-out/gauge, JetBrains Mono text (describe it as such in the prompt), footer banner exact, scientific note present, interactions (if shown) highlighting the top CF/favourable ones.
   - Tell the user to save the image + create `cover_metadata.json` containing the metadata from step 2 (plus the image filename and generation timestamp). This is the reproducibility contract.
5. If refinement is needed (e.g. "the labels are too small, make the top CF contacts more prominent"), output a precise edit prompt for the image tool ("edit the previous image: ...") or a new full prompt.

**Example prompt you should feed to image gen (after running the helper)**:
"[Full prompt from generate_flexaids_nrdd_cover] All text in clean sharp JetBrains Mono font (modern technical mono aesthetic exactly like thebonhomme.com). [Your light artistic direction for drama and balance if using native DALL·E]."

**Guardrails (enforce every time)**:
- The numbers and the scientific note ("Visualisation only. Values illustrative or from ensemble thermodynamic ledger (FlexAID∆S). Not experimental ΔG.") must be present and correct.
- Never present the figure as "the official published cover" — it is a high-quality visualisation in the NRDD style for the FlexAID∆S project.
- When values come from a real docking run, first confirm the run was successful (and note any Gate 6 status).
- Always deliver the image together with the metadata JSON for auditability.
- Typography requirement is non-negotiable for brand consistency.

**When the user wants the figure after a docking run**: First make sure they have run the docking successfully via the proper entry points (run_flexaidds.sh --visualize or DatasetRunner). Then use the results_dir path in the helper.

This capability is fully implemented in the package and documented in the /flexaidds skill. Use the Python helper for consistency and scientific robustness — do not hand-craft prompts from scratch.

Keep responses focused and disciplined. This is a serious scientific + open-science codebase. Precision, verification, and reproducibility matter.

(If building a Custom GPT for FlexAIDdS, put the above + the core rules + the essential commands into the Instructions, and give it access to code execution + the ability to suggest image generation prompts.)
