# Codex / Cursor Instructions — FlexAIDdS (Code-Focused + Figure Prompt Generation)

Optimized for Cursor, Windsurf, old Codex-style completions, or any AI code editor that excels at inline code, refactoring, and generating prompts/artifacts.

**Source of truth**: AGENTS.md (root). Always have the user paste the latest AGENTS.md + CLAUDE.md into your project rules / .cursorrules when starting work on this repo.

## Core Rules (Hardcoded in Your Behavior)
- Verify with actual terminal execution (the editor's terminal or user-run commands) before claiming a change works. Show the output of ctest or pytest.
- For tasks >3 steps, the editor or you should propose a todo list (visible in chat or comments).
- After edits, remind the user (or auto-suggest) to commit + push immediately with conventional prefix.
- Fresh build required after CMake or LIB/ changes.
- 0 failing tests.
- Complete the whole prioritized list.

**Project**: FlexAIDdS entropy-driven docking engine. C++26 core in LIB/ + rich Python package in python/flexaidds/. Apache-2.0 only.

**Quick commands** (have these ready in .cursorrules or as slash commands):
```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build -j 8 && ctest --test-dir build --output-on-failure
cd python && pip install -e . && pytest tests/ -q
# Figure
python -c "from flexaidds.figures import generate_flexaids_nrdd_cover; print(generate_flexaids_nrdd_cover(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92, style='dramatic_faces')['prompt'][:800])"
```

## NRDD Cover Figure Generation (Inline + Prompt Artifact)
The new first-class feature is `generate_flexaids_nrdd_cover` in python/flexaidds/figures.py. It is the canonical way to produce the dramatic Nature Reviews Drug Discovery style covers (faces or gauge) with correct thermodynamic values, E-E index, JetBrains Mono + thebonhomme.com typography, CF-favourable interaction emphasis, full branding, and reproducibility metadata.

**How you (as Codex/Cursor AI) should behave**:
- When the user is working on docking results or wants a cover "in the style of the references with E-E = X", suggest or auto-insert the call to the helper.
- Prefer the results_dir path when a real run exists so values come from the actual ensemble ledger.
- The function returns a ready prompt + metadata. Your job is to:
  1. Help the user execute the Python call (inline or in a temp script in the editor).
  2. Take the output prompt and (if the editor has image features or you can generate artifacts) produce a high-quality image prompt ready for the user's image generator.
  3. Generate or suggest the companion cover_metadata.json content (just dump the metadata).
  4. Offer to create a small helper script (e.g. scripts/make_nrdd_cover.py) that wraps the call with the user's current values.
  5. For the actual pixel generation: output the full polished prompt (with explicit "JetBrains Mono (thebonhomme.com style)" and PLIP-style interaction notes) clearly marked "COPY TO [Grok image / DALL-E / Flux / Ideogram / Midjourney]".
  6. After the user pastes the generated image back, help review it (describe differences from spec) and suggest precise image_edit prompts if the platform supports it.
- In .cursorrules or Composer, you can have a custom command like "/flexaids-cover" that runs the helper with sensible defaults or asks for the three numbers + style, then outputs the prompt + metadata file content + commit suggestion.

**Example inline behavior**:
User: "make the dramatic cover for these values after my last docking"
You: 
- Suggest the exact python -c (with their results_dir if detected).
- Once they run it and paste the output, output the full image prompt + "Also create cover_metadata.json with this content: [pretty json]".
- "After you generate the image, we should commit results/.../figures/cover_xxx.png + cover_metadata.json + any wrapper script. Suggested message: 'Add: NRDD dramatic_faces cover for E-E=0.92 using real ledger values from last run'."

**Guardrails you must embed**:
- Always go through the official Python helper — do not invent prompts.
- Metadata + scientific note are mandatory for every figure.
- Typography spec is mandatory in the prompt you deliver.
- When values are from docking, the figure is a visualisation of that specific run's thermodynamic ledger.

This makes figure generation feel native inside the code editor workflow while keeping everything reproducible and tied to the actual science.

For pure code tasks, fall back to the standard AGENTS rules + the "Common Development Tasks" section in CLAUDE.md.

Keep the editor experience fast, precise, and verification-obsessed.
