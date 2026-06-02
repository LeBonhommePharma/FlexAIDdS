# Claude Instructions — FlexAIDdS (Custom for Claude / Claude Projects / Artifacts)

**Source of truth**: `AGENTS.md` (repo root). Paste the latest AGENTS.md into your project knowledge if needed. This is a Claude-optimized custom instruction set focused on development + the new NRDD/FlexAID∆S journal cover figure generation capability.

## Non-Negotiable Core Workflow (from AGENTS.md — enforce strictly)
- **Verify with actual execution before claiming "done", "fixed", "implemented"**. Run the exact build/test command and show clean passing output (ctest --output-on-failure or pytest -q).
- **Use todo_write for any task with 3+ distinct steps**. Exactly one item `in_progress` at a time. Mark completed immediately. Re-read the list before ending turns with pending work.
- **After any code change, the user (or you via tool if allowed) must commit and push immediately**. Use conventional prefixes (Add:, Fix:, Update:, Polish:). Never batch unrelated changes. If git push hangs, kill stale git processes.
- Fresh `cmake -B build ... && cmake --build build` after any CMakeLists.txt or new .cpp/.h under LIB/.
- **Zero test failures before any push or "ready" claim**.
- Complete **every item** on prioritized lists (P0/P1/P2...) before stopping.
- When the user says "run the command/benchmark/test/figure generation", **actually run it** (use your tools or instruct the user precisely). Do not over-explore first.

**Project**: FlexAIDdS — entropy-driven molecular docking (GA + full statistical mechanics thermodynamic ledger for ΔG estimates, entropy rescue of binding modes, tENCoM vibrational, Shannon configurational with GPU). Primary languages C++26 (LIB/), Python (python/flexaidds/ package + bindings + CLI). Strict Apache-2.0 (no GPL/AGPL deps ever). Lead: Louis-Philippe Morency, NRGlab.

**Key commands** (always run the full relevant suite):
```bash
# C++ (recommended for kernel work)
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
ctest --test-dir build --output-on-failure

# Python package only
cd python && pip install -e . && pytest tests/ -q

# Figure generation (new)
python -c "
from flexaidds.figures import generate_flexaids_nrdd_cover
res = generate_flexaids_nrdd_cover(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92, style='dramatic_faces', results_dir='results/your-run-if-any')
print('PROMPT (copy to image tool):', res['prompt'][:500], '...')
print('METADATA for reproducibility:', res['metadata'])
"
# Then use Claude's image generation / Artifacts / paste prompt to external (Grok image, Ideogram, Flux, Midjourney) or Claude Artifacts if generating code for viz.
```

**Full details**: See CLAUDE.md (this repo's Claude guide) and AGENTS.md. Never assume — read them when starting a session.

## New Capability: NRDD / Journal Cover Figure Generation (FlexAID∆S dramatic style)
The package now has first-class support for generating the exact class of Nature Reviews Drug Discovery cover-style illustrations (dramatic personified Entropy (blue icy face/water) vs Enthalpy (fiery lava/protein), or molecular E-E gauge style, central ligand, precise call-out boxes for TΔS / -TΔS / I_E–E index, ΔG=ΔH−TΔS, FlexAID∆S + Le Bonhomme Pharma branding, **JetBrains Mono + thebonhomme.com** clean technical mono typography for all text, PLIP-style color-coded interactions on key/favourable/CF-dominant contacts, reproducibility metadata).

**How to use (always via the Python helper for consistency + scientific robustness)**:
1. (Recommended) After a real docking run that has produced results/ with ledger values: pass `results_dir=...` so it sources realistic ensemble TΔS/ΔH (or use user-provided illustrative values).
2. Call the helper (as in the command above or in code). It returns:
   - A ready-to-use, highly engineered prompt (includes all branding, fonts, scientific guardrail "Visualisation only. Values illustrative or from ensemble thermodynamic ledger (FlexAID∆S). Not experimental ΔG.", E-E index, most favourable CF contacts emphasis).
   - Full metadata dict (every param, timestamp, source, suggested image_gen + image_edit calls, git sha).
3. **Generate the image**:
   - In Claude: Paste the full prompt into your image generation interface (if available in the Claude session), or use Artifacts to output a refined prompt + any accompanying code (e.g. for a local viz script). Then tell the user "Copy this prompt to Grok / Flux / Ideogram / Midjourney for the cover".
   - Always save the returned image + a `cover_metadata.json` containing the metadata from the helper (for audit/reproducibility, just like the docking reproducibility.json).
4. Post-process if needed: Use image editing capabilities (Claude Artifacts or external) with instructions like "enhance the I_E-E gauge contrast, make JetBrains Mono labels crisper, ensure the most favourable CF contacts have the strongest visual weight (thicker lines/labels)".
5. For the skill/agent flow: If user asks for a cover "like the ones with the blue/red faces and E-E = 0.92", run the helper with those numbers, generate, attach metadata.

**Guardrails specific to figures**:
- Never claim the numbers are "the true experimental ΔG" — they are from the ensemble ledger or illustrative.
- When sourcing from docking, first run the docking via the proper scripts (run_flexaidds.sh or DatasetRunner) and confirm Gate 6 or success if relevant.
- Always include the reproducibility metadata.
- Typography: explicitly require "clean sharp JetBrains Mono (modern technical mono like thebonhomme.com)" in any image prompt.
- If using real docking results, prefer the `results_dir` path so values are real (from load_results + thermo).

**Example full flow (Claude should output this or execute via tools)**:
```bash
python -c '
from flexaidds.figures import generate_flexaids_nrdd_cover
res = generate_flexaids_nrdd_cover(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92, style="dramatic_faces")
print(res["prompt"])
import json; open("cover_metadata.json", "w").write(json.dumps(res["metadata"], indent=2))
'
# Then generate image with the printed prompt (Claude image / user pastes to preferred generator).
# Review the image for: dramatic composition matching references, correct numbers in call-outs, JetBrains Mono text, prominent favourable interactions, correct branding, scientific note present.
# Commit the prompt (or script that generates it), the image, cover_metadata.json, and any code changes.
```

**When user asks for a figure**: Use todo list (e.g. 1. source/validate values 2. call helper 3. generate image via platform tool 4. save metadata + review for accuracy/branding 5. commit). Verify by describing the generated image or (if possible) re-running commands.

Keep responses precise, disciplined, and focused on the scientific codebase. Always start substantial sessions by confirming the core rules and reading AGENTS.md + CLAUDE.md if context is cold.

This custom set is optimized for Claude's strengths in long context, code artifacts, careful reasoning, and iterative design.
