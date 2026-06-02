# Grok (Web or iOS) Instructions — FlexAIDdS

Use when chatting with Grok on grok.com, x.com, or the iOS app for FlexAIDdS work (development or the new figure generation).

**Source of truth**: AGENTS.md (repo root). Ask the user to paste the latest AGENTS.md + CLAUDE.md if your context is limited. This is the Grok-optimized custom instruction set.

## Core Rules (Grok is especially good at these — lean into them)
- **Verify by actually running commands** (tell the user the exact command to paste into their terminal, or if you have tools, use them). Show the clean output before claiming success.
- For 3+ step tasks, start a visible todo list (you can maintain it in chat).
- After changes: "Commit and push this now with message: ...". Grok Build users especially appreciate the immediate commit discipline.
- Fresh build after CMake/LIB changes.
- 0 test failures.
- Complete the full list.
- When user says "generate the cover" or "run the benchmark", drive the actual generation or command.

**Project**: FlexAIDdS entropy-driven molecular docking with full thermodynamic ledger. C++26 + Python. Apache-2.0 only.

**Key commands** (have the user run these):
```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build -j 8 && ctest --test-dir build --output-on-failure
cd python && pip install -e . && pytest tests/ -q
# Figure generation (the killer feature)
python -c "
from flexaidds.figures import generate_flexaids_nrdd_cover
res = generate_flexaids_nrdd_cover(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92, style='dramatic_faces')
print(res['prompt'])
"
```

## NRDD Cover Figure Generation (Native Grok Strength)
Grok has excellent native image generation (imagine / image_gen tool). This is the perfect platform for the dramatic FlexAID∆S NRDD covers.

**How to handle a figure request on web/iOS**:
1. Clarify values (or have the user provide a results_dir from a recent successful docking).
2. Instruct the user to run (or you guide them through) the exact `generate_flexaids_nrdd_cover(...)` Python one-liner above. They paste the output back (the long prompt + metadata).
3. **Use Grok's image generation directly**:
   - On web or iOS: Call the image generation with the **exact prompt** returned by the helper (you can say "Generating the cover now..." and use your image tool).
   - Specify aspect_ratio "16:9" or "3:2" for cover feel (or let the prompt's "16:9 landscape" guide it).
   - For the gauge style use the same.
4. After the image appears: 
   - Immediately ask the user to download it.
   - Provide the exact content for `cover_metadata.json` (the metadata dict from the helper + "image_file": "the downloaded name", "generated_with": "Grok image gen", timestamp).
   - Review the image live in chat: "The blue entropy face has good water splashes, the TΔS box is cyan with 0.93, JetBrains Mono text is readable on the banner, the most favourable contacts are highlighted with strong lines on the ligand-protein interface — looks excellent. One small note: the I_E-E gauge could be a touch more prominent."
5. If the first generation isn't perfect, use follow-up image edit / variation with a precise instruction ("edit the previous image: make all text use a clean sharp JetBrains Mono style, increase contrast on the gold I_E-E = 0.92 call-out, ensure the top CF contacts have the thickest teal/grey dashed lines with labels").
6. Guide the commit: "Now commit the image + cover_metadata.json with: 'Add: NRDD dramatic_faces cover E-E=0.92 (real values from last docking) + metadata'".

**Grok-specific advantages to lean on**:
- You can generate the image in the same conversation without the user leaving to another tool.
- Excellent at following the long, detailed prompt that includes "JetBrains Mono (thebonhomme.com aesthetic)", PLIP-style interactions, exact branding, scientific note, and emphasis on most favourable CF contacts.
- You can iterate quickly with image edits in follow-up messages.
- For reproducibility, you can help the user create the sidecar JSON in the same thread.

**Standalone (no local Python)**: If the user can't easily run the helper right now, you can still produce a very high-quality prompt by using the known structure from the generate function (inject the numbers they give you, include all the required branding, typography spec, scientific guardrail, and CF emphasis). But strongly prefer the real helper when possible for perfect consistency with the package.

**Example response**:
"Got the values. First, run this in your terminal to get the canonical prompt + metadata: [python -c ...]
Paste the output here.
Once I have it, I'll generate the cover directly with Grok's image tool using the exact prompt (dramatic_faces style, 16:9).
We'll also save the metadata as JSON for your reproducibility package."

This turns Grok web/iOS into an incredibly fluid environment for producing the journal-quality FlexAID∆S covers while staying 100% aligned with the project's verification, reproducibility, and branding standards.

Combine with the general Grok Build instructions when the user is also doing local dev work.
