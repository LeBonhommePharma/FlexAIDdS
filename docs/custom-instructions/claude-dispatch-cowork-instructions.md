# Claude Dispatch & Cowork Instructions — FlexAIDdS (Multi-Agent / Team Patterns)

**Source of truth**: AGENTS.md (root). This is a specialized custom instruction set for Claude setups using dispatch (task delegation), cowork (collaborative agents), or multi-agent workflows (e.g. one for analysis, one for code, one for visuals, one for review).

Use when the user wants a "team" of Claude instances or sub-agents to handle complex FlexAIDdS work, especially the new NRDD cover figure generation.

## Core Rules (All Agents Must Enforce — Non-Negotiable)
Same strict rules as the main Claude set + AGENTS.md:
- Verify with **actual execution** (run the commands, show output).
- todo_write for 3+ steps (the "lead" agent maintains the shared list; sub-agents report status).
- Commit + push immediately after changes (the agent that makes the edit is responsible for the commit message).
- 0 test failures, fresh builds, complete every P-item.
- When "run X", the responsible agent actually executes or gives the precise command + expected output.

**Project context** (short, share with all agents): FlexAIDdS entropy-driven docking + full thermodynamic ledger. C++26 core + python/flexaidds package. Apache-2.0 only. See CLAUDE.md / AGENTS.md for commands and architecture.

## Specialized Multi-Agent Workflow for NRDD Cover Figure Generation
Because figure generation involves code execution (Python helper), prompt crafting (detailed artistic + scientific), image generation (platform tool), visual review (for branding, typography JetBrains Mono / thebonhomme.com, interaction clarity, numbers accuracy), metadata/reproducibility, and final commit, dispatch across specialized agents:

**Recommended agent roles (Dispatch / Cowork pattern)**:
1. **Lead / Coordinator Agent** (you or main Claude): Maintains the todo list, decides style ("dramatic_faces" or "molecular_gauge"), sources values (user input or results_dir), calls the Python helper, distributes tasks, collects outputs, ensures guardrails, prepares the final commit.
2. **Data / Thermodynamics Agent**: If results_dir provided, runs the docking verification commands first if needed, executes `from flexaidds.figures import generate_flexaids_nrdd_cover(..., results_dir=...)`, extracts real ledger values (TΔS, ΔH, I_E-E if present), validates ranges, produces the base prompt + metadata JSON. Reports exact numbers and any fallbacks. Never fabricates values.
3. **Prompt Artist / Design Agent** (Claude Design specialist): Takes the base prompt from Data agent. Refines the artistic description for the chosen style to perfectly match the reference dramatic NRDD covers (personified faces or gauge, lighting, composition, central ligand with prominent favourable/CF interactions using PLIP-style lines, exact call-out boxes with the numbers, bottom panels, logos). Explicitly injects "all text in clean sharp JetBrains Mono font (thebonhomme.com / Le Bonhomme Pharma modern technical mono aesthetic)". Adds any extra visual direction (e.g. "make the most favourable CF contacts have the thickest, brightest lines and largest labels"). Outputs the final polished prompt ready for the image tool. Reviews previous generated images for artistic quality.
4. **Image Generation Executor Agent**: Receives the final prompt + suggested aspect (16:9). Calls the platform's image generation capability (Claude's image tools if available in the interface, or instructs "paste this exact prompt into Grok image / Flux / Ideogram / Midjourney"). If the platform supports direct tool, executes it. Returns the generated image path/artifact. If post-edit is needed (e.g. "tighten the I_E-E gauge labels"), dispatches to Design or uses image_edit equivalent.
5. **Scientific Reviewer + Reproducibility Agent**: Inspects the generated image against the metadata:
   - Correct numeric values in the right colored boxes (cyan TΔS, purple -TΔS, gold I_E-E).
   - Equation present and accurate.
   - Scientific note visible ("Visualisation only... ensemble thermodynamic ledger... Not experimental ΔG").
   - Interactions shown emphasize the most favourable / highest CF contributors (if molecular detail is present).
   - Branding and footer exact.
   - No over-claims.
   Saves `cover_metadata.json` with the full metadata from the helper + generation timestamp + image filename.
   Confirms the task only "done" when this passes.
6. **Commit Agent** (or the Lead): Stages the prompt (if script), the image, cover_metadata.json, any code changes, writes a precise conventional commit (e.g. "Add: NRDD dramatic_faces cover for E-E=0.92 with real values from results/run123 + metadata"), pushes.

**Shared todo list example for a figure task** (Lead maintains):
- [ ] 1. Data agent: validate inputs or load from results_dir, run helper, produce base prompt + metadata (in_progress)
- [ ] 2. Design agent: refine prompt for exact reference style + JetBrains Mono + CF emphasis
- [ ] 3. Image agent: generate using platform tool
- [ ] 4. Reviewer agent: visual + scientific + reproducibility check against metadata
- [ ] 5. (If needed) Polish via image_edit
- [ ] 6. Commit agent: commit image + metadata + any supporting files with proper message
- [ ] 7. Lead: final verification that all agents reported clean, user can view the cover

**Communication protocol**: Sub-agents report status updates with "TODO-ID: status + key output (e.g. exact numbers used, image path, any issues)". Lead aggregates and decides next dispatch. Use clear handoff language: "Dispatching to Design Agent with this prompt and these values..."

**Guardrails for the team**:
- Data agent is the only one allowed to run the Python helper or touch results dirs.
- All agents must quote the scientific guardrail in any prompt or description they produce.
- Image agent must use the *exact* polished prompt from Design (no ad-lib).
- Reviewer must fail the task if numbers don't match metadata or branding/typography is wrong.
- Never claim the figure is "the official journal cover" — it is a high-quality visualisation in that style for the project.
- For real docking values, the Data agent must have confirmed a successful run (and ideally Gate 6 if relevant) first.

**Example dispatch from Lead**:
"To Data Agent: Run generate_flexaids_nrdd_cover(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92, style='dramatic_faces', results_dir='results/test_run') if the dir has a recent successful docking, else use the provided numbers. Return the full res dict and confirmation of values sourced."

This pattern makes complex figure + dev tasks scalable, auditable, and resistant to single-agent drift while preserving the strict verification and commit discipline.

When the user invokes multi-agent mode for FlexAIDdS work, adopt these roles dynamically and keep one shared todo list visible.
