# Grok Build (CLI / Terminal) Instructions — FlexAIDdS

This is the native environment for Grok Build (the current TUI/CLI session with tool access including image_gen, run_terminal_command, search_replace, todo_write, etc.).

**Source of truth**: `AGENTS.md` (repo root) + the full project `CLAUDE.md`. The canonical `/flexaidds` skill lives at `.grok/skills/flexaidds/SKILL.md`. For Astex entropy benchmarks, also read `.agents/skills/flexaidds-benchmarking/SKILL.md`. These instructions are the "Grok Build" specialization, focused on leveraging the actual tools you have (terminal, file edit, image generation, todo management) while doing FlexAIDdS work — especially the automated NRDD cover figure generation.

## Repository Hygiene (from AGENTS.md)
- Never commit `.env`, `.env.*`, or `.envrc` files.
- Never add machine-specific absolute paths to committed skills or shared scripts. Use repo-relative paths or `FLEXAIDDS_*` env vars.
- Run `python3 scripts/check_repo_hygiene.py` before pushing skill or agent-instruction changes.

## Non-Negotiable Grok Build Rules (You Must Follow These Exactly)
- **todo_write for any task with 3+ steps**. Start with a clean list (merge: false at the beginning of a new major task). Keep **exactly one** item `in_progress`. Mark completed **immediately** when done — never batch. Re-read the current todo list before ending any turn that still has pending/in-progress work.
- **Verify with actual execution** using your run_terminal_command tool (or by instructing precise commands). Show the real output (ctest --output-on-failure, pytest -q, python -c "...", image paths from generation, etc.) before claiming anything is done.
- **After any code or file change** (search_replace, edit, new file creation), you **must** commit and push immediately using run_terminal_command with proper `git add`, `git commit -m "Conventional: message"`, and `git push`. Use conventional prefixes (Add:, Fix:, Update:, Polish:, etc.). Do not batch multiple logical changes into one commit. If push hangs, `pkill -x git` or `kill $(pgrep -f git)` and retry. Check `git config core.fsmonitor` if issues persist.
- Fresh configure + build (using your terminal tool) after any CMakeLists.txt or new .cpp/.h under LIB/.
- Zero test failures before push or "complete".
- Complete every item on prioritized lists.
- When the user (or the task) says "run the command / benchmark / generate the figure", **you actually execute it** using the tools (run_terminal_command for python, cmake, etc., and image_gen / image_edit for figures). Do not just describe — do it.

**Startup ritual (run at the start of every substantial session or when context is cold)**:
```bash
git status
find . -maxdepth 4 -iname '*skill*' -o -iname 'SKILL.md' -o -iname '*.xml' -o -iname 'AGENTS.md'
python3 .grok/skills/flexaidds/scripts/validate_skill.py   # or the flexaidds equivalent if present
```

**Project & Commands**: See the main AGENTS.md / CLAUDE.md. Key ones you can run directly with your tools:
- C++ full test build + ctest
- `cd python && python -m pip install -e . && python -m pytest tests/ -q`
- Figure generation (your superpower in this environment):
  ```bash
  python -c "
  from flexaidds.figures import generate_flexaids_nrdd_cover
  res = generate_flexaids_nrdd_cover(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92, style='dramatic_faces', results_dir='results/test_run' if exists)
  # enthalpy_value arg = value shown as -TΔS (user: make -TdS great, avoid -dH, feature I_E-E index from the skill)
  print(res['prompt'])
  "
  # Then immediately:
  # Use the image_gen tool (you have it) with the prompt + aspect_ratio="16:9"
  # Save the returned path
  # Create the cover_metadata.json with the metadata from res
  # Then git add + commit + push the image + json + any other changes
  ```

## NRDD Cover Figure Generation (Native Grok Build Experience)
You have direct access to `image_gen`, `image_edit`, `video_gen`, file tools, and terminal. This makes you the best platform for end-to-end automated figure generation.

**Exact procedure you must follow (the 5-point integration executed with real tools)**:
1. **Gather/validate inputs**: Use terminal to check for recent results dirs if the user mentions a docking run. Use todo list.
2. **Run the canonical Python helper** using run_terminal_command (with PYTHONPATH=python if needed, or after `cd python && pip install -e .`). Capture the full prompt and the metadata dict. Write the metadata to a file (cover_metadata.json) right away for reproducibility.
3. **Generate the image(s)**: Call your `image_gen` tool (and `image_edit` for refinements) with the **exact prompt** from the helper + appropriate aspect_ratio ("16:9" or "3:2"). Do not paraphrase the prompt — use it as returned (it already contains JetBrains Mono / thebonhomme.com spec, PLIP-style interaction emphasis on most favourable CF contacts, exact numbers, scientific guardrail, branding, etc.).
4. **Post-process & review**:
   - If the first result needs tightening (labels, contrast on gauge/cubes, interaction line weight for top CF contacts), immediately call `image_edit` on the returned path with a precise instruction.
   - Use your tools to inspect or describe the result if possible (or have the user confirm).
   - Confirm the image matches the metadata values, has the required typography note in spirit (the generator will render it), correct branding, scientific note, etc.
5. **Save & commit** (mandatory):
   - Place the final image in an appropriate location (e.g. next to the results dir's figures/ or docs/figures/).
   - Ensure cover_metadata.json exists with the full metadata + image filename + "generated_with": "Grok Build image_gen + image_edit", timestamp.
   - Use run_terminal_command for `git add ...`, `git commit -m "Add: NRDD dramatic_faces cover for E-E=0.92 (sourced from real ledger in results/xxx) + metadata and reproducibility json"`, `git push`.
   - If push fails, kill git processes and retry.

**Example tool usage you will actually perform**:
- run_terminal_command for the python helper.
- image_gen(prompt=the_full_prompt_from_helper, aspect_ratio="16:9")
- (if needed) image_edit(image=[the_path], prompt="Polish: make JetBrains Mono labels crisper, increase weight on the most favourable CF contacts dashed lines, ensure the gold I_E-E box has high contrast...")
- run_terminal_command for git add / commit / push of the artifacts + metadata.

**When the user says "generate the cover like the ones with the faces"**:
You start a todo list, run the helper with the right style and numbers (0.93/1.4/0.92 or from their run), generate with your image tool, polish with edit if necessary, write the metadata file, commit everything with a proper message, and show the final paths + confirmation that the ritual (including the new figure step) is complete.

**Grok Build advantages**:
- You can do the entire loop (code execution + image generation + file writing + git commit/push) in one coherent session without the user switching tools.
- Perfect for "after my last docking run, make the cover".
- You can also generate the animation variants or interaction-focused details using the same helper + video_gen / further edits.
- Always tie back to the reproducibility package the project loves (the metadata json next to the image is the figure equivalent of reproducibility.json).

For general development work, fall back to the full AGENTS.md + CLAUDE.md rules and use your tools (search_replace for edits, terminal for builds/tests, todo_write, etc.).

You are in the best possible environment to demonstrate the "verify then commit" + high-quality figure generation workflow. Do it visibly with tool calls.

Start every new major user request (especially anything involving figures or changes) with the startup ritual + a fresh todo list.
