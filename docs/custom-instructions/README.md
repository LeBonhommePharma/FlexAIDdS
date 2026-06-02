# Custom Platform Instructions for FlexAIDdS

This directory contains tailored, copy-paste-ready instruction sets / system prompts / custom GPT / project knowledge packs for different AI coding and design environments.

**Goal**: Give users of each platform (Claude, ChatGPT, Cursor/Codex, Grok web/iOS, Grok Build CLI) the exact same high-quality, disciplined, reproducible experience when working on FlexAIDdS — especially the new automated NRDD journal cover figure generation using `generate_flexaids_nrdd_cover` + platform-native image tools, while strictly following the core rules from `AGENTS.md`.

## Files
- `claude-instructions.md` — General Claude (Projects, Artifacts, long context).
- `claude-dispatch-cowork-instructions.md` — Multi-agent / dispatch / coworker patterns (one agent for data/helper, one for design/prompt polish, one for image execution, one for scientific + reproducibility review).
- `claude-design-instructions.md` — Pure design / artistic direction focus (visual language, typography JetBrains Mono + thebonhomme.com, composition, interaction clarity on most favourable CF contacts).
- `chatgpt-instructions.md` — ChatGPT / Custom GPTs / Projects (DALL·E integration + the 5-point handler).
- `codex-cursor-instructions.md` — Cursor, Windsurf, Codex-style editors (inline code + prompt artifact generation, .cursorrules friendly).
- `grok-web-ios-instructions.md` — Grok on the web or iOS app (native image gen + chat iteration).
- `grok-build-cli-instructions.md` — This exact Grok Build CLI / terminal environment (direct tool calls to image_gen, run_terminal_command, search_replace, todo_write, mandatory commit/push via tools, startup ritual).

## How to use
1. Pick the file matching your platform.
2. Paste the entire content into the appropriate place:
   - Claude: Project instructions / custom knowledge / system prompt.
   - ChatGPT: Custom GPT instructions or Project knowledge.
   - Cursor: .cursorrules or Composer rules.
   - Grok: Custom instructions or just reference in chat.
3. For figure work, the instructions tell the AI to run the real `python/flexaidds/figures.py:generate_flexaids_nrdd_cover` helper (for reproducibility + correct values + metadata), then use the platform's image generation (Grok's image_gen tool in Build, DALL·E in ChatGPT, etc.), save metadata, and commit.

All versions include:
- The non-negotiable AGENTS.md workflow (verify with real execution, todo_write, immediate commit+push, 0 failures, complete lists, run when asked).
- Full project context.
- The NRDD cover figure capability (dramatic_faces and molecular_gauge styles matching the reference images, E-E index, JetBrains Mono + thebonhomme.com typography, PLIP-inspired interactions on key/CF-favourable contacts, scientific guardrails, reproducibility via metadata JSON).
- Platform-specific details for calling image generation and handling the prompt + metadata handoff.

## Relationship to other files
- Primary rules: `AGENTS.md` (repo root).
- Detailed Claude guide: `CLAUDE.md`.
- The actual implementation: `python/flexaidds/figures.py` (the handler) + `.grok/skills/flexaidds/SKILL.md` (the /flexaidds skill manifest).
- Previous chatgpt instructions: the short root `chatgpt-instructions.md` (this docs/ version is the expanded, figure-aware one).

When the main AGENTS.md or the figure code changes, these custom versions should be refreshed.

This structure makes the powerful new figure generation (and all the disciplined dev practices) feel native no matter which AI the researcher or developer prefers.
