# ChatGPT Instructions — FlexAIDdS

**This root file is a pointer, not a source of rules.**

- **Source of truth for all workflow rules:** `AGENTS.md` (repo root).
- **Maintained ChatGPT-specific instructions:** `docs/custom-instructions/chatgpt-instructions.md`
  (figure-generation contract, DALL·E usage, reproducibility metadata, etc.).

When starting a ChatGPT session for FlexAIDdS, paste the latest `AGENTS.md`,
`CLAUDE.md`, and `docs/custom-instructions/chatgpt-instructions.md`.

Do not add or restate rules in this file — duplicated rules drift out of sync.
Update `AGENTS.md` first, then propagate the delta into the derived files
(`CLAUDE.md`, `.grok/skills/flexaidds/SKILL.md`, `docs/custom-instructions/`).
