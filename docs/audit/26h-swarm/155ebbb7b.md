# Audit: `155ebbb7b` — Dependabot `actions/checkout` 6.0.2 → 7.0.0

| Field | Value |
|-------|--------|
| **Full SHA** | `155ebbb7b6a002a80e786fde70ac560cdd74905e` |
| **Short** | `155ebbb7b` |
| **Parents** | `8071a5582` (1st / base) · `8d7173cfa` (2nd / PR tip) |
| **Subject** | Merge pull request #255 from LeBonhommePharma/dependabot/github_actions/actions/checkout-7.0.0 |
| **PR** | [#255](https://github.com/LeBonhommePharma/FlexAIDdS/pull/255) — `Build(deps): Bump actions/checkout from 6.0.2 to 7.0.0` |
| **Author / merge** | Dependabot bot (branch) · merged by **LeBonhommePharma** · 2026-07-14 20:20:54 −0400 |
| **Swarm index** | #76 (`docs/audit/26h-swarm/INDEX.md`) |
| **Audit scope** | Dependabot safety for GitHub Actions pin only. **No source edits** in this audit. |
| **Audit date** | 2026-07-15 |
| **Method** | `git show` / merge-tree diff; SHA↔tag resolution against `actions/checkout`; release notes + security PR #2454; workflow trigger/input inventory; PR check rollup. |
| **Verdict** | **SAFE / PASS** — pure immutable SHA pin bump; security-hardening major with no blast radius on this repo’s triggers; residual **comment hygiene** only (`# v6` left stale). |

---

## 1. Executive summary

This merge is a **Dependabot GitHub Actions dependency update** only. It rewrites **20** `uses: actions/checkout@…` pins across **12** workflow files from:

| | Old | New |
|--|-----|-----|
| **Pin** | `de0fac2e4500dabe0009e67214ff5f5447ce83dd` | `9c091bb21b7c1c1d1991bb908d89e4e9dddfe3e0` |
| **Resolved tag** | `v6.0.2` (tag object) | `v7.0.0` **and** floating `v7` (both currently = this commit) |
| **Inline comment (unchanged)** | `# v6` | `# v6` ⚠️ **stale** |

**No application code, no `LIB/`, no CMake, no Python package, no scoring/ranking, no GA budget, no benchmark protocol.** Diff is line-for-line pin replacement (+20 / −20).

**Dependabot safety: ACCEPT.** Prefer keeping the **immutable commit SHA** (already done) over floating `@v7`. The only follow-up worth a tiny hygiene PR is rewriting comments to `# v7` / `# v7.0.0`.

---

## 2. Change surface

| Path | Δ | Checkout sites |
|------|---|----------------|
| `.github/workflows/benchmark-tier1.yml` | 1 line | 1 |
| `.github/workflows/benchmark-tier2.yml` | 3 lines | 3 |
| `.github/workflows/ci.yml` | 5 lines | 5 |
| `.github/workflows/codeql.yml` | 1 line | 1 |
| `.github/workflows/coverage.yml` | 1 line | 1 |
| `.github/workflows/license-scan.yml` | 1 line | 1 |
| `.github/workflows/perf.yml` | 1 line | 1 |
| `.github/workflows/pypi-release.yml` | 2 lines | 2 |
| `.github/workflows/release.yml` | 1 line | 1 |
| `.github/workflows/sanitizers.yml` | 1 line | 1 |
| `.github/workflows/tsan.yml` | 1 line | 1 |
| `.github/workflows/update-site.yml` | 2 lines | 2 |
| **Total** | **12 files, +20 / −20** | **20 pins** |

**Out of scope / untouched at this SHA:**

- No `metal-self-hosted.yml` in the workflow tree at this commit.
- After the merge, **every** `actions/checkout@` reference under `.github/workflows` resolves to the new SHA (no partial update / split-brain pins).

---

## 3. SHA authenticity & supply-chain checks

### 3.1 Tag resolution (authoritative)

Verified via `git ls-remote https://github.com/actions/checkout.git`:

| Ref | Object SHA |
|-----|------------|
| `refs/tags/v6.0.2` | `de0fac2e4500dabe0009e67214ff5f5447ce83dd` ✅ matches **old** pin |
| `refs/tags/v6.0.3` | `9f698171ed81b15d1823a05fc7211befd50c8ae0` (annotated → peel) — **not** used; jump skipped patch |
| `refs/tags/v7.0.0` | `9c091bb21b7c1c1d1991bb908d89e4e9dddfe3e0` ✅ matches **new** pin |
| `refs/tags/v7` (floating major) | `9c091bb21b7c1c1d1991bb908d89e4e9dddfe3e0` ✅ same as v7.0.0 at audit time |

GitHub release page for `v7.0.0` (published 2026-06-18) lists the same commit and marks it **GitHub-verified** (GPG key `B5690EEEBB952194`). Tag `v7.0.0` is a **lightweight** tag pointing at that commit (API `git/tags` 404 is expected for lightweight tags).

### 3.2 Pinning posture (good)

Workflows use **full 40-char commit SHAs**, not floating `@v7` or `@main`. That is the correct Dependabot / supply-chain posture:

1. Dependabot can still propose bumps when a new release appears.
2. A compromised floating major tag cannot retarget historical workflow runs without a new commit to *this* repo.
3. Reproducible CI: checkout action content is fixed for a given workflow revision.

**Residual floating-tag risk:** none in this repo’s workflow YAML (SHA-pinned). Operators who copy-paste `@v7` elsewhere should still prefer SHA+comment.

### 3.3 Author / provenance of the PR

| Check | Result |
|-------|--------|
| PR author | `app/dependabot` (bot) ✅ expected |
| Branch | `dependabot/github_actions/actions/checkout-7.0.0` ✅ Dependabot naming |
| Ecosystem | `github-actions` (matches `.github/dependabot.yml` weekly schedule) ✅ |
| Reviews | **none** recorded on PR #255 |
| Merged by | `LeBonhommePharma` (human) |
| Diff content | SHA substitution only; no unexpected files, no `permissions:` / secret / script changes |

**Finding F1 (INFO):** Unreviewed Dependabot merges are common for pin bumps but reduce the chance a human notices the **stale `# v6` comment** (F2). Process-only; not a security defect given the trivial diff.

---

## 4. What changed in `actions/checkout` v7.0.0 (security relevance)

Upstream release notes (`v6.0.3...v7.0.0`):

| Upstream change | Security impact | Impact on FlexAIDdS |
|-----------------|-----------------|---------------------|
| **Block fork PR checkout for `pull_request_target` and `workflow_run`** ([actions/checkout#2454](https://github.com/actions/checkout/pull/2454)) | **Positive / hardening.** Default-deny unless `allow-unsafe-pr-checkout: true`. Prevents classic “checkout untrusted fork + privileged secrets” footgun. | **None operationally.** Repo workflows at this SHA use only `push`, `pull_request`, `schedule`, `workflow_dispatch`, `release` — **zero** `pull_request_target` / `workflow_run` (confirmed by `rg` over `.github/workflows/`). No job needs the new opt-in. |
| ESM upgrade + dependency bumps (`@actions/core`, `tool-cache`, `js-yaml`, `flatted`, remove `uuid`, etc.) | Routine; reduces known dep surface. | Transparent if `dist/index.js` runs. |
| Runtime `using: node24` | Requires GitHub-hosted runners that ship Node 24 for JS actions (GA on current `ubuntu-latest` / `macos-latest` / `windows-latest` as of 2026). | **No regression vs previous pin:** `v6.0.2` and `v6.0.3` already declared `using: node24`. |
| “update error wording” (#2467) | Cosmetic. | None. |

**Skipped v6.0.3** (SHA-256 repo checkout fixes, merge-commit SHA regex). This monorepo is standard SHA-1 GitHub hosting; those fixes are irrelevant here. Landing on v7 includes the post-6.0.2 line.

---

## 5. Workflow input compatibility

Checkout `with:` blocks used in this repo **at the merge tip**:

| Input | Where used | Still valid on v7? |
|-------|------------|--------------------|
| *(none — defaults)* | Most jobs | ✅ |
| `submodules: false` | `ci.yml` (cxx / python jobs) | ✅ (default is already false) |
| `fetch-depth: 0` | `coverage.yml`, `pypi-release.yml` (sdist + wheel), `update-site.yml` deploy | ✅ full history still supported |

**Not used:** custom `token`, `ssh-key`, `ref`, `repository`, `persist-credentials`, `path`, `sparse-checkout`, `allow-unsafe-pr-checkout`.

No workflow is expected to break solely from the pin bump.

---

## 6. Findings table

| ID | Severity | Finding |
|----|----------|---------|
| **F1** | Info | PR #255 merged without review. Diff is mechanical; acceptable for Dependabot Actions pins if a human still skims the SHA/tag mapping (done in this audit). |
| **F2** | **Low** | **All 20** post-merge lines still comment `# v6` while the pin is **`v7.0.0`**. Dependabot does not rewrite human comments. Misleads auditors and future greps (`rg '# v6'`). **Fix:** mechanical comment rewrite to `# v7.0.0` (or `# v7`) in a follow-up hygiene commit — no behavior change. |
| **F3** | Info | Jump 6.0.2 → 7.0.0 skips 6.0.3. Acceptable; v7 supersedes. |
| **F4** | **Positive** | v7’s fork-PR block for privileged triggers is a net security win for any future `pull_request_target` / `workflow_run` workflow. Current tree does not use those triggers. |
| **F5** | Info | PR check rollup includes unrelated failures (Coverage Analysis, Pure Python results, Windows Python bindings smoke) while **core C++ matrix, CodeQL, license-scan, hygiene, skill package, ASan/UBSan, tsan, tier-1 benchmark** succeeded. Failures are **not attributable** to the checkout pin (no checkout logic in those test scripts). Do not treat this merge as “CI green overall” history — but do not blame checkout either. |
| **F6** | Info | `package.json` for v7 reports `version: 7.0.0`, `engines.node: >=24`. Aligns with action.yml `node24`. |

---

## 7. Science / ranking / reproducibility impact

| Dimension | Impact |
|-----------|--------|
| **Ranking / CF scoring / StatMech / BindingMode** | **NO** — workflows only |
| **GA budget / DoF / EVAL_SCALE** | **NO** |
| **Claim / benchmark result semantics** | **NO** (CI orchestration only) |
| **Local-first / iCloud policy** | **NO** |
| **Reproducibility of CI environment** | **Minor positive** if v7 fixes checkout edge cases; **SHA pin preserves bit-for-bit action identity** for this workflow revision |
| **Tests adequate for this change** | N/A unit tests; validation is SHA authenticity + trigger compatibility (this audit) + post-merge CI green on checkout-using jobs |

---

## 8. Dependabot safety checklist (this PR)

| Check | Status |
|-------|--------|
| Ecosystem allowed by `.github/dependabot.yml` | ✅ `github-actions` weekly |
| Diff limited to declared dependency pins | ✅ 12 workflow files only |
| No unexpected permission / secret / script injection | ✅ |
| Pin is full commit SHA (not floating major alone) | ✅ |
| SHA matches published release tag `v7.0.0` | ✅ |
| No GPL / license matrix change | ✅ (Actions runtime, not product deps) |
| No `LIB/` / science path | ✅ |
| Privileged-trigger footgun introduced | ❌ none; hardening only |
| Required new input for existing workflows | ❌ none |
| Comment / docs aligned with version | ⚠️ F2 stale `# v6` |

---

## 9. Verdict

**SAFE_MERGE / PASS — Dependabot update is legitimate and low risk.**

Keep the merge. Optional hygiene (non-blocking):

1. Rewrite `# v6` → `# v7.0.0` on all `actions/checkout@9c091bb…` lines.
2. Prefer continue-as-SHA pins on future Dependabot PRs; reject any PR that switches to floating `@v7` without a SHA.
3. If a future workflow introduces `pull_request_target` or `workflow_run`, **do not** set `allow-unsafe-pr-checkout: true` without an explicit threat-model review (v7 will refuse fork PR checkout by default — keep that).

**No source edits recommended as part of this audit.** Ranking, scoring, and claim science are unaffected.
