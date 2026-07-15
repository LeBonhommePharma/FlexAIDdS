# Audit: 00a7b6eae — Fix: Formula sha256 for v2.0.3 release tarball

## Summary (2-4 sentences)
Commit `00a7b6eae6216e89e2dafd2121cfeffb398089f6` replaces the intentional all-zero Homebrew `sha256` placeholder in `Formula/flexaidds.rb` with the real digest of the GitHub tag archive for `v2.0.3`. This is the standard second half of the project’s two-step release pattern (release PR ships placeholder + tag; follow-up PR pins the archive checksum so `brew install` integrity checks pass). **Independent re-download and dual-hash verification in this audit confirm an exact match** against the live archive at `https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.3.tar.gz`. No C++/Python engine, GA, scoring, ranking, or dataset code is touched.

## Severity: LOW

Packaging-only integrity fix; verified correct. Residual risk is only the usual GitHub source-archive stability class (not introduced by this digest) and the lack of an automated formula-checksum CI gate.

## Commit facts

| Field | Value |
|-------|--------|
| Full SHA | `00a7b6eae6216e89e2dafd2121cfeffb398089f6` |
| Short | `00a7b6eae` |
| Subject | Fix: Formula sha256 for v2.0.3 release tarball |
| Author | LP \<lp@thebonhomme.com\> |
| AuthorDate | 2026-07-15 00:36:35 -0400 |
| Parent | `388e17336ceb965bdc0a8c8c392b1ac03537a981` (Merge PR #267 — Release v2.0.3) |
| Files | `Formula/flexaidds.rb` only (+1 / −2 lines) |
| PR | [#268](https://github.com/LeBonhommePharma/FlexAIDdS/pull/268) → merge `342d6650de100537a69907d36a9d3468af132efc` |
| `git describe` | `v2.0.3-1-g00a7b6eae` (exactly one commit after the release tag) |

### Diff (entire change)

```diff
--- a/Formula/flexaidds.rb
+++ b/Formula/flexaidds.rb
@@ -3,9 +3,8 @@ class Flexaidds < Formula
   homepage "https://github.com/LeBonhommePharma/FlexAIDdS"
   # v2.0.3 includes flexaid_core Metal OBJCXX membership (PR #260) so stable
   # --with-metal links. Still ships production MC_st0r5.2_6.dat (md5 9dc93717…).
-  # sha256 placeholder until the annotated tag exists; follow-up sets the real digest.
   url "https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.3.tar.gz"
-  sha256 "0000000000000000000000000000000000000000000000000000000000000000"
+  sha256 "6c8442fc672a127db354ff3b6e08a2252e8c921372d902d062ecbf4296aef186"
   license "Apache-2.0"
```

Commit message documents the generation recipe:

```bash
curl -sL https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.3.tar.gz | shasum -a 256
```

PR body claims the digest was verified twice before merge.

## Tag / archive relationship

| Item | Value |
|------|--------|
| Annotated tag object | `8d300cf59c22cb8ccfc258738570003cb793e50a` (`refs/tags/v2.0.3`) |
| Peeled commit (`v2.0.3^{}`) | `388e17336ceb965bdc0a8c8c392b1ac03537a981` |
| Remote `ls-remote` | Tag object + peeled commit match local |
| GitHub Release | [v2.0.3](https://github.com/LeBonhommePharma/FlexAIDdS/releases/tag/v2.0.3) published 2026-07-15T04:35:38Z |
| Archive URL | `https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.3.tar.gz` → codeload redirect |

**Chicken-and-egg is intentional and correct:** the tagged tree still contains the placeholder `sha256 "0000…0"` inside `Formula/flexaidds.rb`. Homebrew does **not** verify that inner string; it verifies the **downloaded tarball bytes** against the formula definition loaded from the tap (monorepo default branch after PR #268). Pinning the real digest must happen *after* the tag exists, so it cannot live inside the tagged tree without a retag. Same pattern as:

- `421fb8256` — v2.0.1 formula sha256  
- `9e923bafa` — v2.0.2 formula sha256  

## SHA256 correctness audit (this session)

### Live archive re-download

```text
URL:     https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.3.tar.gz
size:    66002413 bytes (≈62.9 MiB compressed)
file(1): gzip compressed data, from Unix, original size modulo 2^32 270571520
prefix:  FlexAIDdS-2.0.3/
entries: 2476 paths listed by tar -tzf
```

### Digests (independent tools, same download)

| Tool | Result |
|------|--------|
| Formula claim | `6c8442fc672a127db354ff3b6e08a2252e8c921372d902d062ecbf4296aef186` |
| `shasum -a 256` | `6c8442fc672a127db354ff3b6e08a2252e8c921372d902d062ecbf4296aef186` |
| Python `hashlib.sha256` (1 MiB chunks) | `6c8442fc672a127db354ff3b6e08a2252e8c921372d902d062ecbf4296aef186` |
| **MATCH** | **YES** |

### Format validation

| Check | Result |
|-------|--------|
| Length 64 hex chars | PASS |
| Lowercase `[0-9a-f]` only | PASS |
| Not all-zeros placeholder | PASS |
| URL line still points at `v2.0.3` tag archive | PASS (unchanged) |
| Formula file on current `HEAD` still carries this digest | PASS (no later rewrite of the v2.0.3 sha) |

### Content spot-check of archive

- Top-level directory: `FlexAIDdS-2.0.3/` (GitHub’s standard `repo-tag` layout for tag archives).
- Embedded formula still has placeholder zeros and the “until the annotated tag exists” comment — **expected** for tree at `388e17336`.
- Embedded formula does **not** contain the real digest — correct; retagging was not required.

**Conclusion: the formula `sha256` is correct for the published `v2.0.3` GitHub source archive as of this audit (2026-07-15).**

## Findings

### F1. Real digest matches live GitHub archive (POSITIVE / INFO)
- Evidence: Dual independent hashes of a full re-download equal the formula string exactly (see table above). PR #268 already claimed double verification; this audit is a third independent confirmation.
- Why it matters: A wrong digest would hard-fail every stable `brew install` / `brew reinstall` with a checksum mismatch and block distribution of the Metal-link fix (v2.0.3 purpose).
- Fix recommendation: None. Keep the two-step release ritual documented (placeholder in release PR → sha PR after tag).

### F2. Placeholder-to-real is the correct post-tag follow-up (INFO)
- Evidence: Parent `388e17336` (and thus the tag tree) shipped:
  ```ruby
  sha256 "0000000000000000000000000000000000000000000000000000000000000000"
  ```
  with an explicit comment that a follow-up sets the real digest. This commit is `v2.0.3-1-g00a7b6eae` and only edits that line (+ removes the now-stale comment).
- Why it matters: Auditors sometimes flag “formula inside tarball still has zeros” as a defect. It is not: Homebrew loads the formula from the tap clone (default branch), then fetches and checks the `url` artifact.
- Fix recommendation: Optional docs note in release checklist: “never expect the tag tarball’s Formula to self-describe its own archive hash.”

### F3. Without this commit, stable Homebrew installs were effectively unusable (LOW historical / FIXED)
- Evidence: All-zero is not a valid content hash of any real archive; `brew` would reject the download (or never pass checksum verification). Release PR #267 correctly deferred the real digest until after the annotated tag existed.
- Why it matters: Window between tag publish and merge of #268 is a brief broken-stable interval for monorepo-tap users who `brew update` mid-window.
- Fix recommendation: For future releases, land the sha PR immediately after tag (as done here — ~1 minute after release merge timestamps). Optionally automate with a release workflow that computes the archive sha and opens/merges the formula PR.

### F4. No ranking, scoring, thermodynamics, or GA budget impact (INFO)
- Evidence: Diff is confined to one Ruby formula string and one comment line. No `LIB/`, `python/`, dataset YAML, campaign scripts, or CMake sources.
- Why it matters for science/repro: AGENTS.md ranking/ΔG/DoF rules are out of scope for this commit. Science claims do not change.
- Fix recommendation: N/A.

### F5. Residual: GitHub auto-generated archive bit-stability (LOW residual, pre-existing class)
- Evidence: Formula uses GitHub’s **auto-generated** tag tarball (`archive/refs/tags/…`), not an uploaded release asset with a frozen blob. Historically, rare GitHub archive-generator changes have invalidated Homebrew checksums ecosystem-wide.
- Why it matters: A future GitHub-side regeneration could make today’s correct digest fail without any repo change. Project has used the same pattern successfully for v2.0.1 / v2.0.2.
- Fix recommendation (optional hardening, not required by this commit): Prefer attaching a release asset tarball under `gh release upload` and pointing `url`/`sha256` at that asset; or document a re-pin procedure if brew checksum suddenly fails.

### F6. No automated CI assertion of formula sha vs live archive (LOW)
- Evidence: Change has no test file. PR test plan is manual (“Archive sha256 fetched twice”). This audit ran the same check offline; CI does not appear to gate formula digests on every push.
- Why it matters: Human error on a future release (wrong tag URL, typo’d hex, hash of local git-archive instead of GitHub layout) would only fail at user install time.
- Fix recommendation (optional): Small CI job on formula changes: download `url`, compute sha256, assert equality with the `sha256` line; fail if still all zeros on `main` after a tag exists.

### F7. Security / supply chain (INFO)
- Evidence: Setting a correct sha256 **enables** Homebrew’s integrity check for the stable source path. No secrets, tokens, absolute machine paths, or new network endpoints beyond the already-declared GitHub URL. Apache-2.0 license line unchanged.
- Fix recommendation: None for this commit.

### F8. Hygiene (INFO)
- Evidence: Single file under `Formula/`; no `.env`, no `/Users/…` paths, no GPL introduction. Conventional `Fix:` prefix. Matches monorepo-as-tap model used by `lebonhommepharma/flexaidds`.
- Fix recommendation: N/A.

## Ranking/scoring impact: NO

Zero changes to CF/contact-function scoring, GA operators, clustering, StatMech, BindingMode, or any election metric. Packaging integrity only.

## Reproducibility impact: YES (positive)

- **Positive:** Stable Homebrew consumers can fetch a fixed tree (`v2.0.3` @ `388e17336`) and verify bytes before build — required for trustworthy distribution of the Metal `flexaid_core` link fix.
- **Neutral/expected:** The archive’s *internal* formula still shows zeros; installers must use formula from post-#268 default branch (normal for this release process).
- **Not in scope of this commit:** Whether a full `brew install --build-from-source` succeeds on a given Mac (Metal toolchain, disk space, OpenMP) — that is install-environment, not checksum correctness.

## Tests adequate: N/A (manual verification sufficient for this class)

No unit tests required for a one-line checksum pin. Adequacy criterion is **external hash verification against the published artifact**, which:

1. Author claimed (×2) in PR #268  
2. This audit repeated with `shasum` + `hashlib` on a fresh download → **MATCH**

Optional future automation is nice-to-have (F6), not a merge blocker for this already-landed fix.

## Verdict: APPROVE

Ship / keep as-is. The formula `sha256` for `v2.0.3` is **cryptographically correct** relative to the live GitHub tag archive; the change is minimal, follows established v2.0.1/v2.0.2 practice, and has no scientific ranking surface. No source follow-up is required for correctness of this commit. Optional hardening (release-asset URL, CI formula-sha gate) is process improvement for later releases, not a defect in `00a7b6eae`.

## Audit metadata

| Item | Value |
|------|--------|
| Audit date | 2026-07-15 |
| Repo | `/Users/lp.more/Projects/FlexAIDdS` |
| Method | `git show` full patch; tag peel + `ls-remote`; PR #268 JSON; live `curl` of tag tarball; `shasum -a 256` + Python `hashlib.sha256`; tar listing + extract of embedded formula |
| Source edits in this audit | **None** (report-only) |
| Related commits | Parent release `388e17336` / tag tree; PR merge `342d6650d`; prior formula-sha `9e923bafa` (v2.0.2), `421fb8256` (v2.0.1) |
