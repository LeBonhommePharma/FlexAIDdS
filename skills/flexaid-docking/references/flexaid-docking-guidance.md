# FlexAID Docking — Agent Guidance

This is the long-form reference for the `flexaid-docking` skill. Load it when
the agent needs more depth than `SKILL.md` provides — for example before
writing a multi-step plan, before touching scoring/ranking code, or before
producing user-facing prose that talks about "free energy".

---

## 1. Vocabulary contract

Keep these terms intact and use them as defined below:

- **FlexAID** — the historical genetic-algorithm docking engine; the legacy
  native binary is `FlexAID`.
- **FlexAIDdS** — the current entropy-aware fork; "dS" denotes the entropy
  (∆S) extensions. The repo root is `FlexAIDdS/`.
- **FlexAID∆S** — the canonical Unicode spelling used in publications and
  the Python package docstring (`flexaidds` package, `__init__.py`).
- **docking** — sampling ligand poses against a receptor and ranking them.
- **ensemble analysis** — aggregating a population of poses into binding
  modes and computing distribution-level quantities.
- **thermodynamic ledger** — the running record of partition-function
  contributions, configurational entropy, and vibrational entropy terms
  produced by `statmech.cpp`, `BindingMode.cpp`, `encom.cpp`, and the
  `ShannonThermoStack/` module.
- **CF / contact-function scoring proxy** — the Voronoi contact-function
  energy used as the ranking proxy (`Vcontacts.cpp`, `vcfunction.cpp`). It
  is **not** a true binding free energy; it is a fast surrogate.

When writing user-facing text:

- If the number on screen came out of the CF proxy, say "score" or
  "contact-function score", not "∆G" or "binding free energy".
- If the number came out of `StatMechEngine` or `BindingMode` with
  vibrational correction, qualify it: "configurational free-energy
  estimate", "ensemble ∆F", or similar — and name the estimator.

---

## 2. Repo-aware safety rules

Before any non-trivial edit:

1. `git status` and `git branch --show-current`. Confirm what branch you are
   on and what is uncommitted.
2. Skim `CLAUDE.md` for the relevant section. The build system, test
   command, and licensing rules are authoritative there.
3. If you are about to touch any of:
   - `LIB/gaboom.cpp`
   - `LIB/Vcontacts.cpp` / `LIB/vcfunction.cpp`
   - `LIB/statmech.cpp`
   - `LIB/BindingMode.cpp`
   - `LIB/encom.cpp`
   - `LIB/ShannonThermoStack/`
   - any file under `LIB/tENCoM/` or `LIB/LigandRingFlex/`
   then assume you are touching the scientific stack. Stop and confirm with
   the user before changing scoring weights, ranking rules, or estimator
   formulas. Pure refactors (renames, moving code, fixing warnings) are
   fine but should still come with a build + `ctest` run.

Never use:

- `git push --force` to `master`
- `git reset --hard` on an unstashed working tree
- `git branch -D` on a branch you did not create this session
- `--no-verify` on commit
- `rm -rf` outside `build*/` and `WRK/`

---

## 3. CF-proxy vs thermodynamic claims

A common failure mode: the agent reads a score from a result file and writes
"binding ∆G = -8.3 kcal/mol". This is wrong unless that score actually came
out of `StatMechEngine` with calibrated terms.

Decision table:

| Source field | What it actually is | OK to call it |
|---|---|---|
| `cf_score` / Voronoi CF output | Contact-function proxy | "CF score", "contact-function score" |
| `binding_mode.delta_F` from `BindingMode.cpp` | Ensemble free-energy estimate from partition function | "ensemble ∆F estimate (CF + StatMech)" |
| `encom.S_vib` | Vibrational entropy from ENCoM | "vibrational entropy (ENCoM)" |
| `tencm` differential output | Torsional vibrational ∆S | "torsional vibrational ∆S (tENCoM)" |
| `shannon_thermo` output | Configurational Shannon entropy | "configurational Shannon entropy" |
| Anything from `boltz2.py` | External Boltz-2 predictor result | "Boltz-2 predicted affinity" |

If unsure, name the file and the field rather than guessing the physical
meaning.

---

## 4. Standard tasks and what to do

### 4.1 "Run docking on receptor X with ligand Y"

1. Confirm the build target exists: `build/FlexAID` or `build/flexaids`.
2. If not, build per `CLAUDE.md` (`cmake .. && cmake --build . -j`).
3. Locate the inputs the user gave. If they are SMILES strings or PDB IDs,
   ask before fetching from the network.
4. Run the binary with explicit paths; capture stdout/stderr.
5. Load the results dir with `python -m flexaidds <dir>` to produce a
   structured summary.
6. Report results with the vocabulary contract above. Do not relabel CF
   scores as ∆G.

### 4.2 "Plan the next phase of the thermodynamic roadmap"

Produce a chunked plan. Each chunk should:

- Name the exact files it edits.
- Specify a unit test added or modified.
- State whether ranking can change. If yes, gate behind a CMake/Python
  flag (e.g. `FLEXAIDS_ENABLE_<feature>`) defaulting to OFF.
- Be mergeable independently. No multi-PR dependencies in a single chunk
  unless explicitly requested.

Template:

```
Step N — <title>
  Files: <comma list>
  Test: tests/<name>.cpp or python/tests/<name>.py
  Ranking impact: none | gated behind <flag>
  Build & test:
    cmake --build build -j && ctest --test-dir build
```

### 4.3 "Validate the skill packaging"

```bash
python3 skills/flexaid-docking/scripts/validate_skill.py
```

The script must exit 0. If it does not, fix the underlying file rather than
silencing the check.

### 4.4 "Review a PR / branch"

1. `git log --oneline master..HEAD` — list the commits.
2. `git diff master...HEAD --stat` — summary.
3. For each file in the diff that lives under `LIB/`, check whether it is
   in the scientific stack list (section 2). If yes, look for an
   accompanying test diff.
4. Cross-check commit message prefixes (`Fix:`, `Add:`, `Update:`,
   `Refactor:`, etc.) against the diff contents.
5. Report findings; never `gh pr merge` without explicit confirmation.

---

## 5. Planning template for Codex / Grok Build hand-off

When the user asks for a roadmap to hand off to Codex or Grok Build, use
this skeleton. It is deliberately small per chunk so each step survives
context loss.

```
ROADMAP: <feature name>

Pre-flight
  - Branch: feature/<slug>
  - CI baseline: <commit sha of last green master>
  - Files frozen for this feature: <list> (no other changes allowed in
    these chunks)

Chunk 1 — Tests & flag
  - Add: tests/test_<slug>.cpp (or python/tests/test_<slug>.py)
  - Add: CMake option FLEXAIDS_ENABLE_<SLUG>=OFF
  - Behavior: feature flag exists, default off, test confirms off-path
    matches current behavior.

Chunk 2 — Implementation behind flag
  - Touch: <files>
  - Test added in Chunk 1 now also covers on-path.
  - Off-path ranking unchanged (regression test).

Chunk 3 — Documentation
  - Update: CLAUDE.md (Common Tasks section)
  - Update: docs/IMPLEMENTATION_ROADMAP.md
  - No code edits.

Chunk 4 — Default-on (optional, separate PR)
  - Flip flag default to ON.
  - Update tests that pinned old behavior.
  - Requires explicit reviewer sign-off naming the ranking-change risk.
```

---

## 6. Things that look like bugs but are not

- The `WRK/` and `build*/` directories are scratch space; ignore them in
  diffs and do not commit their contents.
- `build_sources.ignore` is a deliberate filter file (see commit history).
  Don't "fix" it without context.
- The repo has multiple `build-*` dirs (`build-ci-fix`, `build-test`) used
  for CI debugging. Don't delete them assuming they are stale.

---

## 7. When you are blocked

If a file or claim the user references cannot be found, surface that
immediately. Do not invent a path. The user can correct you faster than the
agent can guess.

If the build fails because a `.cpp` was added but not registered in
`CMakeLists.txt`, the fix is in `CMakeLists.txt` — see "Common Build
Pitfalls" in `CLAUDE.md`.
