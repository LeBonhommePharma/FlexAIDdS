# Systemic fix: dataset roster divergence

Four commits, unpushed, on `main` at `/Users/lp.more/Projects/FlexAIDdS`. Each was
verified by **exporting its own committed tree** (`git archive <sha>`, no `.git`) and
running the gates, the test suites and a `c++ -fsyntax-only` compile there — not by
inspecting the working tree.

| sha | item | subject |
|-----|------|---------|
| `5f5728f4` | S1 | generate the engine's PDB code lists from a declared canonical source |
| `3a4fa608` | S2 | record the roster a run ACTUALLY executed and make it checkable |
| `f5341315` | S3 | gate the two live copies of every dataset definition against each other |
| `2f575912` | S5 | every roster must name a source that resolves, or admit it cannot |

S4 has no commit of its own: its deliverable is a measurement plus the CI wiring that
landed in `3a4fa608`.

---

## What is now structurally impossible

**The engine cannot execute a list nobody declared.**
`DatasetRunner::{astex_diverse,hap2,casf2016}_codes()` no longer contain string
literals. They call `flexaids::generated::dataset_codes(slug)`, which throws on a slug
with no declaration. `LIB/generated/dataset_codes.inc` is emitted by
`scripts/gen_dataset_codes.py` from exactly one canonical source per dataset. There is
no second copy of those 427 codes in the C++ to drift.

**A hand-edit of the generated header, or an edit to a canonical source without
regenerating, fails.** `gen_dataset_codes.py --check` regenerates and diffs, printing
both shas. This is not hypothetical: adding the `primary_source:` block to
`benchmarks/datasets/hap2.yaml` during S5 changed that file's sha256, and the export
test of the S5 commit **failed on exactly that** — the committed header recorded
`hap2.yaml sha256=e6bc9d3e8350…` while the generator produced `8dcf34e6a025…` — until
the header was regenerated. The mechanism fired on its author. (The digests quoted in
the merge simulation below, `committed 567733… != generated 0df22d…`, are a *different*
check: they are whole-header shas from the `/tmp/mergesim` run, not this failure.)

**A run that does not record its roster is now detectable as such.** Schema 2 of
`RUN_RECEIPT.json` carries `executed_codes`, a rehashable `executed_codes.txt` sidecar,
its sha256, `declared_codes_count`, and the canonical source path/sha/provenance/primary
the list claims to implement. A pre-schema-2 receipt reports `UNCHECKABLE_SCHEMA` with
zero comparisons — absence of evidence, never a pass.

**A dataset definition cannot stay silent about its source.** A definition declaring a
roster must carry `primary_source:` with either a DOI that is registered in the
CrossRef-backed cache, or `status: UNVERIFIED` **with a reason**. Free prose naming a
paper is explicitly rejected — that is what `published_source: "Gaudreault &
Najmanovich 2015 JCIM Table 2"` was on `astex_diverse.yaml` while its 85 codes were
72/85 of the Hartshorn set.

### The byte-identity proof for the codegen

Acceptance criterion: a codegen that silently changes the executed set is strictly
worse than the hardcoding it replaces. `tests/fixtures/dataset_codes/hardcoded_baseline.json`
captures the three literal lists as they stood at `ed17390b`
(blob `a08ff36dcc92`). Generated vs baseline, in value **and order**:

| list | n | sha256(codes) | order-equal to baseline |
|------|---|---------------|-------------------------|
| `astex_diverse` | 85 | `6ef2d142f2e1f3f5…` | yes |
| `hap2` | 59 | `5c06aff58d6e9e17…` | yes |
| `casf2016` | 283 | `0ecf19b76e24dfb6…` | yes |

Proved twice more at runtime: the **rebuilt** `benchmark_datasets --list-codes` returns
85 / 59 / 283 codes matching the generated lists in order, and a corrupted-`.inc`
control makes the compile fail, so the clean compile discriminates.

---

## What is merely gated

**The two live dataset-definition trees.** `benchmarks/datasets` (harness) and
`python/flexaidds/dataset_runner/datasets` (`DatasetConfig.from_yaml`) are both read at
runtime. Neither was deleted. Rule 6 compares 221 keys across 13 paired definitions and
reports 4 unpaired; 24 live divergences are quarantined with dated reasons. Nothing new
can drift silently, but the 24 remain. See the ledger's S3 blocker for the four
decisions a collapse needs.

**Provenance labels, not provenance.** All 25 declarations added in S5 are
`verified: false` or `status: UNVERIFIED`. The gate proves a source is *named* and
*resolves*; it never proves the roster was extracted from it. Four datasets name a DOI
whose cached CrossRef metadata corroborates the paper the repo already named — hap2
`10.1021/acs.jcim.5b00078` (Gaudreault), casf2016 `10.1021/acs.jcim.8b00545` (Su),
astex_nonnative `10.1021/ci8002254` (Verdonk), dude37 `10.1021/jm300687e` (Mysinger) —
and each says in the file itself that no locator has been recorded for any row. No DOI
was invented.

**Quarantine, not repair.** 131 ledger entries: 85 `LIGAND_CODE_ABSENT`, 18
`CONSTANT_COLUMN`, 9 `MIRROR_VALUE_DIVERGENCE`, 9 `CANONICAL_ONLY_KEY`, 4
`MIRROR_UNPAIRED`, 2 `MIRROR_ONLY_KEY`, 2 `PDB_CODE_ABSENT`, 2
`PRIMARY_SOURCE_UNDECLARED`. Each carries a dated, specific reason; a test fails if any
live divergence lacks one; `DEFECT_LEDGER_DRIFT` fails if an entry stops reproducing.

### Non-vacuity

Every rule prints a live denominator, and the selftest aborts **VOID** if any is zero.
Current sweep: 1703 code checks, 52 column checks, 289 ligand checks, 221 mirror keys
over 13 paired definitions, 29 `primary_source` fields inspected, 10 declared DOIs
checked, 7 roster CSVs checked for per-row provenance. `--selftest` passes known-good
and fails **all** 4 known-bad, 5 mirror and 5 provenance arms, with both must-pass
provenance arms clean. The receipt validator's selftest passes a truthful receipt, fails
7 ways of lying, and performs 36 field comparisons.

Two things were fixed rather than documented around. The identity gate's own known-good
and known-bad CSV fixtures gained a per-row `source_locator` with distinct values, so
each arm fails for the defect it encodes and nothing else. And `CANONICAL.md`'s first
draft described the annotation-prefix allowance on the wrong side of the mirror; the
**code** was changed to match the stated intent (mirror-only always fires; canonical-only
fires unless prefix-declared), which surfaced 9 real findings on `naloxone_ss`, and a
test now asserts the doc and the constant agree.

No gate is reported VOID.

---

## The worktree merge: what protects, what does not

`.claude/worktrees/agent-a8b03e38b85304ca7/` still holds every original defect. Measured
read-only against HEAD; the tree was not edited.

| file | diff lines vs HEAD |
|------|--------------------|
| `LIB/DatasetRunner.cpp` | 2200 — still the literal lists, incl. `"3QGS"` and `"3RP3"` (0 in HEAD) |
| `benchmarks/astex_diverse/astex_diverse_set.csv` | 170 — 85 rows of `LIG,2.0,2.0` |
| `benchmarks/astex_nonnative/astex_non_native_set.csv` | 150 |
| `benchmarks/datasets/casf2016.yaml` | 126 |
| `python/flexaidds/dataset_runner/datasets/casf2016.yaml` | 126 |
| `.github/workflows/ci.yml` | 49 |
| `benchmarks/protocols/astex85_target_manifest.json` | **0** |

**What protects.** The HEAD tree with those six files reverted to the worktree's version
— a merge that takes their side — fails both gates: identity `rc=1` with 79 findings
(73 `LIGAND_CODE_ABSENT`, `PDB_CODE_ABSENT` naming **3QGS**, `UNCACHED` naming **3RP3**,
`CONSTANT_COLUMN`, `ROW_PROVENANCE_UNDECLARED`, 2 `PRIMARY_SOURCE_UNDECLARED`) and
`gen_dataset_codes.py --check` `rc=1`. Both gates run in the `repo_hygiene` job, which
triggers on `push, pull_request, workflow_dispatch` — verified by parsing the workflow,
which invokes all four gate scripts from that job. The nine gate **scripts, tests and
fixtures** are absent on the worktree side, so a merge adds rather than deletes them.

**What does not.** `.github/workflows/ci.yml` is the exception, and it is the one that
matters: it is **present** in the worktree and differs by 49 lines, and the worktree's
copy invokes **none** of the four gate scripts — not the two added here, and not the
pre-existing `check_dataset_identity.py` and `check_dois.py` wiring either. A merge that
takes the worktree's `ci.yml` therefore leaves every gate on disk and wired to nothing,
removing the only automatic protection while the gate files themselves survive and look
intact. The merge simulation above deliberately did not overlay `ci.yml`, so that case
was measured by parsing the two workflows rather than simulated; treat `ci.yml` as a
merge conflict that must be resolved by hand toward HEAD.

Beyond that: a failing check blocks a merge only if `repo_hygiene` is marked **required**
in branch protection — a GitHub setting, not verifiable from inside the repo, and not
verified. A force-push or a local merge that never opens a PR is not covered. And
`astex85_target_manifest.json` is **identical** in both trees: its circular-v1 defect is
live in HEAD today, so a merge reverts nothing there and S4 does not address it.

---

## Remaining switch-over steps for the corrected roster

Nothing was switched: a docking lane was declared to be reading the live definitions, so
`astex_diverse.yaml`, `astex_diverse_set.csv` and `astex85_target_manifest.json` are
byte-identical to their prior state. (The named batch directory's `DONE_BATCH` reads
COMPLETE and no FlexAIDdS process was running, but the freeze was honoured as instructed.)

The live and corrected rosters differ on 13 of 85, measured:

- **in `astex_diverse_hartshorn85.csv` only** — 1GKC 1HVY 1HWI 1HWW 1IG3 1JLA 1LRH 1MMV 1MZC 1OYT 1SQN 1TOW 1XOQ
- **in the live list only** — 1IGJ 1MQ6 1TW6 1Y6R 2BYS 2C3I 2CET 2CGR 2D3U 2GBP 2HB1 2HR7 2J62

To switch over, in one commit:

1. In `scripts/gen_dataset_codes.py`, change the `astex_diverse` registry entry's
   `source` to `benchmarks/astex_diverse/astex_diverse_hartshorn85.csv` and set
   `provenance="ROSTER_CSV"`, `primary_source="10.1021/jm061277y"`.
2. Run `python3 scripts/gen_dataset_codes.py`, then `--check` to confirm.
3. Update `tests/fixtures/dataset_codes/hardcoded_baseline.json` for `astex_diverse`.
   The byte-identity test **will** fail until you do; that is its purpose, and updating
   it is the deliberate act of switching the executed set.
4. Update `targets:` in `astex_diverse.yaml` in **both** trees (rule 6 requires them to
   agree), and add the `primary_source:` block that the freeze prevented.
5. Delete the two `astex_diverse.yaml::PRIMARY_SOURCE_UNDECLARED` quarantine entries and
   re-check the `PDB_CODE_ABSENT`/`CONSTANT_COLUMN` entries that reference the old roster.
6. Rebuild and confirm `benchmark_datasets --benchmark astex --list-codes` returns the
   corrected 85, then run one short arm and check the receipt with
   `scripts/check_run_receipt_codes.py`.

---

## Everything I could not verify

- **Branch protection.** Whether `repo_hygiene` is a *required* check is a GitHub
  setting; not readable from the repo. Without it, CI reports but does not block.
- **That any roster matches its cited paper.** No table extraction was performed for any
  of the four DOIs named in S5. Every one is `verified: false` by design.
- **`casf2016`'s executed list.** 283 codes with `provenance=EXTRACTED_FROM_CODE` /
  `primary_source=UNVERIFIED` — the source CSV was lifted out of the C++ literal, so it
  proves only what the code already said. Both `casf2016.yaml` copies declare a 29-code
  tier-1 subset that overlaps the executed list on only 12 of 29. The prior audit's
  finding that only 115 of the 283 overlap the published CASF-2016 core set is repeated
  from that audit, **not** re-measured here.
- **The generated list against the definition YAML.** No rule compares them. Measured
  today: `astex_diverse` and `hap2` agree as sets with the generated list (hap2 also in
  order); `casf2016` does not (29 vs 283). A rule for this is the obvious next gate.
- **`astex85_target_manifest.json`.** Frozen; its circular-v1 provenance is untouched and
  identical in the worktree.
- **Full CI.** Only the gates and test suites named in the ledger were run. The C++ test
  suite and the rest of the CI matrix were not executed; `benchmark_datasets` was built
  and linked, and `DatasetRunner.cpp`/`RunReceipt.cpp` pass `-fsyntax-only` with the real
  `compile_commands.json` flags.
- **`.orig` and build-output copies.** `python/build/lib/**` is untracked build output and
  was not touched; `astex_diverse_set.csv.orig` is still present and still scanned.

---

## Rules honoured

Nothing under `/Users/lp.more/flexaidds_results` was renamed or removed; the short
verification arm wrote to `/tmp`. The other agent's worktree was read, never written.
DOI resolution went through the committed CrossRef-backed cache, never `doi.org -L`.
Nothing was pushed.
