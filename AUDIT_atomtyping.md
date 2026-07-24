# AUDIT — Theme B: ATOM TYPING

**Auditor:** adversarial code review (OPS/CI discipline)
**Repo:** /Users/lp.more/Projects/FlexAIDdS @ main `d623c45ea`
**Commits:** `3f56c34b1`, `3bdde9932`, `db74575dd`
**Verdict:** SERIOUS_ISSUES

---

## Ground truth established (verified, not trusted)

Matrix `MC_st0r5.2_6.dat` parsed directly: 820 upper-triangular pairs over 40 types, 1-based.
Non-zero partner counts per type (excluding self), the fact that decides every remap:

| type | name | non-zero partners | live? |
|---|---|---|---|
| 1 | C.1 | 10 | **LIVE** |
| 2 | C.2 | 25 | live |
| 6 | N.1 | 5 | **LIVE** |
| 7 | N.2 | 13 | **LIVE** |
| 8 | N.3 | **0** | **DEAD** |
| 10 | N.ar | 20 | live |
| 11 | N.am | 21 | live |
| 25 | Br | 7 | live |
| 26 | I | 3 | **LIVE (sparse)** |

Type table (`LIB/nrgrank_matrix.h:156-179`): there IS a dedicated iodine row (26) and dedicated
C.1 (1), N.1 (6), N.2 (7) rows. **Only N.3 (row 8) is genuinely dead.**

**No env gate.** `Mol2Reader.cpp`, `SdfReader.cpp`, `top.cpp` mapping functions are plain if-chains;
the only `getenv` calls are `FLEXAIDDS_DEBUG_TYPES` debug prints. **Every remap changes the default
scoring path unconditionally** — invariant #1 (DEFAULT-OFF/bit-identity) does not apply as an escape
hatch here, so each remap is acceptable ONLY if demonstrably correct.

**Blast-radius mechanics** (`top.cpp:1398-1443` override): on the SDF path, the ProcessLigand
override applies `sybyl_name_to_canonical_vct(BonMol-perceived-name)` and fires whenever
`canon != old_t`. It only touches C/N/O/S (`continue` for all other elements) — **so iodine is never
touched by the override; the I remap reaches scoring only through `SdfReader::element_to_flexaid_type`
directly.** Benchmark = 85 Astex SDF ligands (no MOL2 inputs).

Astex blast radius, measured by parsing the 85 ligand SDFs:
- **Iodine: 0 ligands.** I→Br has zero benchmark impact; matters only for external iodinated inputs.
- **Triple bonds (sp carbon = C.1): 3 ligands** — 1Z95, 2CGR, 2D3U (all nitriles, C≡N). Their nitrile
  carbon changes row 1→row 2.
- N.2 (imine) blast radius depends on BonMol perception per-ligand (not statically enumerable here),
  but is non-trivial across the set.

---

## Per-commit findings

### `3f56c34b1` — "atom type mapping corrections (N.1, N.2->N.ar, I->Br, C.1->C.2)"
**What it really does (verified from diff + parent):** four unconditional remaps in the SYBYL→VCT
tables of Mol2Reader.cpp, SdfReader.cpp (I only), top.cpp.
**default_behavior_changed: YES** (unconditional; no gate).
**correctness: SUSPECT.** Only ONE of the four targets a dead row. The other three abandon LIVE,
type-exact rows for surrogates on a data-sparsity argument, with no pose-quality validation:

- **N.2→N.ar (10)** — the sharpest problem. Parent top.cpp already mapped **N.2→7 (its own LIVE row,
  13 entries)**. This commit changes it to 10 (N.ar). The commit message itself concedes "Row 7 is
  also live (13 entries), not dead." The stated rationale ("imine is an acceptor; N.am reversed the
  H-bond sign") justifies moving OFF N.am(11) — but the correct destination is N.2's own row **7**,
  not the aromatic-nitrogen row 10. An imine (R₂C=NR) is not aromatic; routing it to the pyridine-type
  row is chemically wrong when the exact row exists and is live. Net effect: the ProcessLigand/SDF
  path, which was ALREADY correct (7), is DEGRADED to 10. This is a regression dressed as a fix.
  `LIB/top.cpp:78`, `LIB/Mol2Reader.cpp:41`.
- **C.1→C.2 (2)** — row 1 (C.1) is LIVE (10 entries), not dead. Abandons the dedicated sp-carbon row
  for the sp2 row on "sp C rare in PDB." Changes default scoring for nitrile carbons in 3 Astex
  targets (1Z95, 2CGR, 2D3U). Unvalidated. `LIB/Mol2Reader.cpp:35`, `LIB/top.cpp:72`.
- **I→Br (25)** — row 26 (I) is LIVE-but-sparse (3 entries), NOT dead. Halogen surrogate onto Br
  (7 entries). Partner profiles are disjoint (I: {C.2,C.3,O.2}; Br: {C.ar,C.cat,N.4,N.ar,N.pl3,S.3,
  SOLVENT}). Defensible as a sparsity mitigation, but it is a modeling choice, not a correctness fix,
  and it is ungated + unvalidated. Zero Astex impact. `LIB/Mol2Reader.cpp:72`, `LIB/SdfReader.cpp:75`.
- **N.1→6 (Mol2Reader)** — the only near-clean one. Old Mol2Reader N.1→11 (N.am, a donor); new →6
  (own live row, 5 entries). top.cpp already mapped N.1→6, so no benchmark change (MOL2 path only;
  Astex has no MOL2). Direction is defensible (own row > donor surrogate). `LIB/Mol2Reader.cpp:40`.

**severity: MEDIUM.** **verdict: NEEDS_CHANGE** — revert N.2 to 7 (its own live row); gate/validate
C.1 and I surrogates or drop them. No positive validation for any of the four (commit message admits
the smoke test "does NOT confirm the fixes work").

### `3bdde9932` — "N.3 -> N.am on ProcessLigand path (top.cpp)"
**What it really does (verified):** `top.cpp:79` N.3 `8→11`, syncing top.cpp with Mol2Reader (which
already did N.3→11) and SdfReader (generic N→11).
**default_behavior_changed: YES**, and this is the largest-blast-radius commit (message claims 20+
Astex targets).
**correctness: SOUND.** This is a genuine bug fix and the ONE remap justified independent of CF
magnitude. Traced: SDF generic N → element type 11 (live). If BonMol perceives it as N.3, pre-fix
`canon=8`, and since `8 != old_t(11)` and old_t∉{4,10,15}, the override at `top.cpp:1420` **actively
overwrote a live row-11 atom onto DEAD row 8, zeroing every contact for that atom.** Post-fix
`canon=11==old_t` → override does not fire → atom stays live. Row 8 = 0 non-zero entries confirmed, so
"zeroing contacts" is an unambiguous bug and restoring them is correct. The reported CF deltas
(1IA1 −43→−71) are consistent with recovering zeroed contacts.
**severity: MEDIUM** (changes default for many targets — but justified by a provably dead row).
**verdict: MAKES_SENSE.** *Caveat:* the offered evidence is CF magnitude; with measured
Spearman(CF,RMSD)≈0, "more negative CF" is NOT evidence of better docking. The dead-row argument
stands on its own; the CF numbers should not be read as an accuracy claim.

### `db74575dd` — "update iodine expected type to 25 (BR row)"
**What it really does (verified):** `tests/test_mol2_sdf_reader.cpp:489` expected I type `26→25`.
**default_behavior_changed: NO** (test-only).
**correctness: SOUND** — matches the code change in 3f56c34b1.
**severity: NONE. verdict: MAKES_SENSE.** Note: it enshrines the I→Br remap as "expected," so if that
remap is later judged wrong the test must move with it. (Observation, not a finding: the pb_vdw parity
test `tests/test_pb_vdw_parity.cpp:143` still carries `{26,"I"}` in its `type_element` reverse map and
iterates "I" at :159 — dormant, since no atom is typed 26 anymore; belongs to the pb_vdw theme, not
this one.)

---

## Theme summary

Three commits, one motivation ("sync the ligand typing paths / stop scoring on bad rows"), but the
correctness splits cleanly on **whether the target row is dead or live**:

- **`3bdde9932` is correct** *because* row 8 (N.3) is provably DEAD — the override was zeroing live
  atoms onto it; the fix restores them. Keep it.
- **`3f56c34b1` is suspect** *because* three of its four remaps (C.1, N.2, I) move OFF LIVE,
  type-exact rows onto surrogates, on a data-sparsity rationale, ungated and with zero pose-quality
  validation. The N.2 case is worse than unvalidated: it **regresses the ProcessLigand/SDF path from a
  previously-correct N.2→7 to N.ar→10**, sending a non-aromatic imine to the aromatic-N row when its
  own live row (7, 13 entries) exists.
- **`db74575dd` is a trivial, correct test sync.**

None of these is env-gated, so all alter the DEFAULT scoring path. That is defensible for the N.3 fix
(dead row) and not defensible for the live-row surrogates without validation. Given
Spearman(CF,RMSD)≈0, the CF deltas offered as evidence cannot support any accuracy claim.

### Single highest-severity finding
**`3f56c34b1`: N.2→N.ar (type 10) regresses a previously-correct default mapping.** The parent already
mapped N.2→7 (its own LIVE energy row, 13 non-zero entries) on the ProcessLigand/SDF path; this commit
overwrites it to the aromatic-nitrogen row (10) for a non-aromatic sp2 imine. It changes the default
score of every SDF/MOL2 ligand whose nitrogen BonMol perceives as N.2, is ungated, unvalidated for
pose quality, and — uniquely among the remaps — makes the SDF path *worse* than before rather than
better. Correct fix: map N.2→7 in both Mol2Reader and top.cpp. `LIB/top.cpp:78`,
`LIB/Mol2Reader.cpp:41`.
