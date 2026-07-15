# Critical Runtime Data — FlexAIDδS Runtime Files

This directory contains (or manages) the precomputed atom-type interaction matrices
(**MC_*.dat**) **and** definition files (**AMINO*.def**, **NUCLEOTIDES*.def**) required
by the FlexAIDδS binary at runtime.

## Why These Files Matter

The core docking engine uses a **Voronoi contact-function (CF) scoring proxy** during
genetic algorithm search. This scoring relies on large precomputed lookup tables
(MC_*.dat files) that encode pairwise interaction energies between atom types.

Without these files the binary will fail early with messages such as:

    ERROR: Could not open file .../MC_st0r5.2_6.dat

## Files Currently Bundled with This Skill

**Interaction matrices (MC_*.dat):**
- `MC_st0r5.2_6.dat` — Primary matrix
- `MC_10p_3.dat`
- `MC_5p_norm_P10_M2_2.dat`

**Definition files (*.def):**

These files provide atom typing, covalent connectivity, and side-chain flexibility
definitions required by the FlexAIDδS binary.

**AMINO*.def (amino acids — 20 standard residues)**
- `AMINO.def` (version 2011.12.08) — Current recommended file. Matches the modern MC matrices.
- `AMINO8.def`, `AMINO12.def`, `AMINO26.def` — Legacy variants (from ~2000) that use entirely different atom type numbering. Using the wrong variant with current matrices will cause incorrect typing/scoring.
- File format (per-residue blocks):
  ```
  RESIDU XXX
  ATMTYP  <serial> <type_code> <name> <r/m> <parent...>
  CONECT  <atom> <bonded...>
  FLEDIH  <bond1> <bond2> ...     # rotatable side-chain dihedrals for GA sampling
  RESEND
  ```
- **ATMTYP columns** (practical interpretation):
  - Column 2 (type_code): Internal numeric code used for radii, VdW, and scoring parameters (e.g. 11=N backbone, 3=CA, 13=O, 12=NHx, 14=OH, 5=CZ, etc.).
  - Column 4: `r` = rigid (usually backbone), `m` = movable (side chain).
  - Later columns: parent atom indices for building the residue tree.
- **FLEDIH lines**: Explicitly list which bonds are treated as flexible during docking. These directly control the conformational search space sampled by the genetic algorithm. Residues with no FLEDIH (e.g. ALA, GLY, PRO) have no side-chain sampling from this file.
- Critical for execution (correct atom typing), configuration (which torsions are active), and analysis (understanding sampled degrees of freedom).

**Key residues with FLEDIH (flexible dihedrals) in the 2011 AMINO.def:**
- ARG (4), LYS (4), GLN/GLU/MET (3), ILE/LEU/PHE/TRP/TYR/HIS/ASN/ASP (2), and several with 1 (CYS, SER, THR, VAL).
- Full list available in the source AMINO.def or by running analysis on the file.

**NUCLEOTIDES*.def**
- Equivalent definitions for RNA/DNA bases and backbone (supports nucleic acid
  docking and the NATURaL module).

These files (along with the MC matrices) must be present in the binary’s base
directory at runtime. They were taken from complete WRK/ installations and are
now bundled in the skill for full self-containment.

## Additional Runtime Files

The following supporting files are also commonly required in a complete runtime
data pack (now bundled in this skill):

- `Lovell_LIB.dat` — Rotamer library used for side-chain sampling.
- `rotobs.lst` — Rotamer observation statistics (largest file).
- `SYBYL_emat.dat` — SYBYL atom type energy matrix.
- `M6_cons_3.dat`, `nrg_mat_BEST_*.dat`, `scr_*.dat` — Scoring, energy, and constraint support matrices.

The `ensure_docking_data.py --info` (or the `inspect_definition_files.py` helper)
will report on the presence and health of the entire set.

Both tools now automatically select the right balance:
- Normal interactive use → rich diagnostics by default.
- CI or constrained environments → automatic lightweight behavior.

You rarely need to specify `--quick` or `--info` manually anymore.

## How the Skill Manages These Files

The recommended tool is:

```bash
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py
```

This script will:
- Look first inside this skill’s own `data/` directory (highest priority).
- Then search common user and system locations.
- Support explicit sourcing from another installation via `--source`.
- Place **both** matrices and all *.def files next to the binary (in its base path).

For the common case of “I have another working FlexAIDδS installation”, use the
convenience wrapper:

```bash
python3 .grok/skills/flexaidds/scripts/copy_docking_data_from_install.py \
    --source /path/to/my/other/working/flexaidds
```

## Keeping the Data Up to Date

If you obtain newer or additional matrix or definition files from a more complete
installation, simply drop the new `*.dat` and `*.def` files into this `data/`
directory (or use the copy script above). The ensure tool will automatically
pick them up on future runs.

## Notes for Advanced Users / CI

- The skill treats these files as **runtime assets**, not source code.
- They are intentionally checked into the skill repository for portability.
- If you are building FlexAIDδS from source, you may generate or obtain these
  matrices as part of your build process and then feed them to the skill via
  `--source`.

For the authoritative list of required matrices and how the binary locates them,
see the FlexAIDδS source (particularly `top.cpp` and related data-loading code).

## Reproducibility Tooling for These Files

The inspector and DatasetRunner wrapper automatically record cryptographic hashes of *every* file listed in this document (plus the matrices and extras) for any serious run:

```bash
# Rich hash table + ready-to-paste block for your lab notebook or paper
python3 .grok/skills/flexaidds/scripts/inspect_definition_files.py --reproducibility

# Full campaign with professional one-pager Validation Summary + manifest
python3 .grok/skills/flexaidds/scripts/dataset_runner.py --all --tier 2 --package
```

The resulting `VALIDATION_SUMMARY.md` (inside the generated zip) contains a clean table of every critical file hash together with environment capture. This is the mechanism that makes FlexAIDδS work attractive to pharma teams and reviewers who demand full computational provenance.