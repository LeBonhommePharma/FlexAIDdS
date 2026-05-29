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
- `AMINO.def` (version 2011.12.08) — Current recommended file.
- `AMINO8.def`, `AMINO12.def`, `AMINO26.def` — Legacy/variant versions that use different
  atom type numbering schemes or cutoff distances.
- Each residue is defined with:
  - `ATMTYP` lines: atom serial, numeric type code, atom name, rigid/movable flag,
    and parent indices for building internal coordinates.
  - `CONECT` lines: explicit covalent bonding information.
  - `FLEDIH` lines: which bonds are treated as rotatable dihedrals (side-chain
    flexibility sampling during the genetic algorithm).
- Critical for protein atom typing and for determining which torsions are sampled.

**NUCLEOTIDES*.def**
- Equivalent definitions for RNA/DNA bases and backbone (supports nucleic acid
  docking and the NATURaL module).

These files (along with the MC matrices) must be present in the binary’s base
directory at runtime. They were taken from complete WRK/ installations and are
now bundled in the skill for full self-containment.

## How the Skill Manages These Files

The recommended tool is:

```bash
python3 .grok/skills/flexaid-docking/scripts/ensure_docking_data.py
```

This script will:
- Look first inside this skill’s own `data/` directory (highest priority).
- Then search common user and system locations.
- Support explicit sourcing from another installation via `--source`.
- Place **both** matrices and all *.def files next to the binary (in its base path).

For the common case of “I have another working FlexAIDδS installation”, use the
convenience wrapper:

```bash
python3 .grok/skills/flexaid-docking/scripts/copy_docking_data_from_install.py \
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