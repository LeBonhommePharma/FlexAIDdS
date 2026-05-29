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
- `AMINO.def`, `AMINO8.def`, `AMINO12.def`, `AMINO26.def`
- `NUCLEOTIDES.def`, `NUCLEOTIDES8.def`, `NUCLEOTIDES12.def`, `NUCLEOTIDES26.def`

These files were located on the original development system (from complete WRK/
installations) and copied into the skill to make it fully self-contained and
portable for both matrices and definition data.

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