# Atom Types in FlexAID∆S: The 40-Type NRGDock System

FlexAID∆S uses a 40-type atom classification system derived from the SYBYL forcefield types. Each heavy atom in a receptor or ligand is assigned one of these 40 types, and that type is used as the row and column index into the 40×40 NRGDock energy matrix (`MC_st0r5.2_6.dat`). Getting the type right is essential: a mis-typed atom scores against the wrong row of the matrix and contributes either noise or the wrong sign to CF.

This document describes all 40 types, explains the biological and chemical rationale for each, and documents the four type assignments that were incorrect in the original FlexAID and have been fixed in FlexAID∆S.

---

## The 40 Canonical Types

The types are numbered 1–40. This numbering is the canonical index used internally in all source files. The table below gives the type number, the SYBYL/chemical name as used in FlexAID∆S, a brief chemical description, and representative examples.

| # | Name | Chemical description | Examples |
|:--|:-----|:---------------------|:---------|
| 1 | C.1 | sp carbon | ≡C– (alkyne terminal), ·C≡N (nitrile carbon) |
| 2 | C.2 | sp2 carbon | Alkene =CH₂, carbonyl C=O, conjugated C |
| 3 | C.3 | sp3 carbon | Aliphatic CH₃, CH₂, CH, C |
| 4 | C.AR | Aromatic carbon | Benzene ring C, pyridine C |
| 5 | C.CAT | Carbocation / guanidinium C | Guanidinium resonance; rare in drug-like ligands |
| 6 | N.1 | sp nitrogen | Nitrile N (≡N), isocyanate N |
| 7 | N.2 | sp2 imine nitrogen | Imine =N– (pyridine-exo, Schiff base, not aromatic) |
| 8 | N.3 | Aliphatic amine (sp3) | –NH₂, –NHR, –NR₂ (non-protonated, non-amide) |
| 9 | N.4 | Quaternary ammonium N | –N⁺R₄ (fully substituted, positive charge) |
| 10 | N.AR | Aromatic nitrogen | Pyridine N, pyrimidine N, imidazole ring N |
| 11 | N.AM | Amide nitrogen | –NH–C=O, –NR–C=O; peptide bond N |
| 12 | N.PL3 | Planar sp2 amine | Aniline (resonance delocalized), guanidine terminal N |
| 13 | O.2 | sp2 oxygen | Carbonyl O (C=O in ketone, ester, amide) |
| 14 | O.3 | sp3 oxygen | Alcohol –OH, ether –O–, phosphate O |
| 15 | O.CO2 | Carboxylate oxygen | –COO⁻ in ionized carboxylic acid |
| 16 | O.AR | Aromatic oxygen | Furan ring O |
| 17 | S.2 | sp2 sulfur | Thioketone >C=S, dithioester |
| 18 | S.3 | sp3 sulfur | Thiol –SH, thioether –S–, disulfide |
| 19 | S.O | Sulfoxide S | R₂S=O |
| 20 | S.O2 | Sulfone S | R₂SO₂ |
| 21 | S.AR | Aromatic sulfur | Thiophene ring S |
| 22 | P.3 | sp3 phosphorus | Phosphate –PO₄²⁻, phosphonate –PO₃H, phosphine |
| 23 | F | Fluorine | –F (organofluorine) |
| 24 | CL | Chlorine | –Cl (organochlorine) |
| 25 | BR | Bromine | –Br (organobromine); **also used for I in FlexAID∆S** |
| 26 | I | Iodine | –I; sparse in PDB training set — remapped to BR in FlexAID∆S |
| 27 | SE | Selenium | Selenomethionine Se, selenocysteine Se |
| 28 | MG | Magnesium | Mg²⁺ cofactor |
| 29 | SR | Strontium | Sr²⁺ (uncommon; crystallographic artifact in some structures) |
| 30 | CU | Copper | Cu²⁺ / Cu⁺ in metalloenzymes |
| 31 | MN | Manganese | Mn²⁺ cofactor |
| 32 | HG | Mercury | Hg²⁺ (organomercury, heavy-atom derivatives) |
| 33 | CD | Cadmium | Cd²⁺ (heavy-atom derivative) |
| 34 | NI | Nickel | Ni²⁺ in Ni-dependent enzymes |
| 35 | ZN | Zinc | Zn²⁺ (very common: zinc-finger, carbonic anhydrase, HDAC) |
| 36 | CA | Calcium | Ca²⁺ cofactor (kinases, EF-hand proteins) |
| 37 | FE | Iron | Fe²⁺/Fe³⁺ (heme, iron-sulfur clusters) |
| 38 | CO.OH | Cobalt / hydroxo metal | Co²⁺ or hydroxo-bridged metal centers |
| 39 | DUMMY | Hydrogen + unknown | H (not scored); anything not matching types 1–38 |
| 40 | SOLVENT | Solvent / water | Structural water O; not scored against ligand in CF.com |

---

## Source of Types: MOL2 vs. SDF vs. Internal Perception

The type assigned to each ligand atom depends on the input format:

**MOL2 files** carry explicit SYBYL type strings (`C.3`, `N.ar`, `O.2`, etc.) in the `@<TRIPOS>ATOM` block. `Mol2Reader.cpp` reads these strings and maps them to canonical VCT indices via `sybyl_to_flexaid_type()`. Because MOL2 carries hybridization information, this is the most accurate path.

**SDF/MOL files** carry only element symbols (from the atoms block) and bond-order information. `SdfReader.cpp` maps element symbols to VCT types via `element_to_flexaid_type()`. Because hybridization is not explicit in SDF format, element fallbacks are used: carbon → C.3 (type 3), nitrogen → N.am (type 11), oxygen → O.3 (type 14), sulfur → S.3 (type 18). When more detailed perception is needed, the `ProcessLigand/BonMol` pipeline perceives full hybridization and aromaticity from the connectivity, then maps the perceived SYBYL type to the canonical VCT index via `sybyl_name_to_canonical_vct()` in `top.cpp`. This is the Tier 2 typing path and is used when BonMol is active.

**Receptor atoms** (from PDB files) are typed by `assign_radii_types.cpp`, which uses the PDB residue name and atom name to assign SYBYL-equivalent types. The same canonical type numbering applies.

All three paths must produce the same canonical integer for a given chemical environment — this is a shared source-of-truth invariant enforced by the identical mapping tables in `Mol2Reader.cpp`, `SdfReader.cpp`, and `top.cpp`.

---

## The Four Type Fixes in FlexAID∆S

Four atom-type assignments were incorrect in the original FlexAID code. Each fix corrects a specific mismatch between the chemical nature of the atom and the matrix row it was scored against. The fixes are applied identically in `Mol2Reader.cpp`, `SdfReader.cpp`, and `top.cpp` to ensure consistency across all input paths.

---

### Fix 1: N.2 → N.AR (type 7 → type 10)

**SYBYL name**: `N.2` (sp2 imine nitrogen, as in pyridine-exocyclic imines, Schiff bases, oximes, hydrazones)

**Old mapping**: → N.am (type 11, amide nitrogen)

**New mapping**: → N.AR (type 10, aromatic/sp2 nitrogen)

**Why this was wrong**: N.am (type 11) is an H-bond *donor* — the nitrogen donates its NH to a receptor carbonyl or charged group. The N.2 imine nitrogen is an H-bond *acceptor* — its lone pair accepts an H-bond from a protein backbone NH or Lys/Arg side chain. These two roles have opposite effects on CF: when a donor type is scored against a receptor H-bond acceptor, the contact appears favorable in the matrix (donor–acceptor = stabilizing). When an acceptor type is scored against the same receptor site, the contact also appears favorable. But if the **wrong** type is used, the calculation scores the interaction against the wrong partner type's statistics. For N.2 → N.am, the acceptor nitrogen was scored as a donor, which reverses the H-bond sign and maps contacts to the wrong quadrant of the energy matrix.

**Why N.AR is the right fix**: The vast majority of sp2 imine nitrogens in PDB drug–receptor complexes are chemically analogous to aromatic nitrogens (they are π-delocalized, acceptors, and geometrically similar to pyridine N). N.ar (type 10) has a well-populated matrix row reflecting exactly these contacts. Using N.ar for N.2 aligns the statistical potential with the correct physical chemistry.

**Chemical examples affected**: Pyridine-exo imines in kinase inhibitors (e.g., imatinib C=N–Ar); hydrazone linkers; oximes; many heterocyclic nitrogen atoms that SDF readers perceive as N.2 when not in a fully aromatic ring.

---

### Fix 2: N.3 → N.AM (type 8 → type 11)

**SYBYL name**: `N.3` (aliphatic sp3 amine: –NH₂, –NHR, –NR₂)

**Old mapping**: → N.3 (type 8)

**New mapping**: → N.AM (type 11)

**Why this was wrong**: Type 8 (N.3) has a dead matrix row — all or nearly all entries are zero or near-zero because the PDB training set classified most aliphatic amine contacts under N.am (type 11, amide/generic amine) rather than the strict N.3 category. The practical consequence is that any ligand nitrogen typed as N.3 contributed nothing to CF.com — it was invisible to the scoring function. This is particularly damaging for amine-containing pharmacophores, where the amine often makes a key hydrogen-bond or ionic interaction with the receptor that should dominate CF for that contact.

**Why N.AM is the right fix**: N.am (type 11) is the "generic nitrogen" in the NRGDock system — amide, amine, and other sp3-like nitrogens all end up in this row in the PDB training set. It has a rich, well-populated set of matrix entries. Mapping N.3 → N.am makes the amine contribute realistically to CF.

**Chemical examples affected**: Primary and secondary amines in drug scaffolds (amfetamines, propranolol, metoprolol); basic nitrogens in piperidines and morpholines when read from MOL2 as N.3; any ligand with an unprotonated sp3 amine.

---

### Fix 3: C.1 → C.2 (type 1 → type 2)

**SYBYL name**: `C.1` (sp carbon: alkyne –C≡C–, allene =C=, nitrile C≡N)

**Old mapping**: → C.1 (type 1)

**New mapping**: → C.2 (type 2, sp2 carbon)

**Why this was wrong**: Type 1 (C.1) has only 10 live matrix entries because sp carbons (alkynes, nitriles) are rare in PDB binding sites — the training set simply did not have enough contact statistics to populate the row. While C.1 is not a dead row in the same way as N.3, its sparse population means CF.com contributions from sp carbons are noisy and unreliable.

**Why C.2 is the right fix**: sp carbon is geometrically and electronically similar to sp2 carbon: both are linear or near-planar, both are relatively non-polar, and both pack in binding sites similarly. C.2 (type 2) has a far richer training set derived from the abundant aromatic, vinyl, and carbonyl carbons in drug-like molecules. The fallback to type 2 improves signal strength without introducing a chemical mismatch.

**Note**: This fix applies specifically when the MOL2 SYBYL string is `C.1`. Nitrile carbons in SDF input are typed as `C` (element fallback → C.3, type 3), so SDF input is unaffected by this particular fix.

---

### Fix 4: I → BR (type 26 → type 25)

**Element**: Iodine (I)

**Old mapping**: → I (type 26)

**New mapping**: → BR (type 25, bromine)

**Why this was wrong**: Type 26 (I) has only 3 live matrix entries — the iodine row in MC_st0r5.2_6.dat is almost entirely empty because iodine is extremely rare in PDB drug–receptor co-crystals in the training set era. The few iodinated compounds in the PDB at training time were mostly heavy-atom derivatives for phasing, not true ligands. As a result, an iodinated ligand scored against the type-26 row would appear to make essentially no contacts — its iodine atoms would be phantom atoms contributing nothing to CF.

**Why BR is the right fix**: Bromine (type 25) is electronically and sterically very similar to iodine — both are large, polarizable halogens in Group 17, both commonly engage in halogen bonding with carbonyl oxygens and aromatic rings, and both have similar van der Waals radii relative to the pocket atoms they contact. The bromine row is well-populated because brominated drug scaffolds are common in PDB binding sites. Using type 25 for iodine gives a physically reasonable score where type 26 gives noise.

**Limitation**: This is an approximation. Iodine is ~25% larger than bromine (vdW radii: Br = 1.85 Å, I = 1.98 Å) and forms stronger halogen bonds due to its larger, more polarizable electron cloud. A future version of the energy matrix trained on the current PDB (which contains substantially more iodinated ligands, particularly iodinated kinase inhibitors) would populate type 26 correctly. Until then, type 25 is the best available approximation.

---

## Verifying Atom Types for a New Ligand

Before docking a ligand with unusual chemistry, it is worth checking that its atom types will be assigned correctly.

### From MOL2 input

The SYBYL type strings in the MOL2 `@<TRIPOS>ATOM` block are used directly. Check that:

1. Aromatic nitrogens in heterocycles are typed `N.ar` (not `N.3` or `N.2`).
2. Imine nitrogens (non-aromatic =N–) are typed `N.2` (which maps to N.ar/type 10 in FlexAID∆S).
3. Amide nitrogens are typed `N.am` (which maps to N.am/type 11 directly).
4. Iodine atoms are typed `I` (which maps to BR/type 25 in FlexAID∆S).

A quick check with OpenBabel or RDKit:

```python
from rdkit import Chem
mol = Chem.MolFromMolFile("ligand.sdf")
for atom in mol.GetAtoms():
    print(atom.GetIdx(), atom.GetSymbol(), atom.GetHybridization())
```

Or to convert to MOL2 with SYBYL types:

```bash
obabel ligand.sdf -O ligand.mol2 --gen3D
```

Then inspect the `@<TRIPOS>ATOM` block — each line's 6th field is the SYBYL type.

### From SDF input

SDF input uses element-level fallbacks. The type mapping is:

| Element | FlexAID∆S type | Notes |
|:--------|:--------------:|:------|
| C | C.3 (type 3) | Generic sp3; perception upgrades to C.2/C.ar if BonMol active |
| N | N.am (type 11) | Generic amine/amide |
| O | O.3 (type 14) | Generic sp3 oxygen |
| S | S.3 (type 18) | Generic thioether/thiol |
| P | P.3 (type 22) | Phosphate/phosphonate |
| F | F (type 23) | |
| Cl | CL (type 24) | |
| Br | BR (type 25) | |
| I | BR (type 25) | Remapped; iodine row too sparse |
| Se | SE (type 27) | |
| H | DUMMY (type 39) | Not scored |
| Other | DUMMY (type 39) | Not scored |

For SDF ligands with aromatic rings or sp2 centers, the BonMol Tier-2 perception path (`ProcessLigand`) is strongly recommended — it perceives full hybridization/aromaticity and assigns the correct C.ar, N.ar, O.2, etc. types. Without it, all carbons default to type 3 (C.3), which underestimates aromatic packing contributions in CF.com.

### Diagnostic Output

FlexAID∆S logs the typed atom count per type when `FLEXAIDDS_BENCHMARK=1` is set. The `[EVAL-BUDGET]` and `[ATOM_TYPES]` output lines report the histogram of type assignments for the current ligand. This is useful for confirming that expected types (e.g., C.ar for a phenyl ring) are actually assigned.

---

## Type Completeness and the Energy Matrix

Not all 40 types are equally well-represented in the energy matrix. The rough tiers are:

**Well-populated (high signal, reliable)**: C.2, C.3, C.AR, N.AR, N.AM, O.2, O.3, S.3, CL, BR, ZN, FE, CA, MG

**Moderately populated (usable)**: C.1 (after fix → C.2), N.1, N.4, N.PL3, O.CO2, O.AR, S.2, S.O, S.O2, S.AR, P.3, F, SE, CU, MN, NI, CO.OH

**Sparse / near-dead (use with caution)**: N.3 (type 8, fixed → N.AM in FlexAID∆S), I (type 26, fixed → BR in FlexAID∆S), SR, HG, CD

**DUMMY / not scored**: H, any unrecognized element, SOLVENT

For any ligand containing atom types in the sparse tier, the CF.com score for those atoms should be interpreted cautiously — the low matrix entries mean the statistical potential provides little discriminatory power for those contacts, and the overall CF may underestimate binding contributions from those functional groups.
