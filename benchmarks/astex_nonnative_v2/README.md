# Astex Non-Native Cross-Docking Set — v2 (rebuilt 2026-09-11)

Replaces `benchmarks/astex_nonnative/`, whose roster was **not a cross-docking
set**: 66 of its 74 pairs paired two *different proteins*, and the 65-family C++
table it delegated to (`astex_nonnative_targets()`, `LIB/DatasetRunner.cpp`)
contains only 24 Astex Diverse entries among its 63 distinct native codes. The
old files are kept in place, marked deprecated, as the provenance record.

## What a row is

One row = one cross-dock. The cognate ligand of `ligand_pdb` (an Astex Diverse
entry) is docked into the receptor conformation of `target_pdb` (a **different
deposition of the same protein**). Single direction; there is no reverse row.

## Construction rule and where it comes from

Gaudreault & Najmanovich 2015, *J. Chem. Inf. Model.* 55:1323–1336,
doi `10.1021/acs.jcim.5b00078`, Methods / "Docking Experiments" (journal page K,
PDF page 11), defines the Astex non-native set as *"1112 structures representing
apo and non-apo forms of 65 proteins from the diverse set"*. Difficulty is
characterised by that paper's own measure, Figure 5 (journal page F, PDF page 6):
**the maximal displacement of any atom of the binding-site between the
apo/non-apo and holo forms**, binned `[0,1) [1,2) [2,3) [3,4) [4,5) [5,dmax]`.

Note what the rule does *not* say: the receptor is **not** required to carry a
cognate ligand. Apo receptors are in the set by definition — 354 of the 1982
pairs here have one.

## What is this build's own rule, not the paper's

Verdonk et al. 2008 (`10.1021/ci8002254`), which supplies the published
per-pair membership, is **closed access** — Unpaywall, Semantic Scholar, PMC,
CrossRef TDM and the publisher all returned no open location on 2026-09-11.
Membership here is therefore generated, not reproduced:

| step | rule |
|:--|:--|
| donors | the 85 Hartshorn Table 2 codes (`benchmarks/datasets/astex_diverse_hartshorn85.yaml`) |
| candidate receptors | every X-ray entry sharing the donor site chain's UniProt accession, **released before 2008-01-01**, minus the donor |
| site definition | protein heavy atoms within 6.0 Å of any heavy atom of the donor's cognate ligand (one representative copy) |
| acceptance | CA RMSD ≤ 3.0 Å, ≥ 80 % of site atoms present, ≥ 90 % sequence identity |

The acceptance thresholds are not decoration: a shared UniProt accession alone
admits a *different domain of the same gene product*. Measured and rejected —
`1t46` (c-Kit kinase domain) against `2e9w` (c-Kit extracellular domain), and
`1jla` (HIV-1 reverse transcriptase) against `1exq` (HIV-1 integrase), which
share P04585 because both are the *pol* polyprotein.

**Same construction, different membership**: 79 donors / 1982 pairs here versus
65 proteins / 1112 structures in the paper. The paper's success rates are **not**
baselines for this set, and none are declared.

## Files

| file | rows | contents |
|:--|--:|:--|
| `astex_nonnative_v2.csv` | 1982 | accepted pairs, 23 columns, per-row provenance |
| `astex_nonnative_v2_rejected.csv` | 187 | every excluded candidate with its measured failure reason |

Definition (two live copies, byte-identical):
`benchmarks/datasets/astex_nonnative_v2.yaml` and
`python/flexaidds/dataset_runner/datasets/astex_nonnative_v2.yaml`.

## Difficulty of the set

| bin (Å) | [0,1) | [1,2) | [2,3) | [3,4) | [4,5) | [5,dmax] |
|:--|--:|--:|--:|--:|--:|--:|
| pairs | 65 | 155 | 673 | 441 | 229 | 419 |

Median maximal site displacement 3.16 Å, q95 8.61 Å, max 23.037 Å.
65 pairs fall below 1 Å and are close to a redock; they are kept and labelled
rather than dropped, so a user can exclude them explicitly.

## Reproducing

Every value is fetched or computed, none recalled. UniProt accessions from the
RCSB data API per member (and `_struct_ref` for the site chain specifically);
ligand codes checked to be real non-polymer components of their own entry
(85/85 donors); displacement from gemmi 0.7.5 CA/P sequence-aligned
superposition with residue correspondence taken through the **SEQRES**
(`label_seq`) alignment rather than author numbering — `1g9v` numbers its
alpha chain 401–541, and author-number matching silently compares the wrong
subunit.

Gate: `python3 scripts/check_dataset_identity.py` (offline) — PASS as committed.
