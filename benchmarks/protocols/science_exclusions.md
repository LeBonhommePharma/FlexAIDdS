# Science Exclusions (Normative)

**Status:** Normative for every roster count in the claim tables.
**Rule:** a target excluded for a *scientific* reason is recorded here with its
evidence. A target that merely failed to run is NOT an exclusion — it is a void
cell and is reported as such.

| Target | Excluded from | Reason | Evidence |
|--------|---------------|--------|----------|
| **2HR7** | Astex-84 pose-prediction roster (85 → 84) | The bound species is a **cryoprotectant polyether (PEG-like)**, not a cognate ligand. A pose-prediction endpoint is undefined for it: there is no biologically meaningful "correct" pose to recover, and its 22 heavy atoms in a free polyether give an enormous conformational space. | Measured cost as well as principle, batch `astex85_full_20260830_212437` (555 cells, production settings, 3 restarts): rigid `wall_s` **11933 / 12034 / 12188 s** across seeds 12345 / 777777 / 999999 (3.31–3.39 h), flexible **11011 / 11104 / 11116 s** (3.06–3.09 h) — **all six cells ~3.1–3.4 h, seed spread < 2.2%**, against that batch's median of **315 s**, i.e. **37.9× the median** and the batch maximum. Every cell produced **150 poses**, so the search does not fail; the endpoint is simply not interpretable. |

## Denominator rule

The canonical roster is **84 targets**, frozen at
`benchmarks/protocols/astex85_target_manifest.json` (85 codes) minus the
exclusions above. Roster file: `state/astex85_codes_84.txt`.

Any table reporting `N/85` is either using the unexcluded manifest or has
silently readmitted an exclusion; both must be corrected. Any table reporting
`N/84` must be reconcilable against this file.

## Reading `wall_s` from a receipt

**Read the field, never a display line.** A first draft of the row above quoted
`11933 / 1203 / 1218 s` and inferred a 10× seed spread. `1203` and `1218` are
**clipped prefixes** of `12034` and `12188`, produced by a `cut -c1-N` on a
receipt line for terminal width. The true spread is under 2.2%. Parse
`wall_s=([0-9]+)` from the whole receipt, and treat any wall time that implies a
large intra-target seed spread as a suspected truncation until re-read.

## What is NOT an exclusion

- **1IGJ** — retained. Its apo file contains no crystallographic waters at all,
  which makes it a *control* for the solvent policy (§2b of the admission
  contract), not a defect.
- **1TW6** — retained. Its ligand is a tetrapeptide written as `ATOM` records, so
  single-residue ligand identification is impossible and the reference-ligand
  gate carries a **declared, verified exception** (chain C minus waters == the
  SDF heavy-atom count). A declared exception is not an exclusion.
- Cells that produced zero poses, timed out, or failed a receipt gate are
  **void**, reported under their own count, and never folded into an exclusion.
