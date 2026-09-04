# Science Exclusions (Normative)

**Status:** Normative for every roster count in the claim tables.
**Rule:** a target excluded for a *scientific* reason is recorded here with its
evidence. A target that merely failed to run is NOT an exclusion — it is a void
cell and is reported as such.

| Target | Excluded from | Reason | Evidence |
|--------|---------------|--------|----------|
| **2HR7** | Astex-84 pose-prediction roster (85 → 84) | The bound species is a **cryoprotectant polyether (PEG-like)**, not a cognate ligand. A pose-prediction endpoint is undefined for it: there is no biologically meaningful "correct" pose to recover, and its 22 heavy atoms in a free polyether give an enormous conformational space. | Measured cost as well as principle: at production settings its rigid cells ran **11933 / 1203 / 1218 s** (the first ~200 min, 3 restarts each) against a campaign median of ~500 s, and produced 150 poses per cell — so the search does not fail, the endpoint is simply not interpretable. |

## Denominator rule

The canonical roster is **84 targets**, frozen at
`benchmarks/protocols/astex85_target_manifest.json` (85 codes) minus the
exclusions above. Roster file: `state/astex85_codes_84.txt`.

Any table reporting `N/85` is either using the unexcluded manifest or has
silently readmitted an exclusion; both must be corrected. Any table reporting
`N/84` must be reconcilable against this file.

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
