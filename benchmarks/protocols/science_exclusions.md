# Science Exclusions (Normative)

**Status:** Normative for every roster count in the claim tables.
**Rule:** a target excluded for a *scientific* reason is recorded here with its
evidence. A target that merely failed to run is NOT an exclusion — it is a void
cell and is reported as such. A hard or ambiguous case that remains in the
roster is a **failure when it misses**, not an exclusion.

## Canonical Astex Diverse denominator: N=85 including 2HR7-as-failure

The Hartshorn 2007 Astex Diverse set is **85** complexes. FlexAIDdS claim
aggregation uses that frozen roster (`benchmarks/protocols/astex85_target_manifest.json`).
Rates are always **k/85**. A target absent from a campaign, dropped in admission,
or RMSD/PoseBusters-failing **counts as a failure**. There is **no**
`astex85_codes_84.txt` roster: shrinking the denominator by dropping 2HR7 is
forbidden.

| Target | Status | Why it stays in N=85 | Evidence |
|--------|--------|----------------------|----------|
| **2HR7** | **In roster; miss = failure** | Bound species is a cryoprotectant polyether (PEG / CCD P33), so a pose-prediction endpoint is scientifically awkward — but the deposit is in Hartshorn 2007 and the cache prepares cleanly. Excluding it would inflate success by shrinking N. | Batch `astex85_full_20260830_212437`: every cell produced poses; wall times were outliers (~3 h vs median ~315 s). That cost is not an exclusion criterion. |

There are currently **no science exclusions** from the Astex-85 claim denominator.

## Denominator rule

The canonical roster is **85 targets**, frozen at
`benchmarks/protocols/astex85_target_manifest.json`.

Any table reporting `N/84` by dropping 2HR7 is incorrect for FlexAIDdS claims
and must be restated as N=85 with 2HR7 counted as failure if it missed. Any
table reporting `N/85` must be reconcilable against the manifest.

`expected_baselines.docking_power_top1: 0.70` is a **CI gate only**
(`expected_baselines_role: ci_gate_only`). It is **not** a published or
receipted FlexAIDdS rate. Do not quote 70% as a result.

## Reading `wall_s` from a receipt

**Read the field, never a display line.** A first draft of a 2HR7 cost row quoted
`11933 / 1203 / 1218 s` and inferred a 10× seed spread. `1203` and `1218` are
**clipped prefixes** of `12034` and `12188`, produced by a `cut -c1-N` on a
receipt line for terminal width. The true spread is under 2.2%. Parse
`wall_s=([0-9]+)` from the whole receipt, and treat any wall time that implies a
large intra-target seed spread as a suspected truncation until re-read.

## What is NOT an exclusion

- **2HR7** — retained in N=85; see table above.
- **1IGJ** — retained. Its apo file contains no crystallographic waters at all,
  which makes it a *control* for the solvent policy (§2b of the admission
  contract), not a defect.
- **1TW6** — retained. Its ligand is a tetrapeptide written as `ATOM` records, so
  single-residue ligand identification is impossible and the reference-ligand
  gate carries a **declared, verified exception** (chain C minus waters == the
  SDF heavy-atom count). A declared exception is not an exclusion.
- Cells that produced zero poses, timed out, or failed a receipt gate are
  **void**, reported under their own count, and never folded into an exclusion.
