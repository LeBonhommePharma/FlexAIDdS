# Entropy ablation design (vibrational channel)

**Status:** the BindingMode vibrational correction is **fail-closed to `0.0`**. This document replaces the retired `atoms[0].eigen` eigenvalue formula. It does not enable tENCoM on the ranking path.

## What ranking actually uses

Live BindingMode election (when `T > 0` and CF-rank emission is not forced) is **configurational soft-β** only:

`G̃ = H̃ − T·S̃` over mode members (`LIB/SoftBetaFreeEnergy.h`).

`BindingMode::compute_vibrational_correction()` is still *called* from `compute_energy()`, but it returns `0.0` on every production path. Adding zero cannot change intra-receptor pose order. Do not claim tENCoM elected the pose.

## Why the old formula is retired

An earlier path read `atoms[0].eigen[m][0]` as if it were the *m*-th **eigenvalue**. It is not.

- `atom::eigen` stores normal-mode **eigenvectors** (x/y/z displacements).
- `assign_eigen.cpp` populates real atoms starting at index 1.
- Atom 0 is a sentinel whose `eigen` pointer is NULL (`read_pdb.cpp`).

Reinterpreting an eigenvector component as an eigenvalue would invent a vibrational entropy the model never computed. Until a real eigenvalue channel is wired from ENCoM/tENCoM into BindingMode, the correction is unavailable and must be zero.

## Surface honesty

PDB output always emits:

```
REMARK Vibrational diagnostic = 0.0000 (fail_closed: no eigenvalue channel; atom::eigen is eigenvectors; proxy_only; inert)
```

Omitting the line when the value was `0.0` made the channel look unwired. Emitting `0.0` is a label, not a new ranking term.

Stderr (once per process): `[TENCOM] BindingMode vib correction disabled`.

## `FLEXAIDDS_NO_TENCOM`

This env hook still zeroes the same correction. Because production already fail-closes to `0.0`, the hook is a **structural NULL** on RMSD success: it can only affect a reported predicted_dG column if a future eigenvalue channel is wired, not which pose is emitted today.

See `docs/classic_entropy_ranking.md` for the election contract.
