# Ensemble pipeline (reproducibility contract)

FlexAIDdS must generate and elect ensembles under **one** physicochemically
consistent contract. That is required for **reproducibility**: the same input
must produce the same support Ω, the same soft-β sampling measure, and the same
rank-0 election objective across machines, restarts, and benchmark modes.

```
[1] Frame chart consistency
      Cartesian ⇄ genes identity gate (CI)

[2] Pocket support Ω_cleft
      ligandable top-K · cleft-centroid confinement · valid spheres/grid

[3] Soft-β CF sampling (SMFREE)
      same β=1/T · niche diversity · restarts × pockets

[4] Classic entropy election
      ACF / BindingMode H−T·S · vib additive · no CF re-sort
```

Layers 1–2 are the **geometry of the integral**. Layer 3 **samples** the CF
measure on that support. Layer 4 is the **estimator** of which basin wins.

Never mix physical-kcal Boltzmann weights into gene search, never elect on a
different objective than you sample, and never sample outside the site support.

| Layer | Role | Failure if broken |
|-------|------|-------------------|
| 1 Frame | Gene chart represents Cartesian poses | Native pose not on manifold |
| 2 Pocket | Integration domain Ω_cleft | Whole-protein / empty grid |
| 3 SMFREE | Soft-β CF niche sampling | Rank-only or physical-kB collapse |
| 4 Election | Soft-β ACF rank-0 | CF false-minimum elected |

## Environment knobs (audit / CI)

| Variable | Effect |
|----------|--------|
| `FLEXAIDDS_FRAME_CHART_STRICT=1` | Abort if native-seed RMSD > 0.1 Å |
| `FLEXAIDDS_SMFREE_REQUIRE_T=1` | Abort if SMFREE runs with T=0 |
| `FLEXAIDDS_FORCE_CF_RANK_EMISSION=1` | Rollback layer 4 to P3b lowest-CF |
| `FLEXAIDDS_CLASSIC_ENTROPY_RANKING=0` | Same as force CF emission |
| `FLEXAIDDS_CLEFT_SPHERE_FILE` | Explicit multi-cleft sphere PDB (Ω_cleft) |
| `FLEXAIDDS_SCORE_NATIVE=1` + `FLEXAIDDS_RMSDST` | Emit `[NATIVE-SEED-RMSD]` / `[FRAME_CHART]` |

## Greppable log tokens

- `[FRAME_CHART] status=ok|warn|fail ...`
- `SITE-CONFINE: cleft-centroid ...` (explicit multi-cleft)
- `[SMFREE] gen=... beta_sel=... T=...`
- `[ENTROPY_RANK] classic FlexAID: rank-0 by ACF ...`

## Tests

```bash
ctest --test-dir build -R 'EnsemblePipeline|ClassicEntropy|read_spheres|CleftCavity' --output-on-failure
```

See also: `docs/classic_entropy_ranking.md` (layer 4 detail).
