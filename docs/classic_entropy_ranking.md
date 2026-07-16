# Classic FlexAID entropy ranking

**Product:** When `temperature > 0`, rank-0 is elected by classic FlexAID **soft-β free energy** (ACF / BindingMode `H − T·S`), not by raw CF and not by physical kcal “binding affinity” ledgers.

**Vibrational entropy stays.** FlexAIDdS still adds the ENCoM/tENCoM correction (`−T·S_vib`) on the BindingMode ranking energy. Classic soft-β configurational ranking is restored *and* vib remains an additive FlexAIDdS term — vib is not stripped for “classic purity.”

**Search** still optimizes the CF/contact-function scoring proxy. **Election** is entropy-aware (configurational soft-β + optional vib).

## Contract (original FlexAID)

| Piece | Classic FlexAID | FlexAIDdS default ranking |
|--------|-----------------|---------------------------|
| β | `1/T` (not `1/(kB T)`) | same for configurational weights |
| Cluster free energy (ACF) | soft-β cluster free energy | same; elects CF-path emission |
| BindingMode F | `G̃ = H̃ − T·S̃` over **mode members** (≡ ACF) | `G̃ + (−T·S_vib) [+ NATURaL]` |
| DatasetRunner S1 Softβ | — | same `G̃` over heads + `.mcf` **only if** `FLEXAIDDS_SOFTBETA_ELECTION=1` (**default OFF**) |
| Shared math | — | `LIB/SoftBetaFreeEnergy.h` |
| Rank-0 (engine) | lowest ACF / lowest F | same product role when T>0 |
| Rank-0 (DatasetRunner) | — | Softβ S1 opt-in; else CF / legacy ZH (no Softβ claim) |

See **`docs/implementation/softbeta_election_policy.md`** for Softβ vs sampling, TEMPER vs kT, and BCR=0 limits.

| Layer | Elects rank-0? |
|-------|----------------|
| Soft-β ACF / classic BindingMode F (+ vib) | **Yes** |
| Physical kB StatMech “affinity” / Shannon CSV / G_bind logs | **No** (diagnostic) |

`compute_vibrational_correction()` is unchanged; classic ranking **calls it**. Disable vib via existing FA normal-modes / weight knobs, not by forcing CF emission.

## Config

```json
"thermodynamics": {
  "temperature": 300,
  "clustering_algorithm": "CF",
  "classic_entropy_ranking": true,
  "force_cf_rank_emission": false
}
```

| Knob | Default | Effect |
|------|---------|--------|
| `temperature` | 300 | `0` disables entropy ranking |
| `classic_entropy_ranking` | **true** | soft-β ACF / BindingMode F elects rank-0 |
| `force_cf_rank_emission` | **false** | if **true**, restores P3b lowest-CF emission |

Env overrides:

- `FLEXAIDDS_FORCE_CF_RANK_EMISSION=1` → CF emission (rollback)
- `FLEXAIDDS_CLASSIC_ENTROPY_RANKING=0` → same rollback

## Easy rollback (do not need git revert for product path)

1. Set `thermodynamics.force_cf_rank_emission: true`, **or**
2. `export FLEXAIDDS_FORCE_CF_RANK_EMISSION=1`, **or**
3. `classic_entropy_ranking: false`

That restores commit `cd9004d` behavior: emit rank-0 by lowest representative CF after optional ACF sort.

Full code revert of this feature: revert the PR / branch that touches:

- `LIB/cluster.cpp` (CF re-sort gate)
- `LIB/BindingMode.cpp` / `.h` (classic F + soft-β Pose weight)
- `LIB/flexaid.h` (`force_cf_rank_emission`)
- `LIB/config_parser.cpp` / `config_defaults.h`
- `tests/test_classic_entropy_ranking.cpp`

## Code map

| File | Change |
|------|--------|
| `LIB/SoftBetaFreeEnergy.h` | Shared `G̃ = H̃ − T·S̃ ≡ E_min − T ln Z` |
| `LIB/cluster.cpp` | ACF via SoftBeta; skip post-ACF CF re-sort unless `force_cf` or `T==0` |
| `LIB/BindingMode.cpp` | Classic SoftBeta G̃ over mode members (+ vib); physical ledger unchanged |
| `LIB/DatasetRunner.cpp` | S1 elect min SoftBeta G̃ (dock T); `LEGACY_ZH` rollback |
| `LIB/flexaid.h` | `force_cf_rank_emission` |

## Success metric

Pose success = RMSD / PoseBusters — **not** lowest CF and not claimed true ΔG.

Live exhibit (pre-fix 1HNN): ACF-best cluster (freq 29) was emitted as rank 3; CF champion was rank 0. Classic ranking puts ACF-best first.

## Ensemble pipeline (layers 1–3)

Classic election is layer 4 of the full reproducibility contract (frame chart →
pocket Ω → soft-β SMFREE → ACF election). See **`docs/ensemble_pipeline.md`**.

## 1HNN ACF-vs-CF ablation (offline)

No re-dock required. The gate is pure election policy over an existing ensemble:

```bash
# Built-in pre-fix 1HNN numbers (CF champion vs dense ACF basin)
python3 scripts/acf_vs_cf_ablation.py --synthetic-1hnn

# Live target directory with <PDB>.cad (+ optional rank PDBs for rep CF)
python3 scripts/acf_vs_cf_ablation.py results/.../1HNN --json

# Unit tests (no C++ binary)
python3 -m pytest tests/test_acf_vs_cf_ablation.py -q
ctest --test-dir build -R ClassicEntropyRankingTests --output-on-failure
```

Expected synthetic verdict: **election flip** — classic elects cluster 3 (ACF≈−263, freq 29); `force_cf` elects cluster 0 (CF≈−189.9).
