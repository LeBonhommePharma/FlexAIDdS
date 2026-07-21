# CMAES_INTEGRATION.md — FlexAIDdS CMA-ES search backend

**License:** Apache-2.0  
**Scope:** Opt-in CMA-ES search operator for FlexAIDdS. Scoring, CF aggregation, and
GA sources stay on their existing paths.  
**On-ramp:** `.swarm/cmaes/chunk3_onramp/artifacts/apply_integration.sh`  
**Acceptance protocol:** [`VALIDATION.md`](VALIDATION.md) items **A–G** (host dispatch).  
**Sandbox proven (do not re-edit):** VALIDATION.md **P1–P7**.

---

## 1. Enabling CMA-ES

```bash
export FLEXAIDDS_SEARCH=cmaes    # also accepts CMAES
# unset or any other value → default GA path (unchanged)
```

The gate lives only in `LIB/top.cpp`, between markers:

```text
// FLEXAIDDS_CMAES_BEGIN
… if (FLEXAIDDS_SEARCH == cmaes|CMAES) → cmaes_run_dock + fill_chromosomes
  else → GA(...)
// FLEXAIDDS_CMAES_END
```

Default behavior with the env unset is **bit-identical to the pre-integration GA path**
(the `else` arm calls `GA(...)`). ParallelDock / campaign / screen modes are
untouched — CMA-ES only replaces the **standard single-GA** call site.

Install wiring from a clean tree:

```bash
bash .swarm/cmaes/chunk3_onramp/artifacts/apply_integration.sh
# or, after landing chunk1 artifacts:
SWARM_ARTIFACTS=.swarm/cmaes/chunk1_adapter/artifacts \
  bash .swarm/cmaes/chunk3_onramp/artifacts/apply_integration.sh
```

Re-running the script is a **no-op** when markers and sources are already present.

---

## 2. Eval budget equivalence (claim A/B)

GA and CMA-ES arms must be **eval-matched**, not “same CLI flags by habit.”

| Axis | GA | CMA-ES |
|------|----|--------|
| Population / offspring | `pop` | `λ` (lambda) |
| Generations / iterations | `gens` | `gens` |
| Total CF evaluations | `pop × gens` | `λ × gens` |

**Claim budget used in VALIDATION item E / KICKOFF:**

\[
\lambda \times \text{gens} \equiv \text{pop} \times \text{gens}
= 1000 \times 2000 = 2 \times 10^{6}
\]

- Keep **generations fixed**; scale population / λ for DoF budget policy
  (`FLEXAIDDS_EVAL_SCALE_DIHEDRAL`, optional `FLEXAIDDS_BUDGET_SCALE`) — see
  `AGENTS.md` scientific guardrails.
- Mode `0` (gen-scale) is legacy; mode `-1` (fixed pop+gen) is oracle-ceiling only.
- Never “match” a campaign by freezing both axes because a log printed `1000×6000`
  — that was the **base**, not the effective budget (`[EVAL-BUDGET]` logs).

Short smoke budgets (e.g. `100×50`) are allowed for wiring/doctor checks only;
they **do not** close VALIDATION **D/E**.

---

## 3. Exactly five seam functions

The adapter (`LIB/cmaes_search.{cpp,h}`) may call **only** these engine symbols
for scoring / gene plumbing (P2 / `nm -uC` footprint):

| Seam | Defined in | Role |
|------|------------|------|
| `set_gene_lim` | `gaboom.cpp` | Copy FA min/max/del/map into `genlim[]` |
| `set_bins` | `gaboom.cpp` | Bin widths / `nbin` from gene limits |
| `eval_chromosome` | `gaboom.cpp` | Build pose + evaluate CF via target (`ic2cf`) |
| `get_cf_evalue` | `ic2cf.cpp` | Full CF aggregator (GA fitness / StatMech energy) |
| `get_apparent_cf_evalue` | `ic2cf.cpp` | Display / ranking aggregator |

Nothing else engine-side may appear as an **undefined** symbol on the adapter
object (P2). Those five are **defined** in the engine so the full binary links
(P3).

**Integration surface is additive only (P4):**

| Touched | How |
|---------|-----|
| `LIB/cmaes_search.cpp`, `LIB/cmaes_search.h` | New TUs only |
| `LIB/CMakeLists.txt` | `cmaes_search.cpp` in `FLEXAID_CORE_SOURCES` **after** `gaboom.cpp` |
| `LIB/top.cpp` | Opt-in ternary + include, marker-guarded |

| **Untouched** | |
|---------------|--|
| `LIB/ic2cf.cpp` | No edits |
| `LIB/gaboom.cpp` | No edits |

Ranking, clustering, and output order stay on the existing post-search path
unless a separate thermodynamic-integration change is explicitly requested
(see `AGENTS.md`).

---

## 4. Entropy trace CSV format

CMA-ES emits a per-generation collapse trace (used by VALIDATION **F** and
`analysis/collapse_fingerprint.py` from chunk4).

**Canonical columns** (header required; comment lines `#…` allowed before header):

```text
gen,H_search,H_energy,F,best_cf,n_evals
```

| Column | Required | Meaning |
|--------|----------|---------|
| `gen` | yes | Generation / iteration index (0-based) |
| `H_search` | yes | Search / distribution Shannon entropy (nats) |
| `H_energy` | yes | Energy-histogram Shannon entropy (nats) |
| `F` | yes | Free-energy / Helmholtz proxy from the search ledger |
| `best_cf` | yes | Best CF (scoring proxy) seen so far |
| `n_evals` | optional | Cumulative CF evaluations |

Flexible aliases accepted by `collapse_fingerprint.py` include
`generation`/`iter` for `gen`, `hsearch` for `H_search`, `free_energy` for `F`,
`cf_best` for `best_cf`, etc. (see that script’s `_COLUMN_ALIASES`).

### Single-basin nuance (P6)

On a **mock / single-basin** objective the energy histogram can collapse to one
occupied bin of width matching the population, so:

- `H_energy → ln(λ)` (e.g. `ln(64) ≈ 4.1589` for λ=64) is expected, not a bug.
- `H_search` is a diversity metric (0.5×softmax CF weights + 0.5×per-dim allele
  histogram Shannon). Rank-only log-weights are intentionally **not** used
  (they are generation-invariant). On a collapsing population, `H_search`
  decreases (e.g. host smoke Δ≈1 nats over a short CMA run). It may drop from
  ~`+0.95` nats toward smaller positive values as the
  search distribution concentrates (`→ −63.5` nats in the sandbox mock).
- `F → −2.4641` was the sandbox free-energy endpoint on that mock.

On the **real rugged CF surface** (VALIDATION **F**), expect a richer
`H_energy` trajectory; do not demand the mock single-basin numbers.

Example mock (chunk4 testdata shape):

```csv
gen,H_search,H_energy,F,best_cf,n_evals
0,0.950000,4.158883,-0.100000,12.500000,64
…
9,0.000020,0.100000,-2.464100,1.000000,640
```

Fingerprint:

```bash
python3 analysis/collapse_fingerprint.py path/to/entropy_trace.csv \
  --out fingerprint.json
```

---

## 5. How to A/B vs GA

Goal: same complex, same eval budget, **only** the search operator differs.

```bash
# Shared budget (example claim budget)
export FLEXAIDDS_GA_POPULATION=1000   # or CLI --ga-population 1000
export FLEXAIDDS_GA_GENERATIONS=2000  # or CLI --ga-generations 2000
# Ensure effective evals ≈ 2e6 (check [EVAL-BUDGET] logs)

COMPLEX=1G9V   # one Astex complex for host D/E

# Arm A — GA (default)
unset FLEXAIDDS_SEARCH
# … run dock → out_ga/  (result.csv / elected pose / RMSD + best CF)

# Arm B — CMA-ES
export FLEXAIDDS_SEARCH=cmaes
# … run dock → out_cmaes/  (same metrics + entropy_trace.csv)
```

Report at minimum (VALIDATION **D** / **E**):

| Metric | Source |
|--------|--------|
| Elected-pose RMSD (Å) vs crystal | PoseBusters / ordered heavy RMSD path |
| Elected CF / best CF | `result.csv` / REMARK CF |
| ΔRMSD, ΔCF (CMA-ES − GA) | Eval-matched pair only |
| Entropy trace fingerprint | `collapse_fingerprint.py` on CMA-ES arm |

**Success for benchmark claims** still requires RMSD ≤ 2.0 Å **and** PoseBusters
pass — RMSD-only is not enough (`AGENTS.md` / benchmarking skill).

---

## 6. Link to VALIDATION.md items A–G

| # | Item | What this package enables |
|---|------|---------------------------|
| **A** | `cmake configure` + `cmake/ValidateSources.cmake` with new TUs | CMake list includes `cmaes_search.cpp`; configure must list the new sources |
| **B** | Full engine + adapter **link** into one `FlexAIDdS` binary (`BUILD_FLEXAIDDS_FAST`) | Adapter TU in `flexaid_core` OBJECT lib |
| **C** | Adapter vs real `ic2cf` (not mock); snapshot fills for real ligand DOF | 5 seam fns resolve at link; doctor symbols for `FLEXAIDDS_SEARCH` / cmaes |
| **D** | Live dock one Astex complex: elected RMSD + best CF, **both** arms | Env gate A/B; same scoring path |
| **E** | GA-vs-CMA-ES A/B, eval-matched at **2e6** | §2 budget rule |
| **F** | Real entropy trace on rugged surface | §4 CSV → fingerprint |
| **G** | Locked-arch `.sif`, in-container dock, collapse fingerprint INVARIANT | Harness + `analysis/collapse_fingerprint.py` (chunk4/5) |

Close A–G **only** with on-disk artifacts under `validation_evidence/` per the
protocol in [`VALIDATION.md`](VALIDATION.md). A truthful OPEN with a real blocker
beats a fabricated CLOSED.

### Unblock order (from VALIDATION bottom line)

1. Land this package: `LIB/cmaes_search.{cpp,h}`, top/CMake wiring,
   `apply_integration.sh`, `CMAES_INTEGRATION.md` (this file).
2. Land `analysis/collapse_fingerprint.py` (chunk4).
3. Apptainer recipe + `.sif` (chunk5 / harness).
4. Re-run KICKOFF workers; fill A–G from new artifacts only.

---

## 7. Terminology (scoring vs thermodynamics)

Use precise language (`AGENTS.md`):

- **CF / contact-function scoring proxy** — what CMA-ES and GA optimize per eval
  (`get_cf_evalue` / VoronoiCF).
- **Ensemble-derived free energy estimate** / **thermodynamic ledger (F, H, −TS, Cv)** —
  StatMechEngine / BindingMode layers on top of the ensemble.
- Do **not** claim “computed true binding free energy ΔG” unless the full partition
  function + vibrational corrections (tENCoM) + solvent/concentration terms are
  active and validated.

---

## 8. File map after apply

```text
LIB/cmaes_search.h          # adapter API (chunk1): cmaes_run_dock / fill / write_trace
LIB/cmaes_search.cpp        # CMA-ES loop + 5-seam plumbing + entropy CSV
LIB/CMakeLists.txt          # + cmaes_search.cpp after gaboom.cpp
LIB/top.cpp                 # + include + FLEXAIDDS_SEARCH branch (markers)
analysis/                   # home for collapse_fingerprint.py (chunk4)
CMAES_INTEGRATION.md        # this document (repo root when installed)
```

### Adapter entry used by `top.cpp`

```cpp
// Chunk1 free functions (via cmaes_search.h)
int  cmaes_run_dock(FA*, GB*, VC*, genlim*, atom*, resid*, gridpoint*,
                    CmaesTargetFn, const CmaesConfig&, CmaesResult*,
                    std::vector<EntropyTraceSample>* = nullptr);
int  cmaes_fill_chromosomes(const CmaesResult&, int num_genes,
                            chromosome*, int max_chrom, gene* gene_storage);
void cmaes_write_trace_csv(const std::string&, const std::vector<EntropyTraceSample>&);
```

Config fields (chunk1): `population` (λ), `max_evals`, `seed`, `sigma0`,
`write_trace`, `enable_entropy_trace`, `archive_size`.  
Result fields: `best_cf`, `n_evals`, `n_gens`, `status`, `archive_*`.

Optional env: `FLEXAIDDS_CMAES_MAX_EVALS` overrides λ×gens when set to a positive integer.

Markers used by the on-ramp (idempotency):

| Marker | File |
|--------|------|
| `// FLEXAIDDS_CMAES_INCLUDE_BEGIN` / `END` | `LIB/top.cpp` |
| `// FLEXAIDDS_CMAES_BEGIN` / `END` | `LIB/top.cpp` |
| `# FLEXAIDDS_CMAES_CMAKE_BEGIN` / `END` | `LIB/CMakeLists.txt` |

On-ramp tools (swarm artifacts only until copied):

| Path | Role |
|------|------|
| `.swarm/cmaes/chunk3_onramp/artifacts/apply_integration.sh` | Installer |
| `.swarm/cmaes/chunk3_onramp/artifacts/_patch_top_cmaes.py` | top.cpp patch helper |
| `.swarm/cmaes/chunk3_onramp/artifacts/CMAES_INTEGRATION.md` | This document |

---

*Chunk 3 deliverable — integration on-ramp only. Adapter implementation is chunk1;
wiring review chunk2; fingerprint chunk4; harness chunk5; tests chunk6.*
