# Softβ election policy (post pilot8)

**Comparative science hub:** [`COMPARATIVE_SCIENCE_README.md`](COMPARATIVE_SCIENCE_README.md)  
**Status:** production DatasetRunner Softβ S1 is **feature-flagged OFF by default**.  
**Shared math:** `LIB/SoftBetaFreeEnergy.h`  
**Parent contracts:** `docs/classic_entropy_ranking.md`, `docs/implementation/3dsig_shannon_ranking.md`

---

## 1. What Softβ is

Softβ ranks **already-clustered binding modes** with a soft free energy on the CF scoring proxy:

\[
\tilde G_i = \tilde H_i - T\,\tilde S_i
\qquad
\tilde H_i = \sum_{j\in i} p_j\,\mathrm{CF}_j
\qquad
\tilde S_i = -\sum_{j\in i} p_j\ln p_j
\qquad
p_j \propto e^{-\mathrm{CF}_j / T}
\]

Identity: \(\tilde G \equiv E_{\min}-T\ln Z_{\mathrm{local}}\) = cluster **ACF**.

| Claim language | Correct? |
|----------------|----------|
| “CF soft-β ranking proxy \(\tilde G\)” | Yes |
| “Mode election by local free energy on CF + Shannon frequencies” | Yes |
| “True thermodynamic binding ΔG with full solvent / concentration” | **No** (unless full ledger validated and claimed separately) |

Softβ is **not** a sampling method. It does not run GA, does not grow the ensemble, and **cannot** create ≤ 2 Å poses if no near-native exists among emitted heads (**BCR = 0**).

---

## 2. When Softβ helps vs when it cannot

| Situation | Softβ effect |
|-----------|----------------|
| Ensemble contains near-natives **and** a dense near-native basin has middling CF while a sparse false minimum has best CF | Softβ can **reorder** rank-0 toward the dense basin (1HNN-class) |
| All modes are singletons (S̃ = 0) | \(\tilde G = \mathrm{CF}\) → Softβ **identical** to CF rank-0 |
| **BCR = 0** (no head RMSD ≤ 2 Å) | Softβ can only permute losers → **cannot** create S1 ≤ 2 Å success |
| pilot8 Softβ re-rank of current heads | **Invalid success claim** — oracle BCR 0/8 means reordering cannot invent natives |

**Production gate (crystal-blind):** feature flag only.

**Offline diagnostic (not production):** `soft_beta::diagnostic_softbeta_can_help_s1(mode_rmsds)` asks whether a near-native exists so ablation tables can separate “sampling failure” from “election failure”. Crystal-gated election is scientifically invalid for blind claims.

---

## 3. TEMPER vs kT honesty

| Symbol | Meaning in FlexAID / Softβ |
|--------|----------------------------|
| \(T\) / `TEMPER` | Soft temperature in the engine CONFIG; \(\beta = 1/T\) over **CF arbitrary units** |
| \(k_B T\) in kcal/mol | **Not** used for Softβ / ACF / BindingMode classic ranking |
| Arm B `TEMPER 21` | Operator-optimized **engine** soft-T for FO clustering + ACF emission order — **not** “physical 21 K” and not kcal \(kT\) |
| DatasetRunner `soft_T` | Same soft-β scale when Softβ S1 is ON: env `FLEXAIDDS_ELECTION_SOFT_T` → dock TEMPER → 298 fallback |

Never report Softβ \(\tilde G\) differences as kcal/mol experimental free-energy differences unless the full StatMech + solvent ledger is active and labeled as such.

---

## 4. Feature flags (DatasetRunner S1)

| Env | Default | Meaning |
|-----|---------|---------|
| **`FLEXAIDDS_SOFTBETA_ELECTION`** | **0 (OFF)** | Preferred name: Softβ S1 over clustered heads |
| `FLEXAIDDS_ELECTION_SHANNON_F` | 0 | Legacy alias — same bit as Softβ S1 |
| `FLEXAIDDS_ELECTION_LEGACY_ZH` | 0 | Force Softβ OFF (legacy ZH / ≈CF path) |
| `FLEXAIDDS_ELECTION_SOFT_T` | 0 → resolve | Soft-β \(T\) override |

Either ON alias enables Softβ; `LEGACY_ZH=1` always forces OFF.

```bash
# Classic pilot / C0 claim harness — Softβ S1 OFF (default)
unset FLEXAIDDS_SOFTBETA_ELECTION FLEXAIDDS_ELECTION_SHANNON_F
# or explicitly:
export FLEXAIDDS_SOFTBETA_ELECTION=0

# Opt in only when reordering clustered modes is intentional:
export FLEXAIDDS_SOFTBETA_ELECTION=1
# Prefer dock TEMPER; override only for re-rank experiments:
# export FLEXAIDDS_ELECTION_SOFT_T=21
```

`scripts/run_C0_claim_clean.sh` defaults Softβ S1 to **0** (explicit opt-in required).

### Engine cluster ACF emission (E1b — Wave 0)

| Env | Default | Meaning |
|-----|---------|---------|
| **`FLEXAIDDS_ACF_STRICT`** | **0 (OFF)** | When `1`, `LIB/cluster.cpp` uses `soft_beta::free_energy_strict` (UniqueGeometry) for `Clus_ACF` instead of legacy `soft_beta::acf`. Exact-CF-duplicate members no longer deepen \(\tilde G\) via \(T\ln N\) multiplicity. |

Default OFF is **bit-identical** to the pre-E1b product path. Opt in only for pilots measuring election-gap targets (e.g. BCR≤2.5 Å with elected RMSD>2). See `docs/implementation/FORWARD_SUCCESS_RATE_PLAN.md` Wave 0–1.

---

## 5. Layers that must not be conflated

| Layer | What it does | Softβ? | Default |
|-------|--------------|--------|---------|
| GA search | Samples CF landscape | No | Always CF search |
| Arm **B** `TEMPER 21` + `CLUSTA FO` | Density modes + engine ACF emission when T>0 | **Engine** Softβ ACF (classic FlexAID product when T>0) | Arm B protocol |
| Arm **B0** `TEMPER 0` + `CLUSTA CF` | CF clustering / CF emission | Engine Softβ off (T=0) | Arm B0 |
| **DatasetRunner S1** Softβ election | Re-elect rank-0 across restarts via \(\tilde G\) | **Only if flag ON** | **OFF** |
| Live pilot B | TEMPER21 + FO — **not** DatasetRunner Softβ rescoring of frozen CF ensembles | Engine path | See red-pair protocol |

**B FO ≠ DatasetRunner Softβ S1.** Three-engine arm B is FO clustering + TEMPER-driven engine ranking. Do not document pilot B as “Softβ rescoring of CF ensembles” unless DatasetRunner Softβ was actually enabled and logged (`[SOFTBETA-ELECT] Softβ S1 ON`).

---

## 6. Gated election algorithm (production)

```
if FLEXAIDDS_SOFTBETA_ELECTION (or SHANNON_F) and not LEGACY_ZH:
    soft_T = ELECTION_SOFT_T or dock TEMPER or 298
    rank modes by Ĝ_i = soft_beta::mode_free_energy(mode_i, soft_T).G
    elect min Ĝ
    log [SOFTBETA-ELECT] Softβ S1 ON … (not sampling; cannot fix BCR=0)
else:
    do not claim Softβ improvement
    DatasetRunner keeps legacy ZH composite (ranking continuity; not Softβ)
    pure CF rank-0 helper: soft_beta::elect_cf_rank0 / elect_gated(..., false)
    log [SOFTBETA-ELECT] Softβ S1 OFF (default)
```

Crystal RMSD is **not** an input.

**Ranking-preservation note:** Default OFF path in DatasetRunner is still the pre-existing
legacy ZH composite (not a forced rewrite to pure min-CF), so enabling this policy
does not flip claim ranking vs prior Softβ-OFF builds. The pure CF helpers in
`SoftBetaFreeEnergy.h` are the clean gated API for tests and future explicit wiring.

Unit gates: `SoftBetaIdentity::*`, `SoftBetaGatedElection::*` in `tests/test_classic_entropy_ranking.cpp`; ProtocolConfig Softβ defaults in `tests/test_protocol_config.cpp`.

---

## 7. Explicit non-goals (post pilot8)

- Do **not** re-rank pilot8 heads and claim Softβ success when BCR = 0  
- Do **not** enable Softβ as default S1 for classic pilot  
- Do **not** dual-launch docks to “add Softβ sampling”  
- Do **not** claim true ΔG from Softβ \(\tilde G\) alone  

---

## 8. Code map

| File | Role |
|------|------|
| `LIB/SoftBetaFreeEnergy.h` | \(\tilde G\)/ACF math + `elect_gated` / diagnostic helpers |
| `LIB/cluster.cpp` | ACF emission when TEMPER>0 |
| `LIB/BindingMode.cpp` | Classic Softβ mode F when T>0 && !force_cf |
| `LIB/DatasetRunner.cpp` | S1 Softβ election if flag ON; `[SOFTBETA-ELECT]` logs |
| `LIB/ProtocolConfig.{h,cpp}` | Flag plumbing (default OFF) |
| `scripts/run_C0_claim_clean.sh` | Softβ S1 default 0 |
| `docs/classic_entropy_ranking.md` | Engine product path |
| `docs/implementation/3dsig_shannon_ranking.md` | 3Dsig ranking contract |
