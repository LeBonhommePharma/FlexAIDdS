# DualAssembly Cotranscriptional/Cotranslational Docking

**Scope.** This document describes how FlexAIDdS solves the moving-target / moving-ligand
problem that arises when a nascent polypeptide chain (or nascent RNA transcript) and a
protofibril simultaneously evolve in vivo. It introduces an entropy-driven target/ligand
discriminator, three parallel simulation modes, and the temporal schedule that maps a
docking trajectory onto human-cell elongation rates. The Aβ42 protofibril (PDB 5OQV or
2NAO) is used as the canonical target for the validation track. The C++ implementation
lives under `LIB/NATURaL/` and is exposed through the `dual_assembly` CLI and
`scripts/run_dual_assembly_cotranslational.sh`.

---

## 1. Motivation: the Target/Ligand assignment problem

Classical FlexAID assumes a static receptor and a flexible ligand. When both the receptor
(a growing protofibril surface) and the ligand (a nascent chain emerging from the
ribosome) are non-stationary in conformation and in chain length, the classical assignment
is undefined. Picking one by convention (e.g. always treat the protofibril as the target)
hides the question rather than answering it: at early checkpoints the nascent chain has no
defined pocket, and at late checkpoints both partners can carry well-defined surfaces.

We need a rule that:

  1. is grounded in physics, not convention;
  2. is testable against a reciprocal control;
  3. is computable from quantities the existing engine already emits.

Shannon configurational entropy over the binned pose distribution satisfies all three. The
codebase already reports entropy in nats (`LIB/ShannonThermoStack/ShannonThermoStack.h:34`)
and ships two named thresholds (`kHSC_soft_nats = 2·ln 2 ≈ 1.386` and
`kHSC_hard_nats = ln 2 ≈ 0.693`).

---

## 2. Shannon-driven pose-collapse role hint

Let an ensemble of `N` GA-sampled poses be binned into `K` mega-clusters with normalised
occupancies `p_i` (Σ p_i = 1, i = 1…K). The Shannon entropy in nats is

  H(X) = −Σᵢ pᵢ ln pᵢ,  0 ≤ H(X) ≤ ln K.

Let `H_A` denote the entropy of Sim A's backend-projected pose-coordinate ensemble at the
current checkpoint and `H_B` that of Sim B. These are **not** intrinsic isolated-chain or
protofibril conformational entropies. The discriminator is therefore a reproducible
pose-collapse role hint, not a proof that one physical object is always "the receptor":

  • **Hard regime** — `min(H_A, H_B) < ln 2`. The lower-entropy system carries a single
    cluster with > 50% probability mass. The lower-entropy simulation direction is the
    primary role hint for that checkpoint.
    Confidence weight `w = 1`.

  • **Deferred regime** — `min(H_A, H_B) ≥ 2·ln 2`. Both ensembles have effective support
    of ≥ 4 clusters. Neither pose ensemble has collapsed; **no assignment**. The schedule either
    extends the GA pass (more generations until one falls below `kHSC_soft_nats`) or
    inherits the prior checkpoint's assignment with the flag `tl_deferred = true`.

  • **Soft regime** — `ln 2 ≤ min(H_A, H_B) < 2·ln 2`. The lower-entropy system is the
    *probable* target, with assignment weight

      w = (kHSC_soft_nats − H_lower) / (kHSC_soft_nats − kHSC_hard_nats) ∈ [0, 1].

    `w` is logged into the trajectory CSV and propagated downstream as a confidence on the
    A-or-B identity at that timepoint.

### 2.1 Justification

  • `H = 0`  ⇔ a single pose-coordinate bin carries all the mass. This is a collapsed
    docking ensemble, not automatically a rigid isolated molecule.

  • `H < ln 2`  ⇒ p_max > 0.5 (since for the worst case with H = ln 2 and a 2-cluster
    distribution we have p_max = 1/2 exactly; lowering H below ln 2 requires p_max > 1/2).
    A single cluster dominates: the simulation direction has a defined pose basin.

  • `H ≥ 2·ln 2 = ln 4` ⇒ the entropy is at least that of a 4-cluster equi-probable
    distribution. The simulated pose ensemble remains diffuse.

  • Units: the engine returns nats. Comparisons use `kHSC_*_nats` constants. Bits are
    *only* emitted at user-facing reporting boundaries, per the header policy.

### 2.2 Falsification — the reciprocal control track

Without a control, the discriminator above is a definition rather than a test. NATURaL
already emits a `DockingRolePolicy::ReciprocalControl` (see
`LIB/NATURaL/NATURaLDualAssembly.h:103`). The runner generates **both** the forward
(`ProtofibrilAsTarget`) and reciprocal control tracks. The reciprocal track now swaps the
actual Sim A target/ligand arguments instead of only changing the label. If the
pose-collapse role hint remains stable under that swap, the target/ligand conclusion is
stronger; if it flips, downstream analysis must treat the checkpoint as direction-sensitive.
The current runner records both tracks but does not yet compute a cross-track
`tl_assignment_inconsistent` aggregate flag.

---

## 3. Three parallel simulation modes per checkpoint

At checkpoint `k` of chain length `L_k`, run the following GA-based docking simulations,
each on a separate OMP thread group:

  • **Sim A — canonical in-vivo.** Target = full protofibril (5OQV/2NAO), fixed. Ligand =
    nascent chain truncated at `L_k` (extended geometry MVP; ESMFold-aware in the
    post-MVP). Produces ΔG_A(L_k), H_A(L_k), best pose seed for `L_{k+1}`.

  • **Sim B — reciprocal control.** Target = nascent chain at `L_k`. Ligand = single Aβ42
    monomer. Run **only if** `H_chain(L_{k-1}) < kHSC_soft_nats` — i.e. only when the
    chain has actually collapsed to a pocket-bearing state at the prior checkpoint.
    Otherwise B is skipped, ΔG_B = NaN, H_B = +∞.

  • **Sim C — fibril elongation oracle.** Target = protofibril (fixed), Ligand = free
    monomer from solution. Runs every `K` checkpoints (default `K = 5`). Output is fed to
    a `FibrilGrowthOracle` which wraps `target::GrandPartitionFunction` and returns

      Ξ = 1 + z·Z,  z = c_monomer / c°,  c° = 1 M (IUPAC convention),
      p(elongation) = z·Z / Ξ = 1 − 1/Ξ,
      ΔG_elong      = −kT ln Z + kT ln(c°/c_monomer).

    Default `c_monomer = 1 µM` (upper bound for physiological cytosolic free Aβ42; see
    Roher 1996, Bjorkdahl 2008). Gates whether the runner advances the protofibril
    state. When a structural advance callback is installed, an accepted gate returns the
    next protofibril PDB path and subsequent checkpoints use that structure. Without that
    callback, the CSV still records the logical state transition but
    `protofibril_structure_updated = 0`, so no coordinate mutation is claimed.

All three simulations write a `StatMechEngine` per Sim; the trajectory CSV aggregates
their `compute().free_energy` and the Shannon entropy over backend-supplied
pose-coordinate samples (`GAResult::pose_rmsds`).

---

## 4. Temporal schedule — human cell

Using `LIB/NATURaL/RibosomeElongation.h` constants:

  • Translation: 5.6 aa/s → 1 aa per 178 ms (`MEAN_EL_RATE_HUMAN`).
  • Transcription: 25 nt/s → 1 nt per 40 ms (`MEAN_NT_RATE_HUMAN`).
  • Ribosomal exit tunnel: 34 aa (`TUNNEL_LENGTH_AA`). No simulation before `L_k ≥ 40`.
  • Chaperone/NAC occlusion: skip docking while `L_k − TUNNEL_LENGTH_AA < 6`.
  • Eukaryotic transcription is nuclear: `direct_encounter_allowed = false` for
    transcription tracks (`LIB/NATURaL/NATURaLDualAssembly.cpp:184..188`). The
    transcription clock is retained as an internal synchronisation reference (a polysome
    may load before splicing completes), but transcription poses do not dock against the
    protofibril directly. The transcription tracks still contribute their `H` to the
    discriminator when the user enables polysome-mode in a future extension; for now their
    Sim A is skipped and only Sim B/C remain, recording an "ineligible-for-encounter"
    flag.

For a 300-aa target sequence with `--checkpoint-interval 10`:
  – first translation checkpoint at `L = 40`, scheduled time
    `t_40 = 1/k_ini + 40/5.6 ≈ 17.14 s` with the default human initiation rate
    (`k_ini = 0.1 s⁻¹`),
  – successive checkpoints every 1.78 s real-time,
  – wall-time is backend-dependent; the shipped CLI is synthetic and must not be used as
    production docking evidence.

---

## 5. C++ class design

All three new types live in `namespace natural` and reuse existing primitives without
modification.

### 5.1 `NascentChainScheduler`

Owns the ordered list of checkpoints derived from
`natural::build_parallel_growth_schedule(tracks)`. Tracks are produced via
`natural::make_human_protofibril_tracks(transcript_nt, peptide_aa,
include_reciprocal_controls=true)`.

```cpp
struct Checkpoint {
    int                 idx;            // 0-based checkpoint index
    int                 L_k;            // chain length (aa or nt)
    double              t_arrival_s;    // real-time since initiation
    GrowthProcess       process;        // Transcription | Translation
    DockingRolePolicy   role_policy;    // ProtofibrilAsTarget | ReciprocalControl
    int                 track_index;    // index into the tracks vector
    std::string         track_name;     // pretty name
    bool                in_tunnel;      // true if L_k ≤ tunnel_length
    bool                chaperone_shielded; // L_k - tunnel < 6
};

struct CheckpointOutcome {
    double H_A_nats     = std::numeric_limits<double>::infinity();
    double H_B_nats     = std::numeric_limits<double>::infinity();
    double dG_A_kcal    = std::numeric_limits<double>::quiet_NaN();
    double dG_B_kcal    = std::numeric_limits<double>::quiet_NaN();
    char   tl_primary   = '?';   // 'A' = chain-as-target, 'B' = protofibril-as-target, '?' = deferred
    double tl_weight    = 0.0;   // confidence ∈ [0,1]
    bool   tl_deferred  = true;
    std::string tl_basis = "pose_entropy_heuristic";
    double p_elong      = std::numeric_limits<double>::quiet_NaN();
    double dG_elong     = std::numeric_limits<double>::quiet_NaN();
    bool   sim_c_evaluated = false;
    bool   sim_c_gated_in = false;
    int    protofibril_state_index = 0;
    bool   protofibril_structure_updated = false;
};

class NascentChainScheduler {
public:
    NascentChainScheduler(std::vector<InVivoAssemblyTrack> tracks, int checkpoint_interval);

    bool has_next() const;
    Checkpoint next();                                  // advance the cursor
    void record(const Checkpoint& ck, CheckpointOutcome out);
    const std::vector<std::pair<Checkpoint, CheckpointOutcome>>& history() const;
};
```

### 5.2 `FibrilGrowthOracle`

Thin wrapper around `target::GrandPartitionFunction`. One instance per runner, reused
across all Sim C invocations.

```cpp
struct ElongationDecision {
    double p_elong;    // ∈ [0, 1)
    double dG_elong;   // kcal/mol (negative = favourable)
    bool   gated_in;   // p_elong ≥ acceptance_threshold
};

class FibrilGrowthOracle {
public:
    explicit FibrilGrowthOracle(double temperature_K = 310.15,
                                double acceptance_threshold = 0.5);

    // Returns the elongation decision for the current monomer ensemble.
    // monomer_engine: StatMechEngine populated by Sim C.
    // c_monomer_M:    free monomer concentration in molar (default 1e-6 M = 1 µM)
    ElongationDecision gate(const statmech::StatMechEngine& monomer_engine,
                            double c_monomer_M = 1.0e-6);

private:
    double T_K_;
    double acceptance_threshold_;
};
```

### 5.3 `DualAssemblyRunner`

Top-level driver. Coordinates Sim A/B/C across NATURaL tracks. Returns a vector of
`(Checkpoint, CheckpointOutcome)` and writes a CSV.

```cpp
struct DualAssemblyConfig {
    std::string   protofibril_pdb;         // required
    std::string   monomer_pdb;             // required if sim_c_enabled
    std::string   sequence_fasta;          // required, residues encoded 1-letter
    int           transcript_nt = 0;       // 0 = derive from sequence × 3 + UTR (post-MVP)
    int           checkpoint_interval = 10;
    int           sim_c_interval = 5;
    double        monomer_conc_M = 1.0e-6;
    double        temperature_K = 310.15;
    int           n_threads = 6;
    bool          include_reciprocal_controls = true;
    bool          sim_c_enabled = true;
    std::string   output_csv = "cotranslational_trajectory.csv";
    std::string   nascent_pdb_dir = ".";   // where per-L_k nascent PDBs are dumped
};

class DualAssemblyRunner {
public:
    DualAssemblyRunner(DualAssemblyConfig cfg,
                       SimAFn sim_a,
                       SimBFn sim_b,
                       SimCFn sim_c,
                       TruncateFn truncate,
                       ProtofibrilAdvanceFn advance_protofibril = nullptr);

    // Run the full trajectory. Returns the per-checkpoint history.
    std::vector<std::pair<Checkpoint, CheckpointOutcome>> run();

private:
    struct TrackState {
        double H_prev_chain_nats;
        char   tl_prev;
        int    checkpoints_since_sim_c;
    };

    std::unordered_map<int, TrackState> track_state_;
    std::string current_protofibril_pdb_;
    int         protofibril_state_index_;

    // Pose-entropy role hint; updates outcome.tl_primary, tl_weight, tl_deferred.
    static void assign_tl(double H_A_nats, double H_B_nats,
                          char prev_tl, CheckpointOutcome& out);
};
```

### 5.4 Pseudocode

```
for ck in scheduler:
    state = track_state[ck.track_index]
    if ck.chaperone_shielded:
        scheduler.record(ck, outcome={ all-NaN, tl_deferred=true })
        continue

    if ck.role_policy == ReciprocalControl:
        (dG_A, H_A) = sim_A(target=nascent_chain, ligand=current_protofibril)
    else:
        (dG_A, H_A) = sim_A(target=current_protofibril, ligand=nascent_chain)

    H_B = +inf; dG_B = NaN
    if state.H_prev_chain < kHSC_soft_nats:
        (dG_B, H_B) = sim_B(target=nascent_chain, ligand=monomer)

    decision = {NaN, NaN, false}
    if (++state.checkpoints_since_sim_c) >= sim_c_interval:
        engine_C  = sim_C(target=current_protofibril, ligand=monomer)
        decision  = oracle.gate(engine_C, cfg.monomer_conc_M)
        if decision.gated_in:
            protofibril_state_index += 1
            if advance_protofibril:
                current_protofibril = advance_protofibril(current_protofibril, monomer)
        state.checkpoints_since_sim_c = 0

    assign_tl(H_A, H_B, state.tl_prev, outcome)

    scheduler.record(ck, outcome={H_A,H_B,dG_A,dG_B,...,decision.p_elong,decision.dG_elong})
    state.H_prev_chain = H_A
    state.tl_prev      = outcome.tl_primary
```

---

## 6. CSV schema

`cotranslational_trajectory.csv` columns:

```
checkpoint_idx, L_k, process, track_name, role_policy, t_arrival_s, in_tunnel,
H_A_nats, H_B_nats, dG_A_kcal, dG_B_kcal, tl_primary, tl_weight, tl_deferred,
tl_basis, p_elong, dG_elong_kcal, sim_c_evaluated, sim_c_gated_in,
protofibril_state_idx, protofibril_structure_updated
```

21 columns. `tl_basis` is currently `pose_entropy_heuristic`. NaN is written as the literal string `NaN`. The CSV is appended-as-you-go so
a long trajectory is partially readable on disk even if the run is interrupted.

---

## 7. Integration with existing engine

  • `LIB/gaboom.cpp:89 calculate_fitness(...)` — black-boxed; one full GA pass per Sim.
  • `LIB/ShannonThermoStack/ShannonThermoStack.h:96 compute_shannon_entropy(values, bins)`
    — invoked on backend-supplied pose-coordinate samples from each GA result.
  • `LIB/statmech.h::StatMechEngine` — one engine per Sim; `compute().free_energy` feeds
    ΔG into the CSV; `log_sum_exp` reused via `GrandPartitionFunction::log_Xi`.
  • `LIB/GrandPartitionFunction.h::GrandPartitionFunction` — one instance per oracle. Two
    entries: the empty site (implicit "1") and `"monomer"` updated each Sim C via
    `add_or_overwrite("monomer", log_Z_monomer, c_monomer_M)`.
  • `LIB/NATURaL/{NATURaLDualAssembly,RibosomeElongation}.{h,cpp}` — unchanged.
  • `LIB/DatasetRunner.{h,cpp}` — not touched. The cotranslational runner is a separate
    top-level driver to avoid cross-cutting changes.

---

## 8. Validation strategy

  1. **Unit tests.** `tests/test_nascent_chain_scheduler.cpp` and
     `tests/test_dual_assembly_runner.cpp`. Cover:
       – schedule ordering against `build_parallel_growth_schedule` (deterministic);
       – T/L discriminator at the three entropy regimes;
       – `FibrilGrowthOracle` against analytic `Ξ = 1 + z·Z` for known z, Z;
       – chaperone-shield gate;
       – tunnel-exclusion (no Sim A while `L_k < 40`).

  2. **Integration smoke test.** A two-layer β-sheet stub PDB shipped under
     `tests/data/fake_protofibril_stub.pdb` plus a 60-aa synthetic sequence; the CSV must
     emit exactly 3 checkpoints (L = 40, 50, 60), 21 columns, no NaN in `H_A_nats`,
     finite `dG_A`.

  3. **Scientific validation (Aβ42).** User downloads PDB 5OQV (or 2NAO), prepares it
     with FlexAID's standard target-prep workflow, then runs the full
     `DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA` (42-aa Aβ42) sequence with
     `--checkpoint-interval 5` (this is a short sequence — finer resolution). The
     falsifiable claim:

       ΔG_A(L_k) becomes monotonically favourable past a critical L_k that matches the
       experimental nucleation length reported by ThT fluorescence kinetics (Knowles 2009;
       Cohen 2013, PNAS 110:9758; Meisl 2014, Nat. Protoc. 11:252).

  4. **T/L falsification.** The Shannon-driven role hint should remain stable when the
     forward and reciprocal-control tracks are flipped. Mismatches are interpreted from the
     paired rows; an explicit aggregate inconsistency flag is not implemented yet.

---

## 9. Known limitations

  • The MVP truncates the nascent chain in extended geometry only. ESMFold integration for
    per-residue secondary-structure-aware truncation is post-MVP.
  • The transcription tracks do not dock directly against the protofibril
    (`direct_encounter_allowed = false`). Polysome-mode (multiple nascent chains at
    staggered lengths sharing one mRNA against a single protofibril) is a future
    extension.
  • The Sim C occupancy probability uses `GrandPartitionFunction` fugacity
    (`z = c/c°`). The reported `dG_elong_kcal` is the apparent concentration-corrected
    value `−kT ln Z + kT ln(c°/c)`, not the intrinsic `F_bound` alone.
  • Debye screening of electrostatics at 150 mM KCl (κ⁻¹ ≈ 7.8 Å) is a flag, not a
    default. FlexAID's `Vcontacts` is short-range and largely unaffected.
  • The fibril-growth oracle assumes a single monomer-binding site per checkpoint and a
    well-mixed monomer pool. Surface-bound oligomeric intermediates (e.g. dodecameric
    Aβ*56) are not represented.

---

## 10. References

  • Hessa et al. 2007 Nature 450:1026 — TM-helix code.
  • Zhao et al. 2011 J. Phys. Chem. B 115:3987 — ribosome master equation.
  • Wohlgemuth et al. 2008 — E. coli elongation rate.
  • Ingolia et al. 2011 — human elongation rate.
  • Pechmann & Frydman 2013 Nat. Struct. Mol. Biol. 20:237 — translational pausing.
  • Knowles et al. 2009 Science 326:1533 — amyloid aggregation kinetics.
  • Cohen et al. 2013 PNAS 110:9758 — Aβ42 secondary nucleation.
  • Meisl et al. 2014 Nat. Protoc. 11:252 — ThT-based mechanistic analysis.
  • Bjorkdahl et al. 2008 — Aβ42 cytosolic concentrations.
