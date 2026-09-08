# DatasetRunner observation and scorer identity

The canonical measurement contract is [METHODOLOGY.md](../METHODOLOGY.md).
These notes describe the BENCH-1/2/3 implementation, not docking-success claims.

`docking_power_top1` and `docking_power_top3` retain their target-level unit.
Their denominator is the declared target roster, including targets with missing
inputs, empty output or failed execution. Reports preserve the target/state
roster and distinguish absent-output observations from engine failures. Mixed
states for the same target are refused before execution; select one state per
run with `run_dataset(..., structural_states=[state])` rather than pool states.
Bootstrap intervals resample declared targets, including empty observations,
and retain each target's full pose election.

For a frozen denominator use the existing manifest explicitly:

```text
--expected-target-manifest benchmarks/protocols/astex85_target_manifest.json
```

The manifest's dataset, N, unique identities and declared sorted-code digest are
validated. Targets absent from execution remain failures in the denominator.
Both manifest file-byte SHA256 and sorted-code SHA256 are recorded. The input
manifest is read without alteration; no second roster is created.

`--resume` restores every schema-2 cached observation, including recorded
failures, and recomputes metrics and gates over fresh plus restored evidence.
Early-termination records, failure exit codes and timeouts survive resume.
Checkpoint compatibility binds the selected roster, dataset configuration,
engine hash, Python scorer source hashes, temperature, objective, relevant
runtime environment and input file hashes. Old, malformed or incompatible
checkpoints fail before dispatch. Baseline changes may reuse compatible
observations but always recompute the verdict. To retry a previously failed
observation or change protocols, use a new results namespace; do not erase or
reinterpret frozen evidence.
Manifest labels use canonical target IDs; checkpoint discovery retains the
scheduled ID spelling, including on case-sensitive filesystems.

The corrected default `--ranking-objective cf_minus_ts` computes the explicitly
named **CF-T*S reranking proxy** using each pose's emitted temperature. It is
neither the ensemble free-energy ledger nor a true binding free energy.
`--ranking-objective legacy_cf_minus_s` is an explicit, dimensionally incorrect
diagnostic objective for comparison. Nonzero S without an emitted finite,
positive T refuses the election. Zero S leaves CF unchanged and requires no
invented temperature.

Pose records preserve CF, raw S, emitted T, T*S, ensemble mean energy, ensemble
free energy, `soft_beta_G`, any emitted `TOTAL_SCORE`, and generator rank as
separate fields. Generator rank follows the final numeric output filename
suffix. `binding_mode` is retained as cluster identity: DensityPeak keeps its
original IDs after sorting the output election. An absent elected rank remains
a generator miss, even if a later output is near-native.
`generator_docking_power_top1` follows the emitted election;
`entropy_reranked_docking_power_top1` and the historical `docking_power_top1`
key use the declared Python proxy objective. Synthetic dry runs suppress all
docking-power endpoints. No RMSD-only endpoint is a validated docking claim.

The repair's regression fixtures establish temperature dependence, a scalar
ranking reversal, unchanged zero-entropy scores, fixed denominator accounting,
identical fresh/full/partial-resume verdicts, and refusal of incompatible caches.
They do not establish a change in real molecular performance.
