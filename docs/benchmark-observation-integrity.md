# DatasetRunner scorer identity

The canonical measurement contract is [METHODOLOGY.md](../METHODOLOGY.md).

The corrected default `--ranking-objective cf_minus_ts` computes the explicitly
named CF-T*S reranking proxy using each pose's emitted temperature. It is
neither the ensemble free-energy ledger nor a true binding free energy.
`--ranking-objective legacy_cf_minus_s` is an explicit, dimensionally incorrect
diagnostic objective for comparison. Nonzero S without an emitted finite,
positive T refuses the election. Zero S leaves CF unchanged and requires no
invented temperature.

Pose records preserve CF, raw S, emitted T, T*S, ensemble mean energy, ensemble
free energy, soft_beta_G, emitted TOTAL_SCORE, generator rank and original
binding-mode ID separately. Generator rank follows the final numeric output
filename suffix; cluster identity may disagree after DensityPeak sorting.
Missing elected output remains a generator miss rather than promoting a later
rank. The generator_docking_power_top1 endpoint follows that election; the
entropy_reranked_docking_power_top1 and historical docking_power_top1 keys use
the declared Python proxy. Synthetic dry runs suppress all docking-power keys.

Regression fixtures establish temperature dependence, a scalar ranking reversal,
unchanged zero-entropy scores, numeric generator election despite disagreeing
cluster IDs, and refusal to substitute a later output for missing rank zero.
These are software checks, not measured molecular-performance changes.
