# Admission and metrics contract

**Canonical methodology:** [METHODOLOGY.md §0.3](../../METHODOLOGY.md#03-admission-identity-missingness-and-repair-evidence).
`AGENTS.md` governs coding conduct. This document maps that methodology to the
CSV interface enforced by `scripts/aggregate_claim_metrics.py`.

## 1. What this report establishes

The aggregator checks **receipt-field completeness and internal consistency**.
Every report carries `evidence_level=validated_receipt_fields` and
`artifacts_verified=false`. Matching SHA-256 strings do not establish that a
validator ran, that the described files still exist, or that those files hash
to the supplied values. Publication still requires the underlying run,
configuration, elected pose, receptor condition, reference, raw validator
outputs, and provenance review.

A stored `claim_ready=1` is a required producer attestation. It is never a
substitute for the measurements and receipts below. Existing historical CSVs
that lack required evidence cannot supply STRICT successes through this tool.
They may still supply explicitly labelled diagnostics when their protocol
eligibility fields are complete.

## 2. Separate metrics and populations

| Metric | Per-observation test | Role |
|---|---|---|
| S1 | Finite nonnegative elected in-place graph-symmetry RMSD at the §0 cutoff; otherwise labelled serial/legacy ordered fallback | RMSD-only diagnostic |
| S2 | S1 and `pb_pass=1`, with a valid elected pose SHA-256 matching `posebusters_pose_sha256` | RMSD/PB receipt diagnostic |
| STRICT | Recomputed serial RMSD, PB, score, validator and protocol receipt conjunction below | Primary **receipt-level** headline |
| S3 | Finite nonnegative `conditional_scanned_pool_ceiling`, else `best_cluster_rmsd` / `rmsd_bcr`, at the §0 cutoff | Conditional scanned-pool diagnostic; never any-pose success |

`rmsd_hungarian` cannot substitute for an elected S1 or STRICT measurement.
`rmsd_to_crystal` retains its serial identity-mapping meaning. Symmetry-corrected
S1 does not change the producer's serial STRICT gate. The ordering arguments in
§0.0 require the same atoms, graph, coordinates, and valid identity mapping;
they do not authenticate a foreign sidecar or an incompatible topology.

All metric numerators count **unique frozen-roster targets** using the declared
observation rule. Off-manifest rows never inflate diagnostic or STRICT rates.
S1/S2/S3 are evaluated on protocol-eligible rows independently of whether STRICT
passes. A tENCoM failure or false `claim_ready` therefore cannot erase a valid
RMSD diagnostic. `N_protocol_eligible` counts those observations; `N_claim` counts
individual STRICT receipt successes, while each metric's `n` counts targets.
Neither count should be confused with the other in multi-seed analyses.

No `success`, `success_s1`, `success_rmsd` or `success_s3` flag can replace a
missing/nonfinite RMSD. For S2, `success_pb` does not override `pb_pass`: it also
contains the producer's serial RMSD gate and is a different predicate from
symmetry-corrected S1 plus PB. Present contradictory PB run/count/failure fields
also defeat S2; independence from STRICT does not permit a diagnostic to ignore
contradictions in its own measurements.

## 3. Required receipt fields

Protocol eligibility requires every one of these:

- Explicit `seed_echo=0` and `native_pose_seeded=0`.
- Explicit `protocol_claim_eligible=1`.
- A syntactically valid row `matrix_md5` equal to the selected expected pin.
  A campaign/default pin alone is not evidence of which matrix a row used.
- If `native_pose_seed_fraction` is supplied, it must be finite zero; a blank fails.

STRICT additionally requires every one of these:

- `claim_ready=1`, `success_rmsd=1`, `success_pb=1`, `pb_pass=1`, and
  `score_pose_consistent=1`.
- Finite nonnegative serial `rmsd_to_crystal` within the §0 cutoff.
- Finite `score_pose_delta` within the producer's absolute tolerance, `1e-4`.
- `pb_backend=bust_cli`, `tencom_status=ok`, and `eigen_status=ok`.
- Syntactically valid 64-hex-character `pose_sha256`, `rmsd_pose_sha256`,
  `posebusters_pose_sha256`, and `tencom_pose_sha256`; the latter three must
  equal the elected pose anchor. Hexadecimal comparison is case insensitive.
- A valid `posebusters_input_sha256`. The PB input may be a converted SDF, so
  its hash is not required to equal the elected PDB hash.
- If `docking_completed` or `docking_exit_code` is present, it must explicitly
  assert completion or zero exit respectively. Missing runtime fields remain
  an identified limitation of older producer schemas; a present blank or
  unknown value cannot be treated as success.

STRICT requires `pb_ran=1`, `pb_n_checks=27`, `pb_n_pass=27`, `pb_n_fail=0`,
positive integer `eigen_n_modes`, and finite `elected_H_vib`. The count 27 is
pinned to the current mandatory `BustCli` check schema; a schema revision must
update the producer, aggregator, and regression controls together. Historical
summaries missing these fields are diagnostics, not full STRICT evidence.

When supplied, `num_poses` must be a positive integer, `pb_failed_keys` empty,
and `rmsd_fail_reason` exactly the success sentinel `none`. Present contradictory
or blank numeric/failure state cannot be overridden by status flags.

Every missing or contradictory STRICT field produces an explicit failure reason.
CLI, `RUN_RECEIPT.json`, and `provenance.json` matrix pins must agree when supplied;
malformed receipts and conflicting pins are errors.

## 4. Identity, denominator, and repeated seeds

A primary rate requires a frozen manifest with schema
`flexaidds.astex.target_manifest/v1`, containing unique uppercase target
codes, matching `N`, and a valid `sha256_of_sorted_codes`. The digest is SHA-256
of the UTF-8/ASCII **comma-joined sorted codes, without a terminal newline**.
The default is `astex85_target_manifest.json`; a different preregistered roster
must be selected explicitly using `--manifest`. An intact digest establishes
internal integrity, not the historical date or independence of preregistration.
Missing/corrupt manifests are errors, never an implicit smaller denominator.

One row per target is labelled `single_observation`. It is not evidence of
majority-of-seeds behavior. More than one observation per target requires
`--expected-seeds` with the prespecified seed list. Per METHODOLOGY §0.3,
strictly more than half of **all expected seeds** must pass; missing seeds do
not pass, and even-seed ties fail. The report lists absent expected target/seed
observations. Every target still contributes at most one numerator unit.

Seed identities are unsigned 64-bit decimal integers, canonicalized before
comparison: `1` and `01` are the same seed. Arbitrary labels and out-of-range
values are rejected. Legacy single-observation rows may omit a seed; repeated
observations may not.

Duplicate target/seed identities are errors, including byte-identical copies.
The tool never clips rates, silently deduplicates, or chooses the best replicate.
Mixed arm or endpoint values, including a mixture of missing and supplied values,
are rejected. `--arm` may explicitly select an arm, and the number of filtered
observations is reported. Separate endpoints must be supplied as separate sources.
Rows lacking all arm/endpoint labels can be analyzed as legacy single-observation
inputs with `unspecified` labels; those labels disclose incomplete protocol
identity and must be resolved against source receipts before publication.

Do not conflate unassessable measurements with physical invalidity. Missing PB
assessments remain named missing observations and cannot enter a STRICT numerator;
they do not shrink its frozen denominator. A separately reported PB-only assessed
subset must state its own assessability rule and count and is a distinct estimand.

## 5. CSV, source, and sidecar contracts

- CSV headers must be unique and nonblank. Every record must have the exact
  header width, with an unambiguous target identity. Malformed input is rejected.
- Per-target `*/result.csv` must contain exactly one observation. No hidden
  second row is silently discarded.
- A directory containing both per-target files and recognized summaries, or
  multiple recognized summaries, is ambiguous. Select a source with `--csv`.
  The tool does not infer source authority from filenames, ordering, or mtimes.
- A usable symmetry sidecar row must have `status=ok`, finite nonnegative RMSD,
  and a valid pose SHA-256. Joins require both target and pose hash. Duplicate
  sidecar target/pose identities and absent/mismatching pose hashes are errors.
  Multiple poses of one target can coexist only under distinct hashes.
- Live files are local-first. CloudDocs must be staged through
  `scripts/icloud_safe_io.py`; this aggregator refuses direct CloudDocs reads.

## 6. Interpretation of validators and election

PoseBusters-derived results are conditional on the receptor and solvent condition
actually supplied to the validator. Claim review must establish the engine-matched
condition from artifacts; the aggregator cannot infer it from `pb_pass`.
The canonical raw PB check set belongs to the versioned `BustCli` schema, not to
a predicate learned from whichever columns happened to be boolean in a corpus.
Blank, missing, or nonboolean check values cannot count as a pass.

Retain raw PB CSV and execution receipts before returning a schema failure.
Keep RMSD reference and PB input provenance separate if their instruments differ.
Historical solvent or aromaticity analyses are observations about their exact
files and invocation; their rates and denominators are not universal admission
constants.

The elected entropy/consensus top-1 and generator CF top-1 are separate estimands.
Preserve separate paths, scores, hashes, and RMSDs. Do not join restarts by rounded
CF value. S3 is the minimum over actually enumerated emitted heads/members; it
must never clear `seed_echo`, rewrite `pose_source`, or become an any-pose claim.

## 7. CLI and validation

```bash
python3 scripts/aggregate_claim_metrics.py --csv results.csv --json metrics.json
python3 scripts/aggregate_claim_metrics.py --csv results.csv --expected-seeds 12345,23456,34567
python3 scripts/aggregate_claim_metrics.py --csv results.csv --arm A --manifest frozen_targets.json
python3 scripts/aggregate_claim_metrics.py --csv results.csv --symmcorr elected_symmcorr.csv --headline s1
python3 scripts/aggregate_claim_metrics.py --csv legacy.csv --legacy-observed-denominator --diagnostic-only --headline s1
```

The explicit legacy denominator mode produces no STRICT rate. `--headline s3`
requires `--diagnostic-only`. Exit 2 means malformed input or contract violation;
exit 1 means no STRICT receipt successes. In diagnostic-only mode, nonempty
observations can return zero even when STRICT has no successes.

Regression gates: `tests/test_aggregate_claim_metrics.py`, the aggregator wiring
in `tests/test_rmsd_symmcorr.py`, and
`tests/p0_claim_contract/test_fixed_denominator.py`.
