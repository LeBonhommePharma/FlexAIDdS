# Bonhomme Fleet Benchmarks

Bonhomme Fleet is the resumable control plane for `benchmark_datasets`. It is
designed so Codex, Grok Build, Claude Code, Claude Science, or a human can
inspect and resume the same campaign without owning the original terminal.

Fleet does not change docking, ranking, PoseBusters, or tENCoM/Eigen. It pins
those executables and the protocol, leases deterministic target shards, runs
on local storage, verifies the archived tree in iCloud, and writes the accepted
chunk manifest last.

## Scientific Contract

- Production Astex redocking uses `defined-cleft-redock`: the cognate GetCleft
  site is supplied, but crystal ligand coordinates are not a GA seed.
- Official PoseBusters (`bust`) and built-in tENCoM/Eigen are mandatory.
- RMSD, PoseBusters, and tENCoM/Eigen must cite the same elected-pose SHA-256.
- `success_rmsd`, `success_pb`, and `claim_ready` are separate. The primary
  strict rate is `success_pb = RMSD <= 2.0 A AND PoseBusters pass`.
- `claim_ready` additionally requires tENCoM/Eigen, exact-pose provenance,
  score/pose consistency, and a claim-eligible no-seed protocol.
- `best_cluster_rmsd_a` is an any-pose sampling ceiling, not top-1 success.
- A completed process is not necessarily a scientific success. Fleet reports
  execution and scientific counters separately.

## M3 Pro Setup

From the repository root:

```bash
ROOT="$(git rev-parse --show-toplevel)"
BUILD="$ROOT/build_fleet"

xcrun -f metal
xcrun metal -v
cmake -S "$ROOT" -B "$BUILD" \
  -DBUILD_TESTING=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DFLEXAIDS_USE_METAL=ON
cmake --build "$BUILD" --target FlexAIDdS benchmark_datasets test_fleet_runner -j 6
ctest --test-dir "$BUILD" -R 'FleetRunnerTests|DatasetRunnerTests' --output-on-failure
```

Resolve the required validator and iCloud destination:

```bash
source "$HOME/.flexaidds_env"
: "${FLEXAIDDS_RESULTS:?set FLEXAIDDS_RESULTS to the iCloud results directory}"

BUST="${FLEXAIDDS_POSEBUSTERS_BIN:-$ROOT/.venv-posebusters/bin/bust}"
test -x "$BUST"
test -d "$ROOT/benchmarks/astex_diverse/astex_diverse"
df -h / "$FLEXAIDDS_RESULTS"
```

The tracked Astex input manifest uses paths relative to its own location. It
therefore works in any checkout and contains no username-specific path.

## Plan One Astex Campaign

Create one deterministic shard per Astex target. This makes a failed target
independently retryable and prevents a partial 85-target process from being
mistaken for a completed campaign.

```bash
CODES="/private/tmp/astex-diverse-85.txt"
sed -n '/^targets:/,/^metrics:/p' "$ROOT/benchmarks/datasets/astex_diverse.yaml" \
  | awk '/^  - / { print toupper($2) }' > "$CODES"
test "$(wc -l < "$CODES" | tr -d ' ')" = 85

RUN_ID="astex-defined-cleft-$(date -u +%Y%m%dT%H%M%SZ)"
CAMPAIGN="$FLEXAIDDS_RESULTS/campaigns/$RUN_ID"

python3 "$ROOT/python/flexaidds/fleet.py" plan "$CAMPAIGN" \
  --campaign-id "$RUN_ID" \
  --runner "$BUILD/benchmark_datasets" \
  --engine "$BUILD/FlexAIDdS" \
  --posebusters-bin "$BUST" \
  --benchmark "crossdock_json:$ROOT/benchmarks/datasets/benchmark_astex_native_85.json" \
  --dataset astex-diverse \
  --mode defined-cleft-redock \
  --codes-file "$CODES" \
  --chunks 85 \
  --threads 1 \
  --omp-threads 4 \
  --ga-population 1000 \
  --ga-generations 6000 \
  --temperature 298 \
  --job-timeout-seconds 10800 \
  --min-free-gb 5 \
  --env FLEXAIDDS_RESTARTS=5 \
  --env FLEXAIDDS_PARALLEL_RESTARTS=0 \
  --env FLEXAIDDS_VCT_R0=4 \
  --env SHARING_ALPHA=4 \
  --env EVAL_SCALE_DIHEDRAL=-1
```

`manifest.json` and `manifest.sha256` are immutable. Fleet also pins the SHA-256
of `benchmark_datasets`, `FlexAIDdS`, and `bust`. If any binary changes, resume
stops instead of silently mixing methods.

## Run, Monitor, and Resume

Use one worker on the 18 GB M3 Pro. The compute tree stays under
`/private/tmp/flexaidds_fleet`; only a verified archive is committed to iCloud.

```bash
caffeinate -i -s python3 "$ROOT/python/flexaidds/fleet.py" run "$CAMPAIGN" \
  --worker-id m3pro-primary
```

Any agent can monitor without claiming work:

```bash
python3 "$ROOT/python/flexaidds/fleet.py" status "$CAMPAIGN"
df -h / "$FLEXAIDDS_RESULTS"
```

If the worker or machine stops, run the same command with `resume`. Live leases
are left alone; stale leases are fenced with a higher epoch and the old attempt
cannot publish an accepted result.

```bash
caffeinate -i -s python3 "$ROOT/python/flexaidds/fleet.py" resume "$CAMPAIGN" \
  --worker-id m3pro-primary
```

Do not use `--force`, delete claims, edit the manifest, or launch another
campaign into the same directory.

## Aggregate and Report

```bash
python3 "$ROOT/python/flexaidds/fleet.py" aggregate "$CAMPAIGN"
cat "$CAMPAIGN/aggregate/summary.json"
```

Long-term claim artifacts are:

- `manifest.json` plus `manifest.sha256`: protocol and executable pins.
- `attempts/<chunk>/<attempt>.json`: immutable terminal attempt history.
- `artifacts/<chunk>/<attempt>/`: verified iCloud copy of logs, poses, validator
  files, reports, and the C++ Fleet chunk result.
- `chunks/<chunk>/result.json`: accepted commit manifest, written last.
- `aggregate/targets.csv`: one row per elected pose.
- `aggregate/summary.json`: strict `success_pb` and `claim_ready` rates.

The dashboard server is local-only by default:

```bash
python3 "$ROOT/benchmarks/m3pro/dashboard/fleet_status_server.py" \
  --host 127.0.0.1 --port 8787
```

Cross-origin browser access requires an explicit exact `--cors-origin`. Fleet
does not currently add application-level encryption; confidentiality is the
responsibility of local disk and iCloud account protection.
