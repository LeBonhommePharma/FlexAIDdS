# Performance swarm baselines

Versioned JSON artifacts for cross-platform regression checks. See
`docs/PERFORMANCE_LINUX_BASELINE_RUNBOOK.md` for population instructions.

## Files

| File | Platform | Contents |
|------|----------|----------|
| `baseline_macos_metal.json` | macOS Metal | 182 harvested dock timings (Wave 1) |
| `baseline_linux_cpu.json` | Linux AVX2+OpenMP | Microbench stub + empty dock harvest |
| `baseline_linux_cuda.json` | Linux CUDA+AVX2 | Microbench stub (GPU runner required) |

## Schema (`schema_version` 1.0.0)

Shared top-level fields:

- `recorded_at` — ISO-8601 UTC timestamp
- `git` — `{commit, dirty, branch}`
- `platform` — label string (`macos_metal`, `linux_cpu`, `linux_cuda`)
- `cmake_flags` — configure cache used for the run
- `notes` — human-readable provenance

### Microbenchmarks (`benchmarks[]`)

Used by `scripts/compare_perf_baseline.py` and `.github/workflows/perf.yml`:

```json
{
  "name": "tencom",
  "reference_n_res": 200,
  "metrics": { "build_ms_full": 12.3, "sample_ms_full": 4.5 }
}
```

```json
{
  "name": "vcfbatch",
  "args": [200, 20],
  "metrics": { "speedup_vs_scalar": 3.2 }
}
```

### Dock timings (`dock_timings_harvested`)

Used by `scripts/compare_dock_timings.py` and `scripts/harvest_perf_baselines.py`:

```json
{
  "count": 182,
  "records": [
    {
      "campaign": "v128",
      "target": "1G9V",
      "job_key": "v128/1G9V",
      "entry_key": "1G9V",
      "log_path": "/path/to/stderr.log",
      "avg_ms_per_gen": 142.5
    }
  ]
}
```

Harvest parses `TIMING SUMMARY: N gens timed, avg X ms/gen` from `stderr.log`.