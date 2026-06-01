# M3 Pro Benchmark Setup — iCloud 2TB Only

Benchmark environment for **MacBook Pro 14" M3 Pro 18GB** with zero local SSD
space. All data lives exclusively on iCloud Drive (2TB). This is the sole storage
for results, logs, and build artifacts.

## Quick Start

```bash
# 1. One-time cloud storage setup (iCloud 2TB only)
chmod +x benchmarks/m3pro/*.sh
./benchmarks/m3pro/setup_cloud_storage.sh

# 2. Build FlexAID for M3 Pro (Metal ON, all benchmarks ON)
./benchmarks/m3pro/build_m3pro.sh

# 3. Run all benchmarks (kernels + tier-1 + tier-2)
./benchmarks/m3pro/run_benchmarks.sh
```

## Scripts

| Script | Purpose |
|--------|---------|
| `setup_cloud_storage.sh` | Create dirs on iCloud, write `~/.flexaidds_env`, add symlinks (iCloud only) |
| `build_m3pro.sh` | CMake configure + build with Metal, OpenMP, Eigen, all benchmarks |
| `run_benchmarks.sh` | Run kernel + tier-1 + tier-2 benchmarks (writes directly to iCloud) |

## Selective Runs

```bash
./benchmarks/m3pro/run_benchmarks.sh --kernels-only  # dispatch, vcfbatch, tencom
./benchmarks/m3pro/run_benchmarks.sh --tier1-only     # CASF-2016, 5 targets
./benchmarks/m3pro/run_benchmarks.sh --tier2-only     # all 10 datasets
```

## Storage Architecture (iCloud 2TB only)

All writes go directly to iCloud Drive. No other cloud storage is used.

```
iCloud 2TB (PRIMARY + ONLY)              
  FlexAIDdS/
  ├── build/       ← NOT synced to local SSD (symlink target on iCloud)
  ├── benchmark_data/
  ├── results/
  │   ├── kernels/
  │   ├── tier1/
  │   └── tier2/
  └── logs/
```

- Writes go to iCloud first (lowest latency on macOS)
- `build/` excluded from any local sync (rebuild is cheap)
- `mirror_to_gdrive.sh` is now a no-op (Google Drive support removed)

## Memory Budget (18GB Unified)

| Component | Allocation |
|-----------|-----------|
| macOS + cloud sync | 3 GB |
| Metal GPU buffers | 4 GB |
| Tier-1 workers (×4) | 2.5 GB each |
| Tier-2 workers (×2) | 4.5 GB each |

Tier-2 datasets run sequentially (one at a time) to prevent memory pressure.

## Configuration

Hardware profile and all parameters are declared in `m3pro_profile.yaml`.
Environment variables are stored in `~/.flexaidds_env` (auto-sourced from `.zshrc`).

The robust Python failsafe runner (`failsafe_campaign.py`) supports the same environment variables plus `--remote-base`, `--no-remote-sync`, `--lock-dir`, and auto repo detection for portable / Codex / CI use. See its `--help` and module docstring.

For day-to-day use on this exact M3 Pro iCloud-only machine, use the convenience wrapper:
```bash
./benchmarks/m3pro/grok_master_launcher.sh full
```
It provides safe preflight/launch/sync/analyze modes with the recommended healthy settings for Astex + HAP2 campaigns (everything durable on iCloud, local APFS hot paths, max 18GB hardware utilization). See the script header for details.

## Mirror Script (Deprecated)

`mirror_to_gdrive.sh` is kept only as a placeholder and now exits immediately with a message that Google Drive support has been removed. All data stays on your 2TB iCloud Drive.
