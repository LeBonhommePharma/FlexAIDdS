# M3 Pro Benchmark Setup — iCloud + Google Drive

Benchmark environment for **MacBook Pro 14" M3 Pro 18GB** with zero local SSD
space. All data lives on iCloud Drive (2TB, primary) and Google Drive AI Pro
(5TB, async mirror).

## Quick Start

```bash
# 1. One-time cloud storage setup
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
| `setup_cloud_storage.sh` | Create dirs on iCloud/GDrive, write `~/.flexaidds_env`, add symlinks |
| `build_m3pro.sh` | CMake configure + build with Metal, OpenMP, Eigen, all benchmarks |
| `run_benchmarks.sh` | Run kernel + tier-1 + tier-2 benchmarks with cloud sync |
| `mirror_to_gdrive.sh` | Async rsync from iCloud to Google Drive (called automatically) |

## Selective Runs

```bash
./benchmarks/m3pro/run_benchmarks.sh --kernels-only  # dispatch, vcfbatch, tencom
./benchmarks/m3pro/run_benchmarks.sh --tier1-only     # CASF-2016, 5 targets
./benchmarks/m3pro/run_benchmarks.sh --tier2-only     # all 10 datasets
```

## Storage Architecture

```
iCloud 2TB (PRIMARY)              Google Drive 5TB (MIRROR)
  FlexAIDdS/                        FlexAIDdS/
  ├── build/       ← NOT synced     ├── benchmark_data/
  ├── benchmark_data/                ├── results/
  ├── results/                       │   ├── kernels/
  │   ├── kernels/                   │   ├── tier1/
  │   ├── tier1/                     │   └── tier2/
  │   └── tier2/                     └── logs/
  └── logs/
```

- Writes go to iCloud first (lowest latency on macOS)
- `rsync -avz --delete` mirrors to Google Drive after each benchmark phase
- `build/` excluded (rebuild is cheaper than syncing .o files)

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

## Manual Mirror

```bash
./benchmarks/m3pro/mirror_to_gdrive.sh          # foreground
nohup ./benchmarks/m3pro/mirror_to_gdrive.sh &  # background
```
