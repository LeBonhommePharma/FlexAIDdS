# Reproducing the Astex 85 Benchmark on Windows

Native Windows is **not supported** for the FlexAIDdS benchmark runner.
The C++ `DatasetRunner` requires POSIX process management (`fork`/`exec`/`waitpid`/`killpg`),
and MSVC cannot yet compile the C++26 core engine (`/std:c++26` is absent in
VS 2022 17.14 / MSVC 19.44 — see CMakeLists.txt line ~68).

The supported Windows path is **WSL2 with Ubuntu 22.04 LTS**, which runs the
identical Linux binary inside a full Linux kernel on your machine.
Performance is close to native (WSL2 uses hardware virtualisation, not emulation).

---

## One-time WSL2 setup (~5 minutes)

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

This installs WSL2 and the Ubuntu 22.04 LTS distro in one step.
Reboot when prompted, then launch **Ubuntu** from the Start menu.

> **Already have WSL1?** Upgrade with:
> ```powershell
> wsl --set-default-version 2
> wsl --set-version Ubuntu 2
> ```

---

## Inside the Ubuntu terminal

### 1. Install build dependencies

```bash
sudo apt update && sudo apt install -y \
    git cmake build-essential \
    libboost-all-dev libeigen3-dev \
    libssl-dev curl python3 python3-pip
```

GCC 14 is required (Ubuntu 22.04 ships GCC 11 by default):

```bash
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt install -y gcc-14 g++-14
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100
gcc --version   # should print gcc (Ubuntu 14.x) 14.x.x
```

### 2. Clone the repo

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git
cd FlexAIDdS
```

### 3. Run the reproduce script

```bash
bash scripts/reproduce_astex85.sh
```

This is exactly the same script used on macOS and Linux.
Results will appear in `~/FlexAIDdS_reviewer_benchmark/`.

---

## Expected runtime in WSL2

| Hardware | Workers | Estimated wall time |
|---|---|---|
| Modern laptop (8–16 cores) | 4 | ~35–55 minutes |
| Desktop workstation (32+ cores) | 8 | ~20–30 minutes |

The bottleneck is a single long-running target (1OF6, ~34 min sequential).
With 4 workers the total wall time is dominated by this target running in parallel
with 3 others.

To use more workers:

```bash
FLEXAIDDS_BENCH_THREADS=8 FLEXAIDDS_OMP_THREADS=2 bash scripts/reproduce_astex85.sh
```

---

## Accessing results from Windows Explorer

Your WSL2 Ubuntu home directory is accessible at:

```
\\wsl.localhost\Ubuntu\home\<your-username>\FlexAIDdS_reviewer_benchmark\
```

Paste that path into Windows Explorer's address bar to browse the output files.

---

## Why not native Windows?

Three independent blockers prevent a native Windows port today:

1. **MSVC lacks C++26 support** — `/std:c++26` is absent in VS 2022 17.14
   (MSVC 19.44). The engine uses C++26 features (structured bindings in ranges,
   `std::format`, `[[nodiscard]]` on constructors) that don't compile under C++20.

2. **POSIX process model** — `DatasetRunner.cpp` uses `fork()`, `execl()`,
   `waitpid()`, `killpg()`, and `setpgid()` for subprocess management and
   per-job timeout/kill. The Windows fallback (`system()`) loses all timeout
   and kill-on-shutdown functionality; partial results from a crashed run
   cannot be safely cleaned up.

3. **Signal handling** — job cancellation and clean shutdown use `SIGTERM`/`SIGKILL`
   via process groups. Windows has no equivalent without a significant rewrite
   using `CreateProcess` + `TerminateProcess` + `Job Objects`.

WSL2 avoids all three issues cleanly: it runs the genuine Linux binary under a
real Linux kernel, so results are numerically identical to Linux builds.

---

## Verified environment

| Component | Version |
|---|---|
| WSL2 kernel | ≥ 5.15 (default in Windows 11 / Windows 10 22H2) |
| Ubuntu | 22.04 LTS |
| GCC | 14.2 |
| CMake | ≥ 3.28 |
| Eigen3 | ≥ 3.4 |
| Boost | ≥ 1.74 |
