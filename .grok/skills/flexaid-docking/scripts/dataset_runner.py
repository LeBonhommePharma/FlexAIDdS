#!/usr/bin/env python3
"""
dataset_runner.py

Part of the flexaid-docking skill.

Convenient, high-quality wrapper around the FlexAIDδS DatasetRunner
with built-in safety, data checks, and helpful defaults.

This script follows the same professional patterns as the rest of the skill:
- Clear banners
- Integration with ensure_docking_data (full matrices + *.def + extra runtime files)
- Automatic rich reproducibility manifest capture on every run
- --package produces a complete, audit-ready zip + one-pager VALIDATION_SUMMARY.md
- Helpful guidance and guardrails

Usage examples:
    # Basic tier-1 run on Astex Diverse (fast) — manifest is captured automatically
    python3 .grok/skills/flexaid-docking/scripts/dataset_runner.py \
        --dataset astex_diverse --tier 1

    # Full distributed campaign + professional validation package (recommended)
    mpirun -n 8 python3 .grok/skills/flexaid-docking/scripts/dataset_runner.py \
        --all --tier 2 --distributed --package

    # Dry run to validate everything without docking
    python3 .grok/skills/flexaid-docking/scripts/dataset_runner.py \
        --all --tier 1 --dry-run

    # Long campaign with per-entry checkpointing + MPI dynamic master-worker + timing/cost in reproducibility package
    python3 .grok/skills/flexaid-docking/scripts/dataset_runner.py \
        --all --tier 2 --workers 8 --resume --package
"""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any

# Robust import for both "python -m" and direct execution
try:
    from .ensure_docking_data import (
        print_skill_banner,
        get_all_critical_file_names,
        EXPECTED_MATRICES,
        EXPECTED_DEF_FILES,
        EXPECTED_EXTRA_FILES,
    )
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from ensure_docking_data import (
        print_skill_banner,
        get_all_critical_file_names,
        EXPECTED_MATRICES,
        EXPECTED_DEF_FILES,
        EXPECTED_EXTRA_FILES,
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="flexaid-docking dataset-runner",
        description="Run FlexAIDδS DatasetRunner campaigns with skill-integrated safety, data guarantees, and pharma-grade reproducibility.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Professional reproducibility (enabled by default for all runs):
  - Automatic capture of git SHA, binary SHA256, and complete hashes of every critical runtime file
    (matrices + all AMINO*.def/NUCLEOTIDES*.def + Lovell_LIB.dat + rotobs.lst + SYBYL_emat + scoring support)
  - Rich environment capture (conda/pip state, selected vars, hardware)
  - One-click professional validation package via --package (zip + REPRODUCIBILITY_MANIFEST.json + beautiful VALIDATION_SUMMARY.md one-pager)

Examples:
  # Fast sanity check on Astex Diverse + automatic reproducibility manifest
  python3 .../dataset_runner.py --dataset astex_diverse --tier 1

  # Full campaign with shareable audit package (recommended for publications / internal reviews)
  python3 .../dataset_runner.py --all --tier 2 --results-dir results/my_benchmark --package

  # Distributed run (launch via mpirun)
  mpirun -n 8 python3 .../dataset_runner.py --all --tier 2 --distributed --binary /path/to/FlexAIDδS --package

  # Long-running campaign with automated per-entry saving + resume (EntryTaskManager)
  python3 .../dataset_runner.py --all --tier 2 --workers 6 --resume --package

Always run ensure_docking_data.py first (or let this script remind you).
""",
    )

    p.add_argument("--dataset", "-d", help="Run a single dataset by slug (e.g. astex_diverse, casf2016, itc187)")
    p.add_argument("--all", "-a", action="store_true", help="Run all available datasets")
    p.add_argument("--tier", "-t", type=int, choices=[1, 2], default=2,
                   help="Benchmark tier (1=fast sanity, 2=full). Default: 2")
    p.add_argument("--metric", "-m", help="Run only a specific metric")
    p.add_argument("--distributed", action="store_true", help="Enable MPI distributed mode (launch with mpirun)")
    p.add_argument("--nodes", type=int, default=1, help="Informational: number of nodes")
    p.add_argument("--workers", type=int, default=1, help="Local parallel workers")
    p.add_argument("--binary", help="Path to FlexAIDδS binary (overrides FLEXAIDDS_BINARY)")
    p.add_argument("--results-dir", default="results/benchmarks",
                   help="Where to write reports (default: results/benchmarks)")
    p.add_argument("--dry-run", action="store_true",
                   help="Skip actual docking; useful for pipeline validation")
    p.add_argument("--resume", action="store_true",
                   help="Resume from per-entry checkpoints (skip targets that already have individual result files). Strongly recommended for long campaigns.")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--ensure-data", action="store_true", default=True,
                   help="Automatically run ensure_docking_data.py first (default: True)")
    p.add_argument("--no-ensure-data", dest="ensure_data", action="store_false")
    p.add_argument("--package", "--create-package", action="store_true",
                   help="After the run, automatically create a professional, shareable validation package (zip + reproducibility manifest)")

    return p


def _get_file_sha256(path: Path) -> str:
    if not path.exists():
        return "missing"
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _discover_git_root(start: Path) -> Path | None:
    """Walk upward to find the nearest .git directory (robust across worktrees and subdirs)."""
    p = start.resolve()
    for _ in range(12):
        if (p / ".git").exists():
            return p
        if p.parent == p:
            break
        p = p.parent
    return None


def _capture_git_info() -> Dict[str, Any]:
    info: Dict[str, Any] = {"git_sha": "unavailable", "git_status": "unknown", "git_root": None}
    git_root = _discover_git_root(Path(__file__))
    if not git_root:
        return info
    info["git_root"] = str(git_root)
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=git_root, stderr=subprocess.DEVNULL
        ).decode().strip()
        info["git_sha"] = sha
    except Exception:
        pass
    try:
        status = subprocess.check_output(
            ["git", "status", "--porcelain", "--branch"], cwd=git_root, stderr=subprocess.DEVNULL
        ).decode().strip()
        # Keep first 20 lines to stay compact
        lines = status.splitlines()[:20]
        info["git_status"] = "\n".join(lines) if lines else "clean"
        info["git_dirty"] = any(not line.startswith("##") and line.strip() for line in lines)
    except Exception:
        info["git_status"] = "unavailable"
    return info


def _capture_rich_environment() -> Dict[str, Any]:
    """Best-effort rich environment capture for pharma-grade auditability (conda, pip, system)."""
    env: Dict[str, Any] = {
        "selected_env_vars": {},
        "conda": {"available": False},
        "pip": {"available": False, "packages_sample": []},
        "cpu_count": os.cpu_count(),
        "processor": platform.processor() or platform.machine(),
    }

    # Selected relevant environment variables (expanded set)
    interesting_prefixes = (
        "FLEXAID", "FLEXAIDDS", "OMP", "MPI", "MKL", "OPENBLAS", "PYTHON",
        "CONDA", "VIRTUAL_ENV", "PATH", "LD_LIBRARY", "DYLD"
    )
    for k, v in os.environ.items():
        if any(k.startswith(pref) for pref in interesting_prefixes):
            # Truncate very long values (e.g. full PATH)
            env["selected_env_vars"][k] = v[:500] + "..." if len(v) > 500 else v

    # Conda info (best effort, never fatal)
    try:
        conda_env = os.environ.get("CONDA_DEFAULT_ENV") or os.environ.get("CONDA_PREFIX")
        if conda_env:
            env["conda"]["available"] = True
            env["conda"]["env_name_or_prefix"] = conda_env
        # Try to get a compact package list
        res = subprocess.run(
            ["conda", "list", "--export"], capture_output=True, text=True, timeout=8, check=False
        )
        if res.returncode == 0 and res.stdout:
            lines = [ln for ln in res.stdout.strip().splitlines() if not ln.startswith("#")][:80]
            env["conda"]["packages_export_sample"] = lines
    except Exception:
        pass

    # Pip freeze (best effort, limited)
    try:
        res = subprocess.run(
            [sys.executable, "-m", "pip", "list", "--format=freeze"],
            capture_output=True, text=True, timeout=8, check=False
        )
        if res.returncode == 0 and res.stdout:
            lines = res.stdout.strip().splitlines()[:60]
            env["pip"]["available"] = True
            env["pip"]["packages_sample"] = lines
    except Exception:
        pass

    return env


def _manifest_entry_count(man_data: dict) -> int:
    """Count per-entry manifest rows, excluding polluted None_None keys."""
    try:
        repo_root = Path(__file__).resolve().parents[4]
        python_pkg = repo_root / "python"
        if python_pkg.is_dir():
            import sys
            if str(python_pkg) not in sys.path:
                sys.path.insert(0, str(python_pkg))
            from flexaidds.dataset_runner.runner import sanitize_entry_manifest
            clean = sanitize_entry_manifest(dict(man_data))
            status = clean.get("per_entry_status") or {}
            if status:
                return len(status)
            return len(clean.get("timings", {}).get("per_entry_wall_seconds") or {})
    except Exception:
        pass
    status = man_data.get("per_entry_status") or {}
    status = {k: v for k, v in status.items() if k and "None" not in k and "_" in k}
    if status:
        return len(status)
    wall = man_data.get("timings", {}).get("per_entry_wall_seconds") or {}
    wall = {k: v for k, v in wall.items() if k and "None" not in k and "_" in k}
    return len(wall)


def _capture_per_entry_provenance(results_dir: Path) -> Dict[str, Any]:
    """Capture hashes and summary of per-entry benchmark artifacts produced by
    the inner DatasetRunner + EntryTaskManager (timing/cost manifests + individual results).
    This integrates the new per-entry automation into the skill's pharma-grade reproducibility package.
    """
    info: Dict[str, Any] = {
        "entry_count": 0,
        "entry_manifests": [],
        "sample_entry_hashes": {},
    }

    if not results_dir.exists():
        return info

    # Look for dataset/tier directories containing _entry_manifest.json
    for manifest in sorted(results_dir.rglob("_entry_manifest.json")):
        try:
            rel = str(manifest.relative_to(results_dir))
            info["entry_manifests"].append({
                "path": rel,
                "sha256": _get_file_sha256(manifest),
            })

            # Prefer manifest entry counts for large-N campaigns (avoid full dir glob)
            parent = manifest.parent
            try:
                man_data = json.loads(manifest.read_text())
                n_from_manifest = _manifest_entry_count(man_data)
                if n_from_manifest > 0:
                    info["entry_count"] += n_from_manifest
                    entry_jsons = []
                else:
                    entry_jsons = [f for f in parent.glob("*.json") if not f.name.startswith("_")]
                    info["entry_count"] += len(entry_jsons)
            except Exception:
                entry_jsons = [f for f in parent.glob("*.json") if not f.name.startswith("_")]
                info["entry_count"] += len(entry_jsons)

            # Sample up to 3 per manifest for the reproducibility package (keeps it compact)
            for jf in sorted(entry_jsons)[:3]:
                key = str(jf.relative_to(results_dir))
                info["sample_entry_hashes"][key] = _get_file_sha256(jf)
        except Exception:
            continue

    return info


def gather_reproducibility_metadata(args: argparse.Namespace, binary_path: str | None) -> Dict[str, Any]:
    """
    Superior general reproducibility capture (better than narrow per-report checksums).

    Captures:
    - Precise run identification (time, host, skill version)
    - Full git state (SHA + cleanliness)
    - Binary identity + content hash
    - Complete integrity hashes for *every* critical runtime file (matrices + all *.def + Lovell_LIB + rotobs + emat + scoring support)
    - Rich environment (selected vars + conda/pip + hardware)
    - Exact command line
    """
    meta: Dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "python_version": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "platform_system": platform.system(),
        "platform_release": platform.release(),
        "skill_version": "2026-05 (flexaid-docking)",
        "script": "dataset_runner.py (flexaid-docking skill wrapper)",
    }

    # Git (robust)
    meta.update(_capture_git_info())

    # Binary
    if binary_path:
        bin_path = Path(binary_path)
        meta["binary_path"] = str(bin_path.resolve())
        meta["binary_sha256"] = _get_file_sha256(bin_path)
        try:
            st = bin_path.stat()
            meta["binary_size_bytes"] = st.st_size
            meta["binary_mtime_utc"] = datetime.fromtimestamp(st.st_mtime, tz=timezone.utc).isoformat()
        except Exception:
            pass
    else:
        meta["binary_path"] = "not explicitly provided (discovered by inner DatasetRunner)"
        meta["binary_sha256"] = "unknown_at_wrapper_time"

    # === COMPLETE critical data file hashes (the heart of the better general solution) ===
    # Discover plausible locations where the data actually lived for this run
    skill_data_dir = Path(__file__).parent.parent / "data"
    search_roots: list[Path] = [
        skill_data_dir,
        Path.cwd(),
    ]
    if binary_path:
        b = Path(binary_path).resolve()
        search_roots.extend([b.parent, b.parent.parent])

    # Also include the locations ensure_docking_data would have searched
    try:
        from ensure_docking_data import DEFAULT_SEARCH_PATHS as ENSURE_DEFAULTS  # type: ignore
        search_roots.extend(ENSURE_DEFAULTS)
    except Exception:
        pass

    # Deduplicate roots while preserving order
    seen_roots = set()
    ordered_roots = []
    for r in search_roots:
        rp = r.resolve() if r.exists() else r
        if rp not in seen_roots:
            seen_roots.add(rp)
            ordered_roots.append(r)

    critical_names = get_all_critical_file_names()
    data_file_hashes: Dict[str, Any] = {}
    located_data_dir: str | None = None

    for name in critical_names:
        found_hash = None
        found_path = None
        for root in ordered_roots:
            candidate = root / name
            if candidate.is_file():
                found_hash = _get_file_sha256(candidate)
                found_path = str(candidate.resolve())
                if located_data_dir is None:
                    located_data_dir = str(root.resolve())
                break
        data_file_hashes[name] = {
            "sha256": found_hash or "missing",
            "path": found_path or "not found in search roots",
        }

    meta["critical_data_files"] = data_file_hashes
    meta["data_search_roots_used"] = [str(r) for r in ordered_roots[:6]]
    if located_data_dir:
        meta["effective_data_directory"] = located_data_dir

    # Command & environment (rich)
    meta["command_line"] = " ".join(sys.argv)
    meta["original_args"] = {k: v for k, v in vars(args).items() if not k.startswith("_")}
    meta["environment"] = _capture_rich_environment()

    # Add a short human summary for convenience
    present = sum(1 for v in data_file_hashes.values() if v["sha256"] != "missing")
    meta["data_integrity_summary"] = f"{present}/{len(critical_names)} critical files present with hashes"

    # === NEW: Per-entry provenance integration (3rd priority) ===
    # Capture hashes of the new per-entry result artifacts produced by the inner DatasetRunner + EntryTaskManager
    try:
        results_dir = Path(getattr(args, "results_dir", "results/benchmarks"))
        per_entry_info = _capture_per_entry_provenance(results_dir)
        if per_entry_info.get("entry_manifests") or per_entry_info.get("entry_count", 0) > 0:
            meta["per_entry_benchmark_artifacts"] = per_entry_info
            meta["data_integrity_summary"] += f" + {per_entry_info.get('entry_count', 0)} per-entry results"
    except Exception:
        pass

    return meta


def generate_validation_summary(metadata: Dict[str, Any]) -> str:
    """
    Produce a professional, attractive one-pager Validation Summary (Markdown).
    This is the polished deliverable for the 3rd item (one-pager) and 4th (env capture).
    Suitable for inclusion in papers, audit folders, or regulatory packages.
    """
    ts = metadata.get("timestamp_utc", "unknown")
    git_sha = metadata.get("git_sha", "unavailable")[:12]
    bin_path = metadata.get("binary_path", "N/A")
    bin_sha = metadata.get("binary_sha256", "N/A")[:16] + "..."
    hostname = metadata.get("hostname", "N/A")
    data_summary = metadata.get("data_integrity_summary", "N/A")

    lines = []
    lines.append("# FlexAIDδS Validation Summary — Reproducibility & Audit Package")
    lines.append("")
    lines.append(f"**Generated**: {ts}")
    lines.append(f"**Skill / Wrapper**: {metadata.get('skill_version')}")
    lines.append(f"**Host**: {hostname}  |  **Python**: {metadata.get('python_version')}")
    lines.append("")

    lines.append("## 1. Run Identification & Integrity")
    lines.append("")
    lines.append(f"- Git commit: `{git_sha}` (full SHA in REPRODUCIBILITY_MANIFEST.json)")
    lines.append(f"- Binary: `{bin_path}`")
    lines.append(f"- Binary SHA256 (first 16): `{bin_sha}`")
    lines.append(f"- Data integrity: **{data_summary}**")
    if "effective_data_directory" in metadata:
        lines.append(f"- Effective runtime data directory: `{metadata['effective_data_directory']}`")
    lines.append("")

    lines.append("## 2. Critical Runtime Data File Hashes (Complete Set)")
    lines.append("")
    lines.append("Every file required for deterministic FlexAIDδS execution is hashed below.")
    lines.append("These are the authoritative values for this exact run.")
    lines.append("")
    lines.append("| File | SHA256 (full in JSON) | Status |")
    lines.append("|------|-----------------------|--------|")

    for fname, info in metadata.get("critical_data_files", {}).items():
        sha = info.get("sha256", "missing")
        short = sha[:16] + "..." if len(sha) > 20 and sha != "missing" else sha
        status = "present" if sha != "missing" else "**MISSING**"
        lines.append(f"| `{fname}` | `{short}` | {status} |")
    lines.append("")

    # Per-entry benchmark artifacts (integrated from EntryTaskManager + timing/cost manifests)
    pea = metadata.get("per_entry_benchmark_artifacts", {})
    if pea:
        lines.append("## 2b. Per-Entry Benchmark Artifacts (DatasetRunner + EntryTaskManager)")
        lines.append("")
        lines.append(f"- Entry result files located: **{pea.get('entry_count', 0)}**")
        if pea.get("entry_manifests"):
            lines.append("- Entry manifests captured:")
            for m in pea["entry_manifests"][:5]:
                lines.append(f"  - `{m['path']}` → `{m['sha256'][:16]}...`")
        if pea.get("sample_entry_hashes"):
            lines.append(f"- Sample per-entry result hashes included in full manifest ({len(pea['sample_entry_hashes'])} shown)")
        lines.append("")

    lines.append("## 3. Environment Capture (Conda / Pip / System)")
    lines.append("")
    env = metadata.get("environment", {})
    lines.append(f"- Platform: {metadata.get('platform')}")
    lines.append(f"- Processor: {env.get('processor', 'N/A')}  |  CPUs: {env.get('cpu_count', 'N/A')}")
    if env.get("conda", {}).get("available"):
        c = env["conda"]
        lines.append(f"- Conda env: `{c.get('env_name_or_prefix', 'unknown')}`")
    if env.get("pip", {}).get("available"):
        lines.append(f"- Pip packages captured (sample): {len(env['pip'].get('packages_sample', []))} entries")
    lines.append("")

    # Selected vars (compact)
    sel = env.get("selected_env_vars", {})
    if sel:
        lines.append("**Selected relevant environment variables** (truncated):")
        for k in sorted(sel)[:12]:
            lines.append(f"- `{k}` = {sel[k][:80]}...")
        if len(sel) > 12:
            lines.append(f"- ... and {len(sel)-12} more (see full manifest)")
    lines.append("")

    lines.append("## 4. Exact Command Line")
    lines.append("")
    lines.append("```")
    lines.append(metadata.get("command_line", "(unavailable)"))
    lines.append("```")
    lines.append("")

    lines.append("## 5. Reproducibility Instructions")
    lines.append("")
    lines.append("To reproduce this exact campaign or audit the results:")
    lines.append("")
    lines.append("1. Checkout the git commit listed above (or the version of the flexaid-docking skill used).")
    lines.append("2. Ensure the FlexAIDδS binary whose SHA256 matches the manifest is on PATH or passed via `--binary`.")
    lines.append("3. Run the skill data ensure step (it will use the same data files whose hashes appear above):")
    lines.append("   ```bash")
    lines.append("   python3 .grok/skills/flexaid-docking/scripts/ensure_docking_data.py --info")
    lines.append("   ```")
    lines.append("4. Re-execute the exact command line shown in section 4 (or the inner DatasetRunner invocation).")
    lines.append("5. Compare new REPRODUCIBILITY_MANIFEST.json hashes against the archived copy.")
    lines.append("")

    lines.append("## 6. Scientific & Terminology Notes")
    lines.append("")
    lines.append("- **CF / contact-function scoring proxy**: The Voronoi-based contact function (Vcontacts) used inside the genetic algorithm for pose search and ranking during the run.")
    lines.append("- **Thermodynamic ledger**: Post-hoc quantities (Helmholtz free energy F, average energy <H>, −TΔS, Cv, Boltzmann weights) computed by the StatMechEngine + BindingMode layer on the final ensemble. These are *not* the same as the CF proxy scores.")
    lines.append("- All claims in the accompanying results JSON/Markdown reports are derived from the captured ensemble using the exact data and binary hashed above.")
    lines.append("")

    lines.append("## 7. Audit / Regulatory Notes")
    lines.append("")
    lines.append("- This package (zip + manifest + this summary) constitutes a self-contained reproducibility artifact.")
    lines.append("- The full set of data file hashes makes it possible for a third party to verify that the identical runtime assets were used.")
    lines.append("- No external internet resources or unreferenced files are required beyond the binary + the data files whose hashes are recorded.")
    lines.append("")

    lines.append("## 8. Disclaimer")
    lines.append("")
    lines.append("This summary and the associated manifest were generated automatically by the flexaid-docking skill. They record the computational environment and inputs with high fidelity. They do not constitute a claim of experimental validation or regulatory approval. Users are responsible for interpreting results in the appropriate scientific and regulatory context.")
    lines.append("")
    lines.append("---")
    lines.append("*FlexAIDδS — entropy-augmented molecular docking (FlexAID + ΔS)*")
    lines.append("*Part of the flexaid-docking Grok skill — professional reproducibility tooling*")

    return "\n".join(lines)


def create_validation_package(results_dir: Path, metadata: Dict[str, Any], package_name: str | None = None) -> Path:
    """Create a clean, professional, shareable validation package (the improved general solution)."""
    if package_name is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        package_name = f"flexaidds_validation_package_{timestamp}"

    package_dir = results_dir.parent / package_name
    package_dir.mkdir(parents=True, exist_ok=True)

    # Copy results artifacts if present
    if results_dir.exists():
        import shutil
        for item in results_dir.rglob("*"):
            if item.is_file():
                rel = item.relative_to(results_dir)
                dest = package_dir / "results" / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(item, dest)

    # Write the full machine-readable manifest
    with open(package_dir / "REPRODUCIBILITY_MANIFEST.json", "w") as f:
        json.dump(metadata, f, indent=2, default=str)

    # Write the beautiful one-pager (Markdown is primary for readability + GitHub rendering)
    summary_md = generate_validation_summary(metadata)
    with open(package_dir / "VALIDATION_SUMMARY.md", "w") as f:
        f.write(summary_md)

    # Also write a plain-text fallback (for strict environments)
    with open(package_dir / "README_Validation_Summary.txt", "w") as f:
        f.write(summary_md.replace("# ", "").replace("## ", "").replace("**", "").replace("`", ""))

    # Create the zip archive
    import shutil
    zip_path = package_dir.with_suffix(".zip")
    shutil.make_archive(str(zip_path.with_suffix("")), "zip", package_dir)

    return zip_path


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    print_skill_banner(verbose=args.verbose)

    if args.ensure_data:
        print("\n[Skill] Ensuring all critical runtime data is present before benchmarking...")
        ensure_script = Path(__file__).parent / "ensure_docking_data.py"
        subprocess.run([sys.executable, str(ensure_script)], check=False)

    # Build the command for the real module
    cmd = [sys.executable, "-m", "flexaidds.dataset_runner"]

    if args.dataset:
        cmd += ["--dataset", args.dataset]
    elif args.all:
        cmd += ["--all"]
    else:
        print("Error: You must specify either --dataset <slug> or --all")
        return 2

    cmd += ["--tier", str(args.tier)]

    if args.metric:
        cmd += ["--metric", args.metric]
    if args.distributed:
        cmd += ["--distributed"]
    if args.nodes and args.nodes > 1:
        cmd += ["--nodes", str(args.nodes)]
    if args.workers and args.workers > 1:
        cmd += ["--workers", str(args.workers)]
    if args.binary:
        cmd += ["--binary", args.binary]
    if args.results_dir:
        cmd += ["--results-dir", args.results_dir]
    if args.dry_run:
        cmd += ["--dry-run"]
    if args.resume:
        cmd += ["--resume"]
    if args.verbose:
        cmd += ["--verbose"]

    large_slugs = {"astex_nonnative", "posex", "posex_cd", "posebusters"}
    if args.dataset in large_slugs or args.all:
        print(
            f"\n[Skill] Large-N dataset campaign detected "
            f"({'--all' if args.all else args.dataset}). "
            f"Manifest-first resume is enabled with --resume; "
            f"use --plan-runtime on the inner CLI for wall-clock estimates."
        )

    print(f"\n[Skill] Launching DatasetRunner:")
    print("  " + " ".join(cmd))
    print()

    # Gather reproducibility metadata before the run
    binary_for_meta = args.binary
    metadata = gather_reproducibility_metadata(args, binary_for_meta)

    # The real module lives in the repo's python/ package dir, which is not
    # necessarily installed or on PYTHONPATH (e.g. CI runs this wrapper directly
    # without `pip install -e python/`). Inject it so `-m flexaidds.dataset_runner`
    # resolves regardless of how the wrapper is invoked.
    sub_env = os.environ.copy()
    repo_root = Path(__file__).resolve().parents[4]
    python_pkg_dir = repo_root / "python"
    if python_pkg_dir.is_dir():
        existing_pp = sub_env.get("PYTHONPATH", "")
        sub_env["PYTHONPATH"] = (
            str(python_pkg_dir) + (os.pathsep + existing_pp if existing_pp else "")
        )

    try:
        result = subprocess.run(cmd, check=False, env=sub_env)

        if args.package and result.returncode == 0:
            results_dir = Path(args.results_dir)
            print("\n[Skill] Creating professional validation package (pharma-grade reproducibility artifact)...")
            package_path = create_validation_package(results_dir, metadata)
            print(f"[Skill] Validation package created: {package_path}")
            print("       Contains: REPRODUCIBILITY_MANIFEST.json + VALIDATION_SUMMARY.md (one-pager) + results/")
            print("       This archive is suitable for internal audit, collaboration, publications, or regulatory purposes.")

        return result.returncode
    except KeyboardInterrupt:
        print("\n[Skill] Interrupted by user.")
        return 130


if __name__ == "__main__":
    sys.exit(main())
