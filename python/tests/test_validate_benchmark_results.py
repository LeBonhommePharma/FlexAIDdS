"""Tests for the benchmark validation script."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


def _write_summary_csv(path: Path) -> None:
    rows = [
        ("1SQ5", 1, 0.80, 1.0, 1.0, 1.0, 1.00, 12, 101.0, 1),
        ("2HB1", 1, 0.90, 1.0, 1.0, 1.0, 1.10, 12, 102.0, 1),
        ("1R1H", 1, 1.00, 1.0, 1.0, 1.0, 1.20, 12, 103.0, 1),
        ("1T46", 1, 1.10, 1.0, 1.0, 1.0, 1.30, 12, 104.0, 1),
        ("2C69", 1, 1.20, 1.0, 1.0, 1.0, 1.40, 12, 105.0, 1),
        ("1G9V", 1, 0.85, 1.0, 1.0, 1.0, 1.50, 12, 106.0, 1),
        ("1GM8", 1, 1.25, 1.0, 1.0, 1.0, 1.60, 12, 107.0, 1),
        ("1GPK", 1, 0.95, 1.0, 1.0, 1.0, 1.70, 12, 108.0, 1),
        ("1HP0", 1, 1.05, 1.0, 1.0, 1.0, 1.80, 12, 109.0, 1),
        ("1HNN", 1, 1.15, 1.0, 1.0, 1.0, 1.90, 12, 110.0, 1),
    ]
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "pdb_id",
                "best_score",
                "rmsd_to_crystal",
                "predicted_dG",
                "predicted_dH",
                "predicted_TdS",
                "shannon_entropy",
                "num_poses",
                "wall_time_s",
                "success",
            ]
        )
        for row in rows:
            writer.writerow(row)


def _write_manifest(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "slug: astex_diverse_test",
                "name: Astex Diverse Test",
                "docking_mode: self_docking",
                "description: Synthetic manifest for validator regression testing",
                "tier: 1",
                "tier1_subset_size: 5",
                "metrics:",
                "  - docking_power_top1",
                "  - mean_rmsd",
                "  - median_rmsd",
                "expected_baselines:",
                "  docking_power_top1: 0.90",
                "  mean_rmsd: 1.20",
                "  median_rmsd: 1.20",
                "published_baselines:",
                "  docking_power_top1: 0.60",
                "published_source: JCIM 2015 comparator",
                "baseline_tolerance: 0.05",
                "",
            ]
        )
    )


def test_validate_benchmark_results_uses_manifest_comparator(tmp_path):
    summary_csv = tmp_path / "summary.csv"
    manifest = tmp_path / "manifest.yaml"
    out_dir = tmp_path / "figures"
    _write_summary_csv(summary_csv)
    _write_manifest(manifest)

    script = Path(__file__).resolve().parents[2] / "scripts" / "validate_benchmark_results.py"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            str(summary_csv),
            "--manifest",
            str(manifest),
            "--out-dir",
            str(out_dir),
            "--no-plots",
        ],
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parents[2],
        check=False,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert "Published source: JCIM 2015 comparator" in result.stdout
    assert "Exact Binomial Test vs. Published Comparator" in result.stdout
    assert "Significant (p < 0.05): YES" in result.stdout
    assert "OVERALL: PASS" in result.stdout
