#!/usr/bin/env python3
"""Guard: conditional_scanned_pool_ceiling is NOT a full-pop / any-pose ceiling.

Prevents re-deriving the tautology BCR vs pool as 'retained-pool extraction' evidence.
Drives the shipped header + admission contract text on the real repo paths.
"""
from __future__ import annotations

import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
HDR = REPO / "LIB" / "DatasetRunner.h"
CONTRACT = REPO / "benchmarks" / "protocols" / "admission_metrics_contract.md"


class TestPoolCeilingIsNotAnyPose(unittest.TestCase):
    def test_header_aliases_bcr_to_scanned_pool(self) -> None:
        text = HDR.read_text(encoding="utf-8", errors="replace")
        self.assertIn("conditional_scanned_pool_ceiling", text)
        self.assertIn("best_cluster_rmsd", text)
        # Legacy column is explicitly equal to the scanned-pool ceiling.
        self.assertRegex(
            text,
            r"best_cluster_rmsd.*==.*conditional_scanned_pool_ceiling"
            r"|== conditional_scanned_pool_ceiling \(legacy column\)",
        )
        # Explicit not any-pose / not full emission census.
        self.assertRegex(
            text,
            r"not guaranteed any-pose|not full\s+emission census|NOT any-pose",
            msg="DatasetRunner.h must state pool ceiling is not full-pop any-pose",
        )

    def test_admission_contract_marks_s3_emission_only(self) -> None:
        text = CONTRACT.read_text(encoding="utf-8", errors="replace")
        self.assertIn("conditional_scanned_pool_ceiling", text)
        self.assertRegex(
            text,
            r"not any-pose|scanned emission pool only|not.*full emission",
            msg="admission contract must demote pool ceiling from any-pose claim",
        )
        # S3 diagnostic only
        self.assertRegex(text, r"S3.*[Dd]iagnostic|Diagnostic only")


if __name__ == "__main__":
    unittest.main()
