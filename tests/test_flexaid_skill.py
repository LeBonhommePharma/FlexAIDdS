#!/usr/bin/env python3
"""Regression test for the FlexAIDdS agent skill package."""
from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = REPO_ROOT / "skills" / "flexaid-docking" / "scripts" / "validate_skill.py"


class FlexAidSkillPackagingTest(unittest.TestCase):
    def test_skill_validator_passes(self) -> None:
        result = subprocess.run(
            [sys.executable, str(VALIDATOR)],
            cwd=REPO_ROOT,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, output)
        self.assertIn("Validation passed.", output)


if __name__ == "__main__":
    unittest.main()
