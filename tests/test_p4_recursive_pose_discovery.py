#!/usr/bin/env python3
"""Recursive pose discovery for any-persisted ceiling must include r*/ trees."""
from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "p4_best_of_n_diagnostic.py"


def _load():
    spec = importlib.util.spec_from_file_location("p4_best_of_n_diagnostic", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


class TestP4RecursivePoseDiscovery(unittest.TestCase):
    def test_recursive_includes_restart_unique_poses(self) -> None:
        p4 = _load()
        self.assertTrue(hasattr(p4, "discover_pose_pdbs"))
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            code = "1ABC"
            # top-level only pose
            (root / f"{code}_0.pdb").write_text("TOP0\n")
            (root / f"{code}_INI.pdb").write_text("SEED\n")
            # restart unique pose missed by top-level glob semantics
            r3 = root / "r3"
            r3.mkdir()
            (r3 / f"{code}_17.pdb").write_text("RESTART17\n")
            # duplicate content should dedupe
            (r3 / f"{code}_0.pdb").write_text("TOP0\n")

            top = p4.discover_pose_pdbs(str(root), code, recursive=False)
            rec = p4.discover_pose_pdbs(str(root), code, recursive=True)
            self.assertTrue(any(Path(p).name == f"{code}_0.pdb" for p in top))
            self.assertFalse(any("INI" in Path(p).name for p in rec))
            names = {Path(p).name for p in rec}
            self.assertIn(f"{code}_17.pdb", names)
            # unique by content: TOP0 once + RESTART17
            self.assertEqual(len(rec), 2)

    def test_hungarian_rmsd_function_is_shipped(self) -> None:
        p4 = _load()
        import numpy as np

        a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        b = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        r = p4.hungarian_rmsd(a, ["C", "C"], b, ["C", "C"])
        self.assertIsNotNone(r)
        self.assertLess(float(r), 1e-9)


if __name__ == "__main__":
    unittest.main()
