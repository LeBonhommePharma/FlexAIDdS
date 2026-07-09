from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmarks.astex_entropy.entropy import limit_poses_per_target
from benchmarks.astex_entropy.models import PoseRecord


def _pose(target_id: str, pose_id: str, score: str) -> PoseRecord:
    return PoseRecord(
        target_id=target_id,
        mode="native",
        tool="flexaidds",
        pose_id=pose_id,
        pose_sdf=f"/tmp/{pose_id}.sdf",
        receptor_pdb="/tmp/receptor.pdb",
        reference_sdf="/tmp/reference.sdf",
        raw_score=score,
        score_direction="lower",
    )


class AstexEntropyRescoreTests(unittest.TestCase):
    def test_limit_poses_per_target_keeps_best_raw_scores(self) -> None:
        records = [
            _pose("1G9V", "p1", "-10"),
            _pose("1G9V", "p2", "-5"),
            _pose("1G9V", "p3", "-20"),
            _pose("1HP0", "q1", "1"),
            _pose("1HP0", "q2", "0"),
        ]
        limited = limit_poses_per_target(records, 2)
        by_target: dict[str, list[str]] = {}
        for record in limited:
            by_target.setdefault(record.target_id, []).append(record.pose_id)
        self.assertEqual(len(limited), 4)
        self.assertEqual(by_target["1G9V"], ["p3", "p1"])
        self.assertEqual(by_target["1HP0"], ["q2", "q1"])

    def test_limit_zero_returns_all(self) -> None:
        records = [_pose("1G9V", "p1", "-1")]
        self.assertEqual(limit_poses_per_target(records, 0), records)


class FlexaiddsLigandExtractTests(unittest.TestCase):
    def test_extract_keeps_all_high_serial_atoms_with_digit_resname(self) -> None:
        """Regression: RQ3-style resnames and C 10 atom names must not truncate."""
        from benchmarks.astex_entropy.tools import _extract_flexaidds_ligand_pdb

        lines = ["REMARK test"]
        for i in range(25):
            serial = 90001 + i
            name = f"C {i}"
            lines.append(
                f"HETATM{serial} {name} RQ3     1      "
                f"{6.0 + i * 0.1:8.3f}{26.0:8.3f}{35.0:8.3f}  1.00  0.00           C  "
            )
        lines.append("END")
        with tempfile.TemporaryDirectory() as tmp:
            complex_pdb = Path(tmp) / "complex.pdb"
            ligand_pdb = Path(tmp) / "ligand.pdb"
            complex_pdb.write_text("\n".join(lines) + "\n")
            self.assertTrue(_extract_flexaidds_ligand_pdb(complex_pdb, ligand_pdb))
            het = [line for line in ligand_pdb.read_text().splitlines() if line.startswith("HETATM")]
            self.assertEqual(len(het), 25)
            resnames = {line[17:20].strip() for line in het}
            self.assertEqual(resnames, {"RQ3"})


if __name__ == "__main__":
    unittest.main()
