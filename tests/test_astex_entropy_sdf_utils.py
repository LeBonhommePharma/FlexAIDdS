from __future__ import annotations

import tempfile
import unittest
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmarks.astex_entropy.sdf_utils import read_first_sdf_mol, smiles_from_sdf
from benchmarks.astex_entropy.tools import _box_from_ligand
from benchmarks.astex_entropy.validation import rmsd_to_reference


MALFORMED_PEPTIDE_SDF = """ALA
  FlexAIDdS DatasetRunner
  Extracted from structure HETATM/peptide records | FLEXAIDDS_LIGAND_EXTRACTOR_V4
  3  2  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.2000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  2  0  0  0  0
M  END
$$$$
"""


class AstexEntropySdfUtilsTests(unittest.TestCase):
    def test_malformed_counts_line_is_read_for_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sdf = Path(tmp) / "ligand.sdf"
            sdf.write_text(MALFORMED_PEPTIDE_SDF)
            mol = read_first_sdf_mol(sdf, allow_unsanitized=True)
            self.assertEqual(mol.GetNumAtoms(), 3)

            center, size = _box_from_ligand(sdf, padding=4.0)
            self.assertAlmostEqual(center[0], 0.6)
            self.assertEqual(size, (12.0, 12.0, 12.0))

    def test_malformed_counts_line_supports_ordered_rmsd_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sdf = Path(tmp) / "ligand.sdf"
            sdf.write_text(MALFORMED_PEPTIDE_SDF)
            self.assertEqual(rmsd_to_reference(sdf, sdf), 0.0)

    def test_smiles_has_openbabel_fallback(self) -> None:
        obabel = shutil.which("obabel")
        if not obabel:
            homebrew_obabel = Path("/opt/homebrew/bin/obabel")
            obabel = str(homebrew_obabel) if homebrew_obabel.exists() else None
        if not obabel:
            self.skipTest("OpenBabel is not installed")
        with tempfile.TemporaryDirectory() as tmp:
            sdf = Path(tmp) / "ligand.sdf"
            sdf.write_text(MALFORMED_PEPTIDE_SDF)
            self.assertTrue(smiles_from_sdf(sdf, obabel))


if __name__ == "__main__":
    unittest.main()
