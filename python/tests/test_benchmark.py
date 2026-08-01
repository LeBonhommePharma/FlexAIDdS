"""Tests for the comparative benchmark module."""

import json
import math
import os
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from flexaidds.benchmark import (
    BenchmarkResult,
    BenchmarkSummary,
    BenchmarkSystem,
    MethodResult,
    SystemBenchmarkResult,
    compute_rmsd,
    compute_rmsd_superposed,
    enrichment_factor,
    extract_ligand_coords_from_mmcif,
    extract_ligand_coords_from_pdb,
    ic50_to_dg,
    kendall_tau,
    ki_to_dg,
    load_benchmark_dataset,
    pic50_to_dg,
    r_squared,
    roc_auc,
    run_benchmark,
    save_benchmark_dataset,
    spearman_rho,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_system(tmp_path):
    """BenchmarkSystem with minimal synthetic files."""
    pdb_content = textwrap.dedent("""\
        HETATM    1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00           C
        HETATM    2  N1  LIG A   1       4.000   5.000   6.000  1.00  0.00           N
        HETATM    3  O1  LIG A   1       7.000   8.000   9.000  1.00  0.00           O
        END
    """)
    protein_pdb = tmp_path / "receptor.pdb"
    protein_pdb.write_text("ATOM      1  CA  ALA A   1       0.0  0.0  0.0\nEND\n")
    ligand_mol2 = tmp_path / "ligand.mol2"
    ligand_mol2.write_text("@<TRIPOS>MOLECULE\nlig\n1 0\nSMALL\n\n@<TRIPOS>ATOM\n")
    ref_pdb = tmp_path / "reference.pdb"
    ref_pdb.write_text(pdb_content)

    return BenchmarkSystem(
        system_id="test_sys",
        protein_pdb_path=protein_pdb,
        protein_sequence="MKTAYIAKQ",
        ligand_mol2_path=ligand_mol2,
        ligand_smiles="CCO",
        reference_pose_pdb_path=ref_pdb,
        experimental_dg_kcal_mol=-8.5,
        is_active=True,
        pocket_residues=(10, 25),
    )


@pytest.fixture
def sample_mmcif():
    """Synthetic mmCIF string mimicking Boltz-2 output."""
    return textwrap.dedent("""\
        data_boltz2_prediction
        #
        loop_
        _atom_site.group_PDB
        _atom_site.type_symbol
        _atom_site.label_atom_id
        _atom_site.label_comp_id
        _atom_site.label_asym_id
        _atom_site.Cartn_x
        _atom_site.Cartn_y
        _atom_site.Cartn_z
        ATOM   C  CA  ALA  A  10.0  20.0  30.0
        ATOM   C  CB  ALA  A  11.0  21.0  31.0
        HETATM C  C1  LIG  B   1.000   2.000   3.000
        HETATM N  N1  LIG  B   4.000   5.000   6.000
        HETATM O  O1  LIG  B   7.000   8.000   9.000
        HETATM H  H1  LIG  B   1.500   2.500   3.500
        #
    """)


@pytest.fixture
def three_system_result():
    """Pre-built BenchmarkResult with 3 systems for aggregate tests."""
    systems = []
    for i, (fa_rmsd, b2_rmsd, fa_dg, b2_dg, exp_dg, active) in enumerate([
        (1.5, 2.5, -9.0, -8.0, -8.5, True),
        (3.0, 1.0, -7.0, -10.0, -9.5, False),
        (0.8, 1.8, -11.0, -9.5, -10.0, True),
    ]):
        sys = BenchmarkSystem(
            system_id=f"sys_{i}",
            protein_pdb_path=Path("."),
            protein_sequence="M",
            ligand_mol2_path=Path("."),
            ligand_smiles="C",
            reference_pose_pdb_path=Path("."),
            experimental_dg_kcal_mol=exp_dg,
            is_active=active,
        )
        fa = MethodResult(
            method="flexaidds", system_id=f"sys_{i}",
            best_pose_rmsd_angstrom=fa_rmsd,
            predicted_dg_kcal_mol=fa_dg,
            predicted_score=-fa_dg,  # CF: lower=better, so negate
            wall_time_seconds=10.0,
        )
        b2 = MethodResult(
            method="boltz2", system_id=f"sys_{i}",
            best_pose_rmsd_angstrom=b2_rmsd,
            predicted_dg_kcal_mol=b2_dg,
            predicted_score=-b2_dg,  # pIC50: higher=better
            wall_time_seconds=5.0,
        )
        systems.append(SystemBenchmarkResult(system=sys, flexaidds_result=fa, boltz2_result=b2))

    return BenchmarkResult(systems=tuple(systems))


# ---------------------------------------------------------------------------
# Affinity conversions
# ---------------------------------------------------------------------------


class TestAffinityConversions:
    def test_ki_to_dg(self):
        # Ki = 1 nM → ΔG ≈ -12.3 kcal/mol at 298.15 K
        dg = ki_to_dg(1.0)
        assert -13.0 < dg < -11.5

    def test_ki_to_dg_known_value(self):
        # Ki = 1 M (1e9 nM) → ΔG = RT ln(1) = 0
        dg = ki_to_dg(1e9)
        assert abs(dg) < 0.01

    def test_ic50_to_dg(self):
        # IC50 = 2 nM → Ki ≈ 1 nM → same as ki_to_dg(1)
        dg = ic50_to_dg(2.0)
        assert abs(dg - ki_to_dg(1.0)) < 0.01

    def test_pic50_to_dg(self):
        # pIC50 = 9 → IC50 = 1 nM → Ki ≈ 0.5 nM
        dg = pic50_to_dg(9.0)
        assert dg < -12.0

    def test_pic50_round_trip(self):
        # pIC50 = 6 → IC50 = 1 µM
        dg = pic50_to_dg(6.0)
        assert -9.0 < dg < -7.0


# ---------------------------------------------------------------------------
# Coordinate extraction
# ---------------------------------------------------------------------------


class TestExtractPDB:
    def test_basic(self, sample_system):
        coords = extract_ligand_coords_from_pdb(sample_system.reference_pose_pdb_path)
        assert coords.shape == (3, 3)
        # File order (not name-sorted): C1, N1, O1 as written
        np.testing.assert_allclose(coords[0], [1.0, 2.0, 3.0])
        np.testing.assert_allclose(coords[1], [4.0, 5.0, 6.0])
        np.testing.assert_allclose(coords[2], [7.0, 8.0, 9.0])

    def test_preserves_file_order_not_lexicographic_atom_names(self, tmp_path):
        """FlexAID-style names must stay in write order, not C/O/N lex sort."""
        # File order: O then C then N. Lex sort by name would be C, N, O.
        pdb = tmp_path / "flexaid_style.pdb"
        pdb.write_text(
            "HETATM90001 O  1 HUP     1       1.000   0.000   0.000  1.00  0.00           O\n"
            "HETATM90002 C  0 HUP     1       2.000   0.000   0.000  1.00  0.00           C\n"
            "HETATM90003 N  2 HUP     1       3.000   0.000   0.000  1.00  0.00           N\n"
            "END\n"
        )
        coords = extract_ligand_coords_from_pdb(pdb)
        np.testing.assert_allclose(coords[0], [1.0, 0.0, 0.0])
        np.testing.assert_allclose(coords[1], [2.0, 0.0, 0.0])
        np.testing.assert_allclose(coords[2], [3.0, 0.0, 0.0])

    def test_empty_raises(self, tmp_path):
        empty = tmp_path / "empty.pdb"
        empty.write_text("ATOM      1  CA  ALA A   1       0.0  0.0  0.0\nEND\n")
        with pytest.raises(ValueError, match="No ligand"):
            extract_ligand_coords_from_pdb(empty)


class TestExtractMmcif:
    def test_basic(self, sample_mmcif):
        coords = extract_ligand_coords_from_mmcif(sample_mmcif)
        # HETATM only, excluding H: C1, N1, O1
        assert coords.shape == (3, 3)
        np.testing.assert_allclose(coords[0], [1.0, 2.0, 3.0])

    def test_no_atom_site_raises(self):
        with pytest.raises(ValueError, match="No _atom_site"):
            extract_ligand_coords_from_mmcif("data_empty\n#\n")

    def test_no_ligand_raises(self):
        mmcif = textwrap.dedent("""\
            data_test
            loop_
            _atom_site.group_PDB
            _atom_site.type_symbol
            _atom_site.label_atom_id
            _atom_site.label_comp_id
            _atom_site.label_asym_id
            _atom_site.Cartn_x
            _atom_site.Cartn_y
            _atom_site.Cartn_z
            ATOM C CA ALA A 1.0 2.0 3.0
            #
        """)
        with pytest.raises(ValueError, match="No ligand"):
            extract_ligand_coords_from_mmcif(mmcif)


# ---------------------------------------------------------------------------
# RMSD
# ---------------------------------------------------------------------------


class TestComputeRmsd:
    """Docking RMSD is in-place (no superposition). Placement is the signal."""

    def test_identity(self):
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        assert compute_rmsd(coords, coords) == pytest.approx(0.0, abs=1e-10)

    def test_translation_is_the_error(self):
        ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        # Pure +5 Å translation: docking error is 5 Å, not 0.
        pred = ref + np.array([5.0, 0.0, 0.0])
        assert compute_rmsd(pred, ref) == pytest.approx(5.0, abs=1e-8)
        # Superposed form still reports ~0 (shape match) — different question.
        assert compute_rmsd_superposed(pred, ref) == pytest.approx(0.0, abs=1e-8)

    def test_rotation_is_the_error(self):
        ref = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        pred = ref @ R.T
        # In-place: non-zero. Superposed: ~0.
        assert compute_rmsd(pred, ref) > 0.5
        assert compute_rmsd_superposed(pred, ref) == pytest.approx(0.0, abs=1e-8)

    def test_known_single_atom(self):
        ref = np.array([[0.0, 0.0, 0.0]])
        pred = np.array([[1.0, 0.0, 0.0]])
        assert compute_rmsd(pred, ref) == pytest.approx(1.0, abs=1e-10)

    def test_two_points_known(self):
        ref = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        pred = np.array([[0.0, 0.0, 0.0], [2.0, 1.0, 0.0]])  # 1Å offset on one atom
        # In-place: sqrt((0 + 1)/2) = sqrt(0.5)
        assert compute_rmsd(pred, ref) == pytest.approx(math.sqrt(0.5), abs=1e-10)

    def test_kabsch_would_hide_pocket_miss(self):
        """Regression: superposed RMSD collapsed every 1gpk pose to ~3.156 Å.

        A rigid translation of a multi-atom ligand must stay visible to
        ``compute_rmsd`` and invisible to ``compute_rmsd_superposed``.
        """
        ref = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        miss = ref + np.array([6.5, 0.0, 0.0])  # ligand left the pocket
        hit = ref + np.array([0.1, 0.0, 0.0])  # near-native
        assert compute_rmsd(miss, ref) == pytest.approx(6.5, abs=1e-8)
        assert compute_rmsd(hit, ref) == pytest.approx(0.1, abs=1e-8)
        # Superposed: both look like "same shape" → both ~0 (the old bug).
        assert compute_rmsd_superposed(miss, ref) == pytest.approx(0.0, abs=1e-8)
        assert compute_rmsd_superposed(hit, ref) == pytest.approx(0.0, abs=1e-8)
        # The discriminator the gate needs: miss >> hit.
        assert compute_rmsd(miss, ref) > 10 * compute_rmsd(hit, ref)

    def test_shape_mismatch(self):
        a = np.array([[1.0, 2.0, 3.0]])
        b = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        with pytest.raises(ValueError, match="Shape mismatch"):
            compute_rmsd(a, b)

    def test_empty(self):
        a = np.empty((0, 3))
        assert compute_rmsd(a, a) == 0.0


# ---------------------------------------------------------------------------
# Statistical metrics
# ---------------------------------------------------------------------------


class TestSpearmanRho:
    def test_perfect_positive(self):
        x = [1.0, 2.0, 3.0, 4.0]
        y = [10.0, 20.0, 30.0, 40.0]
        assert spearman_rho(x, y) == pytest.approx(1.0, abs=1e-10)

    def test_perfect_negative(self):
        x = [1.0, 2.0, 3.0, 4.0]
        y = [40.0, 30.0, 20.0, 10.0]
        assert spearman_rho(x, y) == pytest.approx(-1.0, abs=1e-10)

    def test_uncorrelated(self):
        x = [1.0, 2.0, 3.0, 4.0]
        y = [1.0, 3.0, 2.0, 4.0]
        rho = spearman_rho(x, y)
        assert -1.0 <= rho <= 1.0

    def test_too_few(self):
        assert spearman_rho([1.0], [2.0]) == 0.0


class TestKendallTau:
    def test_perfect(self):
        x = [1.0, 2.0, 3.0]
        y = [1.0, 2.0, 3.0]
        assert kendall_tau(x, y) == pytest.approx(1.0)

    def test_inverse(self):
        x = [1.0, 2.0, 3.0]
        y = [3.0, 2.0, 1.0]
        assert kendall_tau(x, y) == pytest.approx(-1.0)


class TestRSquared:
    def test_perfect(self):
        x = [1.0, 2.0, 3.0]
        y = [2.0, 4.0, 6.0]
        assert r_squared(x, y) == pytest.approx(1.0, abs=1e-10)

    def test_weak_correlation(self):
        x = [1.0, 2.0, 3.0, 4.0]
        y = [1.0, -1.0, 1.0, -1.0]
        assert r_squared(x, y) < 0.5


class TestRocAuc:
    def test_perfect(self):
        scores = [10.0, 9.0, 1.0, 0.0]
        labels = [True, True, False, False]
        assert roc_auc(scores, labels) == pytest.approx(1.0)

    def test_random(self):
        scores = [1.0, 2.0, 3.0, 4.0]
        labels = [True, False, True, False]
        auc = roc_auc(scores, labels)
        assert 0.0 <= auc <= 1.0

    def test_inverted(self):
        scores = [0.0, 1.0, 2.0, 3.0]
        labels = [True, True, False, False]
        auc = roc_auc(scores, labels, higher_is_better=True)
        assert auc < 0.5

    def test_empty(self):
        assert roc_auc([], []) == 0.0


class TestEnrichmentFactor:
    def test_perfect(self):
        scores = [10.0, 9.0, 8.0, 1.0, 0.5, 0.1, 0.0, 0.0, 0.0, 0.0]
        labels = [True, True, True, False, False, False, False, False, False, False]
        ef = enrichment_factor(scores, labels, fraction=0.1)
        # Top 10% = 1 item → should be active → EF = 1/(3/10) = 3.33
        assert ef > 3.0

    def test_empty(self):
        assert enrichment_factor([], [], 0.01) == 0.0


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


class TestBenchmarkSystem:
    def test_construction(self, sample_system):
        assert sample_system.system_id == "test_sys"
        assert sample_system.experimental_dg_kcal_mol == -8.5
        assert sample_system.pocket_residues == (10, 25)

    def test_frozen(self, sample_system):
        with pytest.raises(AttributeError):
            sample_system.system_id = "other"

    def test_to_dict_round_trip(self, sample_system):
        d = sample_system.to_dict()
        assert d["system_id"] == "test_sys"
        assert d["pocket_residues"] == [10, 25]
        restored = BenchmarkSystem.from_dict(d)
        assert restored.system_id == sample_system.system_id


class TestMethodResult:
    def test_construction(self):
        mr = MethodResult(
            method="flexaidds", system_id="sys_0",
            best_pose_rmsd_angstrom=1.5, predicted_dg_kcal_mol=-9.0,
        )
        assert mr.method == "flexaidds"
        assert mr.best_pose_rmsd_angstrom == 1.5


# ---------------------------------------------------------------------------
# BenchmarkResult
# ---------------------------------------------------------------------------


class TestBenchmarkResult:
    def test_n_systems(self, three_system_result):
        assert three_system_result.n_systems == 3

    def test_to_records(self, three_system_result):
        records = three_system_result.to_records()
        assert len(records) == 3
        assert records[0]["system_id"] == "sys_0"
        assert records[0]["flexaidds_rmsd"] == 1.5
        assert records[0]["boltz2_rmsd"] == 2.5

    def test_to_json_round_trip(self, three_system_result, tmp_path):
        json_path = tmp_path / "results.json"
        three_system_result.to_json(json_path)
        loaded = BenchmarkResult.from_json(json_path)
        assert loaded.n_systems == 3
        assert loaded.systems[0].flexaidds_result.best_pose_rmsd_angstrom == 1.5

    def test_to_json_string(self, three_system_result):
        text = three_system_result.to_json()
        assert text is not None
        data = json.loads(text)
        assert data["n_systems"] == 3

    def test_to_csv(self, three_system_result):
        csv_text = three_system_result.to_csv()
        assert csv_text is not None
        assert "sys_0" in csv_text
        lines = csv_text.strip().split("\n")
        assert len(lines) == 4  # header + 3 data rows

    def test_summary(self, three_system_result):
        summary = three_system_result.summary()
        assert summary.n_systems == 3
        # FlexAIDdS: rmsds [1.5, 3.0, 0.8], success (< 2.0) = 2/3
        assert summary.flexaidds_success_rate == pytest.approx(2 / 3)
        # Boltz-2: rmsds [2.5, 1.0, 1.8], success (< 2.0) = 2/3
        assert summary.boltz2_success_rate == pytest.approx(2 / 3)
        # Both methods should have timing
        assert summary.flexaidds_mean_time_seconds == pytest.approx(10.0)
        assert summary.boltz2_mean_time_seconds == pytest.approx(5.0)
        # Correlations should be computable
        assert summary.rank_correlation_spearman is not None
        assert summary.rank_correlation_kendall is not None
        assert summary.flexaidds_dg_r_squared is not None
        assert summary.boltz2_dg_r_squared is not None


# ---------------------------------------------------------------------------
# Dataset I/O
# ---------------------------------------------------------------------------


class TestDatasetIO:
    def test_round_trip(self, sample_system, tmp_path):
        json_path = tmp_path / "dataset.json"
        save_benchmark_dataset([sample_system], json_path)
        loaded = load_benchmark_dataset(json_path)
        assert len(loaded) == 1
        assert loaded[0].system_id == "test_sys"
        assert loaded[0].experimental_dg_kcal_mol == -8.5

    def test_relative_paths(self, tmp_path):
        dataset_dir = tmp_path / "benchmark"
        dataset_dir.mkdir()
        (dataset_dir / "proteins").mkdir()
        (dataset_dir / "ligands").mkdir()
        (dataset_dir / "refs").mkdir()

        protein = dataset_dir / "proteins" / "rec.pdb"
        protein.write_text("ATOM 1\n")
        ligand = dataset_dir / "ligands" / "lig.mol2"
        ligand.write_text("@<TRIPOS>\n")
        ref = dataset_dir / "refs" / "ref.pdb"
        ref.write_text("HETATM 1\n")

        data = {
            "systems": [{
                "system_id": "sys1",
                "protein_pdb_path": "proteins/rec.pdb",
                "protein_sequence": "MKT",
                "ligand_mol2_path": "ligands/lig.mol2",
                "ligand_smiles": "C",
                "reference_pose_pdb_path": "refs/ref.pdb",
            }]
        }
        json_path = dataset_dir / "dataset.json"
        with open(json_path, "w") as f:
            json.dump(data, f)

        loaded = load_benchmark_dataset(json_path)
        assert loaded[0].protein_pdb_path == dataset_dir / "proteins" / "rec.pdb"


# ---------------------------------------------------------------------------
# Runner functions (mocked)
# ---------------------------------------------------------------------------


class TestRunBenchmark:
    @patch("flexaidds.benchmark.run_flexaidds")
    @patch("flexaidds.benchmark.run_boltz2")
    def test_both_methods(self, mock_boltz2, mock_flexaidds, sample_system):
        mock_flexaidds.return_value = MethodResult(
            method="flexaidds", system_id="test_sys",
            best_pose_rmsd_angstrom=1.2, predicted_dg_kcal_mol=-9.0,
            predicted_score=-15.0, n_poses=50, wall_time_seconds=120.0,
        )
        mock_boltz2.return_value = MethodResult(
            method="boltz2", system_id="test_sys",
            best_pose_rmsd_angstrom=1.8, predicted_dg_kcal_mol=-8.5,
            predicted_score=7.2, n_poses=5, wall_time_seconds=30.0,
        )

        result = run_benchmark([sample_system])
        assert result.n_systems == 1
        assert result.systems[0].flexaidds_result.best_pose_rmsd_angstrom == 1.2
        assert result.systems[0].boltz2_result.predicted_score == 7.2

    @patch("flexaidds.benchmark.run_flexaidds")
    @patch("flexaidds.benchmark.run_boltz2")
    def test_on_error_skip(self, mock_boltz2, mock_flexaidds, sample_system):
        mock_flexaidds.side_effect = RuntimeError("FlexAID crashed")
        mock_boltz2.return_value = MethodResult(
            method="boltz2", system_id="test_sys",
        )

        result = run_benchmark([sample_system], on_error="skip")
        assert result.systems[0].flexaidds_result is None
        assert result.systems[0].boltz2_result is not None

    @patch("flexaidds.benchmark.run_flexaidds")
    @patch("flexaidds.benchmark.run_boltz2")
    def test_on_error_raise(self, mock_boltz2, mock_flexaidds, sample_system):
        mock_flexaidds.side_effect = RuntimeError("FlexAID crashed")
        with pytest.raises(RuntimeError, match="FlexAID crashed"):
            run_benchmark([sample_system], on_error="raise")

    @patch("flexaidds.benchmark.run_flexaidds")
    @patch("flexaidds.benchmark.run_boltz2")
    def test_progress_callback(self, mock_boltz2, mock_flexaidds, sample_system):
        mock_flexaidds.return_value = MethodResult(method="flexaidds", system_id="test_sys")
        mock_boltz2.return_value = MethodResult(method="boltz2", system_id="test_sys")

        calls = []
        result = run_benchmark(
            [sample_system],
            progress_callback=lambda sid, idx, total: calls.append((sid, idx, total)),
        )
        assert calls == [("test_sys", 0, 1)]

    @patch("flexaidds.benchmark.run_flexaidds")
    @patch("flexaidds.benchmark.run_boltz2")
    def test_single_method(self, mock_boltz2, mock_flexaidds, sample_system):
        mock_boltz2.return_value = MethodResult(method="boltz2", system_id="test_sys")

        result = run_benchmark([sample_system], methods=("boltz2",))
        assert result.systems[0].flexaidds_result is None
        assert result.systems[0].boltz2_result is not None
        mock_flexaidds.assert_not_called()


# ---------------------------------------------------------------------------
# Import test
# ---------------------------------------------------------------------------


class TestImport:
    def test_import_from_package(self):
        from flexaidds import BenchmarkSystem, BenchmarkResult, run_benchmark
        assert BenchmarkSystem is not None
        assert BenchmarkResult is not None
        assert run_benchmark is not None
