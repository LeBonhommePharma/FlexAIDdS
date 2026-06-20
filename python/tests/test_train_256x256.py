"""Tests for the 256×256 soft contact matrix training pipeline.

Tests cover: data structures, parsers, contact enumeration,
frequency matrix construction, inverse Boltzmann, ridge regression,
L-BFGS refinement, CASF validation, and 256→40 projection.

All tests use synthetic data — no PDBbind download required.
"""

import math
import os
import tempfile
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

from flexaidds.train_256x256 import (
    Atom,
    Complex,
    ContactPair,
    TrainingConfig,
    _is_hbond_acceptor_base,
    _is_hbond_donor_base,
    build_contact_matrix,
    build_reference_matrix,
    enumerate_contacts,
    inverse_boltzmann,
    lbfgs_refine,
    parse_mol2_atoms,
    parse_pdb_atoms,
    parse_pdbbind_index,
    ridge_fit,
    validate_casf,
    validate_projection,
    CONTACT_CUTOFF,
    kB_kcal,
    TEMPERATURE,
)
from flexaidds.energy_matrix import (
    EnergyMatrix,
    encode_256_type,
    decode_256_type,
)


# ─── fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_dir():
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture
def simple_pdb(tmp_dir):
    """Minimal PDB with 3 atoms in close proximity."""
    path = tmp_dir / "protein.pdb"
    path.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  N   ALA A   1       1.500   0.000   0.000  1.00  0.00           N  \n"
        "ATOM      3  O   ALA A   1       0.000   1.500   0.000  1.00  0.00           O  \n"
        "END\n"
    )
    return str(path)


@pytest.fixture
def simple_mol2(tmp_dir):
    """Minimal MOL2 with 2 atoms near origin."""
    path = tmp_dir / "ligand.mol2"
    path.write_text(
        "@<TRIPOS>MOLECULE\ntest\n2 0 0 0 0\nSMALL\nNO_CHARGES\n\n"
        "@<TRIPOS>ATOM\n"
        "  1 C1     2.000   0.000   0.000   C.3       1 LIG   0.050\n"
        "  2 N1     0.000   0.000   2.000   N.3       1 LIG  -0.300\n"
        "@<TRIPOS>BOND\n"
    )
    return str(path)


@pytest.fixture
def synthetic_complexes():
    """5 synthetic complexes with known contacts and affinities."""
    rng = np.random.RandomState(42)
    complexes = []
    for i in range(5):
        prot_atoms = [
            Atom(index=k, name=f"CA{k}", element="C",
                 x=rng.uniform(-5, 5), y=rng.uniform(-5, 5),
                 z=rng.uniform(-5, 5), charge=0.0, base_type=2,
                 type_256=encode_256_type(2, False, False))
            for k in range(20)
        ]
        lig_atoms = [
            Atom(index=k, name=f"L{k}", element="N",
                 x=rng.uniform(-2, 2), y=rng.uniform(-2, 2),
                 z=rng.uniform(-2, 2), charge=-0.2, base_type=7,
                 type_256=encode_256_type(7, True, True))
            for k in range(5)
        ]
        contacts = enumerate_contacts(prot_atoms, lig_atoms, cutoff=10.0)
        pkd = 5.0 + i * 0.5
        dg = -kB_kcal * TEMPERATURE * math.log(10) * pkd

        complexes.append(Complex(
            pdb_code=f"test{i:04d}",
            protein_atoms=prot_atoms,
            ligand_atoms=lig_atoms,
            contacts=contacts,
            pKd=pkd,
            deltaG=dg,
        ))
    return complexes


# ─── Atom dataclass tests ─────────────────────────────────────────────────────

class TestAtom:
    def test_coords_property(self):
        a = Atom(1, "CA", "C", 1.0, 2.0, 3.0)
        np.testing.assert_array_equal(a.coords, [1.0, 2.0, 3.0])

    def test_default_type(self):
        a = Atom(0, "X", "X", 0, 0, 0)
        assert a.base_type == 31
        assert a.type_256 == 0


# ─── H-bond role classification ─────────────────────────────────────────────

class TestHBondRoleClassification:
    def test_n_sp3_is_donor_acceptor(self):
        assert _is_hbond_donor_base(7, 0.0) is True
        assert _is_hbond_acceptor_base(7, 0.0) is True

    def test_nitrile_n_is_acceptor_only(self):
        assert _is_hbond_donor_base(5, 0.0) is False
        assert _is_hbond_acceptor_base(5, 0.0) is True

    def test_quaternary_n_is_donor_only(self):
        assert _is_hbond_donor_base(8, 0.5) is True
        assert _is_hbond_acceptor_base(8, 0.5) is False

    def test_carbonyl_o_is_acceptor_only(self):
        assert _is_hbond_donor_base(12, -0.5) is False
        assert _is_hbond_acceptor_base(12, -0.5) is True

    def test_o_sp3_is_donor_acceptor(self):
        assert _is_hbond_donor_base(13, 0.0) is True
        assert _is_hbond_acceptor_base(13, 0.0) is True

    def test_carbon_is_neither_role(self):
        assert _is_hbond_donor_base(2, 0.0) is False
        assert _is_hbond_acceptor_base(2, 0.0) is False


# ─── PDB parser tests ────────────────────────────────────────────────────────

class TestParsePDB:
    def test_basic_parsing(self, simple_pdb):
        atoms = parse_pdb_atoms(simple_pdb)
        assert len(atoms) == 3

    def test_element_types(self, simple_pdb):
        atoms = parse_pdb_atoms(simple_pdb)
        elements = {a.element for a in atoms}
        assert elements == {"C", "N", "O"}

    def test_coordinates(self, simple_pdb):
        atoms = parse_pdb_atoms(simple_pdb)
        ca = atoms[0]
        assert ca.x == pytest.approx(0.0)
        assert ca.y == pytest.approx(0.0)
        assert ca.z == pytest.approx(0.0)

    def test_256_type_assigned(self, simple_pdb):
        atoms = parse_pdb_atoms(simple_pdb)
        for a in atoms:
            assert 0 <= a.type_256 <= 255

    def test_empty_file(self, tmp_dir):
        path = tmp_dir / "empty.pdb"
        path.write_text("END\n")
        atoms = parse_pdb_atoms(str(path))
        assert len(atoms) == 0

    def test_hetatm(self, tmp_dir):
        path = tmp_dir / "het.pdb"
        path.write_text(
            "HETATM    1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00           C  \n"
            "END\n"
        )
        atoms = parse_pdb_atoms(str(path))
        assert len(atoms) == 1


# ─── MOL2 parser tests ───────────────────────────────────────────────────────

class TestParseMOL2:
    def test_basic_parsing(self, simple_mol2):
        atoms = parse_mol2_atoms(simple_mol2)
        assert len(atoms) == 2

    def test_sybyl_type_mapping(self, simple_mol2):
        atoms = parse_mol2_atoms(simple_mol2)
        c_atom = [a for a in atoms if a.element == "C"][0]
        assert c_atom.base_type == 2  # C.3 → base 2

    def test_charge_parsing(self, simple_mol2):
        atoms = parse_mol2_atoms(simple_mol2)
        n_atom = [a for a in atoms if a.element == "N"][0]
        assert n_atom.charge == pytest.approx(-0.3)

    def test_hbond_assignment(self, simple_mol2):
        atoms = parse_mol2_atoms(simple_mol2)
        n_atom = [a for a in atoms if a.element == "N"][0]
        _, donor, acceptor = decode_256_type(n_atom.type_256)
        assert donor is True
        assert acceptor is True

    def test_empty_mol2(self, tmp_dir):
        path = tmp_dir / "empty.mol2"
        path.write_text("@<TRIPOS>MOLECULE\ntest\n0 0 0\nSMALL\n\n@<TRIPOS>ATOM\n@<TRIPOS>BOND\n")
        atoms = parse_mol2_atoms(str(path))
        assert len(atoms) == 0


# ─── PDBbind index parser ────────────────────────────────────────────────────

class TestParsePDBbindIndex:
    def test_basic_parsing(self, tmp_dir):
        path = tmp_dir / "INDEX_general_PL_data.2020"
        path.write_text(
            "# comment line\n"
            "1a0q  2020  2.30  5.40  Kd=0.0040uM\n"
            "1b0r  2019  1.80  6.20  Ki=0.63nM\n"
        )
        affinities = parse_pdbbind_index(str(path))
        assert len(affinities) == 2
        assert affinities["1a0q"] == pytest.approx(5.40)
        assert affinities["1b0r"] == pytest.approx(6.20)

    def test_comment_skipping(self, tmp_dir):
        path = tmp_dir / "index.txt"
        path.write_text("# header\n# more comments\n1abc  2020  1.0  3.5  Kd\n")
        aff = parse_pdbbind_index(str(path))
        assert len(aff) == 1

    def test_empty_file(self, tmp_dir):
        path = tmp_dir / "empty.txt"
        path.write_text("# nothing\n")
        assert parse_pdbbind_index(str(path)) == {}


# ─── contact enumeration ─────────────────────────────────────────────────────

class TestEnumerateContacts:
    def test_close_atoms_detected(self, simple_pdb, simple_mol2):
        prot = parse_pdb_atoms(simple_pdb)
        lig = parse_mol2_atoms(simple_mol2)
        contacts = enumerate_contacts(prot, lig, cutoff=5.0)
        assert len(contacts) > 0

    def test_distant_atoms_excluded(self):
        prot = [Atom(0, "CA", "C", 0, 0, 0, 0, 2, encode_256_type(2, False, False))]
        lig = [Atom(0, "N1", "N", 100, 100, 100, 0, 7, encode_256_type(7, True, True))]
        contacts = enumerate_contacts(prot, lig, cutoff=4.5)
        assert len(contacts) == 0

    def test_contact_distance_accuracy(self):
        prot = [Atom(0, "CA", "C", 0, 0, 0, 0, 2, encode_256_type(2, False, False))]
        lig = [Atom(0, "N1", "N", 3, 0, 0, 0, 7, encode_256_type(7, True, True))]
        contacts = enumerate_contacts(prot, lig, cutoff=5.0)
        assert len(contacts) == 1
        assert contacts[0].distance == pytest.approx(3.0)

    def test_empty_atom_lists(self):
        assert enumerate_contacts([], [], cutoff=5.0) == []

    def test_contact_types_assigned(self):
        t_a = encode_256_type(2, False, False)
        t_b = encode_256_type(7, True, True)
        prot = [Atom(0, "CA", "C", 0, 0, 0, 0, 2, t_a)]
        lig = [Atom(0, "N1", "N", 1, 0, 0, 0, 7, t_b)]
        contacts = enumerate_contacts(prot, lig, cutoff=5.0)
        assert contacts[0].type_a == t_a
        assert contacts[0].type_b == t_b

    def test_brute_force_fallback(self):
        """Test that brute-force fallback gives same results as KD-tree."""
        from flexaidds.train_256x256 import _enumerate_contacts_brute
        prot = [
            Atom(0, "CA", "C", 0, 0, 0, 0, 2, encode_256_type(2, False, False)),
            Atom(1, "CB", "C", 5, 0, 0, 0, 2, encode_256_type(2, False, False)),
        ]
        lig = [
            Atom(0, "N1", "N", 1, 0, 0, 0, 7, encode_256_type(7, True, True)),
        ]
        brute = _enumerate_contacts_brute(prot, lig, cutoff=4.5)
        kdtree = enumerate_contacts(prot, lig, cutoff=4.5)
        assert len(brute) == len(kdtree)


# ─── frequency matrix ────────────────────────────────────────────────────────

class TestBuildContactMatrix:
    def test_shape(self, synthetic_complexes):
        freq = build_contact_matrix(synthetic_complexes)
        assert freq.shape == (256, 256)

    def test_non_negative(self, synthetic_complexes):
        freq = build_contact_matrix(synthetic_complexes)
        assert np.all(freq >= 0)

    def test_symmetric(self, synthetic_complexes):
        freq = build_contact_matrix(synthetic_complexes)
        np.testing.assert_array_equal(freq, freq.T)

    def test_empty(self):
        freq = build_contact_matrix([])
        assert freq.shape == (256, 256)
        assert np.all(freq == 0)


class TestBuildReferenceMatrix:
    def test_shape(self, synthetic_complexes):
        ref = build_reference_matrix(synthetic_complexes)
        assert ref.shape == (256, 256)

    def test_positive(self, synthetic_complexes):
        ref = build_reference_matrix(synthetic_complexes)
        assert np.all(ref > 0)  # Laplace smoothing ensures > 0

    def test_symmetric(self, synthetic_complexes):
        ref = build_reference_matrix(synthetic_complexes)
        np.testing.assert_allclose(ref, ref.T, atol=1e-10)

    def test_empty_returns_ones(self):
        ref = build_reference_matrix([])
        np.testing.assert_array_equal(ref, np.ones((256, 256)))


# ─── inverse Boltzmann ───────────────────────────────────────────────────────

class TestInverseBoltzmann:
    def test_shape(self):
        freq = np.ones((256, 256))
        ref = np.ones((256, 256))
        ib = inverse_boltzmann(freq, ref)
        assert ib.shape == (256, 256)

    def test_equal_frequencies_zero(self):
        """When observed == reference, potential should be ~0."""
        freq = np.ones((256, 256)) * 10.0
        ref = np.ones((256, 256)) * 10.0
        ib = inverse_boltzmann(freq, ref)
        np.testing.assert_allclose(ib, 0.0, atol=1e-10)

    def test_enriched_negative(self):
        """Enriched contacts should be favourable (negative)."""
        freq = np.ones((256, 256))
        ref = np.ones((256, 256))
        freq[10, 20] = 100.0  # enriched
        ib = inverse_boltzmann(freq, ref)
        assert ib[10, 20] < 0

    def test_depleted_positive(self):
        """Depleted contacts should be unfavourable (positive)."""
        freq = np.ones((256, 256)) * 10.0
        ref = np.ones((256, 256)) * 10.0
        freq[10, 20] = 0.5  # below Laplace floor
        ib = inverse_boltzmann(freq, ref)
        assert ib[10, 20] > 0


# ─── ridge regression ────────────────────────────────────────────────────────

class TestRidgeFit:
    def test_shape(self, synthetic_complexes):
        mat = ridge_fit(synthetic_complexes)
        assert mat.shape == (256, 256)

    def test_symmetric(self, synthetic_complexes):
        mat = ridge_fit(synthetic_complexes)
        np.testing.assert_allclose(mat, mat.T, atol=1e-10)

    def test_empty(self):
        mat = ridge_fit([])
        assert mat.shape == (256, 256)
        assert np.all(mat == 0)

    def test_regularisation_effect(self, synthetic_complexes):
        """Higher alpha should produce smaller values."""
        mat_low = ridge_fit(synthetic_complexes, alpha=0.1)
        mat_high = ridge_fit(synthetic_complexes, alpha=100.0)
        assert np.abs(mat_high).mean() <= np.abs(mat_low).mean()


# ─── L-BFGS refinement ──────────────────────────────────────────────────────

class TestLBFGSRefine:
    @pytest.mark.skipif(
        not __import__('flexaidds.train_256x256', fromlist=['HAS_SCIPY']).HAS_SCIPY,
        reason="scipy not available"
    )
    def test_shape(self, synthetic_complexes):
        init = np.zeros((256, 256))
        refined = lbfgs_refine(init, synthetic_complexes, max_iter=5)
        assert refined.shape == (256, 256)

    @pytest.mark.skipif(
        not __import__('flexaidds.train_256x256', fromlist=['HAS_SCIPY']).HAS_SCIPY,
        reason="scipy not available"
    )
    def test_symmetric(self, synthetic_complexes):
        init = np.zeros((256, 256))
        refined = lbfgs_refine(init, synthetic_complexes, max_iter=5)
        np.testing.assert_allclose(refined, refined.T, atol=1e-10)

    def test_too_few_complexes(self):
        """Should return input matrix when < 3 complexes."""
        init = np.ones((256, 256)) * 3.14
        result = lbfgs_refine(init, [], max_iter=10)
        np.testing.assert_array_equal(result, init)


# ─── CASF validation ─────────────────────────────────────────────────────────

class TestValidateCASF:
    def test_empty(self):
        result = validate_casf(np.zeros((256, 256)), [])
        assert result["pearson_r"] == 0.0

    def test_returns_metrics(self, synthetic_complexes):
        mat = np.random.RandomState(42).randn(256, 256)
        result = validate_casf(mat, synthetic_complexes)
        assert "pearson_r" in result
        assert "rmse" in result
        assert "n_complexes" in result

    @pytest.mark.skipif(
        not __import__('flexaidds.train_256x256', fromlist=['HAS_SCIPY']).HAS_SCIPY,
        reason="scipy not available"
    )
    def test_perfect_prediction(self):
        """With a trivially constructed matrix, should get high correlation."""
        t_a = encode_256_type(2, False, False)
        t_b = encode_256_type(7, True, True)
        complexes = []
        for i in range(10):
            dg = -5.0 - i * 0.5
            contacts = [ContactPair(t_a, t_b, 3.0)] * (i + 1)
            complexes.append(Complex(
                pdb_code=f"t{i}", protein_atoms=[], ligand_atoms=[],
                contacts=contacts, pKd=0, deltaG=dg,
            ))
        # Matrix where mat[t_a, t_b] = predicted_per_contact
        # score = mat[t_a,t_b] * n_contacts → should correlate with dG
        mat = np.zeros((256, 256))
        mat[t_a, t_b] = -0.5  # each contact contributes -0.5
        result = validate_casf(mat, complexes)
        # More contacts → more negative score → should correlate with
        # more negative dG
        assert result["pearson_r"] > 0.5


# ─── 256→40 projection validation ────────────────────────────────────────────

class TestValidateProjection:
    def test_with_synthetic_reference(self, tmp_dir):
        """Create a synthetic 10-type reference and validate projection."""
        # Create a small reference .dat
        ref_mat = np.random.RandomState(99).randn(10, 10)
        ref_mat = (ref_mat + ref_mat.T) / 2
        ref = EnergyMatrix(10, ref_mat)
        ref_path = str(tmp_dir / "ref.dat")
        ref.to_dat_file(ref_path)

        # Create 256 matrix (all zeros)
        mat256 = np.zeros((256, 256))
        result = validate_projection(mat256, ref_path)
        assert "projection_r" in result
        assert "projection_rmse" in result


# ─── ContactPair tests ───────────────────────────────────────────────────────

class TestContactPair:
    def test_default_area(self):
        c = ContactPair(10, 20, 3.5)
        assert c.area == 1.0

    def test_attributes(self):
        c = ContactPair(5, 10, 2.5, 1.5)
        assert c.type_a == 5
        assert c.type_b == 10
        assert c.distance == 2.5
        assert c.area == 1.5


# ─── Complex tests ───────────────────────────────────────────────────────────

class TestComplex:
    def test_default_contacts(self):
        c = Complex("test", [], [])
        assert c.contacts == []
        assert c.pKd == 0.0
        assert c.deltaG == 0.0


# ─── full pipeline (synthetic, no real data) ─────────────────────────────────

class TestTrainingPipeline:
    def test_sippl_then_ridge(self, synthetic_complexes):
        """Test the Sippl + ridge combination."""
        freq = build_contact_matrix(synthetic_complexes)
        ref = build_reference_matrix(synthetic_complexes)
        sippl = inverse_boltzmann(freq, ref)
        ridge = ridge_fit(synthetic_complexes)

        combined = 0.7 * sippl + 0.3 * ridge
        combined = (combined + combined.T) / 2.0

        assert combined.shape == (256, 256)
        np.testing.assert_allclose(combined, combined.T, atol=1e-10)

    def test_binary_roundtrip(self, synthetic_complexes, tmp_dir):
        """Test that trained matrix survives binary save/load."""
        freq = build_contact_matrix(synthetic_complexes)
        ref = build_reference_matrix(synthetic_complexes)
        mat = inverse_boltzmann(freq, ref)

        em = EnergyMatrix(256, mat)
        path = str(tmp_dir / "test_trained.bin")
        em.to_binary(path)

        loaded = EnergyMatrix.from_binary(path)
        # float64 → float32 → float64 loses precision
        np.testing.assert_allclose(
            loaded.matrix, mat.astype(np.float32).astype(np.float64),
            atol=1e-6
        )

    def test_projection_after_training(self, synthetic_complexes):
        """Trained 256→40 projection should produce a valid 40×40 matrix."""
        freq = build_contact_matrix(synthetic_complexes)
        ref = build_reference_matrix(synthetic_complexes)
        mat = inverse_boltzmann(freq, ref)

        em = EnergyMatrix(256, mat)
        proj = em.project_to_40()

        assert proj.ntypes == 40
        assert proj.matrix.shape == (40, 40)
        assert proj.is_symmetric
