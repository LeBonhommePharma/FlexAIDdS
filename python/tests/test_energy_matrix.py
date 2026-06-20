"""Tests for energy_matrix.py — pure Python, no C++ bindings needed."""

import math
import os
import struct
import tempfile

import numpy as np
import pytest

from flexaidds.energy_matrix import (
    MATRIX_256_SIZE,
    SHNN_MAGIC,
    SHNN_VERSION,
    SYBYL_RADII,
    SYBYL_TYPE_NAMES,
    DensityPoint,
    EnergyMatrix,
    MatrixEntry,
    base_to_sybyl,
    decode_256_type,
    encode_256_type,
    parse_dat_file,
    sybyl_to_base,
    write_dat_file,
)


# ── fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_dir():
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def sample_10type_dat(tmp_dir):
    """Create a minimal 10-type scalar energy matrix .dat file."""
    path = os.path.join(tmp_dir, "test_10.dat")
    lines = []
    roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
    val = 1.0
    for i in range(10):
        for j in range(i, 10):
            label = f"{roman[i]}-{roman[j]}".rjust(10)
            lines.append(f"{label} = {val:10.4f}\n")
            val += 0.5
    with open(path, "w") as fh:
        fh.writelines(lines)
    return path


@pytest.fixture
def sample_density_dat(tmp_dir):
    """Create a 3-type matrix with density function entries."""
    path = os.path.join(tmp_dir, "test_density.dat")
    lines = [
        "       1-1 =     1.5000\n",                           # scalar
        "       1-2 = 0.1 -2.0 0.3 0.0 0.5 3.0\n",            # density
        "       1-3 =     0.7500\n",                           # scalar
        "       2-2 = 0.0 1.0 0.2 2.0 0.4 -1.0 0.8 0.5\n",   # density
        "       2-3 =    -1.2500\n",                           # scalar
        "       3-3 =     2.0000\n",                           # scalar
    ]
    with open(path, "w") as fh:
        fh.writelines(lines)
    return path


@pytest.fixture
def sample_256_matrix():
    """Create a simple 256×256 matrix with known values."""
    rng = np.random.RandomState(42)
    matrix = rng.randn(256, 256).astype(np.float64)
    # Symmetrise
    matrix = (matrix + matrix.T) / 2.0
    return EnergyMatrix(256, matrix)


# ── 256-type encoding tests ─────────────────────────────────────────────────

class TestEncode256Type:
    def test_basic_encoding(self):
        code = encode_256_type(base_type=3, donor=True, acceptor=True)
        base, donor, acceptor = decode_256_type(code)
        assert base == 3
        assert donor is True
        assert acceptor is True

    def test_zero_encoding(self):
        code = encode_256_type(0, donor=False, acceptor=False)
        assert code == 0
        base, donor, acceptor = decode_256_type(code)
        assert base == 0
        assert donor is False
        assert acceptor is False

    def test_max_encoding(self):
        code = encode_256_type(63, donor=True, acceptor=True)
        assert code == 255
        base, donor, acceptor = decode_256_type(code)
        assert base == 63
        assert donor is True
        assert acceptor is True

    def test_all_256_codes_roundtrip(self):
        for code in range(256):
            base, donor, acceptor = decode_256_type(code)
            reconstructed = encode_256_type(base, donor, acceptor)
            assert reconstructed == code

    def test_clamping(self):
        code = encode_256_type(100, donor=10, acceptor=True)  # out of range
        base, donor, acceptor = decode_256_type(code)
        assert 0 <= base <= 63
        assert donor is True
        assert acceptor is True


# ── SYBYL projection tests ──────────────────────────────────────────────────

class TestSYBYLProjection:
    def test_sybyl_type_names_length(self):
        assert len(SYBYL_TYPE_NAMES) == 40

    def test_sybyl_radii_keys(self):
        assert set(SYBYL_RADII.keys()) == set(range(1, 41))

    def test_base_to_sybyl_range(self):
        for base in range(64):
            sybyl = base_to_sybyl(base)
            assert 1 <= sybyl <= 40, f"base={base} → sybyl={sybyl} out of range"

    def test_sybyl_to_base_range(self):
        for sybyl in range(1, 41):
            base = sybyl_to_base(sybyl)
            assert 0 <= base <= 63

    def test_roundtrip_identity_for_canonical_types(self):
        # First 22 base types should round-trip through sybyl
        for base in range(22):
            sybyl = base_to_sybyl(base)
            recovered = sybyl_to_base(sybyl)
            assert recovered == base, (
                f"base={base} → sybyl={sybyl} → base={recovered}")

    def test_carbon_ar_hetadj(self):
        # Base type 26 (C_ar_hetadj) maps to SYBYL C.AR (4)
        assert base_to_sybyl(26) == 4

    def test_carbon_pi_bridging(self):
        # Base type 27 (C_pi_bridging) maps to SYBYL C.2 (2)
        assert base_to_sybyl(27) == 2


# ── MatrixEntry evaluate tests ──────────────────────────────────────────────

class TestMatrixEntry:
    def test_scalar_evaluate(self):
        entry = MatrixEntry(type1=0, type2=1, is_scalar=True,
                            scalar_value=2.5)
        assert entry.evaluate(0.0) == 2.5
        assert entry.evaluate(0.5) == 2.5
        assert entry.evaluate(1.0) == 2.5

    def test_density_evaluate_interpolation(self):
        entry = MatrixEntry(type1=0, type2=0, is_scalar=False,
                            density_points=[
                                DensityPoint(0.1, -2.0),
                                DensityPoint(0.3, 0.0),
                                DensityPoint(0.5, 3.0),
                            ])
        # Below first point
        assert entry.evaluate(0.05) == 0.0
        # At first point
        assert entry.evaluate(0.1) == pytest.approx(-2.0, abs=1e-6)
        # Midpoint between first two
        assert entry.evaluate(0.2) == pytest.approx(-1.0, abs=1e-6)
        # At second point
        assert entry.evaluate(0.3) == pytest.approx(0.0, abs=1e-6)
        # Beyond last point
        assert entry.evaluate(0.8) == pytest.approx(3.0, abs=1e-6)

    def test_density_empty(self):
        entry = MatrixEntry(type1=0, type2=0, is_scalar=False,
                            density_points=[])
        assert entry.evaluate(0.5) == 0.0


# ── Legacy .dat I/O tests ───────────────────────────────────────────────────

class TestLegacyDatIO:
    def test_parse_scalar_dat(self, sample_10type_dat):
        mat = EnergyMatrix.from_dat_file(sample_10type_dat)
        assert mat.ntypes == 10
        assert mat.matrix.shape == (10, 10)
        assert mat.is_symmetric

    def test_roundtrip_scalar(self, sample_10type_dat, tmp_dir):
        mat = EnergyMatrix.from_dat_file(sample_10type_dat)
        out_path = os.path.join(tmp_dir, "roundtrip.dat")
        mat.to_dat_file(out_path)
        mat2 = EnergyMatrix.from_dat_file(out_path)
        np.testing.assert_allclose(mat.matrix, mat2.matrix, atol=1e-3)

    def test_parse_density_dat(self, sample_density_dat):
        mat = EnergyMatrix.from_dat_file(sample_density_dat)
        assert mat.ntypes == 3
        # Check scalar entry
        assert mat.entries[(0, 0)].is_scalar
        assert mat.entries[(0, 0)].scalar_value == pytest.approx(1.5)
        # Check density entry
        assert not mat.entries[(0, 1)].is_scalar
        assert len(mat.entries[(0, 1)].density_points) == 3
        # Symmetry
        assert mat.entries[(1, 0)].density_points == mat.entries[(0, 1)].density_points

    def test_evaluate_density_entry(self, sample_density_dat):
        mat = EnergyMatrix.from_dat_file(sample_density_dat)
        # Evaluate density function at known point
        val = mat.evaluate(0, 1, 0.2)  # midpoint of first density
        assert val == pytest.approx(-1.0, abs=1e-6)

    def test_convenience_functions(self, sample_10type_dat, tmp_dir):
        ntypes, matrix = parse_dat_file(sample_10type_dat)
        assert ntypes == 10
        assert matrix.shape == (10, 10)
        out = os.path.join(tmp_dir, "conv.dat")
        write_dat_file(out, ntypes, matrix)
        nt2, m2 = parse_dat_file(out)
        np.testing.assert_allclose(matrix, m2, atol=1e-3)

    def test_invalid_line_count(self, tmp_dir):
        path = os.path.join(tmp_dir, "bad.dat")
        with open(path, "w") as fh:
            fh.write("1.0\n2.0\n3.0\n4.0\n")  # 4 lines, not N*(N+1)/2
        with pytest.raises(ValueError, match="Invalid line count"):
            EnergyMatrix.from_dat_file(path)


# ── 256×256 binary I/O tests ────────────────────────────────────────────────

class TestBinaryIO:
    def test_roundtrip(self, sample_256_matrix, tmp_dir):
        path = os.path.join(tmp_dir, "test.shnn")
        sample_256_matrix.to_binary(path)
        loaded = EnergyMatrix.from_binary(path)
        assert loaded.ntypes == 256
        # float32 precision loss expected
        np.testing.assert_allclose(
            loaded.matrix, sample_256_matrix.matrix, atol=1e-5)

    def test_invalid_magic(self, tmp_dir):
        path = os.path.join(tmp_dir, "bad.shnn")
        with open(path, "wb") as fh:
            fh.write(b"BAAD")
            fh.write(struct.pack("<II", 1, 256))
            fh.write(np.zeros(256 * 256, dtype=np.float32).tobytes())
        with pytest.raises(ValueError, match="Invalid magic"):
            EnergyMatrix.from_binary(path)

    def test_invalid_dimension(self, tmp_dir):
        path = os.path.join(tmp_dir, "bad_dim.shnn")
        with open(path, "wb") as fh:
            fh.write(SHNN_MAGIC)
            fh.write(struct.pack("<II", 1, 128))  # wrong dim
            fh.write(np.zeros(128 * 128, dtype=np.float32).tobytes())
        with pytest.raises(ValueError, match="Expected 256"):
            EnergyMatrix.from_binary(path)

    def test_wrong_ntypes_for_binary(self):
        mat = EnergyMatrix(10, np.zeros((10, 10)))
        with pytest.raises(ValueError, match="Binary format requires"):
            mat.to_binary("/tmp/should_not_exist.shnn")

    def test_binary_header_format(self, sample_256_matrix, tmp_dir):
        path = os.path.join(tmp_dir, "header_check.shnn")
        sample_256_matrix.to_binary(path)
        with open(path, "rb") as fh:
            magic = fh.read(4)
            version, dim = struct.unpack("<II", fh.read(8))
        assert magic == SHNN_MAGIC
        assert version == SHNN_VERSION
        assert dim == 256

    def test_file_size(self, sample_256_matrix, tmp_dir):
        path = os.path.join(tmp_dir, "size_check.shnn")
        sample_256_matrix.to_binary(path)
        expected = 4 + 4 + 4 + 256 * 256 * 4  # magic + version + dim + data
        assert os.path.getsize(path) == expected


# ── projection tests ─────────────────────────────────────────────────────────

class TestProjection:
    def test_project_to_40_shape(self, sample_256_matrix):
        proj = sample_256_matrix.project_to_40()
        assert proj.ntypes == 40
        assert proj.matrix.shape == (40, 40)

    def test_project_to_40_symmetry(self, sample_256_matrix):
        proj = sample_256_matrix.project_to_40()
        assert proj.is_symmetric

    def test_project_requires_256(self):
        mat = EnergyMatrix(10, np.zeros((10, 10)))
        with pytest.raises(ValueError, match="project_to_40 requires"):
            mat.project_to_40()

    def test_project_preserves_sign_pattern(self):
        # All positive 256×256 → all positive 40×40
        mat = EnergyMatrix(256, np.ones((256, 256)))
        proj = mat.project_to_40()
        assert np.all(proj.matrix >= 0)

    def test_projection_averages_correctly(self):
        # Set all entries for base_type=0 (C_sp, SYBYL=1) to 5.0
        mat = np.zeros((256, 256))
        for code_i in range(256):
            base_i = code_i & 0x3F
            if base_i == 0:
                for code_j in range(256):
                    base_j = code_j & 0x3F
                    if base_j == 0:
                        mat[code_i, code_j] = 5.0
        em = EnergyMatrix(256, mat)
        proj = em.project_to_40()
        # SYBYL type 1 (C.1) = base_to_sybyl(0) = 1, 0-indexed = 0
        assert proj.matrix[0, 0] == pytest.approx(5.0, abs=1e-6)


# ── lookup / evaluate tests ─────────────────────────────────────────────────

class TestLookup:
    def test_lookup_matches_numpy(self, sample_256_matrix):
        mat = sample_256_matrix
        for i in range(0, 256, 17):
            for j in range(0, 256, 19):
                assert mat.lookup(i, j) == mat.matrix[i, j]

    def test_evaluate_scalar(self, sample_10type_dat):
        mat = EnergyMatrix.from_dat_file(sample_10type_dat)
        val = mat.evaluate(0, 0, 0.5)
        assert val == pytest.approx(mat.matrix[0, 0], abs=1e-3)


# ── symmetry and utilities ───────────────────────────────────────────────────

class TestUtilities:
    def test_symmetrise(self):
        mat = EnergyMatrix(3, np.array([[1, 2, 3],
                                         [4, 5, 6],
                                         [7, 8, 9]], dtype=np.float64))
        assert not mat.is_symmetric
        mat.symmetrise()
        assert mat.is_symmetric
        assert mat.matrix[0, 1] == pytest.approx(3.0)  # (2+4)/2
        assert mat.matrix[1, 0] == pytest.approx(3.0)

    def test_repr(self, sample_256_matrix):
        r = repr(sample_256_matrix)
        assert "256" in r
        assert "symmetric" in r
