"""Tests for flexaidds.docking — exercising Python-only paths."""

import math
import tempfile
import textwrap
from pathlib import Path

import pytest

from flexaidds.docking import Docking, BindingMode, BindingPopulation, Pose


# ── Fixtures ──────────────────────────────────────────────────────────────────

MINIMAL_CONFIG = textwrap.dedent("""\
    PDBNAM receptor.pdb
    INPLIG ligand.mol2
    TEMPER 300
    METOPT GA
    NRGOUT 3
""")

REMARK_PDB_CONTENT = textwrap.dedent("""\
    REMARK Binding Mode:1 Best CF in Binding Mode:-11.23 Binding Mode Frequency:5
    REMARK 0.45678 RMSD to ref. structure
    ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00           N
    END
""")


def _write_config(path: Path, lines: list) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def config_file(tmp_path):
    cfg = tmp_path / "test.inp"
    cfg.write_text(MINIMAL_CONFIG)
    return cfg


@pytest.fixture
def docking(config_file):
    return Docking(str(config_file))


# ── Docking.__init__ / _parse_config ─────────────────────────────────────────

class TestDockingInit:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            Docking(str(tmp_path / "nonexistent.inp"))

    def test_receptor_property(self, docking):
        assert docking.receptor == "receptor.pdb"

    def test_ligand_property(self, docking):
        assert docking.ligand == "ligand.mol2"

    def test_temperature_property(self, docking):
        assert docking.temperature == 300

    def test_optimization_method(self, docking):
        assert docking.optimization_method == "GA"

    def test_repr(self, docking):
        assert "test.inp" in repr(docking)

    def test_list_keys_always_lists(self, config_file):
        """OPTIMZ and FLEXSC should always be lists even when absent."""
        d = Docking(str(config_file))
        assert isinstance(d._config["OPTIMZ"], list)
        assert isinstance(d._config["FLEXSC"], list)

    def test_multiple_optimz(self, tmp_path):
        cfg = tmp_path / "multi.inp"
        cfg.write_text("PDBNAM x.pdb\nOPTIMZ res1\nOPTIMZ res2\n")
        d = Docking(str(cfg))
        assert d._config["OPTIMZ"] == ["res1", "res2"]

    def test_raises_when_config_missing(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Config file not found"):
            Docking(str(tmp_path / "nonexistent.inp"))


# ── _parse_remark_pdb ─────────────────────────────────────────────────────────

class TestParseRemarkPdb:
    def test_valid_remark(self, tmp_path):
        pdb = tmp_path / "pose.pdb"
        pdb.write_text(REMARK_PDB_CONTENT)
        result = Docking._parse_remark_pdb(pdb, temperature=300.0)
        assert result is not None
        mode_idx, pose = result
        assert mode_idx == 1
        assert abs(pose.energy - (-11.23)) < 1e-6
        assert abs(pose.rmsd - 0.45678) < 1e-6

    def test_missing_remark_returns_none(self, tmp_path):
        pdb = tmp_path / "plain.pdb"
        pdb.write_text("ATOM      1  N   GLY A   1       1.0   2.0   3.0\nEND\n")
        assert Docking._parse_remark_pdb(pdb, temperature=300.0) is None

    def test_boltzmann_weight_positive(self, tmp_path):
        pdb = tmp_path / "pose.pdb"
        pdb.write_text(REMARK_PDB_CONTENT)
        _, pose = Docking._parse_remark_pdb(pdb, temperature=300.0)
        assert pose.boltzmann_weight > 0.0

    def test_boltzmann_weight_formula(self, tmp_path):
        pdb = tmp_path / "pose.pdb"
        pdb.write_text(REMARK_PDB_CONTENT)
        _, pose = Docking._parse_remark_pdb(pdb, temperature=300.0)
        kB = 0.001987206
        beta = 1.0 / (kB * 300.0)
        expected = math.exp(-beta * (-11.23))
        assert abs(pose.boltzmann_weight - expected) < 1e-10


# ── Pose ──────────────────────────────────────────────────────────────────────

class TestPose:
    def test_to_dict_keys(self):
        p = Pose(index=1, energy=-10.5, rmsd=0.3, boltzmann_weight=0.4)
        d = p.to_dict()
        assert "index" in d
        assert "energy_kcal_mol" in d
        assert "rmsd_angstrom" in d

    def test_to_dict_values(self):
        pose = Pose(index=3, energy=-8.0, rmsd=0.5, boltzmann_weight=0.25)
        d = pose.to_dict()
        assert d["index"] == 3
        assert d["energy_kcal_mol"] == pytest.approx(-8.0)
        assert d["rmsd_angstrom"] == pytest.approx(0.5)
        assert d["boltzmann_weight"] == pytest.approx(0.25)

    def test_optional_fields_default_none(self):
        p = Pose(index=0, energy=-5.0)
        assert p.rmsd is None
        assert p.coordinates is None

    def test_default_boltzmann_weight_is_zero(self):
        pose = Pose(index=0, energy=-5.0)
        assert pose.boltzmann_weight == 0.0


# ── BindingMode ───────────────────────────────────────────────────────────────

class TestBindingMode:
    def test_free_energy_no_cpp(self):
        mode = BindingMode(cpp_binding_mode=None)
        mode._poses = [Pose(0, -10.5), Pose(1, -9.8), Pose(2, -11.0)]
        # Free energy from StatMechEngine: F = -kT ln Z ≤ min(E)
        assert mode.free_energy < -11.0  # ensemble F is lower than best pose
        assert math.isfinite(mode.free_energy)

    def test_free_energy_single_pose(self):
        mode = BindingMode(cpp_binding_mode=None)
        mode._poses = [Pose(0, -10.5)]
        # Single pose: F ≈ E (no ensemble spread)
        assert abs(mode.free_energy - (-10.5)) < 1e-6

    def test_free_energy_empty_is_inf(self):
        mode = BindingMode()
        assert mode.free_energy == float("inf")

    def test_n_poses_python_path(self):
        mode = BindingMode()
        mode._poses = [Pose(0, -10.0), Pose(1, -9.0)]
        assert mode.n_poses == 2

    def test_len_equals_n_poses(self):
        mode = BindingMode()
        mode._poses = [Pose(i, float(-i)) for i in range(5)]
        assert len(mode) == 5

    def test_repr_contains_n_poses(self):
        mode = BindingMode()
        mode._poses = [Pose(0, -10.0)]
        assert "n_poses=1" in repr(mode)

    def test_enthalpy_computed_without_cpp(self):
        mode = BindingMode(cpp_binding_mode=None)
        mode._poses = [Pose(0, -10.0), Pose(1, -9.0)]
        # Enthalpy is the Boltzmann-weighted average, not inf
        assert math.isfinite(mode.enthalpy)
        assert -10.0 <= mode.enthalpy <= -9.0  # between min and max

    def test_entropy_computed_without_cpp(self):
        mode = BindingMode(cpp_binding_mode=None)
        mode._poses = [Pose(0, -10.0), Pose(1, -9.0)]
        # With 2 different energies, entropy should be > 0
        assert mode.entropy > 0.0

    def test_entropy_zero_for_empty_mode(self):
        mode = BindingMode(cpp_binding_mode=None)
        assert mode.entropy == 0.0

    def test_get_thermodynamics_works_without_cpp(self):
        mode = BindingMode(cpp_binding_mode=None)
        mode._poses = [Pose(0, -10.0), Pose(1, -9.5)]
        thermo = mode.get_thermodynamics()
        assert math.isfinite(thermo.free_energy)
        assert math.isfinite(thermo.entropy)
        assert thermo.temperature == 300.0

    def test_get_thermodynamics_empty_raises(self):
        mode = BindingMode(cpp_binding_mode=None)
        # Empty mode returns inf free energy but doesn't crash
        thermo = mode.get_thermodynamics()
        assert thermo.free_energy == float("inf")


# ── BindingPopulation ─────────────────────────────────────────────────────────

class TestBindingPopulation:
    def _make_mode(self, energy: float) -> BindingMode:
        m = BindingMode()
        m._poses = [Pose(0, energy)]
        return m

    def test_rank_by_free_energy_sorted(self):
        pop = BindingPopulation()
        for e in [-8.0, -11.0, -9.5]:
            pop.add_mode(self._make_mode(e))
        ranked = pop.rank_by_free_energy()
        energies = [m.free_energy for m in ranked]
        assert energies == sorted(energies)

    def test_n_modes(self):
        pop = BindingPopulation()
        assert pop.n_modes == 0
        pop.add_mode(self._make_mode(-10.0))
        assert pop.n_modes == 1

    def test_len_equals_n_modes(self):
        modes = [self._make_mode(-10.0 - i) for i in range(3)]
        pop = BindingPopulation(modes=modes)
        assert len(pop) == 3

    def test_getitem(self):
        modes = [self._make_mode(-10.0 - i) for i in range(3)]
        pop = BindingPopulation(modes=modes)
        assert pop[0] is modes[0]

    def test_iter(self):
        modes = [self._make_mode(-10.0 - i) for i in range(3)]
        pop = BindingPopulation(modes=modes)
        assert list(pop) == modes

    def test_repr(self):
        pop = BindingPopulation(temperature=310.0)
        assert "310.0K" in repr(pop)

    def test_shannon_entropy_empty(self):
        pop = BindingPopulation()
        assert pop.get_shannon_entropy() == 0.0

    def test_shannon_entropy_positive(self):
        modes = [self._make_mode(-10.0 - i * 2.0) for i in range(3)]
        pop = BindingPopulation(modes=modes)
        S = pop.get_shannon_entropy()
        assert S > 0.0

    def test_shannon_entropy_uniform_higher(self):
        """A uniform distribution should have higher Shannon S than a skewed one."""
        # Uniform-ish: similar energies
        uniform_modes = [self._make_mode(-10.0 - i * 0.01) for i in range(5)]
        pop_uniform = BindingPopulation(modes=uniform_modes)

        # Skewed: one dominant low energy
        skewed_modes = [self._make_mode(-50.0)]
        for i in range(4):
            skewed_modes.append(self._make_mode(-5.0 - i))
        pop_skewed = BindingPopulation(modes=skewed_modes)

        assert pop_uniform.get_shannon_entropy() > pop_skewed.get_shannon_entropy()

    def test_deltaG_matrix_empty(self):
        pop = BindingPopulation()
        matrix = pop.get_deltaG_matrix()
        assert matrix == []

    def test_deltaG_matrix_dimensions(self):
        modes = [self._make_mode(-10.0 - i) for i in range(3)]
        pop = BindingPopulation(modes=modes)
        matrix = pop.get_deltaG_matrix()
        assert len(matrix) == 3
        for row in matrix:
            assert len(row) == 3

    def test_deltaG_matrix_antisymmetric(self):
        modes = [self._make_mode(-10.0 - i * 3.0) for i in range(3)]
        pop = BindingPopulation(modes=modes)
        matrix = pop.get_deltaG_matrix()
        for i in range(3):
            assert abs(matrix[i][i]) < 1e-10
            for j in range(i + 1, 3):
                assert abs(matrix[i][j] + matrix[j][i]) < 1e-10


# ── Docking._parse_config – keyword type dispatch ────────────────────────────

class TestDockingParseConfig:
    def _docking(self, tmp_path: Path, lines: list) -> Docking:
        cfg = tmp_path / "test.inp"
        _write_config(cfg, lines)
        return Docking(str(cfg))

    # --- string keywords ----------------------------------------------------

    def test_string_keyword_pdbnam(self, tmp_path):
        d = self._docking(tmp_path, ["PDBNAM receptor.pdb"])
        assert d.receptor == "receptor.pdb"

    def test_string_keyword_inplig(self, tmp_path):
        d = self._docking(tmp_path, ["INPLIG ligand.mol2"])
        assert d.ligand == "ligand.mol2"

    def test_string_keyword_metopt(self, tmp_path):
        d = self._docking(tmp_path, ["METOPT GA"])
        assert d.optimization_method == "GA"

    # --- float keywords -----------------------------------------------------

    def test_float_keyword_clrmsd(self, tmp_path):
        d = self._docking(tmp_path, ["CLRMSD 1.5"])
        assert d._config["CLRMSD"] == pytest.approx(1.5)
        assert isinstance(d._config["CLRMSD"], float)

    def test_float_keyword_extra_tokens_ignored(self, tmp_path):
        """Only the first token after the keyword is used for float keys."""
        d = self._docking(tmp_path, ["CLRMSD 2.0 # comment"])
        assert d._config["CLRMSD"] == pytest.approx(2.0)

    def test_float_keyword_bad_value_stored_as_string(self, tmp_path):
        d = self._docking(tmp_path, ["CLRMSD NOT_A_FLOAT"])
        assert isinstance(d._config["CLRMSD"], str)

    # --- integer keywords ---------------------------------------------------

    def test_int_keyword_temper(self, tmp_path):
        d = self._docking(tmp_path, ["TEMPER 310"])
        assert d.temperature == 310
        assert isinstance(d._config["TEMPER"], int)

    def test_int_keyword_bad_value_stored_as_string(self, tmp_path):
        d = self._docking(tmp_path, ["TEMPER NOT_INT"])
        assert isinstance(d._config["TEMPER"], str)

    # --- list keywords (OPTIMZ / FLEXSC) ------------------------------------

    def test_list_keyword_multiple_entries(self, tmp_path):
        d = self._docking(tmp_path, [
            "OPTIMZ option1",
            "OPTIMZ option2",
            "OPTIMZ option3",
        ])
        assert d._config["OPTIMZ"] == ["option1", "option2", "option3"]

    def test_list_keyword_initialised_empty(self, tmp_path):
        d = self._docking(tmp_path, ["PDBNAM receptor.pdb"])
        assert d._config["OPTIMZ"] == []
        assert d._config["FLEXSC"] == []

    def test_flexsc_collected(self, tmp_path):
        d = self._docking(tmp_path, ["FLEXSC chainA", "FLEXSC chainB"])
        assert d._config["FLEXSC"] == ["chainA", "chainB"]

    # --- flag keywords ------------------------------------------------------

    def test_flag_keyword_exchet(self, tmp_path):
        d = self._docking(tmp_path, ["EXCHET"])
        assert d._config["EXCHET"] is True

    def test_flag_keyword_rotobs(self, tmp_path):
        d = self._docking(tmp_path, ["ROTOBS"])
        assert d._config["ROTOBS"] is True

    def test_flag_keyword_usedee(self, tmp_path):
        d = self._docking(tmp_path, ["USEDEE"])
        assert d._config["USEDEE"] is True

    def test_absent_flag_not_in_config(self, tmp_path):
        d = self._docking(tmp_path, ["PDBNAM receptor.pdb"])
        assert "EXCHET" not in d._config

    # --- comment and blank line handling ------------------------------------

    def test_hash_comment_lines_ignored(self, tmp_path):
        d = self._docking(tmp_path, [
            "# This is a comment",
            "PDBNAM receptor.pdb",
        ])
        assert d._config["PDBNAM"] == "receptor.pdb"

    def test_blank_lines_ignored(self, tmp_path):
        d = self._docking(tmp_path, ["", "   ", "PDBNAM receptor.pdb", ""])
        assert d._config["PDBNAM"] == "receptor.pdb"

    def test_line_shorter_than_6_chars_ignored(self, tmp_path):
        d = self._docking(tmp_path, ["AB", "PDBNAM receptor.pdb"])
        assert d._config["PDBNAM"] == "receptor.pdb"

    # --- unknown keywords ---------------------------------------------------

    def test_unknown_keyword_stored_as_raw_string(self, tmp_path):
        d = self._docking(tmp_path, ["UNKNWN some_value"])
        assert d._config["UNKNWN"] == "some_value"

    def test_unknown_keyword_no_value(self, tmp_path):
        d = self._docking(tmp_path, ["UNKNWN"])
        assert d._config["UNKNWN"] == ""

    # --- full realistic config ----------------------------------------------

    def test_realistic_config(self, tmp_path):
        lines = [
            "# FlexAID config",
            "PDBNAM /data/receptor.pdb",
            "INPLIG /data/ligand.mol2",
            "METOPT GA",
            "TEMPER 300",
            "CLRMSD 2.0",
            "OPTIMZ ROTRAN",
            "OPTIMZ DIHANG",
            "EXCHET",
            "USEDEE",
        ]
        d = self._docking(tmp_path, lines)
        assert d.receptor == "/data/receptor.pdb"
        assert d.ligand == "/data/ligand.mol2"
        assert d.optimization_method == "GA"
        assert d.temperature == 300
        assert d._config["CLRMSD"] == pytest.approx(2.0)
        assert d._config["OPTIMZ"] == ["ROTRAN", "DIHANG"]
        assert d._config["EXCHET"] is True
        assert d._config["USEDEE"] is True

    # --- property accessors when keys absent --------------------------------

    def test_receptor_none_when_absent(self, tmp_path):
        d = self._docking(tmp_path, ["METOPT GA"])
        assert d.receptor is None

    def test_ligand_none_when_absent(self, tmp_path):
        d = self._docking(tmp_path, ["PDBNAM receptor.pdb"])
        assert d.ligand is None

    def test_temperature_none_when_absent(self, tmp_path):
        d = self._docking(tmp_path, ["PDBNAM receptor.pdb"])
        assert d.temperature is None

    def test_optimization_method_none_when_absent(self, tmp_path):
        d = self._docking(tmp_path, ["PDBNAM receptor.pdb"])
        assert d.optimization_method is None


# ── Docking.run – raises FileNotFoundError when binary not found ──────────────

class TestDockingRun:
    def test_run_raises_file_not_found_when_no_binary(self, tmp_path):
        cfg = tmp_path / "test.inp"
        _write_config(cfg, ["PDBNAM receptor.pdb"])
        d = Docking(str(cfg))
        with pytest.raises(FileNotFoundError):
            d.run()


# ── BindingPopulation.compute_global_thermodynamics – needs C++ core ─────────

_core_available: bool
try:
    from flexaidds.thermodynamics import StatMechEngine as _SME
    _SME(300.0)
    _core_available = True
except Exception:
    _core_available = False


class TestComputeGlobalThermodynamics:
    """compute_global_thermodynamics() aggregates all mode pose energies via
    StatMechEngine and returns a Thermodynamics object."""

    def _pop_with_real_modes(self) -> BindingPopulation:
        pop = BindingPopulation()
        for energy in (-10.0, -9.0):
            mode = BindingMode(cpp_binding_mode=None)
            mode._poses.append(Pose(0, energy))
            pop.add_mode(mode)
        return pop

    def test_returns_thermodynamics_object(self):
        from flexaidds.thermodynamics import Thermodynamics
        pop = self._pop_with_real_modes()
        result = pop.compute_global_thermodynamics()
        assert isinstance(result, Thermodynamics)

    def test_free_energy_is_finite(self):
        pop = self._pop_with_real_modes()
        result = pop.compute_global_thermodynamics()
        assert math.isfinite(result.free_energy)

    def test_free_energy_lower_than_best_pose(self):
        pop = self._pop_with_real_modes()
        result = pop.compute_global_thermodynamics()
        # Ensemble F ≤ min(E) for canonical ensemble
        assert result.free_energy <= -9.0


# ── DiFT torsional integration into BindingMode ────────────────────────────────

class TestBindingModeTorsional:
    """The DiFT torsional free-energy / entropy contribution folds into a
    BindingMode additively and opt-in: absent attached potentials the mode is a
    pure configurational ensemble; attaching them adds a first-principles ΔS."""

    @staticmethod
    def _cosine_profile(n=128, k=1, amp=3.0, phase=math.pi):
        """A single-cosine torsional profile over [0, 2π): a confined bond."""
        return [amp * math.cos(k * (2.0 * math.pi * j / n) - phase)
                for j in range(n)]

    def _confined_mode(self):
        mode = BindingMode(cpp_binding_mode=None)
        mode._poses = [Pose(0, -10.0), Pose(1, -9.0)]
        return mode

    def test_no_torsional_by_default(self):
        mode = self._confined_mode()
        assert mode.torsional_free_energy == 0.0
        assert mode.torsional_entropy == 0.0
        assert mode.torsional_score.n_bonds == 0

    def test_free_energy_unchanged_without_potentials(self):
        mode = self._confined_mode()
        f_config = mode.free_energy
        # Attaching an empty list keeps it a pure configurational ensemble.
        mode.set_torsional_potentials([])
        assert mode.free_energy == f_config

    def test_confinement_penalty_is_positive(self):
        mode = self._confined_mode()
        f_config = mode.free_energy
        mode.set_torsional_profiles([self._cosine_profile()])
        # A confined bond loses entropy → −TΔS > 0 → free energy rises.
        assert mode.torsional_free_energy > 0.0
        assert mode.torsional_entropy < 0.0
        assert mode.free_energy > f_config

    def test_energy_term_needs_dihedral_state(self):
        mode = self._confined_mode()
        # Without dihedral angles: only the state-independent −TΔS penalty.
        mode.set_torsional_profiles([self._cosine_profile()])
        assert mode.torsional_score.energy == 0.0
        assert mode.torsional_score.minus_TS > 0.0

    def test_energy_rises_away_from_well_minimum(self):
        # Profile amp·cos(φ − π) = −amp·cos φ: minimum at φ = 0, maximum at φ = π.
        prof = self._cosine_profile(amp=4.0, phase=math.pi)
        m_min = self._confined_mode()
        m_min.set_torsional_profiles([prof], dihedral_angles_rad=[0.0])
        m_top = self._confined_mode()
        m_top.set_torsional_profiles([prof], dihedral_angles_rad=[math.pi])
        assert m_top.torsional_score.energy > m_min.torsional_score.energy + 1.0

    def test_get_thermodynamics_includes_torsional(self):
        mode = self._confined_mode()
        base = mode.get_thermodynamics()
        mode.set_torsional_profiles([self._cosine_profile()])
        withtors = mode.get_thermodynamics()
        assert withtors.free_energy > base.free_energy
        assert withtors.entropy < base.entropy  # torsional entropy is a loss
        assert withtors.temperature == base.temperature

    def test_free_rotor_adds_nothing(self):
        mode = self._confined_mode()
        f_config = mode.free_energy
        # A flat profile (V ≡ const) is a free rotor: zero excess entropy.
        mode.set_torsional_profiles([[0.0] * 64])
        assert abs(mode.torsional_free_energy) < 1e-9
        assert abs(mode.free_energy - f_config) < 1e-9

    def test_two_bonds_penalize_more_than_one(self):
        prof = self._cosine_profile()
        one = self._confined_mode()
        one.set_torsional_profiles([prof])
        two = self._confined_mode()
        two.set_torsional_profiles([prof, prof])
        assert two.torsional_free_energy > one.torsional_free_energy
        assert two.torsional_score.n_bonds == 2

    def test_set_potentials_directly(self):
        from flexaidds.dift import make_bond_torsion
        mode = self._confined_mode()
        rbt = make_bond_torsion(self._cosine_profile(), gene_index=0)
        mode.set_torsional_potentials([rbt], dihedral_angles_rad=[0.0])
        assert mode.torsional_score.n_bonds == 1
        assert mode.torsional_score.minus_TS > 0.0

    def test_result_is_cached_until_reset(self):
        mode = self._confined_mode()
        mode.set_torsional_profiles([self._cosine_profile()])
        first = mode.torsional_score
        assert mode.torsional_score is first  # cached identity
        mode.set_torsional_profiles([self._cosine_profile(), self._cosine_profile()])
        assert mode.torsional_score is not first  # reset on re-attach
