"""Tests for flexaidds.io – core PDB REMARK parsing and ID inference.

All tests here are pure-Python: no C++ extension needed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from flexaidds.io import (
    _coerce_value,
    _normalize_key,
    infer_mode_id,
    infer_pose_rank,
    parse_remark_map,
    parse_pose_result,
    is_ion,
    Atom,
)


# ===========================================================================
# _coerce_value
# ===========================================================================


class TestCoerceValue:
    def test_integer_whole_float(self):
        """'1.0' should become the int 1, not the float 1.0."""
        assert _coerce_value("1.0") == 1
        assert isinstance(_coerce_value("1.0"), int)

    def test_plain_integer(self):
        assert _coerce_value("42") == 42
        assert isinstance(_coerce_value("42"), int)

    def test_negative_integer(self):
        assert _coerce_value("-7") == -7

    def test_real_float(self):
        result = _coerce_value("-10.5")
        assert result == pytest.approx(-10.5)
        assert isinstance(result, float)

    def test_scientific_notation(self):
        result = _coerce_value("1.38e-23")
        assert result == pytest.approx(1.38e-23)

    def test_boolean_true(self):
        assert _coerce_value("true") is True
        assert _coerce_value("True") is True
        assert _coerce_value("TRUE") is True

    def test_boolean_false(self):
        assert _coerce_value("false") is False
        assert _coerce_value("False") is False

    def test_plain_string(self):
        assert _coerce_value("some_value") == "some_value"

    def test_empty_string(self):
        assert _coerce_value("") == ""

    def test_strips_trailing_semicolons_and_commas(self):
        assert _coerce_value("300;") == 300
        assert _coerce_value("-5.0,") == pytest.approx(-5.0)

    def test_whitespace_is_stripped(self):
        assert _coerce_value("  42  ") == 42


# ===========================================================================
# _normalize_key
# ===========================================================================


class TestNormalizeKey:
    """Every alias defined in _normalize_key must map to its canonical name."""

    @pytest.mark.parametrize("raw, expected", [
        ("bindingmode",      "binding_mode"),
        ("binding_mode_id",  "binding_mode"),
        ("modeid",           "binding_mode"),
        ("clusterid",        "cluster_id"),
        ("cf_app",           "cf_app"),
        ("cfapp",            "cf_app"),
        ("app_cf",           "cf_app"),
        ("delta_g",          "free_energy"),
        ("dg",               "free_energy"),
        ("freeenergy",       "free_energy"),
        ("enthalpy_like",    "enthalpy"),
        ("energy_std",       "std_energy"),
        ("sigma_energy",     "std_energy"),
        ("cv",               "heat_capacity"),
        ("temp",             "temperature"),
    ])
    def test_alias(self, raw, expected):
        assert _normalize_key(raw) == expected

    def test_unknown_key_passthrough(self):
        assert _normalize_key("my_custom_key") == "my_custom_key"

    def test_strips_leading_trailing_whitespace(self):
        assert _normalize_key("  cf_app  ") == "cf_app"

    def test_uppercased_input_is_lowercased(self):
        # Keys are normalised to lower-case before alias lookup
        assert _normalize_key("TEMP") == "temperature"

    def test_non_alphanumeric_collapsed_to_underscore(self):
        # e.g. "free energy" → "free_energy"
        result = _normalize_key("free energy")
        assert result == "free_energy"

    def test_consecutive_separators_collapsed(self):
        result = _normalize_key("free__energy")
        assert result == "free_energy"

    def test_leading_trailing_underscores_stripped(self):
        result = _normalize_key("_cf_app_")
        assert result == "cf_app"


# ===========================================================================
# parse_remark_map
# ===========================================================================


class TestParseRemarkMap:
    def test_equals_delimiter(self):
        lines = ["REMARK binding_mode = 2"]
        result = parse_remark_map(lines)
        assert result["binding_mode"] == 2

    def test_colon_delimiter(self):
        lines = ["REMARK CF: -10.5"]
        result = parse_remark_map(lines)
        # "CF" normalises to "cf" (lowercase, no alias)
        assert result.get("cf") == pytest.approx(-10.5)

    def test_space_delimiter(self):
        lines = ["REMARK TEMPERATURE 300.0"]
        result = parse_remark_map(lines)
        assert result.get("temperature") == pytest.approx(300.0)

    def test_non_remark_lines_ignored(self):
        lines = [
            "ATOM      1  C   LIG A   1       0.000   0.000   0.000",
            "REMARK CF = -7.0",
            "END",
        ]
        result = parse_remark_map(lines)
        assert list(result.keys()) == ["cf"]

    def test_blank_remark_payload_ignored(self):
        lines = ["REMARK", "REMARK   ", "REMARK CF = -5.0"]
        result = parse_remark_map(lines)
        assert "cf" in result
        assert len(result) == 1

    def test_multiple_keys(self):
        lines = [
            "REMARK binding_mode = 3",
            "REMARK pose_rank = 1",
            "REMARK CF = -15.2",
            "REMARK free_energy = -14.8",
            "REMARK temperature = 300.0",
        ]
        result = parse_remark_map(lines)
        assert result["binding_mode"] == 3
        assert result["pose_rank"] == 1
        assert result["cf"] == pytest.approx(-15.2)
        assert result["free_energy"] == pytest.approx(-14.8)
        assert result["temperature"] == pytest.approx(300.0)

    def test_alias_applied_during_parse(self):
        """Keys like 'delta_g' should be stored under their canonical name."""
        lines = ["REMARK delta_g = -9.5"]
        result = parse_remark_map(lines)
        assert "free_energy" in result
        assert result["free_energy"] == pytest.approx(-9.5)

    def test_first_occurrence_wins_for_space_delimiter(self):
        """Space-delimited entries: first value for a key is kept."""
        lines = [
            "REMARK CF -5.0",
            "REMARK CF -6.0",
        ]
        result = parse_remark_map(lines)
        assert result["cf"] == pytest.approx(-5.0)

    def test_equals_delimiter_overwrites_earlier_value(self):
        """Equals-delimiter entries are always stored (last wins)."""
        lines = [
            "REMARK CF = -5.0",
            "REMARK CF = -6.0",
        ]
        result = parse_remark_map(lines)
        assert result["cf"] == pytest.approx(-6.0)

    def test_flexaid_compound_binding_mode_remark(self):
        """Classic FlexAID multi-field REMARK extracts mode, CF, frequency."""
        lines = [
            "REMARK Binding Mode:1 Best CF in Binding Mode:-5.500 "
            "Binding Mode Frequency:2"
        ]
        result = parse_remark_map(lines)
        assert result["binding_mode"] == 1
        assert result["cf"] == pytest.approx(-5.5)
        assert result["frequency"] == 2

    def test_integer_coercion(self):
        lines = ["REMARK binding_mode = 4.0"]
        result = parse_remark_map(lines)
        assert result["binding_mode"] == 4
        assert isinstance(result["binding_mode"], int)

    def test_empty_input(self):
        assert parse_remark_map([]) == {}

    def test_no_remark_lines(self):
        lines = ["ATOM      1  C   LIG A   1       0.000   0.000   0.000"]
        assert parse_remark_map(lines) == {}


# ===========================================================================
# infer_mode_id
# ===========================================================================


class TestInferModeId:
    """REMARK-based inference takes priority; filename patterns are the fallback."""

    def test_remark_binding_mode_key(self):
        assert infer_mode_id(Path("x.pdb"), {"binding_mode": 5}) == 5

    def test_remark_mode_key(self):
        assert infer_mode_id(Path("x.pdb"), {"mode": 3}) == 3

    def test_remark_cluster_id_key(self):
        assert infer_mode_id(Path("x.pdb"), {"cluster_id": 7}) == 7

    def test_remark_float_coerced_to_int(self):
        assert infer_mode_id(Path("x.pdb"), {"binding_mode": 2.0}) == 2

    def test_filename_binding_mode_pattern(self):
        assert infer_mode_id(Path("binding_mode_4_pose_1.pdb"), {}) == 4

    def test_filename_mode_pattern(self):
        assert infer_mode_id(Path("mode_7.pdb"), {}) == 7

    def test_filename_cluster_pattern(self):
        assert infer_mode_id(Path("cluster_12.pdb"), {}) == 12

    def test_filename_bm_pattern(self):
        assert infer_mode_id(Path("bm_2_pose_1.pdb"), {}) == 2

    def test_filename_case_insensitive(self):
        assert infer_mode_id(Path("Mode_9_pose_1.pdb"), {}) == 9

    def test_default_when_nothing_matches(self):
        assert infer_mode_id(Path("ligand_docked.pdb"), {}) == 1

    def test_remark_takes_priority_over_filename(self):
        """Even if filename says mode 3, REMARK binding_mode=8 wins."""
        assert infer_mode_id(Path("mode_3.pdb"), {"binding_mode": 8}) == 8


# ===========================================================================
# infer_pose_rank
# ===========================================================================


class TestInferPoseRank:
    def test_remark_pose_rank_key(self):
        assert infer_pose_rank(Path("x.pdb"), {"pose_rank": 2}) == 2

    def test_remark_rank_key(self):
        assert infer_pose_rank(Path("x.pdb"), {"rank": 4}) == 4

    def test_remark_pose_key(self):
        assert infer_pose_rank(Path("x.pdb"), {"pose": 3}) == 3

    def test_remark_model_key(self):
        assert infer_pose_rank(Path("x.pdb"), {"model": 6}) == 6

    def test_remark_float_coerced_to_int(self):
        assert infer_pose_rank(Path("x.pdb"), {"pose_rank": 5.0}) == 5

    def test_filename_pose_pattern(self):
        assert infer_pose_rank(Path("mode_1_pose_3.pdb"), {}) == 3

    def test_filename_conformer_pattern(self):
        assert infer_pose_rank(Path("conformer_7.pdb"), {}) == 7

    def test_filename_model_pattern(self):
        assert infer_pose_rank(Path("model_2.pdb"), {}) == 2

    def test_filename_case_insensitive(self):
        assert infer_pose_rank(Path("Pose_5.pdb"), {}) == 5

    def test_default_when_nothing_matches(self):
        assert infer_pose_rank(Path("ligand_docked.pdb"), {}) == 1

    def test_remark_takes_priority_over_filename(self):
        assert infer_pose_rank(Path("pose_9.pdb"), {"pose_rank": 1}) == 1


# ===========================================================================
# parse_pose_result  (integration of the above through a real temp file)
# ===========================================================================


def _write_pdb(path: Path, remarks: list[str]) -> None:
    lines = [f"REMARK {r}\n" for r in remarks]
    lines += ["ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n", "END\n"]
    path.write_text("".join(lines), encoding="utf-8")


class TestParsePoseResult:
    def test_full_remark_block(self, tmp_path):
        p = tmp_path / "binding_mode_2_pose_1.pdb"
        _write_pdb(p, [
            "binding_mode = 2",
            "pose_rank = 1",
            "CF = -42.5",
            "CF_APP = -41.0",
            "free_energy = -40.8",
            "enthalpy = -39.5",
            "entropy = 0.0035",
            "heat_capacity = 0.12",
            "std_energy = 1.1",
            "temperature = 300.0",
            "rmsd = 1.23",
            "rmsd_sym = 0.95",
        ])
        pose = parse_pose_result(p)
        assert pose.mode_id == 2
        assert pose.pose_rank == 1
        assert pose.cf == pytest.approx(-42.5)
        assert pose.cf_app == pytest.approx(-41.0)
        assert pose.free_energy == pytest.approx(-40.8)
        assert pose.enthalpy == pytest.approx(-39.5)
        assert pose.entropy == pytest.approx(0.0035)
        assert pose.heat_capacity == pytest.approx(0.12)
        assert pose.std_energy == pytest.approx(1.1)
        assert pose.temperature == pytest.approx(300.0)
        assert pose.rmsd_raw == pytest.approx(1.23)
        assert pose.rmsd_sym == pytest.approx(0.95)

    def test_falls_back_to_filename_when_no_remarks(self, tmp_path):
        p = tmp_path / "cluster_5_model_2.pdb"
        _write_pdb(p, ["CF = -8.0"])
        pose = parse_pose_result(p)
        assert pose.mode_id == 5
        assert pose.pose_rank == 2

    def test_cf_falls_back_to_cf_app_when_cf_missing(self, tmp_path):
        """If CF is absent, cf field should be populated from CF_APP."""
        p = tmp_path / "mode_1_pose_1.pdb"
        _write_pdb(p, ["CF_APP = -7.5"])
        pose = parse_pose_result(p)
        assert pose.cf == pytest.approx(-7.5)
        assert pose.cf_app == pytest.approx(-7.5)

    def test_optional_fields_are_none_when_absent(self, tmp_path):
        p = tmp_path / "mode_1_pose_1.pdb"
        _write_pdb(p, ["CF = -5.0"])
        pose = parse_pose_result(p)
        assert pose.free_energy is None
        assert pose.enthalpy is None
        assert pose.entropy is None
        assert pose.rmsd_raw is None
        assert pose.rmsd_sym is None
        assert pose.temperature is None

    def test_remarks_dict_populated(self, tmp_path):
        p = tmp_path / "mode_1_pose_1.pdb"
        _write_pdb(p, ["CF = -5.0", "temperature = 300.0"])
        pose = parse_pose_result(p)
        assert "cf" in pose.remarks
        assert "temperature" in pose.remarks


# ===========================================================================
# is_ion
# ===========================================================================


def _make_atom(record: str, resname: str, name: str = " X  ") -> Atom:
    """Minimal Atom helper for is_ion tests."""
    return Atom(
        serial=1, name=name, altloc="", resname=resname,
        chainid="A", resseq=1, icode="",
        x=0.0, y=0.0, z=0.0,
        occupancy=1.0, bfactor=0.0, element="X",
        record=record,
    )


class TestIsIon:
    def test_mg_hetatm_is_ion(self):
        assert is_ion(_make_atom("HETATM", "MG ")) is True

    def test_zn_hetatm_is_ion(self):
        assert is_ion(_make_atom("HETATM", "ZN ")) is True

    def test_ca_hetatm_is_ion(self):
        assert is_ion(_make_atom("HETATM", "CA ")) is True

    def test_na_hetatm_is_ion(self):
        assert is_ion(_make_atom("HETATM", "NA ")) is True

    def test_fe2_hetatm_is_ion(self):
        assert is_ion(_make_atom("HETATM", "FE2")) is True

    def test_cl_hetatm_is_ion(self):
        assert is_ion(_make_atom("HETATM", "CL ")) is True

    def test_iod_hetatm_is_ion(self):
        assert is_ion(_make_atom("HETATM", "IOD")) is True

    def test_protein_ca_atom_record_not_ion(self):
        """ATOM record CA (Cα) must not be classified as calcium."""
        assert is_ion(_make_atom("ATOM", "ALA", " CA ")) is False

    def test_organic_cofactor_not_ion(self):
        """Heme (HEM) is HETATM but not a single-atom ion."""
        assert is_ion(_make_atom("HETATM", "HEM")) is False

    def test_water_not_ion(self):
        assert is_ion(_make_atom("HETATM", "HOH")) is False

    def test_stripped_resname_matches(self):
        """Residue names may have trailing spaces — strip() must handle them."""
        assert is_ion(_make_atom("HETATM", "MG  ")) is True  # 4-char with spaces

    def test_atom_record_with_ion_resname_not_ion(self):
        """ATOM record even with 'MG' resname → not classified as ion."""
        assert is_ion(_make_atom("ATOM", "MG ")) is False
