"""Tests for from_dict / from_json round-trip deserialization on model classes.

Priority 4 coverage.  No C++ extension needed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flexaidds.models import BindingModeResult, DockingResult, PoseResult


# ===========================================================================
# PoseResult.from_dict
# ===========================================================================

class TestPoseResultFromDict:
    def test_round_trip_all_fields(self):
        original = PoseResult(
            path=Path("/data/pose_1.pdb"),
            mode_id=2,
            pose_rank=3,
            cf=-42.5,
            cf_app=-41.0,
            rmsd_raw=1.2,
            rmsd_sym=0.9,
            free_energy=-40.0,
            enthalpy=-38.0,
            entropy=0.0033,
            heat_capacity=0.05,
            std_energy=0.3,
            temperature=300.0,
            remarks={"run_id": 7},
        )
        data = {
            "path": str(original.path),
            "mode_id": original.mode_id,
            "pose_rank": original.pose_rank,
            "cf": original.cf,
            "cf_app": original.cf_app,
            "rmsd_raw": original.rmsd_raw,
            "rmsd_sym": original.rmsd_sym,
            "free_energy": original.free_energy,
            "enthalpy": original.enthalpy,
            "entropy": original.entropy,
            "heat_capacity": original.heat_capacity,
            "std_energy": original.std_energy,
            "temperature": original.temperature,
            "remarks": original.remarks,
        }
        restored = PoseResult.from_dict(data)
        assert restored.path == original.path
        assert restored.mode_id == original.mode_id
        assert restored.pose_rank == original.pose_rank
        assert restored.cf == pytest.approx(original.cf)
        assert restored.cf_app == pytest.approx(original.cf_app)
        assert restored.free_energy == pytest.approx(original.free_energy)
        assert restored.temperature == pytest.approx(original.temperature)
        assert restored.remarks == original.remarks

    def test_defaults_for_missing_keys(self):
        restored = PoseResult.from_dict({})
        assert restored.path == Path("")
        assert restored.mode_id == 0
        assert restored.pose_rank == 0
        assert restored.cf is None
        assert restored.free_energy is None
        assert restored.remarks == {}

    def test_accepts_best_pose_path_alias(self):
        restored = PoseResult.from_dict({"best_pose_path": "/data/pose.pdb"})
        assert restored.path == Path("/data/pose.pdb")

    def test_path_field_takes_precedence_over_alias(self):
        restored = PoseResult.from_dict({
            "path": "/real.pdb",
            "best_pose_path": "/alias.pdb",
        })
        assert restored.path == Path("/real.pdb")

    def test_accepts_path_object(self):
        restored = PoseResult.from_dict({"path": Path("/data/pose.pdb")})
        assert restored.path == Path("/data/pose.pdb")


# ===========================================================================
# BindingModeResult.from_dict
# ===========================================================================

class TestBindingModeResultFromDict:
    def test_round_trip_with_poses(self):
        pose_data = {
            "path": "/data/pose_1.pdb",
            "mode_id": 1,
            "pose_rank": 1,
            "cf": -10.0,
        }
        data = {
            "mode_id": 1,
            "rank": 1,
            "poses": [pose_data],
            "free_energy": -9.5,
            "enthalpy": -8.0,
            "entropy": 0.005,
            "best_cf": -10.0,
            "frequency": 15,
            "temperature": 300.0,
            "metadata": {"run_id": 42},
        }
        restored = BindingModeResult.from_dict(data)
        assert restored.mode_id == 1
        assert restored.rank == 1
        assert restored.n_poses == 1
        assert restored.poses[0].cf == pytest.approx(-10.0)
        assert restored.free_energy == pytest.approx(-9.5)
        assert restored.frequency == 15
        assert restored.metadata == {"run_id": 42}

    def test_defaults_for_missing_keys(self):
        restored = BindingModeResult.from_dict({})
        assert restored.mode_id == 0
        assert restored.rank == 0
        assert restored.poses == []
        assert restored.free_energy is None
        assert restored.metadata == {}

    def test_multiple_poses_preserved(self):
        data = {
            "mode_id": 1,
            "rank": 1,
            "poses": [
                {"path": "a.pdb", "mode_id": 1, "pose_rank": 1, "cf": -12.0},
                {"path": "b.pdb", "mode_id": 1, "pose_rank": 2, "cf": -10.0},
            ],
        }
        restored = BindingModeResult.from_dict(data)
        assert restored.n_poses == 2
        assert restored.poses[0].cf == pytest.approx(-12.0)
        assert restored.poses[1].cf == pytest.approx(-10.0)


# ===========================================================================
# DockingResult.from_dict
# ===========================================================================

class TestDockingResultFromDict:
    def test_round_trip_with_nested_modes(self):
        data = {
            "source_dir": "/data/output",
            "temperature": 300.0,
            "metadata": {"n_pose_files": 4},
            "binding_modes": [
                {
                    "mode_id": 1,
                    "rank": 1,
                    "poses": [
                        {"path": "p1.pdb", "mode_id": 1, "pose_rank": 1, "cf": -10.0},
                    ],
                    "free_energy": -9.5,
                },
                {
                    "mode_id": 2,
                    "rank": 2,
                    "poses": [],
                    "free_energy": -7.0,
                },
            ],
        }
        restored = DockingResult.from_dict(data)
        assert restored.source_dir == Path("/data/output")
        assert restored.temperature == pytest.approx(300.0)
        assert restored.n_modes == 2
        assert restored.binding_modes[0].free_energy == pytest.approx(-9.5)
        assert restored.binding_modes[0].n_poses == 1
        assert restored.metadata == {"n_pose_files": 4}

    def test_from_flat_records(self):
        """Accepts flat records from to_records() (no poses key)."""
        data = {
            "source_dir": "/out",
            "binding_modes": [
                {"mode_id": 1, "rank": 1, "free_energy": -10.0, "best_cf": -12.0},
                {"mode_id": 2, "rank": 2, "free_energy": -8.0},
            ],
        }
        restored = DockingResult.from_dict(data)
        assert restored.n_modes == 2
        assert restored.binding_modes[0].free_energy == pytest.approx(-10.0)
        assert restored.binding_modes[0].poses == []
        assert restored.binding_modes[1].free_energy == pytest.approx(-8.0)

    def test_defaults_for_missing_keys(self):
        restored = DockingResult.from_dict({})
        assert restored.source_dir == Path(".")
        assert restored.binding_modes == []
        assert restored.temperature is None
        assert restored.metadata == {}


# ===========================================================================
# DockingResult.from_json – file-based
# ===========================================================================

class TestDockingResultFromJson:
    def _make_result(self) -> DockingResult:
        pose = PoseResult(
            path=Path("/data/pose_1.pdb"),
            mode_id=1,
            pose_rank=1,
            cf=-42.5,
            free_energy=-41.0,
            temperature=300.0,
        )
        mode = BindingModeResult(
            mode_id=1,
            rank=1,
            poses=[pose],
            free_energy=-41.0,
            best_cf=-42.5,
            temperature=300.0,
        )
        return DockingResult(
            source_dir=Path("/data/output"),
            binding_modes=[mode],
            temperature=300.0,
            metadata={"n_pose_files": 1},
        )

    def test_round_trip_via_file(self, tmp_path):
        original = self._make_result()
        json_path = tmp_path / "results.json"
        original.to_json(json_path)

        restored = DockingResult.from_json(json_path)
        assert restored.n_modes == original.n_modes
        assert restored.temperature == pytest.approx(original.temperature)
        assert restored.binding_modes[0].free_energy == pytest.approx(-41.0)
        assert restored.metadata == {"n_pose_files": 1}

    def test_round_trip_via_string(self):
        original = self._make_result()
        json_text = original.to_json()

        restored = DockingResult.from_json(json_text)
        assert restored.n_modes == original.n_modes
        assert restored.temperature == pytest.approx(original.temperature)

    def test_from_json_invalid_string_raises(self):
        with pytest.raises((json.JSONDecodeError, ValueError)):
            DockingResult.from_json("not valid json {{{")

    def test_round_trip_preserves_mode_count(self, tmp_path):
        pose1 = PoseResult(path=Path("a.pdb"), mode_id=1, pose_rank=1, cf=-10.0)
        pose2 = PoseResult(path=Path("b.pdb"), mode_id=2, pose_rank=1, cf=-8.0)
        mode1 = BindingModeResult(mode_id=1, rank=1, poses=[pose1], free_energy=-9.5)
        mode2 = BindingModeResult(mode_id=2, rank=2, poses=[pose2], free_energy=-7.5)
        original = DockingResult(
            source_dir=tmp_path,
            binding_modes=[mode1, mode2],
            temperature=310.0,
        )
        json_path = tmp_path / "out.json"
        original.to_json(json_path)

        restored = DockingResult.from_json(json_path)
        assert restored.n_modes == 2
        assert restored.binding_modes[0].mode_id == 1
        assert restored.binding_modes[1].mode_id == 2


# ===========================================================================
# DockingResult.from_csv
# ===========================================================================

class TestDockingResultFromCsv:
    def _make_result(self) -> DockingResult:
        pose1 = PoseResult(path=Path("a.pdb"), mode_id=1, pose_rank=1, cf=-42.5)
        pose2 = PoseResult(path=Path("b.pdb"), mode_id=2, pose_rank=1, cf=-35.0)
        mode1 = BindingModeResult(
            mode_id=1, rank=1, poses=[pose1],
            free_energy=-41.0, enthalpy=-40.0, entropy=0.0033,
            best_cf=-42.5, temperature=300.0,
        )
        mode2 = BindingModeResult(
            mode_id=2, rank=2, poses=[pose2],
            free_energy=-34.0, best_cf=-35.0, temperature=300.0,
        )
        return DockingResult(
            source_dir=Path("/data/output"),
            binding_modes=[mode1, mode2],
            temperature=300.0,
        )

    def test_round_trip_via_file(self, tmp_path):
        original = self._make_result()
        csv_path = tmp_path / "results.csv"
        original.to_csv(csv_path)

        restored = DockingResult.from_csv(csv_path)
        assert restored.n_modes == 2
        assert restored.binding_modes[0].mode_id == 1
        assert restored.binding_modes[0].free_energy == pytest.approx(-41.0)
        assert restored.binding_modes[0].best_cf == pytest.approx(-42.5)
        assert restored.binding_modes[1].mode_id == 2

    def test_round_trip_via_string(self):
        original = self._make_result()
        csv_text = original.to_csv()

        restored = DockingResult.from_csv(csv_text)
        assert restored.n_modes == 2
        assert restored.binding_modes[0].free_energy == pytest.approx(-41.0)

    def test_none_values_preserved(self):
        mode = BindingModeResult(mode_id=1, rank=1, poses=[], free_energy=-5.0)
        original = DockingResult(source_dir=Path("."), binding_modes=[mode])
        csv_text = original.to_csv()

        restored = DockingResult.from_csv(csv_text)
        assert restored.binding_modes[0].free_energy == pytest.approx(-5.0)
        assert restored.binding_modes[0].enthalpy is None
        assert restored.binding_modes[0].entropy is None

    def test_int_fields_coerced(self):
        original = self._make_result()
        csv_text = original.to_csv()

        restored = DockingResult.from_csv(csv_text)
        assert isinstance(restored.binding_modes[0].mode_id, int)
        assert isinstance(restored.binding_modes[0].rank, int)

    def test_float_fields_coerced(self):
        original = self._make_result()
        csv_text = original.to_csv()

        restored = DockingResult.from_csv(csv_text)
        assert isinstance(restored.binding_modes[0].free_energy, float)
        assert isinstance(restored.binding_modes[0].best_cf, float)

    def test_empty_csv_yields_no_modes(self):
        mode = BindingModeResult(mode_id=1, rank=1, poses=[])
        original = DockingResult(source_dir=Path("."), binding_modes=[])
        csv_text = original.to_csv()

        restored = DockingResult.from_csv(csv_text)
        assert restored.n_modes == 0


# ===========================================================================
# Chunk 0 — claim firewall survives every serialization round trip
# ===========================================================================

from flexaidds.thermodynamics import (  # noqa: E402
    ClaimValidity,
    EnergyDomain,
    EnsembleMeasure,
    ReferenceState,
    ScientificProvenance,
    ThermodynamicBreakdown,
)

_ENERGY_RECEIPT = "sha256:" + "1a2b3c4d5e6f7089" * 4
_MEASURE_RECEIPT = "sha256:" + "9f8e7d6c5b4a3021" * 4
_REFERENCE_RECEIPT = "sha256:" + "0badc0de1234abcd" * 4


def _physical_provenance() -> ScientificProvenance:
    return ScientificProvenance(
        schema_version=2,
        energy_domain=EnergyDomain.CALIBRATED_KCAL_PER_MOL,
        ensemble_measure=EnsembleMeasure.ENUMERATED_MICROSTATES,
        reference_state=ReferenceState.MATCHED_ASSOCIATION_CYCLE,
        energy_provenance=_ENERGY_RECEIPT,
        measure_provenance=_MEASURE_RECEIPT,
        reference_provenance=_REFERENCE_RECEIPT,
    )


def _result_with(provenance: ScientificProvenance) -> DockingResult:
    mode = BindingModeResult(
        mode_id=1,
        rank=1,
        poses=[PoseResult(path=Path("p1.pdb"), mode_id=1, pose_rank=1, cf=-10.0)],
        free_energy=-9.8,
        proxy_free_energy=-9.8,
        soft_beta_G=-3.25,
        best_cf=-10.0,
        temperature=300.0,
        scientific_provenance=provenance,
    )
    return DockingResult(source_dir=Path("/tmp"), binding_modes=[mode])


class TestFirewallFieldsSurviveRoundTrips:
    def test_json_round_trip_keeps_election_and_provenance(self):
        original = _result_with(_physical_provenance())

        restored = DockingResult.from_json(original.to_json())
        mode = restored.binding_modes[0]

        assert mode.soft_beta_G == pytest.approx(-3.25)
        assert mode.proxy_free_energy == pytest.approx(-9.8)
        assert mode.scientific_provenance == _physical_provenance()
        assert mode.claim_validity is ClaimValidity.BINDING_PHYSICAL
        # The placeholder pose inherits the same evidence, not a blank default.
        assert mode.poses[0].soft_beta_G == pytest.approx(-3.25)
        assert mode.poses[0].claim_validity is ClaimValidity.BINDING_PHYSICAL

    def test_csv_round_trip_keeps_election_and_provenance(self):
        original = _result_with(_physical_provenance())

        restored = DockingResult.from_csv(original.to_csv())
        mode = restored.binding_modes[0]

        assert mode.soft_beta_G == pytest.approx(-3.25)
        assert mode.proxy_free_energy == pytest.approx(-9.8)
        assert mode.scientific_provenance == _physical_provenance()

    def test_from_dict_flat_record_keeps_election_and_provenance(self):
        original = _result_with(_physical_provenance())
        payload = {
            "source_dir": "/tmp",
            "binding_modes": original.to_records(),
        }

        restored = DockingResult.from_dict(payload)
        mode = restored.binding_modes[0]

        assert mode.soft_beta_G == pytest.approx(-3.25)
        assert mode.scientific_provenance == _physical_provenance()

    def test_records_are_json_serialisable(self):
        records = _result_with(_physical_provenance()).to_records()
        # Would raise TypeError if a dataclass object leaked into a record.
        json.dumps(records)


class TestHostileDeserialization:
    @pytest.mark.parametrize(
        "provenance_payload",
        [
            # Serialized validity trying to self-authorize with no evidence.
            {"schema_version": 2, "claim_validity": "binding_physical"},
            # Right vocabulary, wrong schema version.
            {
                "schema_version": 1,
                "energy_domain": "calibrated_kcal_per_mol",
                "ensemble_measure": "enumerated_microstates",
                "reference_state": "matched_association_cycle",
                "energy_provenance": _ENERGY_RECEIPT,
                "measure_provenance": _MEASURE_RECEIPT,
                "reference_provenance": _REFERENCE_RECEIPT,
                "claim_validity": "binding_physical",
            },
            # Non-string receipts.
            {
                "schema_version": 2,
                "energy_domain": "calibrated_kcal_per_mol",
                "ensemble_measure": "enumerated_microstates",
                "reference_state": "matched_association_cycle",
                "energy_provenance": 1,
                "measure_provenance": True,
                "reference_provenance": [_REFERENCE_RECEIPT],
                "claim_validity": "binding_physical",
            },
            # Known-bad digests (empty content + historical filler).
            {
                "schema_version": 2,
                "energy_domain": "calibrated_kcal_per_mol",
                "ensemble_measure": "enumerated_microstates",
                "reference_state": "matched_association_cycle",
                "energy_provenance": "sha256:e3b0c44298fc1c149afbf4c8996fb924"
                "27ae41e4649b934ca495991b7852b855",
                "measure_provenance": "sha256:3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c"
                "5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a",
                "reference_provenance": _REFERENCE_RECEIPT,
            },
            # Low-entropy filler digest.
            {
                "schema_version": 2,
                "energy_domain": "calibrated_kcal_per_mol",
                "ensemble_measure": "weighted_quadrature",
                "reference_state": "matched_association_cycle",
                "energy_provenance": "sha256:" + "a" * 64,
                "measure_provenance": _MEASURE_RECEIPT,
                "reference_provenance": _REFERENCE_RECEIPT,
            },
            # Whitespace-padded receipt: no trimming is performed.
            {
                "schema_version": 2,
                "energy_domain": "calibrated_kcal_per_mol",
                "ensemble_measure": "enumerated_microstates",
                "reference_state": "matched_association_cycle",
                "energy_provenance": " " + _ENERGY_RECEIPT,
                "measure_provenance": _MEASURE_RECEIPT + "\n",
                "reference_provenance": _REFERENCE_RECEIPT,
            },
            # Upper-case scheme prefix is not the literal "sha256:".
            {
                "schema_version": 2,
                "energy_domain": "calibrated_kcal_per_mol",
                "ensemble_measure": "enumerated_microstates",
                "reference_state": "matched_association_cycle",
                "energy_provenance": _ENERGY_RECEIPT.replace("sha256:", "SHA256:"),
                "measure_provenance": _MEASURE_RECEIPT,
                "reference_provenance": _REFERENCE_RECEIPT,
            },
            # Hostile availability-shaped values in place of the whole block.
            {"schema_version": 2, "energy_domain": True},
            {"schema_version": 2, "energy_domain": 1},
            {"schema_version": 2, "energy_domain": None},
            {"schema_version": 2, "energy_domain": []},
            {"schema_version": "2", "energy_domain": "calibrated_kcal_per_mol"},
        ],
    )
    def test_json_payload_cannot_self_authorize(self, provenance_payload):
        payload = {
            "source_dir": "/tmp",
            "binding_modes": [
                {
                    "mode_id": 1,
                    "rank": 1,
                    "free_energy": -9.8,
                    "best_pose_path": "p1.pdb",
                    "scientific_provenance": provenance_payload,
                }
            ],
        }

        restored = DockingResult.from_json(json.dumps(payload))
        mode = restored.binding_modes[0]

        assert mode.claim_validity is ClaimValidity.PROXY_ONLY
        assert not mode.scientific_provenance.allows_canonical_claims()
        assert not mode.scientific_provenance.allows_binding_claims()
        assert mode.poses[0].claim_validity is ClaimValidity.PROXY_ONLY

    @pytest.mark.parametrize(
        "raw", [True, False, 1, 0, None, [], "true", "not-a-dict", ""]
    )
    def test_non_mapping_provenance_falls_closed(self, raw):
        """Availability-shaped junk (true/1/null/[]/"true") is never a witness."""
        payload = {
            "source_dir": "/tmp",
            "binding_modes": [
                {
                    "mode_id": 1,
                    "rank": 1,
                    "best_pose_path": "p1.pdb",
                    "scientific_provenance": raw,
                }
            ],
        }

        mode = DockingResult.from_json(json.dumps(payload)).binding_modes[0]

        assert mode.scientific_provenance == ScientificProvenance()
        assert mode.claim_validity is ClaimValidity.PROXY_ONLY

    def test_csv_claim_validity_column_cannot_self_authorize(self):
        original = _result_with(ScientificProvenance())
        csv_text = original.to_csv().replace("proxy_only", "binding_physical")

        mode = DockingResult.from_csv(csv_text).binding_modes[0]

        assert mode.claim_validity is ClaimValidity.PROXY_ONLY

    def test_provenance_is_immutable_after_construction(self):
        provenance = ScientificProvenance()

        with pytest.raises(Exception):
            provenance.energy_domain = EnergyDomain.CALIBRATED_KCAL_PER_MOL
        with pytest.raises(Exception):
            provenance.energy_provenance = _ENERGY_RECEIPT

        assert provenance.claim_validity is ClaimValidity.PROXY_ONLY

    def test_breakdown_correction_downgrades_ledger_like_cpp(self):
        """Mirror of C++ provenance_for_breakdown: corrections force proxy-only."""
        clean = ThermodynamicBreakdown(
            G_config_kcal_mol=-9.8, provenance=_physical_provenance()
        )
        assert clean.claim_validity is ClaimValidity.BINDING_PHYSICAL

        for kwargs in (
            {"has_vib": True},
            {"has_natural": True},
            {"has_other": True},
            {"G_vib_kcal_mol": -0.5},
            {"G_natural_kcal_mol": 0.25},
            {"G_other_kcal_mol": 1e-9},
        ):
            corrected = ThermodynamicBreakdown(
                G_config_kcal_mol=-9.8,
                provenance=_physical_provenance(),
                **kwargs,
            )
            assert corrected.claim_validity is ClaimValidity.PROXY_ONLY
            # Numerics are untouched by the downgrade.
            assert corrected.G_config_kcal_mol == pytest.approx(-9.8)
            for key, value in kwargs.items():
                assert getattr(corrected, key) == value
