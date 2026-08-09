import csv
import json
from pathlib import Path

import pytest

from flexaidds.models import DockingResult
from flexaidds.results import load_results


def _write_pdb(path: Path, remarks: list[str]) -> None:
    lines = [f"REMARK {line}\n" for line in remarks]
    lines.extend(
        [
            "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n",
            "END\n",
        ]
    )
    path.write_text("".join(lines), encoding="utf-8")


def test_load_results_groups_poses_by_mode(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            "free_energy = -41.0",
            "enthalpy = -40.0",
            "entropy = 0.0033",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_1_pose_2.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 2",
            "CF = -39.0",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        [
            "binding_mode = 2",
            "pose_rank = 1",
            "CF = -35.0",
            "free_energy = -34.2",
            "temperature = 300.0",
        ],
    )

    result = load_results(tmp_path)

    assert result.n_modes == 2
    assert result.temperature == 300.0
    assert result.binding_modes[0].mode_id == 1
    assert result.binding_modes[0].n_poses == 2
    assert result.binding_modes[0].best_cf == -42.5
    assert result.binding_modes[0].free_energy == -41.0
    assert result.binding_modes[1].mode_id == 2
    assert result.binding_modes[1].best_cf == -35.0


def test_load_results_cluster_remark_mode_identity(tmp_path: Path) -> None:
    """Real production Cluster N REMARK lines must not collapse to mode_id=1."""
    for mid, cf in ((0, -86.1), (1, -69.5), (2, -2.9)):
        _write_pdb(
            tmp_path / f"complex_{mid}.pdb",
            [
                f"CF = {cf}",
                f"Cluster {mid}: Rank (top):{mid} Average CF:{cf - 10} Frequency:10",
            ],
        )
    result = load_results(tmp_path)
    mode_ids = sorted(m.mode_id for m in result.binding_modes)
    assert mode_ids == [0, 1, 2]
    # CF-only files should recompute F (labelled)
    for mode in result.binding_modes:
        assert mode.free_energy is not None
        assert mode.metadata.get("ledger_source") == "ensemble_estimate_from_cf"


def test_load_results_engine_ledger_not_recomputed(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            "free_energy = -41.0",
            "enthalpy = -40.0",
            "entropy = 0.0033",
            "temperature = 300.0",
        ],
    )
    result = load_results(tmp_path)
    mode = result.binding_modes[0]
    assert mode.free_energy == pytest.approx(-41.0)
    assert mode.metadata.get("ledger_source") == "engine_remark"


def test_load_results_uses_filename_heuristics_when_remarks_missing(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "mode_7_pose_3.pdb",
        [
            "CF = -11.5",
            "temperature = 298.15",
        ],
    )
    _write_pdb(
        tmp_path / "cluster-12_conformer-4.pdb",
        [
            "CF = -9.25",
            "temperature = 298.15",
        ],
    )

    result = load_results(tmp_path)
    first_mode = result.binding_modes[0]
    first_pose = first_mode.poses[0]
    second_mode = result.binding_modes[1]
    second_pose = second_mode.poses[0]

    assert first_mode.mode_id == 7
    assert first_pose.pose_rank == 3
    assert first_mode.best_cf == -11.5
    assert second_mode.mode_id == 12
    assert second_pose.pose_rank == 4
    assert second_mode.best_cf == -9.25


def test_from_json_round_trip_string(tmp_path: Path) -> None:
    """to_json() → from_json() round-trips mode-level scalars via string."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            "free_energy = -41.0",
            "enthalpy = -40.0",
            "entropy = 0.0033",
            "heat_capacity = 0.012",
            "std_energy = 1.5",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        [
            "binding_mode = 2",
            "pose_rank = 1",
            "CF = -35.0",
            "free_energy = -34.2",
            "temperature = 300.0",
        ],
    )

    original = load_results(tmp_path)
    json_text = original.to_json()
    restored = DockingResult.from_json(json_text)

    assert restored.n_modes == original.n_modes
    assert restored.temperature == original.temperature
    assert str(restored.source_dir) == str(original.source_dir)

    for orig_mode, rest_mode in zip(original.binding_modes, restored.binding_modes):
        assert rest_mode.mode_id == orig_mode.mode_id
        assert rest_mode.rank == orig_mode.rank
        assert rest_mode.free_energy == orig_mode.free_energy
        assert rest_mode.enthalpy == orig_mode.enthalpy
        assert rest_mode.entropy == orig_mode.entropy
        assert rest_mode.heat_capacity == orig_mode.heat_capacity
        assert rest_mode.std_energy == orig_mode.std_energy
        assert rest_mode.best_cf == orig_mode.best_cf
        assert rest_mode.temperature == orig_mode.temperature


def test_from_json_round_trip_file(tmp_path: Path) -> None:
    """to_json(path) → from_json(path) round-trips via a file on disk."""
    _write_pdb(
        tmp_path / "mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -10.0",
            "free_energy = -9.5",
            "temperature = 310.0",
        ],
    )

    original = load_results(tmp_path)
    json_path = tmp_path / "results.json"
    original.to_json(json_path)

    restored = DockingResult.from_json(json_path)

    assert restored.n_modes == 1
    assert restored.temperature == 310.0
    assert restored.binding_modes[0].free_energy == -9.5
    assert restored.binding_modes[0].best_cf == -10.0


def test_from_json_empty_modes() -> None:
    """from_json handles results with no binding modes."""
    payload = json.dumps({
        "source_dir": "/tmp/empty",
        "temperature": 300.0,
        "n_modes": 0,
        "metadata": {},
        "binding_modes": [],
    })
    restored = DockingResult.from_json(payload)
    assert restored.n_modes == 0
    assert restored.temperature == 300.0


# ---------------------------------------------------------------------------
# Backward compatibility: old REMARK formats
# ---------------------------------------------------------------------------

def test_old_uppercase_remark_keys(tmp_path: Path) -> None:
    """Upper-case REMARK keys (older engine versions) are parsed correctly."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "BINDING_MODE = 1",
            "POSE_RANK = 1",
            "CF = -30.0",
            "FREE_ENERGY = -29.5",
            "TEMPERATURE = 300.0",
        ],
    )
    result = load_results(tmp_path)
    assert result.n_modes == 1
    mode = result.binding_modes[0]
    assert mode.best_cf == -30.0


def test_missing_temperature_remark_falls_back(tmp_path: Path) -> None:
    """When the REMARK Temperature line is absent, temperature is None or falls back."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -22.0",
            # NO temperature line
        ],
    )
    result = load_results(tmp_path)
    # Should not raise; temperature may be None when not provided
    assert result.n_modes == 1
    mode = result.binding_modes[0]
    assert mode.best_cf == -22.0


def test_mixed_remark_versions_same_directory(tmp_path: Path) -> None:
    """Directory mixing old-style and new-style REMARK formats loads all modes."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -20.0",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        [
            "BINDING_MODE = 2",  # old-style upper-case
            "POSE_RANK = 1",
            "CF = -18.0",
            "TEMPERATURE = 300.0",
        ],
    )
    result = load_results(tmp_path)
    assert result.n_modes == 2
    cf_values = {m.best_cf for m in result.binding_modes}
    assert -20.0 in cf_values
    assert -18.0 in cf_values


def test_no_pdb_files_returns_empty_result(tmp_path: Path) -> None:
    """A directory with no PDB files returns an empty DockingResult gracefully."""
    result = load_results(tmp_path)
    assert result.n_modes == 0


def test_partial_remarks_only_cf(tmp_path: Path) -> None:
    """Files with only CF REMARK still load; F is recomputed from CF (labelled)."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -11.0",
            "temperature = 298.15",
        ],
    )
    result = load_results(tmp_path)
    assert result.n_modes == 1
    mode = result.binding_modes[0]
    assert mode.best_cf == -11.0
    # Smoothing: recompute ensemble estimate from CF when ledger REMARKs missing
    assert mode.free_energy is not None
    assert mode.free_energy == pytest.approx(-11.0)  # single sample → F ≈ E
    assert mode.metadata.get("ledger_source") == "ensemble_estimate_from_cf"


def test_multiple_poses_per_mode_aggregated(tmp_path: Path) -> None:
    """Multiple poses from one mode return the best CF as best_cf."""
    for rank, cf in [(1, -25.0), (2, -18.0), (3, -12.0)]:
        _write_pdb(
            tmp_path / f"binding_mode_1_pose_{rank}.pdb",
            [
                "binding_mode = 1",
                f"pose_rank = {rank}",
                f"CF = {cf}",
                "temperature = 300.0",
            ],
        )
    result = load_results(tmp_path)
    assert result.n_modes == 1
    assert result.binding_modes[0].best_cf == -25.0
    assert result.binding_modes[0].n_poses == 3


def test_to_csv_round_trip(tmp_path: Path) -> None:
    """to_csv() produces valid CSV that round-trips mode IDs and CF values."""
    for mode_id, cf in [(1, -30.0), (2, -25.5)]:
        _write_pdb(
            tmp_path / f"binding_mode_{mode_id}_pose_1.pdb",
            [
                f"binding_mode = {mode_id}",
                "pose_rank = 1",
                f"CF = {cf}",
                "temperature = 300.0",
            ],
        )
    result = load_results(tmp_path)
    csv_text = result.to_csv()

    assert csv_text.strip()
    lines = [l for l in csv_text.splitlines() if l.strip()]
    # Header + 2 data rows
    assert len(lines) == 3

    reader = csv.DictReader(csv_text.splitlines())
    rows = list(reader)
    assert len(rows) == 2
    cfs = {float(r["best_cf"]) for r in rows}
    assert -30.0 in cfs
    assert -25.5 in cfs


def test_subdirectory_pdb_files_loaded(tmp_path: Path) -> None:
    """PDB files in sub-directories are discovered recursively."""
    subdir = tmp_path / "run1"
    subdir.mkdir()
    _write_pdb(
        subdir / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -17.0",
            "temperature = 300.0",
        ],
    )
    result = load_results(tmp_path)
    assert result.n_modes == 1
    assert result.binding_modes[0].best_cf == -17.0


# ===========================================================================
# Chunk 0 — statistical-mechanics claim firewall (loader side)
#
# C++ is the single source of truth for this vocabulary (LIB/statmech.h,
# kScientificProvenanceSchemaVersion = 2).  These tests pin the Python loader
# to it: the emitted election field must drive mode order, and no PDB may
# self-authorize a physical claim.
# ===========================================================================

from flexaidds.thermodynamics import (  # noqa: E402
    ClaimValidity,
    EnergyDomain,
    EnsembleMeasure,
    ReferenceState,
)

# Structurally valid, non-filler receipts (>= 3 distinct nibbles, 64 hex).
_ENERGY_RECEIPT = "sha256:" + "1a2b3c4d5e6f7089" * 4
_MEASURE_RECEIPT = "sha256:" + "9f8e7d6c5b4a3021" * 4
_REFERENCE_RECEIPT = "sha256:" + "0badc0de1234abcd" * 4


def _proxy_firewall_remarks() -> list[str]:
    """The exact block the C++ engine emits today (always proxy-only)."""
    return [
        "thermo_schema_version = 2",
        "thermo_claim_validity = proxy_only",
        "thermo_energy_domain = cf_arbitrary_units",
        "thermo_ensemble_measure = optimizer_samples",
        "thermo_reference_state = bound_only",
    ]


def test_firewall_remarks_are_lifted_into_typed_provenance(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            *_proxy_firewall_remarks(),
            "proxy_free_energy = -41.0",
            "free_energy = -41.0",
            "soft_beta_G = -12.5",
            "enthalpy = -40.0",
            "temperature = 300.0",
        ],
    )

    mode = load_results(tmp_path).binding_modes[0]
    prov = mode.scientific_provenance

    assert prov.schema_version == 2
    assert prov.energy_domain is EnergyDomain.CF_ARBITRARY_UNITS
    assert prov.ensemble_measure is EnsembleMeasure.OPTIMIZER_SAMPLES
    assert prov.reference_state is ReferenceState.BOUND_ONLY
    assert mode.claim_validity is ClaimValidity.PROXY_ONLY
    assert mode.proxy_free_energy == pytest.approx(-41.0)
    assert mode.soft_beta_G == pytest.approx(-12.5)
    # Pose-level fields carry the same evidence.
    assert mode.poses[0].soft_beta_G == pytest.approx(-12.5)
    assert mode.poses[0].claim_validity is ClaimValidity.PROXY_ONLY


def test_emitted_soft_beta_G_order_beats_legacy_free_energy_order(
    tmp_path: Path,
) -> None:
    """Election order wins when soft_beta_G disagrees with free_energy.

    Mode 1 has the better (lower) legacy free_energy; mode 2 has the better
    (lower) emitted election objective.  The loader must reproduce the engine's
    election, so mode 2 must be rank 1.
    """
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            *_proxy_firewall_remarks(),
            "proxy_free_energy = -90.0",
            "free_energy = -90.0",   # legacy says mode 1 is best
            "soft_beta_G = -1.0",    # election says mode 1 is worst
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        [
            "binding_mode = 2",
            "pose_rank = 1",
            "CF = -35.0",
            *_proxy_firewall_remarks(),
            "proxy_free_energy = -10.0",
            "free_energy = -10.0",
            "soft_beta_G = -99.0",   # election winner
            "temperature = 300.0",
        ],
    )

    result = load_results(tmp_path)
    ranked = {m.rank: m.mode_id for m in result.binding_modes}

    assert ranked[1] == 2
    assert ranked[2] == 1
    assert result.top_mode().mode_id == 2
    assert result.top_mode().soft_beta_G == pytest.approx(-99.0)


def test_legacy_free_energy_order_is_unchanged_without_election_field(
    tmp_path: Path,
) -> None:
    """No soft_beta_G anywhere -> byte-for-byte the historical ordering."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        ["binding_mode = 1", "pose_rank = 1", "CF = -42.5", "free_energy = -10.0"],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        ["binding_mode = 2", "pose_rank = 1", "CF = -35.0", "free_energy = -90.0"],
    )

    result = load_results(tmp_path)
    ranked = {m.rank: m.mode_id for m in result.binding_modes}

    assert ranked[1] == 2  # lowest legacy free_energy
    assert ranked[2] == 1
    assert all(m.soft_beta_G is None for m in result.binding_modes)


@pytest.mark.parametrize(
    "hostile_remarks",
    [
        # A file that simply declares itself physical.
        ["thermo_claim_validity = binding_physical"],
        # Full physical vocabulary but wrong schema version.
        [
            "thermo_schema_version = 1",
            "thermo_claim_validity = binding_physical",
            "thermo_energy_domain = calibrated_kcal_per_mol",
            "thermo_ensemble_measure = enumerated_microstates",
            "thermo_reference_state = matched_association_cycle",
            f"thermo_energy_provenance = {_ENERGY_RECEIPT}",
            f"thermo_measure_provenance = {_MEASURE_RECEIPT}",
            f"thermo_reference_provenance = {_REFERENCE_RECEIPT}",
        ],
        # Correct schema, but prose instead of receipts.
        [
            "thermo_schema_version = 2",
            "thermo_claim_validity = binding_physical",
            "thermo_energy_domain = calibrated_kcal_per_mol",
            "thermo_ensemble_measure = enumerated_microstates",
            "thermo_reference_state = matched_association_cycle",
            "thermo_energy_provenance = calibrated by the authors",
            "thermo_measure_provenance = exhaustive enumeration",
            "thermo_reference_provenance = matched cycle",
        ],
        # Known-bad digests: empty-content SHA-256 and the historical filler.
        [
            "thermo_schema_version = 2",
            "thermo_energy_domain = calibrated_kcal_per_mol",
            "thermo_ensemble_measure = enumerated_microstates",
            "thermo_reference_state = matched_association_cycle",
            "thermo_energy_provenance = sha256:"
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
            "thermo_measure_provenance = sha256:"
            "3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a",
            f"thermo_reference_provenance = {_REFERENCE_RECEIPT}",
        ],
        # Single-nibble filler digest (fewer than 3 distinct hex values).
        [
            "thermo_schema_version = 2",
            "thermo_energy_domain = calibrated_kcal_per_mol",
            "thermo_ensemble_measure = weighted_quadrature",
            "thermo_reference_state = matched_association_cycle",
            "thermo_energy_provenance = sha256:" + "a" * 64,
            f"thermo_measure_provenance = {_MEASURE_RECEIPT}",
            f"thermo_reference_provenance = {_REFERENCE_RECEIPT}",
        ],
        # Booleans / numbers where evidence strings belong.
        [
            "thermo_schema_version = true",
            "thermo_energy_domain = calibrated_kcal_per_mol",
            "thermo_ensemble_measure = enumerated_microstates",
            "thermo_reference_state = matched_association_cycle",
            "thermo_energy_provenance = 1",
            "thermo_measure_provenance = true",
            "thermo_reference_provenance = 0",
        ],
        # Unknown vocabulary must not be treated as an unlock.
        [
            "thermo_schema_version = 2",
            "thermo_energy_domain = experimental_free_energy",
            "thermo_ensemble_measure = magic_sampler",
            "thermo_reference_state = assumed_standard_state",
            f"thermo_energy_provenance = {_ENERGY_RECEIPT}",
            f"thermo_measure_provenance = {_MEASURE_RECEIPT}",
            f"thermo_reference_provenance = {_REFERENCE_RECEIPT}",
        ],
    ],
)
def test_hostile_pdb_cannot_self_authorize_a_physical_claim(
    tmp_path: Path, hostile_remarks: list[str]
) -> None:
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        ["binding_mode = 1", "pose_rank = 1", "CF = -42.5", *hostile_remarks],
    )

    mode = load_results(tmp_path).binding_modes[0]

    assert mode.claim_validity is ClaimValidity.PROXY_ONLY
    assert not mode.scientific_provenance.allows_canonical_claims()
    assert not mode.scientific_provenance.allows_binding_claims()
    assert mode.poses[0].claim_validity is ClaimValidity.PROXY_ONLY


def test_absent_firewall_block_fails_closed(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        ["binding_mode = 1", "pose_rank = 1", "CF = -42.5", "free_energy = -41.0"],
    )

    mode = load_results(tmp_path).binding_modes[0]

    assert mode.scientific_provenance.schema_version == 2  # dataclass default
    assert mode.scientific_provenance.energy_domain is EnergyDomain.UNCLASSIFIED
    assert mode.claim_validity is ClaimValidity.PROXY_ONLY
    assert mode.proxy_free_energy is None
    assert mode.soft_beta_G is None


def test_mode_provenance_is_the_weakest_of_its_poses(tmp_path: Path) -> None:
    """Disagreeing pose evidence collapses the mode to the proxy-only default."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        ["binding_mode = 1", "pose_rank = 1", "CF = -42.5", *_proxy_firewall_remarks()],
    )
    _write_pdb(
        tmp_path / "binding_mode_1_pose_2.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 2",
            "CF = -39.0",
            "thermo_schema_version = 2",
            "thermo_energy_domain = model_scale",
            "thermo_ensemble_measure = optimizer_samples",
            "thermo_reference_state = none",
        ],
    )

    mode = load_results(tmp_path).binding_modes[0]

    assert mode.scientific_provenance.energy_domain is EnergyDomain.UNCLASSIFIED
    assert mode.claim_validity is ClaimValidity.PROXY_ONLY


@pytest.mark.parametrize("hostile_value", ["true", "false", "null", "[]", "None"])
def test_boolean_or_junk_energy_remark_is_not_a_measurement(
    tmp_path: Path, hostile_value: str
) -> None:
    """``free_energy = true`` must stay absent, never become 1.0."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            f"free_energy = {hostile_value}",
            f"soft_beta_G = {hostile_value}",
            f"proxy_free_energy = {hostile_value}",
        ],
    )

    pose = load_results(tmp_path).binding_modes[0].poses[0]

    assert pose.free_energy is None
    assert pose.soft_beta_G is None
    assert pose.proxy_free_energy is None
