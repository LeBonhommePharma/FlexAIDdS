"""Tests for DatasetRunner benchmark path resolution and FlexAID output parsing."""

import json
import tempfile
from pathlib import Path

import pytest

from flexaidds.dataset_runner.astex_targets import (
    ASTEX_DIVERSE_CODES,
    ASTEX_NONNATIVE_BY_NAME,
    parse_crossdock_entry_id,
)
from flexaidds.dataset_runner.data_paths import (
    find_ligand_file,
    find_structure_pdb,
    resolve_astex_nonnative_paths,
    resolve_benchmark_paths,
)
from flexaidds.dataset_runner.runner import (
    DatasetConfig,
    DatasetRunner,
    KNOWN_LARGE_DATASETS,
    load_large_dataset_catalog,
)


@pytest.fixture
def astex_diverse_tree(tmp_path: Path) -> Path:
    root = tmp_path / "astex_diverse"
    for code in ("1GPK", "1MQ6"):
        tdir = root / code
        tdir.mkdir(parents=True)
        (tdir / f"{code}_holo.pdb").write_text("ATOM\n")
        (tdir / f"{code}_ligand.sdf").write_text(
            "\n".join(
                [
                    "ligand",
                    "  FlexAIDdS",
                    "",
                    "  2  0  0  0  0  0  0  0  0  0999 V2000",
                    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
                    "    1.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0",
                    "M  END",
                ]
            )
        )
    return root


@pytest.fixture
def astex_nonnative_tree(tmp_path: Path) -> Path:
    root = tmp_path / "astex_nonnative"
    for code in ("1G9V", "2ACE", "1EVE"):
        tdir = root / code
        tdir.mkdir(parents=True)
        (tdir / f"{code}.pdb").write_text("ATOM\n")
        (tdir / f"{code}_ligand.sdf").write_text(
            "\n".join(
                [
                    "ligand",
                    "  FlexAIDdS",
                    "",
                    "  1  0  0  0  0  0  0  0  0  0999 V2000",
                    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
                    "M  END",
                ]
            )
        )
    return root


def test_astex_diverse_codes_count():
    assert len(ASTEX_DIVERSE_CODES) == 85


def test_parse_crossdock_entry_id():
    assert parse_crossdock_entry_id("1G9V_1EVE") == ("1G9V", "1EVE")
    assert parse_crossdock_entry_id("ACE") is None


def test_find_structure_pdb_astex_diverse(astex_diverse_tree: Path):
    rec = find_structure_pdb(astex_diverse_tree, "1gpk", state="holo")
    assert rec is not None
    assert rec.name == "1GPK_holo.pdb"


def test_resolve_astex_nonnative_family_holo(astex_nonnative_tree: Path):
    rec, ligs = resolve_astex_nonnative_paths(astex_nonnative_tree, "ACE", "holo")
    assert rec is not None
    assert rec.name == "1G9V.pdb"
    assert len(ligs) == 1
    assert ligs[0].name == "1G9V_ligand.sdf"


def test_resolve_astex_nonnative_crossdock(astex_nonnative_tree: Path):
    rec, ligs = resolve_astex_nonnative_paths(astex_nonnative_tree, "1G9V_1EVE", "crossdock")
    assert rec is not None
    assert rec.name == "1EVE.pdb"
    assert ligs[0].name == "1G9V_ligand.sdf"


def test_resolve_astex_nonnative_alternative(astex_nonnative_tree: Path):
    rec, ligs = resolve_astex_nonnative_paths(astex_nonnative_tree, "ACE", "alternative")
    assert rec is not None
    assert rec.parent.name in {"2ACE", "1EVE"}
    assert ligs[0].name == "1G9V_ligand.sdf"


def test_scheduled_work_items_yaml_nonnative():
    cfg = DatasetConfig(
        slug="astex_nonnative",
        name="t",
        description="",
        targets=["ACE", "CA2"],
        structural_states=["holo", "apo", "alternative"],
    )
    items = cfg.scheduled_work_items(tier=2)
    assert len(items) == 6
    assert ("ACE", "holo") in items


def test_parse_flexaid_output_cf_and_rmsd(tmp_path: Path):
    ref = tmp_path / "ref.sdf"
    ref.write_text(
        "\n".join(
            [
                "ligand",
                "  FlexAIDdS",
                "",
                "  2  0  0  0  0  0  0  0  0  0999 V2000",
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
                "    1.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0",
                "M  END",
            ]
        )
    )
    pose = tmp_path / "flexaid_0.pdb"
    pose.write_text(
        "\n".join(
            [
                "REMARK CF=-12.50000",
                "HETATM    1  C   LIG A   1       0.100   0.000   0.000",
                "HETATM    2  N   LIG A   1       1.400   0.000   0.000",
                "END",
            ]
        )
    )
    poses = DatasetRunner._parse_flexaid_output(
        tmp_path, "1gpk", "lig", "holo", reference_ligand=ref
    )
    assert len(poses) == 1
    assert poses[0].enthalpy_score == pytest.approx(-12.5)
    assert poses[0].rmsd >= 0.0


def test_load_large_dataset_catalog_has_astex_nonnative():
    catalog = load_large_dataset_catalog("astex_nonnative")
    if not catalog:
        pytest.skip("large_dataset_entry_catalogs.json not present")
    assert len(catalog) >= KNOWN_LARGE_DATASETS["astex_nonnative"]


def test_dry_run_astex_diverse_yaml_targets():
    from flexaidds.dataset_runner.runner import DatasetRunner as DR

    datasets_dir = Path(__file__).resolve().parents[1] / "flexaidds" / "dataset_runner" / "datasets"
    cfg = DR(dry_run=True).load_dataset_config(datasets_dir / "astex_diverse.yaml")
    assert len(cfg.targets) == 85
    assert cfg.targets[0] == "1gpk"