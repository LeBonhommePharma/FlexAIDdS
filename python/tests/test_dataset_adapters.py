"""Unit tests for flexaidds.dataset_adapters (affinity normalization + registry)."""

from __future__ import annotations

import math

import pytest

from flexaidds.dataset_adapters import (
    ADAPTER_REGISTRY,
    TIER_WEIGHTS,
    BindingMOADAdapter,
    DatasetMetadata,
    PDBbindAdapter,
    create_adapter,
    normalize_affinity,
)
from flexaidds.train_256x256 import TEMPERATURE, kB_kcal


R = kB_kcal
LN10 = math.log(10)


def test_normalize_affinity_deltag_passthrough():
    assert normalize_affinity(-8.5, "deltaG") == -8.5
    assert normalize_affinity(-8.5, "Delta_G") == -8.5


def test_normalize_affinity_kd_ki():
    kd = 1e-6  # 1 µM
    expected = R * TEMPERATURE * math.log(kd)
    assert normalize_affinity(kd, "Kd") == pytest.approx(expected)
    assert normalize_affinity(kd, "Ki") == pytest.approx(expected)


def test_normalize_affinity_ic50_cheng_prusoff():
    ic50 = 2e-6
    # Ki ≈ IC50 / 2
    expected = R * TEMPERATURE * math.log(ic50 / 2.0)
    assert normalize_affinity(ic50, "IC50") == pytest.approx(expected)


def test_normalize_affinity_pkd_pki():
    pkd = 9.0  # nM-scale
    expected = -R * TEMPERATURE * LN10 * pkd
    assert normalize_affinity(pkd, "pKd") == pytest.approx(expected)
    assert normalize_affinity(pkd, "pKi") == pytest.approx(expected)


def test_normalize_affinity_pic50():
    pic50 = 8.0
    pkd_approx = pic50 - math.log10(2)
    expected = -R * TEMPERATURE * LN10 * pkd_approx
    assert normalize_affinity(pic50, "pIC50") == pytest.approx(expected)


def test_normalize_affinity_rejects_nonpositive():
    with pytest.raises(ValueError, match="Non-positive"):
        normalize_affinity(0.0, "Kd")
    with pytest.raises(ValueError, match="Non-positive"):
        normalize_affinity(-1.0, "IC50")


def test_normalize_affinity_unknown_unit():
    with pytest.raises(ValueError, match="Unknown affinity unit"):
        normalize_affinity(1.0, "kcal")


def test_normalize_affinity_custom_temperature():
    kd = 1e-9
    t = 310.0
    expected = R * t * math.log(kd)
    assert normalize_affinity(kd, "Kd", temperature=t) == pytest.approx(expected)


def test_tier_weights_ordering():
    assert TIER_WEIGHTS[1] > TIER_WEIGHTS[2] > TIER_WEIGHTS[3] > TIER_WEIGHTS[4]


def test_create_adapter_registry():
    assert create_adapter("itc_187").name() == "itc_187"
    assert create_adapter("binding_moad").name() == "binding_moad"
    assert create_adapter("bindingdb").name() == "bindingdb"
    assert create_adapter("chembl").name() == "chembl"
    assert create_adapter("dude").name() == "dude"
    assert create_adapter("dekois2").name() == "dekois2"

    core = create_adapter("pdbbind_core")
    assert isinstance(core, PDBbindAdapter)
    assert core.name() == "pdbbind_core"

    refined = create_adapter("pdbbind_refined")
    assert refined.name() == "pdbbind_refined"

    general = create_adapter("pdbbind_general")
    assert general.name() == "pdbbind_general"


def test_create_adapter_unknown():
    with pytest.raises(ValueError, match="Unknown dataset adapter"):
        create_adapter("not_a_real_dataset")


def test_adapter_registry_keys_stable():
    expected = {
        "pdbbind_core",
        "pdbbind_refined",
        "pdbbind_general",
        "itc_187",
        "binding_moad",
        "bindingdb",
        "chembl",
        "dude",
        "dekois2",
    }
    assert expected.issubset(set(ADAPTER_REGISTRY.keys()))


def test_dataset_metadata_to_dict():
    meta = DatasetMetadata(
        name="itc_187",
        version="v1.0",
        n_complexes=10,
        reliability_tier=1,
        weight=1.0,
    )
    d = meta.to_dict()
    assert d["name"] == "itc_187"
    assert d["reliability_tier"] == 1
    assert d["n_complexes"] == 10


def test_moad_parse_binding_string():
    adapter = BindingMOADAdapter()
    parsed = adapter._parse_binding_string("Kd = 12.0 nM")
    assert parsed is not None
    measure, value_m = parsed
    assert measure.lower() == "kd"
    assert value_m == pytest.approx(12.0e-9)

    parsed_u = adapter._parse_binding_string("Ki~1.5 uM")
    assert parsed_u is not None
    assert parsed_u[0].lower() == "ki"
    assert parsed_u[1] == pytest.approx(1.5e-6)

    assert adapter._parse_binding_string("no affinity here") is None


def test_dude_is_validation_only():
    dude = create_adapter("dude")
    assert dude.is_training_dataset is False
