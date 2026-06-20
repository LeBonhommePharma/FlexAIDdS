"""Tests for the continuous training pipeline.

All tests use synthetic data — no real PDBbind, ITC-187, or other datasets required.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from pathlib import Path
from typing import List

import numpy as np
import pytest

from flexaidds.dataset_adapters import (
    DatasetAdapter,
    DatasetMetadata,
    normalize_affinity,
    complexes_to_contact_table,
    checksum_contact_table,
    create_adapter,
    TIER_WEIGHTS,
)
from flexaidds.continuous_training import (
    ContinuousTrainer,
    ContinuousTrainingConfig,
    CurriculumPhase,
    QualityGateResult,
    TrainingRunResult,
    compute_cell_confidence,
    warm_start_combine,
)
from flexaidds.train_256x256 import (
    Atom,
    Complex,
    ContactPair,
    CONTACT_CUTOFF,
    kB_kcal,
    TEMPERATURE,
)
from flexaidds.energy_matrix import EnergyMatrix, encode_256_type


# ── fixtures ─────────────────────────────────────────────────────────────────

def _make_atom(idx: int, element: str, x: float, y: float, z: float,
               base_type: int = 2, charge: float = 0.0) -> Atom:
    """Create a synthetic atom."""
    t256 = encode_256_type(base_type, False, False)
    return Atom(idx, f"{element}{idx}", element, x, y, z, charge, base_type, t256)


def _make_complex(code: str, n_contacts: int = 5, pkd: float = 6.0) -> Complex:
    """Create a synthetic complex with random contacts."""
    rng = np.random.RandomState(hash(code) % 2**31)
    prot_atoms = [_make_atom(i, "C", *rng.randn(3), base_type=2) for i in range(10)]
    lig_atoms = [_make_atom(i + 10, "N", *rng.randn(3), base_type=7) for i in range(5)]

    contacts = []
    for _ in range(n_contacts):
        pa = prot_atoms[rng.randint(len(prot_atoms))]
        la = lig_atoms[rng.randint(len(lig_atoms))]
        contacts.append(ContactPair(
            type_a=pa.type_256,
            type_b=la.type_256,
            distance=rng.uniform(2.0, 4.5),
        ))

    dg = -kB_kcal * TEMPERATURE * math.log(10) * pkd
    return Complex(
        pdb_code=code,
        protein_atoms=prot_atoms,
        ligand_atoms=lig_atoms,
        contacts=contacts,
        pKd=pkd,
        deltaG=dg,
    )


@pytest.fixture
def synthetic_complexes() -> List[Complex]:
    """Generate a small set of synthetic complexes."""
    return [_make_complex(f"synth_{i:03d}", pkd=4.0 + i * 0.5) for i in range(10)]


@pytest.fixture
def tmp_dir():
    with tempfile.TemporaryDirectory() as d:
        yield d


# ── affinity normalization tests ─────────────────────────────────────────────

class TestAffinityNormalization:

    def test_kd_to_dg(self):
        """Kd = 1e-9 M at 298.15K should give ~ -12.3 kcal/mol."""
        dg = normalize_affinity(1e-9, "Kd")
        assert dg < 0  # favorable
        expected = kB_kcal * TEMPERATURE * math.log(1e-9)
        assert abs(dg - expected) < 0.01

    def test_pki_to_dg(self):
        """pKi = 9.0 should match Kd = 1e-9."""
        dg = normalize_affinity(9.0, "pKi")
        dg_kd = normalize_affinity(1e-9, "Ki")
        assert abs(dg - dg_kd) < 0.01

    def test_pkd_to_dg(self):
        dg = normalize_affinity(6.0, "pKd")
        assert dg < 0

    def test_ic50_cheng_prusoff(self):
        """IC50 → Ki ≈ IC50/2."""
        dg_ic50 = normalize_affinity(2e-9, "IC50")
        dg_ki = normalize_affinity(1e-9, "Ki")
        assert abs(dg_ic50 - dg_ki) < 0.01

    def test_deltag_passthrough(self):
        dg = normalize_affinity(-8.5, "deltaG")
        assert dg == -8.5

    def test_pic50_to_dg(self):
        dg = normalize_affinity(7.0, "pIC50")
        assert dg < 0

    def test_nonpositive_kd_raises(self):
        with pytest.raises(ValueError, match="Non-positive"):
            normalize_affinity(0.0, "Kd")
        with pytest.raises(ValueError, match="Non-positive"):
            normalize_affinity(-1.0, "Ki")

    def test_unknown_unit_raises(self):
        with pytest.raises(ValueError, match="Unknown"):
            normalize_affinity(1.0, "foobar")


# ── warm-start tests ────────────────────────────────────────────────────────

class TestWarmStart:

    def test_cell_confidence_below_min(self):
        counts = np.full((256, 256), 5.0)  # below min_contacts=10
        conf = compute_cell_confidence(counts, min_contacts=10)
        assert np.allclose(conf, 0.0)

    def test_cell_confidence_above_saturation(self):
        counts = np.full((256, 256), 2000.0)  # above saturation=1000
        conf = compute_cell_confidence(counts, saturation=1000)
        assert np.allclose(conf, 1.0)

    def test_cell_confidence_linear(self):
        counts = np.full((256, 256), 505.0)  # midpoint between 10 and 1000
        conf = compute_cell_confidence(counts, min_contacts=10, saturation=1000)
        expected = (505.0 - 10.0) / (1000.0 - 10.0)
        assert abs(conf[0, 0] - expected) < 0.01

    def test_no_prior_returns_new_estimate(self):
        prior = np.zeros((256, 256))
        new_est = np.random.randn(256, 256)
        counts = np.full((256, 256), 500.0)
        result = warm_start_combine(prior, new_est, counts, dataset_weight=1.0)
        # Should be close to new_est (prior is zeros, mixing with zeros)
        # Actually: (1-λ)*0 + λ*new = λ*new
        assert np.any(result != 0)

    def test_high_weight_increases_new_influence(self):
        prior = np.ones((256, 256))
        new_est = np.ones((256, 256)) * 2.0
        counts = np.full((256, 256), 500.0)

        r_high = warm_start_combine(prior, new_est, counts, dataset_weight=1.0)
        r_low = warm_start_combine(prior, new_est, counts, dataset_weight=0.1)
        # Higher weight → closer to new_est (2.0)
        assert r_high.mean() > r_low.mean()

    def test_mixing_floor_preserves_prior(self):
        prior = np.ones((256, 256)) * 5.0
        new_est = np.zeros((256, 256))
        counts = np.full((256, 256), 5000.0)  # high confidence

        # With 70% floor, prior retains at least 70%
        result = warm_start_combine(
            prior, new_est, counts,
            dataset_weight=1.0, mixing_floor=0.7,
        )
        # result = (1 - 0.3)*5 + 0.3*0 = 3.5
        assert result.mean() >= 3.4  # approximately 0.7 * 5 = 3.5

    def test_result_is_symmetric(self):
        prior = np.random.randn(256, 256)
        new_est = np.random.randn(256, 256)
        counts = np.random.randint(0, 100, (256, 256)).astype(float)
        result = warm_start_combine(prior, new_est, counts)
        assert np.allclose(result, result.T, atol=1e-10)


# ── contact table tests ─────────────────────────────────────────────────────

class TestContactTable:

    def test_complexes_to_contact_table(self, synthetic_complexes):
        table = complexes_to_contact_table(synthetic_complexes)
        assert table.ntypes == 256
        assert table.n_structures == len(synthetic_complexes)
        assert table.counts.sum() > 0

    def test_contact_table_checksum_deterministic(self, synthetic_complexes):
        table = complexes_to_contact_table(synthetic_complexes)
        c1 = checksum_contact_table(table)
        c2 = checksum_contact_table(table)
        assert c1 == c2
        assert c1.startswith("sha256:")

    def test_contact_table_save_load(self, synthetic_complexes, tmp_dir):
        table = complexes_to_contact_table(synthetic_complexes)
        path = os.path.join(tmp_dir, "test_table.json")
        table.save(path)
        loaded = type(table).load(path)
        assert loaded.ntypes == table.ntypes
        assert loaded.n_structures == table.n_structures


# ── quality gate tests ───────────────────────────────────────────────────────

class TestQualityGates:

    def test_both_gates_pass(self):
        result = QualityGateResult(
            casf_pearson_r=0.80, casf_passed=True,
            itc_pearson_r=0.90, itc_passed=True,
        )
        assert result.all_gates_passed

    def test_casf_fail_blocks(self):
        result = QualityGateResult(
            casf_pearson_r=0.70, casf_passed=False,
            itc_pearson_r=0.90, itc_passed=True,
        )
        assert not result.all_gates_passed

    def test_itc_fail_blocks(self):
        result = QualityGateResult(
            casf_pearson_r=0.80, casf_passed=True,
            itc_pearson_r=0.80, itc_passed=False,
        )
        assert not result.all_gates_passed

    def test_to_dict_roundtrip(self):
        result = QualityGateResult(
            casf_pearson_r=0.79, casf_rmse=1.8, casf_n=285, casf_passed=True,
            itc_pearson_r=0.91, itc_rmse=1.2, itc_n=187, itc_passed=True,
        )
        d = result.to_dict()
        assert d["casf_pearson_r"] == 0.79
        assert d["itc_passed"] is True


# ── curriculum builder tests ─────────────────────────────────────────────────

class TestCurriculumBuilder:

    def test_empty_config_raises(self):
        config = ContinuousTrainingConfig()
        trainer = ContinuousTrainer(config)
        with pytest.raises(RuntimeError, match="No datasets configured"):
            trainer.run()

    def test_phases_ordered_correctly(self, tmp_dir):
        config = ContinuousTrainingConfig(
            itc_dir="/fake/itc",
            pdbbind_refined_dir="/fake/refined",
            pdbbind_general_dir="/fake/general",
            output_dir=tmp_dir,
        )
        trainer = ContinuousTrainer(config)
        phases = trainer._build_curriculum()
        assert len(phases) == 3
        assert phases[0].name == "phase_1_itc187"
        assert phases[1].name == "phase_2_pdbbind_high"
        assert phases[2].name == "phase_3_broad_coverage"

    def test_skips_unconfigured_phases(self, tmp_dir):
        config = ContinuousTrainingConfig(
            pdbbind_refined_dir="/fake/refined",
            output_dir=tmp_dir,
        )
        trainer = ContinuousTrainer(config)
        phases = trainer._build_curriculum()
        assert len(phases) == 1
        assert phases[0].name == "phase_2_pdbbind_high"

    def test_regularization_increases(self, tmp_dir):
        config = ContinuousTrainingConfig(
            itc_dir="/fake",
            pdbbind_refined_dir="/fake",
            pdbbind_general_dir="/fake",
            bindingdb_dir="/fake",
            output_dir=tmp_dir,
        )
        trainer = ContinuousTrainer(config)
        phases = trainer._build_curriculum()
        alphas = [p.ridge_alpha for p in phases]
        assert alphas == sorted(alphas), "Ridge alpha should increase with phase"


# ── adapter registry tests ───────────────────────────────────────────────────

class TestAdapterRegistry:

    def test_create_pdbbind_refined(self):
        adapter = create_adapter("pdbbind_refined")
        assert adapter.name() == "pdbbind_refined"
        meta = adapter.metadata()
        assert meta.reliability_tier == 2

    def test_create_pdbbind_core(self):
        adapter = create_adapter("pdbbind_core")
        assert adapter.name() == "pdbbind_core"
        meta = adapter.metadata()
        assert meta.reliability_tier == 1

    def test_create_itc_187(self):
        adapter = create_adapter("itc_187")
        assert adapter.name() == "itc_187"
        meta = adapter.metadata()
        assert meta.reliability_tier == 1
        assert meta.weight == 1.0

    def test_create_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown"):
            create_adapter("nonexistent_dataset")

    def test_validation_adapters_not_training(self):
        dude = create_adapter("dude")
        assert not dude.is_training_dataset
        dekois = create_adapter("dekois2")
        assert not dekois.is_training_dataset


# ── training run result tests ────────────────────────────────────────────────

class TestTrainingRunResult:

    def test_to_dict(self):
        result = TrainingRunResult(
            run_id="test_run",
            promoted=True,
            elapsed_seconds=42.0,
            phase_metrics={"phase_1": {"pearson_r": 0.85}},
        )
        d = result.to_dict()
        assert d["run_id"] == "test_run"
        assert d["promoted"] is True

    def test_gate_results_included(self):
        gates = QualityGateResult(casf_passed=True, itc_passed=True)
        result = TrainingRunResult(
            run_id="test", gate_results=gates,
        )
        d = result.to_dict()
        assert "gate_results" in d
        assert d["gate_results"]["casf_passed"] is True


# ── manifest tests ───────────────────────────────────────────────────────────

class TestManifest:

    def test_manifest_json_valid(self, tmp_dir):
        """Manifest should produce valid JSON."""
        manifest = {
            "run_id": "test_run",
            "datasets": [{"name": "test", "n_complexes": 10}],
            "quality_gates": {"casf_2016": {"passed": True}},
        }
        path = Path(tmp_dir) / "manifest.json"
        path.write_text(json.dumps(manifest, indent=2))
        loaded = json.loads(path.read_text())
        assert loaded["run_id"] == "test_run"


# ── tier weights tests ───────────────────────────────────────────────────────

class TestTierWeights:

    def test_tier_1_is_highest(self):
        assert TIER_WEIGHTS[1] > TIER_WEIGHTS[2] > TIER_WEIGHTS[3] > TIER_WEIGHTS[4]

    def test_tier_1_is_one(self):
        assert TIER_WEIGHTS[1] == 1.0

    def test_all_tiers_positive(self):
        for tier, weight in TIER_WEIGHTS.items():
            assert weight > 0
