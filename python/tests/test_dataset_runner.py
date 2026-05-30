"""Unit tests for new DatasetRunner / EntryTaskManager features (per-entry, cost-aware, manager)."""

import json
import tempfile
from pathlib import Path

import pytest

from flexaidds.dataset_runner.runner import EntryTaskManager


def test_entry_manager_basic_local():
    work = [("t1", "holo"), ("t2", "holo")]
    mgr = EntryTaskManager(work, n_workers=2)

    def fake(item):
        return (*item, [], 0.1, "")

    res = mgr.run(fake)
    assert len(res) == 2


def test_entry_manager_cost_hints_sorting():
    work = [("expensive", "holo"), ("cheap", "holo")]
    hints = {"cheap_holo": 1.0, "expensive_holo": 100.0}
    mgr = EntryTaskManager(work, cost_hints=hints)
    # Cheapest should be first after __init__ sorting
    assert mgr.work_items[0] == ("cheap", "holo")


def test_entry_manager_load_cost_hints_from_manifest():
    with tempfile.TemporaryDirectory() as tmp:
        man = Path(tmp) / "_entry_manifest.json"
        man.write_text(json.dumps({
            "timings": {
                "per_entry_cost_cpu_seconds": {"1a30_holo": 4.2}
            }
        }))
        hints = EntryTaskManager.load_cost_hints_from_manifest(man)
        assert hints["1a30_holo"] == 4.2


def test_entry_manager_hybrid_pool():
    mgr = EntryTaskManager([("t", "holo")], n_workers=2)

    def fake(item):
        return (*item, [], 0.05, "")

    res = mgr.run(fake)
    assert len(res) == 1
