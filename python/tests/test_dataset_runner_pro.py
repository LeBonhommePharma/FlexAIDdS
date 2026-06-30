"""Tests for the standalone dataset_runner_pro.py helper script."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "dataset_runner_pro.py"

if not SCRIPT_PATH.exists():
    pytest.skip("dataset_runner_pro.py missing (optional)", allow_module_level=True)


def load_module():
    if not SCRIPT_PATH.exists():
        pytest.skip(f"dataset_runner_pro.py not present at {SCRIPT_PATH} (pro script optional)")
    spec = importlib.util.spec_from_file_location("dataset_runner_pro", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class FakeFuture:
    def __init__(self, result=None, error=None):
        self._result = result
        self._error = error

    def result(self):
        if self._error is not None:
            raise self._error
        return self._result


class FakeExecutor:
    def __init__(self, max_workers=None):
        self.max_workers = max_workers
        self.submissions = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def submit(self, func, *args):
        self.submissions.append((func, args))
        try:
            return FakeFuture(result=func(*args))
        except Exception as exc:  # pragma: no cover - exercised via .result()
            return FakeFuture(error=exc)


def test_discover_benchmarks(tmp_path, monkeypatch):
    module = load_module()
    fake_root = tmp_path / "repo"
    (fake_root / "tests" / "benchmarks" / "casf2016").mkdir(parents=True)
    (fake_root / "tests" / "benchmarks" / "crossdock").mkdir(parents=True)
    (fake_root / "tests" / "benchmarks" / "casf2016" / "run_casf2016.py").write_text("")
    (fake_root / "tests" / "benchmarks" / "crossdock" / "run_crossdock.py").write_text("")
    monkeypatch.setattr(module, "__file__", str(fake_root / "dataset_runner_pro.py"))
    assert module.discover_benchmarks() == ["casf2016", "crossdock"]


def test_list_benchmarks_does_not_require_paths(monkeypatch, capsys):
    module = load_module()
    monkeypatch.setattr(module, "discover_benchmarks", lambda: ["casf2016", "litpcba"])
    module.main(["--list-benchmarks"])
    out = capsys.readouterr().out
    assert "casf2016" in out
    assert "litpcba" in out


def test_run_script_raises_on_called_process_error(monkeypatch):
    module = load_module()

    def fake_run(*args, **kwargs):
        raise module.subprocess.CalledProcessError(returncode=2, cmd=["python", "x.py"])

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    with pytest.raises(module.subprocess.CalledProcessError):
        module.run_script(["python", "x.py"], timeout=10)


def test_run_script_raises_on_timeout(monkeypatch):
    module = load_module()

    def fake_run(*args, **kwargs):
        raise module.subprocess.TimeoutExpired(cmd=["python", "x.py"], timeout=10)

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    with pytest.raises(module.subprocess.TimeoutExpired):
        module.run_script(["python", "x.py"], timeout=10)


def test_main_writes_all_summary_formats(tmp_path, monkeypatch):
    module = load_module()

    results_root = tmp_path / "results"
    data_root = tmp_path / "data"
    output_root = tmp_path / "out"
    (results_root / "casf2016").mkdir(parents=True)
    (data_root / "CASF-2016").mkdir(parents=True)

    def fake_run(results_dir, data_dir, output_dir, timeout):
        return {
            "benchmark": "casf2016",
            "status": "ok",
            "returncode": 0,
            "duration_s": 0.01,
            "report": str(output_dir / "casf2016_report.json"),
            "error": "",
        }

    captured = {}

    class CapturingExecutor(FakeExecutor):
        def __init__(self, max_workers=None):
            captured["max_workers"] = max_workers
            super().__init__(max_workers=max_workers)

    monkeypatch.setitem(module.BENCHMARKS, "casf2016", (fake_run, 0))
    monkeypatch.setattr(module, "ProcessPoolExecutor", CapturingExecutor)
    monkeypatch.setattr(module, "as_completed", lambda futures: list(futures))
    monkeypatch.setitem(sys.modules, "yaml", SimpleNamespace(safe_dump=lambda data, fp, sort_keys=False: fp.write("ok: true\n")))

    module.main([
        "--results", str(results_root),
        "--data-root", str(data_root),
        "--benchmarks", "casf2016",
        "--output", str(output_root),
        "--summary-format", "all",
        "--max-workers", "1",
    ])

    assert (output_root / "summary.json").is_file()
    assert (output_root / "summary.yaml").is_file()
    assert (output_root / "summary.html").is_file()
    assert captured["max_workers"] == 1


def test_main_logs_missing_results(caplog):
    module = load_module()
    with caplog.at_level("ERROR"):
        module.main(["--data-root", "/tmp/does-not-matter"])
    assert "--results is required" in caplog.text
