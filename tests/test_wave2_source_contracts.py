"""Wave 2 source contracts: live register_result, quoted PoseBusters cache path."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASET_RUNNER = ROOT / "LIB" / "DatasetRunner.cpp"


def _fetch_posebusters_body(text: str) -> str:
    start = text.find("DatasetRunner::fetch_posebusters")
    end = text.find("DatasetRunner::fetch_bindingdb_itc")
    assert start != -1 and end != -1 and end > start
    return text[start:end]


def test_register_result_call_is_not_commented_out():
    text = DATASET_RUNNER.read_text(encoding="utf-8")
    live = [
        line.strip()
        for line in text.splitlines()
        if "register_result(sess)" in line and not line.strip().startswith("//")
    ]
    assert live, "register_result(sess) must actually run in DatasetRunner"
    assert any("ts_it2" in line for line in live)


def test_fetch_posebusters_quotes_cache_path_for_execve():
    body = _fetch_posebusters_body(DATASET_RUNNER.read_text(encoding="utf-8"))
    assert "shell_quote(repo_dir)" in body
    assert "+ repo_dir +" not in body
