"""Source contract for startup-vs-docking failure diagnosis.

WHAT THIS PREVENTS
------------------
Until 2026-09-20 every runtime failure in LIB/benchmark_datasets.cpp left the
process through one line::

    ERROR: Incomplete docking; inspect per-target runtime fields and child logs

`grep -rn "Incomplete docking" LIB/` returned exactly one emission site, and
benchmark_runtime_exit_code() (LIB/DatasetRunner.h) folds six distinct
conditions into one nonzero return. So an engine that aborted at startup having
read no atoms, an engine that ran fine and found no pose, a missing input and a
timeout all produced the same sentence and the same exit code 2.

That is what made the ghost-TMPDIR incident expensive. Nine cells of an
85-target Astex campaign died because $TMPDIR named a directory an agent
harness had swept; the parent reported it as incomplete docking, and it read
like nine hard targets for six hours.

WHY A SOURCE TEST AND NOT ONLY A UNIT TEST
------------------------------------------
tests/test_startup_diagnosis.cpp proves the classifier classifies. It cannot
prove three properties that live between files and therefore need a tree-wide
assertion:

1. the classifier is actually WIRED into both exit paths of main() — a correct
   classifier nobody calls is what the old code effectively had;
2. the exit-code table in the banner comment still agrees with
   exit_code_for() — documentation drift here hands a shell driver a number
   that means something else;
3. the substrings used to detect the engine's terminal handler still match the
   format strings LIB/top.cpp actually prints. That coupling crosses a file
   boundary, is invisible to the compiler, and silently degrades
   engine_startup_abort back into the generic class if top.cpp is reworded.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "LIB" / "benchmark_datasets.cpp"
TOP = REPO / "LIB" / "top.cpp"

# Frozen v1 column prefix of <dataset>_failure_diagnosis.csv. APPEND-ONLY: a
# new column goes after child_stderr_log/schema_version, so this string must
# remain a PREFIX of the header forever. Inserting a column in the middle
# silently shifts every later named field in every downstream reader.
CSV_V1_PREFIX = (
    "pdb_id,failure_reason,failure_exit_code,docking_exit_code,num_poses,"
    "docking_completed,stuck,wall_time_s,child_stderr_log"
)

# Classes this site must NOT invent. Each is unreachable or unknowable here;
# the rationale is in the startup_diag banner in LIB/benchmark_datasets.cpp.
FORBIDDEN_CLASSES = ("poses_written_but_all_filtered", "input_missing")


def _runner_text() -> str:
    return RUNNER.read_text(encoding="utf-8", errors="replace")


def _vocabulary(text: str) -> set[str]:
    """Terms returned by startup_diag::to_string() for an enumerated class.

    Only `case Class::X: return "term";` counts. The trailing
    `return "unclassified";` is the fall-through for a value outside the enum,
    paired with exit_code_for()'s fall-through to the legacy code 2; it is not
    a member of the vocabulary and must not be required to have a code of its
    own. Matching bare returns here conflated the two.
    """
    block = re.search(r"inline const char\* to_string\(Class c\)\s*\{(.*?)\n\}",
                      text, re.S)
    assert block, "to_string(Class) not found — has the classifier been removed?"
    return set(re.findall(r'case Class::\w+:\s*return\s+"([a-z_]+)";', block.group(1)))


def _exit_codes(text: str) -> dict[str, int]:
    """case Class::X: return N;  ->  {to_string(X): N} via enumerator order."""
    to_str = re.search(r"inline const char\* to_string\(Class c\)\s*\{(.*?)\n\}",
                       text, re.S)
    ec = re.search(r"inline int exit_code_for\(Class c\)\s*\{(.*?)\n\}", text, re.S)
    assert to_str and ec, "to_string / exit_code_for not both present"
    names = dict(re.findall(r"case Class::(\w+):\s*return\s+\"([a-z_]+)\";",
                            to_str.group(1)))
    codes = dict(re.findall(r"case Class::(\w+):\s*return\s+(\d+);", ec.group(1)))
    return {names[k]: int(v) for k, v in codes.items() if k in names}


def _documented_codes(text: str, vocab: set[str]) -> dict[str, int]:
    """The banner table: `//    20   engine_startup_abort  ...`."""
    out: dict[str, int] = {}
    for num, term in re.findall(r"^//\s+(\d+)\s+([a-z_]+)", text, re.M):
        if term in vocab:
            out[term] = int(num)
    return out


# ── 1. the classifier is wired into every exit path ──────────────────────

def test_the_single_emission_site_now_names_a_class():
    """The old one-size message must not survive as the only output.

    The assertion is deliberately per-LINE. An earlier version of this guard
    read `assert "failure_reason" in ln or "failure_reason" in text`, which is
    unfalsifiable: the sidecar CSV header literal alone contains the token, so
    the second disjunct is always true and a reintroduced bare emission would
    have passed. Caught in review; the escape hatch is gone.
    """
    text = _runner_text()
    emissions = [ln.strip() for ln in text.splitlines()
                 if "ERROR: Incomplete docking" in ln
                 and not ln.strip().startswith("//")]
    assert emissions, (
        "the literal was dropped entirely; an existing log grep now misses")
    bare = [ln for ln in emissions if "failure_reason=" not in ln]
    assert not bare, (
        "emission prints the bare historical sentence with no class on the same "
        f"line: {bare}")


def test_both_exit_paths_diagnose():
    """`--benchmark all` and the single-benchmark path must both classify.

    A classifier wired into only one of them would leave `--benchmark all`
    exactly as opaque as before.
    """
    text = _runner_text()
    assert text.count("diagnose_and_publish(") >= 3, (
        "expected the definition plus a call in each of the two branches of main")
    assert "print_human_summary(worst" in text, "run summary never printed"


def test_legacy_predicate_still_decides_pass_fail():
    """Renaming failures must not reclassify any run as passing.

    The granular code is applied only inside a `runtime_exit_code != 0` /
    `legacy_rc != 0` guard, so a class can never zero a nonzero outcome.
    """
    text = _runner_text()
    assert "benchmark_runtime_exit_code" in text
    # Count first. A for-loop over a regex that matches nothing passes
    # silently, which is the same vacuity defect as the disjunct removed from
    # test_the_single_emission_site_now_names_a_class above: a rename of the
    # assignment would retire this guard without failing it.
    # Search CODE, not prose: a rationale comment above an assignment can be
    # longer than the window and push the enclosing guard out of view, which
    # would fail an otherwise-correct file. Drop whole-line comments only —
    # stripping to end-of-line everywhere could cut inside a string literal.
    code = "\n".join(ln for ln in text.splitlines()
                     if not ln.strip().startswith("//"))
    sites = list(re.finditer(r"worst\.process_exit_code", code))
    assert len(sites) >= 2, (
        f"expected the granular assignment in both branches of main, found "
        f"{len(sites)} — guard would be vacuous")
    for m in sites:
        window = code[max(0, m.start() - 400):m.end() + 60]
        assert "!= 0" in window, (
            "granular exit code assigned outside a failure guard — this could "
            f"turn a failing run into a passing one; context: ...{window[-200:]}")


# ── 2. vocabulary is closed, total, and honest ───────────────────────────

def test_every_class_has_a_code_and_a_severity():
    text = _runner_text()
    vocab = _vocabulary(text)
    assert len(vocab) >= 10, f"vocabulary implausibly small: {sorted(vocab)}"
    codes = _exit_codes(text)
    missing = vocab - set(codes)
    assert not missing, f"classes with no exit code: {sorted(missing)}"
    sev = re.search(r"inline int severity_rank\(Class c\)\s*\{(.*?)\n\}", text, re.S)
    assert sev, "severity_rank not found"
    ranked = set(re.findall(r"case Class::(\w+):", sev.group(1)))
    enum_names = set(re.findall(r"case Class::(\w+):",
                                re.search(r"to_string\(Class c\)\s*\{(.*?)\n\}",
                                          text, re.S).group(1)))
    assert enum_names <= ranked, (
        f"classes with no severity rank: {sorted(enum_names - ranked)}")


def test_exit_codes_are_distinct_and_do_not_reuse_the_legacy_four():
    codes = _exit_codes(_runner_text())
    dupes = {c for c in codes.values() if list(codes.values()).count(c) > 1}
    assert not dupes, f"exit codes shared between classes: {sorted(dupes)}"
    for term, code in codes.items():
        if term == "none":
            assert code == 0, "none must remain exit 0"
        else:
            assert code > 3, (
                f"{term} uses {code}, colliding with a pre-existing code "
                "(1 usage, 2 generic/fleet-empty, 3 fleet publish)")


def test_no_invented_classes():
    """Categories the code cannot actually distinguish must stay absent."""
    vocab = _vocabulary(_runner_text())
    present = [c for c in FORBIDDEN_CLASSES if c in vocab]
    assert not present, (
        f"{present} is in the vocabulary but this site cannot produce or read "
        "it; see the startup_diag banner for why")


# ── 3. the documented table matches the implementation ───────────────────

def test_banner_exit_code_table_matches_the_code():
    text = _runner_text()
    vocab = _vocabulary(text)
    documented = _documented_codes(text, vocab)
    implemented = _exit_codes(text)
    assert documented, "no exit-code table found in the banner comment"
    drift = {t: (documented[t], implemented.get(t))
             for t in documented if documented[t] != implemented.get(t)}
    assert not drift, f"banner/code drift (documented, implemented): {drift}"
    undocumented = {t for t, c in implemented.items() if c > 3 and t not in documented}
    assert not undocumented, f"classes with no banner entry: {sorted(undocumented)}"


# ── 4. the cross-file coupling to the engine's own handlers ──────────────

def test_terminal_handler_markers_match_top_cpp():
    """The detected substrings must be ones LIB/top.cpp actually prints.

    This is the only place the classifier reads prose, and it reads prose the
    engine emits from a format string in this repository. If top.cpp is
    reworded and this is not, engine_startup_abort silently degrades to the
    generic class and the incident becomes invisible again.
    """
    text = _runner_text()
    fn = re.search(r"line_has_engine_terminal_handler\(const std::string& line\)"
                   r"\s*\{(.*?)\n\}", text, re.S)
    assert fn, "marker predicate not found"
    markers = re.findall(r'line\.find\("([^"]+)"\)', fn.group(1))
    assert markers, "predicate matches nothing"
    top = TOP.read_text(encoding="utf-8", errors="replace")
    for m in markers:
        assert m in top, (
            f"marker {m!r} no longer appears in LIB/top.cpp; the classifier is "
            "looking for a string the engine does not print")


def test_the_marker_guard_can_actually_fail():
    """Non-vacuity: a marker top.cpp does not print must be rejected.

    A guard that has never been observed failing is not a guard. This runs the
    same containment check the test above runs, with a string chosen to be
    absent, and asserts it fails.
    """
    top = TOP.read_text(encoding="utf-8", errors="replace")
    bogus = "Catastrophic engine meltdown: "
    assert bogus not in top, "fixture string unexpectedly present in top.cpp"
    # Same predicate, opposite expectation — proves the check discriminates.
    assert not all(m in top for m in ["Fatal error: ", bogus])


# ── 5. the sidecar schema is append-only ─────────────────────────────────

def test_csv_header_keeps_the_v1_prefix():
    text = _runner_text()
    hdr = re.search(r"inline const char\* diagnosis_csv_header\(\)\s*\{(.*?)\n\}",
                    text, re.S)
    assert hdr, "diagnosis_csv_header not found"
    literal = "".join(re.findall(r'"([^"]*)"', hdr.group(1)))
    assert literal.startswith(CSV_V1_PREFIX), (
        "column inserted or renamed rather than appended.\n"
        f"  expected prefix: {CSV_V1_PREFIX}\n"
        f"  actual header:   {literal}")
    assert literal.split(",")[0] == "pdb_id", "pdb_id must remain the join key"


def test_header_and_row_emit_the_same_number_of_fields():
    text = _runner_text()
    hdr = re.search(r"diagnosis_csv_header\(\)\s*\{(.*?)\n\}", text, re.S)
    row = re.search(r"inline std::string diagnosis_csv_row\(.*?\)\s*\{(.*?)\n\}",
                    text, re.S)
    assert hdr and row
    n_cols = len("".join(re.findall(r'"([^"]*)"', hdr.group(1))).split(","))
    # Each emitted field is followed by << ',' except the last one.
    n_seps = len(re.findall(r"<<\s*','", row.group(1)))
    assert n_seps == n_cols - 1, (
        f"header declares {n_cols} columns, row emits {n_seps + 1}")


def test_existing_result_csv_writers_were_not_touched():
    """The sidecar is a NEW file. No pre-existing CSV may be rewritten here.

    <dataset>_results.csv and per-target result.csv are written by
    DatasetRunner::write_report and the per-complex writer. A post-pass in this
    file that reopened either of them would be a read-modify-write of a file a
    live campaign is appending to.
    """
    text = _runner_text()
    for forbidden in ("_results.csv", "result.csv"):
        offenders = [ln.strip() for ln in text.splitlines()
                     if forbidden in ln and not ln.strip().startswith("//")]
        assert not offenders, (
            f"{forbidden} referenced in non-comment code: {offenders[:3]}")
    assert "_failure_diagnosis" in text, "sidecar basename missing"


# ── 6. the rendered JSON is really JSON ──────────────────────────────────

@pytest.mark.skipif(shutil.which("clang++") is None and shutil.which("g++") is None,
                    reason="no C++ compiler available")
def test_rendered_json_parses():
    """Compile the classifier alone and json.loads its output.

    The unit test asserts escaping with substring checks; only a real parser
    proves the document is well formed.
    """
    cxx = shutil.which("clang++") or shutil.which("g++")
    driver = r"""
#define FLEXAIDS_BENCHMARK_DATASETS_CLASSIFIER_ONLY 1
#include "LIB/benchmark_datasets.cpp"
#include <iostream>
int main() {
    startup_diag::RunDiagnosis run;
    run.dataset_name = "Astex Diverse";
    run.dominant = startup_diag::Class::EngineStartupAbort;
    run.process_exit_code = startup_diag::exit_code_for(run.dominant);
    run.n_targets = 2;
    run.startup_abort_targets = {"1G9V"};
    run.first_abort_pdb_id = "1G9V";
    run.first_abort_log_path = "/out/1G9V/stderr.log";
    run.first_abort_stderr_tail = {"Fatal error: in temp_directory_path: \"/tmp/x\"\ttab"};
    startup_diag::TargetDiagnosis a, b;
    a.facts.pdb_id = "1G9V"; a.cls = startup_diag::Class::EngineStartupAbort;
    b.facts.pdb_id = "1HWI"; b.cls = startup_diag::Class::None;
    run.targets = {a, b};
    std::cout << startup_diag::render_diagnosis_json(run);
    return 0;
}
"""
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "json_probe.cpp"
        src.write_text(driver, encoding="utf-8")
        exe = Path(td) / "json_probe"
        build = subprocess.run(
            [cxx, "-std=c++20", "-O0", f"-I{REPO}", str(src), "-o", str(exe)],
            capture_output=True, text=True, cwd=REPO)
        assert build.returncode == 0, build.stderr[-2000:]
        out = subprocess.run([str(exe)], capture_output=True, text=True)
        assert out.returncode == 0, out.stderr[-2000:]
        doc = json.loads(out.stdout)

    assert doc["run_failure_reason"] == "engine_startup_abort"
    assert doc["process_exit_code"] > 3
    assert doc["counts"] == {"engine_startup_abort": 1, "none": 1}
    assert doc["startup_abort_targets"] == ["1G9V"]
    assert "temp_directory_path" in doc["first_startup_abort"]["stderr_tail"][0]
    assert doc["schema_version"] >= 1
