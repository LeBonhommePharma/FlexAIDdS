#!/usr/bin/env python3
"""Fail-closed gate on dataset IDENTITY - pytest face of scripts/check_dataset_identity.py.

The 2026-09-10 dataset identity audit found that the Astex Diverse roster differs
from the published set on 13 of 85 targets, that ligand_id was the literal string
``LIG`` on all 85 rows, that resolution_A was ``2.0`` on all 85 rows and wrong
against RCSB on 67, and that the CASF-2016 C++ list named two PDB codes that do
not exist. Nothing caught it because nothing asked the PDB.

These tests are offline: they read the committed cache at
``tests/fixtures/dataset_identity/pdb_entry_cache.json``. Refresh it explicitly
with ``python3 scripts/check_dataset_identity.py --refresh``.

NON-VACUITY is asserted directly. The DOI gate that preceded this one shipped a
rule that matched ZERO fields on the live tree while passing its own fixture -- a
green run that checked nothing. Every test below either asserts a non-zero
denominator or proves a rule fires on a known-bad fixture.

Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
SPDX-License-Identifier: Apache-2.0
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "check_dataset_identity", ROOT / "scripts" / "check_dataset_identity.py")
check_dataset_identity = importlib.util.module_from_spec(_SPEC)
sys.modules["check_dataset_identity"] = check_dataset_identity
_SPEC.loader.exec_module(check_dataset_identity)

CDI = check_dataset_identity


@pytest.fixture(scope="module")
def cache():
    doc = json.loads(CDI.CACHE_PATH.read_text(encoding="utf-8"))
    return doc["entries"]


@pytest.fixture(scope="module")
def defects():
    return json.loads(CDI.DEFECTS_PATH.read_text(encoding="utf-8"))["quarantined"]


@pytest.fixture(scope="module")
def sweep(cache, defects):
    return CDI.check(ROOT, cache, defects)


# ---------------------------------------------------------------------------
# Non-vacuity: every rule must inspect a non-zero number of live fields
# ---------------------------------------------------------------------------

def test_sweep_reaches_the_dataset_declarations(sweep):
    _, stats = sweep
    assert stats["declarations"] >= 10, stats
    assert stats["csv_rosters"] >= 4, stats
    assert stats["yaml_definitions"] >= 10, stats
    assert stats["cpp_lists"] >= 1, "the hardcoded C++ code lists are not being scanned"


@pytest.mark.parametrize("rule,key,floor", [
    ("PDB_CODE_ABSENT", "code_checks", 400),
    ("CONSTANT_COLUMN", "column_checks", 20),
    ("LIGAND_CODE_ABSENT", "ligand_checks", 50),
])
def test_rule_is_not_vacuous(sweep, rule, key, floor):
    """A rule that inspects nothing passes for free. This is the trap the DOI
    gate fell into, so each rule's live denominator is asserted, not printed."""
    _, stats = sweep
    assert stats[key] >= floor, (
        "rule %s inspected %d fields (floor %d) -- it is vacuous or the scanner "
        "stopped reaching the roster files" % (rule, stats[key], floor))


def test_cache_is_populated_and_has_both_outcomes(cache):
    assert len(cache) >= 400, len(cache)
    present = [c for c, v in cache.items() if v.get("exists")]
    absent = [c for c, v in cache.items() if not v.get("exists")]
    assert present, "no resolvable entry in the cache"
    assert absent, ("no absent code in the cache, so the PDB_CODE_ABSENT rule "
                    "cannot be exercised offline")
    with_ligand = [c for c in present if cache[c].get("nonpolymer_comp_ids")]
    assert len(with_ligand) >= 100, len(with_ligand)


def test_cache_controls_were_verified_at_refresh():
    doc = json.loads(CDI.CACHE_PATH.read_text(encoding="utf-8"))
    ctrl = doc["controls"]
    assert ctrl["positive_ok"] is True
    assert ctrl["negative_resolved"] is False
    assert cache_entry_exists(doc, CDI.POSITIVE_CONTROL) is True
    assert cache_entry_exists(doc, CDI.NEGATIVE_CONTROL) is False


def cache_entry_exists(doc, code):
    return bool(doc["entries"].get(code, {}).get("exists"))


def test_positive_control_carries_its_known_ligand(cache):
    comps = cache[CDI.POSITIVE_CONTROL]["nonpolymer_comp_ids"]
    assert CDI.POSITIVE_CONTROL_LIGAND in comps, comps


# ---------------------------------------------------------------------------
# The gate itself
# ---------------------------------------------------------------------------

def test_no_unquarantined_identity_defects(sweep):
    findings, _ = sweep
    assert not findings, "\n".join("[%s] %s %s" % (f.kind, f.path, f.detail) for f in findings)


def test_quarantine_ledger_has_no_drift(sweep):
    """Every quarantined defect must still reproduce, or the ledger is stale and
    is silently tolerating something that is no longer there."""
    findings, _ = sweep
    drift = [f for f in findings if f.kind == "DEFECT_LEDGER_DRIFT"]
    assert not drift, "\n".join(f.detail for f in drift)


def test_every_quarantine_entry_carries_a_reason(defects):
    assert defects, "the ledger is empty: either nothing is open, or it is not being read"
    missing = [k for k, v in defects.items() if not (v.get("reason") or "").strip()]
    assert not missing, missing


def test_quarantine_is_not_a_blanket(defects):
    """A quarantine must name a file, a rule and a specific field. A wildcard
    entry would disarm the gate for a whole file."""
    for key in defects:
        parts = key.split("::")
        assert len(parts) == 3, key
        assert all(p.strip() for p in parts), key
        assert "*" not in key, key


# ---------------------------------------------------------------------------
# Rule-level pass/fail proofs against synthetic fixtures
# ---------------------------------------------------------------------------

def _write(tmp_path: Path, name: str, body: str) -> Path:
    (tmp_path / name).write_text(body, encoding="utf-8")
    return tmp_path


def test_known_good_fixture_passes(tmp_path, cache):
    """A clean roster CSV: existing codes, real ligands, no constant column, and a
    per-row provenance locator (rule 9 fires on a roster that records none, so a
    fixture without one would fail for a reason this test is not about)."""
    good = [c for c, v in cache.items() if v.get("exists") and v.get("nonpolymer_comp_ids")][:3]
    body = "pdb_id,ligand_id,resolution_A,source_locator\n" + "".join(
        "%s,%s,%s,row %d\n" % (c, cache[c]["nonpolymer_comp_ids"][0], 1.4 + 0.1 * i, i)
        for i, c in enumerate(good))
    _write(tmp_path, "good_set.csv", body)
    findings, stats = CDI.check(tmp_path, cache, {}, files=["good_set.csv"])
    assert stats["code_checks"] == len(good)
    assert stats["ligand_checks"] == len(good)
    assert not findings, findings


def test_absent_pdb_code_fails(tmp_path, cache):
    absent = next(c for c, v in cache.items() if not v.get("exists"))
    ok = next(c for c, v in cache.items() if v.get("exists") and v.get("nonpolymer_comp_ids"))
    body = ("pdb_id,ligand_id,resolution_A\n%s,%s,1.5\n%s,%s,1.9\n"
            % (absent, cache[ok]["nonpolymer_comp_ids"][0], ok,
               cache[ok]["nonpolymer_comp_ids"][0]))
    _write(tmp_path, "bad_set.csv", body)
    findings, _ = CDI.check(tmp_path, cache, {}, files=["bad_set.csv"])
    assert any(f.kind == "PDB_CODE_ABSENT" and absent in f.detail for f in findings), findings


def test_constant_column_fails(tmp_path, cache):
    good = [c for c, v in cache.items() if v.get("exists") and v.get("nonpolymer_comp_ids")][:3]
    body = "pdb_id,ligand_id,resolution_A\n" + "".join(
        "%s,%s,2.0\n" % (c, cache[c]["nonpolymer_comp_ids"][0]) for c in good)
    _write(tmp_path, "const_set.csv", body)
    findings, _ = CDI.check(tmp_path, cache, {}, files=["const_set.csv"])
    assert any(f.kind == "CONSTANT_COLUMN" and "resolution_A" in f.detail for f in findings), findings


def test_absent_ligand_code_fails(tmp_path, cache):
    """The exact live defect: a ligand code that is not in its own entry."""
    good = [c for c, v in cache.items() if v.get("exists") and v.get("nonpolymer_comp_ids")][:3]
    fake = "LIG"
    assert all(fake not in cache[c]["nonpolymer_comp_ids"] for c in good)
    body = "pdb_id,ligand_id,resolution_A\n" + "".join(
        "%s,%s,%s\n" % (c, fake, 1.4 + 0.1 * i) for i, c in enumerate(good))
    _write(tmp_path, "lig_set.csv", body)
    findings, _ = CDI.check(tmp_path, cache, {}, files=["lig_set.csv"])
    hits = [f for f in findings if f.kind == "LIGAND_CODE_ABSENT"]
    assert len(hits) == len(good), findings


def test_declared_constant_is_accepted(tmp_path, cache):
    """A constant that is declared with a justification must NOT fail -- otherwise
    the rule cannot distinguish a measured constant from a placeholder."""
    good = [c for c, v in cache.items() if v.get("exists") and v.get("nonpolymer_comp_ids")][:3]
    (tmp_path / "decl_set.csv").write_text(
        "pdb_id,ligand_id,experimental_method\n" + "".join(
            "%s,%s,X-ray\n" % (c, cache[c]["nonpolymer_comp_ids"][0]) for c in good),
        encoding="utf-8")
    (tmp_path / "decl.yaml").write_text(
        "slug: t\n"
        "targets:\n" + "".join("  - %s\n" % c.lower() for c in good) +
        "roster_csv: decl_set.csv\n"
        "constant_columns:\n"
        "  - column: experimental_method\n"
        "    justification: 'measured; all entries are X-ray per RCSB'\n",
        encoding="utf-8")
    findings, _ = CDI.check(tmp_path, cache, {}, files=["decl_set.csv", "decl.yaml"])
    assert not [f for f in findings if f.kind == "CONSTANT_COLUMN"], findings
    # ... and the same file without the declaration must fail, or the exemption
    # is doing nothing.
    findings2, _ = CDI.check(tmp_path, cache, {}, files=["decl_set.csv"])
    assert [f for f in findings2 if f.kind == "CONSTANT_COLUMN"], findings2


def test_quarantine_suppresses_exactly_one_finding(tmp_path, cache):
    good = [c for c, v in cache.items() if v.get("exists") and v.get("nonpolymer_comp_ids")][:3]
    body = "pdb_id,ligand_id,resolution_A\n" + "".join(
        "%s,%s,2.0\n" % (c, cache[c]["nonpolymer_comp_ids"][0]) for c in good)
    _write(tmp_path, "q_set.csv", body)
    before, _ = CDI.check(tmp_path, cache, {}, files=["q_set.csv"])
    key = "q_set.csv::CONSTANT_COLUMN::resolution_A"
    after, _ = CDI.check(tmp_path, cache, {key: {"reason": "test"}}, files=["q_set.csv"])
    assert len(after) == len(before) - 1, (before, after)


def test_pdb_code_pattern_is_anchored():
    """An unanchored code pattern matches years and page numbers. The Hartshorn
    paper contains 2007, 1000 and 5000 as plain text."""
    for good in ("1t46", "2BSM", "3qgs"):
        assert CDI._PDB_CODE.match(good), good
    for bad in ("2007", "1000", "5000", "abcd", "1t4", "1t466", ""):
        assert not CDI._PDB_CODE.match(bad), bad


def test_selftest_entrypoint_reports_non_vacuous():
    assert CDI.selftest(ROOT) == 0

# ---------------------------------------------------------------------------
# MIRROR_DIVERGENCE - the two live copies of every dataset definition
# ---------------------------------------------------------------------------

def _mirror_tree(root, canon: dict, mirror: dict) -> None:
    import yaml
    (root / CDI.CANONICAL_DATASET_DIR).mkdir(parents=True, exist_ok=True)
    (root / CDI.MIRROR_DATASET_DIR).mkdir(parents=True, exist_ok=True)
    for name, doc in canon.items():
        (root / CDI.CANONICAL_DATASET_DIR / name).write_text(yaml.safe_dump(doc))
    for name, doc in mirror.items():
        (root / CDI.MIRROR_DATASET_DIR / name).write_text(yaml.safe_dump(doc))


def test_mirror_rule_inspects_a_nonzero_number_of_live_keys(cache):
    """The C1c defect was a rule whose regex matched ZERO fields while passing its
    own fixture. This asserts the mirror rule is not that on the real tree."""
    findings, stats = CDI.check(ROOT, cache, {})
    assert stats["mirror_pairs"] >= 10, stats
    assert stats["mirror_keys_compared"] >= 100, (
        "mirror rule compared only %s keys on the real tree" % stats.get("mirror_keys_compared"))


def test_agreeing_copies_produce_no_mirror_finding(tmp_path):
    doc = {"name": "probe", "slug": "probe", "targets": ["1T46"]}
    _mirror_tree(tmp_path, {"probe.yaml": doc}, {"probe.yaml": dict(doc)})
    findings, stats = CDI.compare_mirror(tmp_path)
    assert not findings, findings
    assert stats["mirror_keys_compared"] == 3, stats


def test_a_diverging_value_fails(tmp_path):
    doc = {"name": "probe", "slug": "probe", "targets": ["1T46"]}
    _mirror_tree(tmp_path, {"probe.yaml": doc},
                 {"probe.yaml": dict(doc, targets=["1T46", "2BSM"])})
    findings, _ = CDI.compare_mirror(tmp_path)
    assert [f for f in findings if f.kind == "MIRROR_VALUE_DIVERGENCE"], findings


def test_a_key_only_in_the_mirror_fails(tmp_path):
    """A mirror carrying data the canonical copy lacks is a second source of truth."""
    doc = {"name": "probe", "slug": "probe"}
    _mirror_tree(tmp_path, {"probe.yaml": doc}, {"probe.yaml": dict(doc, targets=["1T46"])})
    findings, _ = CDI.compare_mirror(tmp_path)
    assert [f for f in findings if f.kind == "MIRROR_ONLY_KEY"], findings


def test_an_annotation_key_only_in_the_canonical_copy_is_allowed(tmp_path):
    """The asymmetry is deliberate and applies to the CANONICAL side only."""
    doc = {"name": "probe", "slug": "probe"}
    _mirror_tree(tmp_path,
                 {"probe.yaml": dict(doc, claude_reevaluated_baselines={"a": 1})},
                 {"probe.yaml": dict(doc)})
    findings, _ = CDI.compare_mirror(tmp_path)
    assert not findings, findings


def test_a_non_annotation_key_only_in_the_canonical_copy_fails(tmp_path):
    """The Python runner resolves a dataclass default for it, silently. Leaving this
    direction unchecked let naloxone_ss drop nine declared keys from the mirror,
    three of which change the resolved config."""
    doc = {"name": "probe", "slug": "probe"}
    _mirror_tree(tmp_path, {"probe.yaml": dict(doc, baseline_tolerance=0.1)},
                 {"probe.yaml": dict(doc)})
    findings, _ = CDI.compare_mirror(tmp_path)
    assert [f for f in findings if f.kind == "CANONICAL_ONLY_KEY"], findings


def test_an_annotation_key_only_in_the_MIRROR_is_not_excused(tmp_path):
    """The allowance is one-directional: the prefixes excuse nothing in the mirror."""
    doc = {"name": "probe", "slug": "probe"}
    _mirror_tree(tmp_path, {"probe.yaml": dict(doc)},
                 {"probe.yaml": dict(doc, claude_reevaluated_by="probe")})
    findings, _ = CDI.compare_mirror(tmp_path)
    assert [f for f in findings if f.kind == "MIRROR_ONLY_KEY"], findings


def test_the_documented_allowance_matches_the_code():
    """CANONICAL.md described this allowance in the opposite direction to the code
    once; a reader trusting the doc believed the gate covered a case it skipped."""
    doc = (ROOT / "benchmarks" / "datasets" / "CANONICAL.md").read_text(encoding="utf-8")
    assert "CANONICAL_ONLY_KEY" in doc, "CANONICAL.md does not document the rule"
    for pfx in CDI._CANONICAL_ONLY_PREFIXES:
        assert pfx.rstrip("_") in doc, "CANONICAL.md does not name the prefix %r" % pfx
    canon_idx = doc.index("CANONICAL_ONLY_KEY")
    mirror_idx = doc.index("MIRROR_ONLY_KEY")
    # The exemption sentence must sit under CANONICAL_ONLY_KEY, not MIRROR_ONLY_KEY.
    exempt_idx = doc.index("Exempt **only** for the declared")
    assert canon_idx < exempt_idx, "the exemption is not documented under CANONICAL_ONLY_KEY"
    assert "unconditionally" in doc[mirror_idx:canon_idx], (
        "MIRROR_ONLY_KEY is not documented as firing unconditionally")


def test_an_unpaired_definition_fails(tmp_path):
    _mirror_tree(tmp_path, {"probe.yaml": {"name": "probe"}}, {})
    findings, stats = CDI.compare_mirror(tmp_path)
    assert [f for f in findings if f.kind == "MIRROR_UNPAIRED"], findings
    assert stats["mirror_unpaired"] == 1, stats


def test_every_live_mirror_divergence_is_quarantined_with_a_reason(cache):
    """Nothing about the two live trees may be tolerated silently. Each open
    divergence must carry a dated reason in the ledger."""
    import json
    raw = CDI.compare_mirror(ROOT)[0]
    assert raw, "the mirror rule found nothing on the real tree -- it would be vacuous"
    ledger = json.loads(CDI.DEFECTS_PATH.read_text())["quarantined"]
    unexplained = []
    for f in raw:
        key_part = f.detail.split("'")[1] if "'" in f.detail else "file"
        key = "%s::%s::%s" % (f.path, f.kind, key_part)
        if key not in ledger or not ledger[key].get("reason"):
            unexplained.append(key)
    assert not unexplained, (
        "live mirror divergences with no quarantine reason: %s" % unexplained)


# ---------------------------------------------------------------------------
# S5 - provenance, not just syntax (rules 7-9)
# ---------------------------------------------------------------------------

def _registry():
    return CDI.load_doi_registry()


def test_doi_registry_loads_and_has_both_polarities():
    """The resolution arms need a registered AND an unregistered DOI to be real."""
    reg = _registry()
    assert reg, "no DOI registry loaded from %s" % CDI.DOI_CACHE_PATH
    assert any(reg.values()), "no registered DOI in the cache"
    assert any(not v for v in reg.values()), "no unregistered DOI in the cache"


def test_provenance_rule_inspects_a_nonzero_number_of_live_fields():
    _, stats = CDI.check_provenance(ROOT, _registry())
    assert stats["definitions_with_a_roster"] >= 20, stats
    assert stats["primary_sources_inspected"] >= 20, stats
    assert stats["declared_dois_checked"] >= 1, (
        "no declared primary_source DOI was checked on the real tree: the "
        "PRIMARY_SOURCE_UNRESOLVED rule is vacuous here")
    assert stats["row_provenance_checks"] >= 5, stats


def _prov(tmp_path, doc=None, csv_body=None):
    import yaml
    files = []
    if doc is not None:
        (tmp_path / "probe.yaml").write_text(yaml.safe_dump(doc))
        files.append("probe.yaml")
    if csv_body is not None:
        (tmp_path / "probe_set.csv").write_text(csv_body)
        files.append("probe_set.csv")
    return CDI.check_provenance(tmp_path, _registry(), files)[0]


def test_a_roster_with_no_primary_source_fails(tmp_path):
    f = _prov(tmp_path, {"name": "p", "slug": "p", "targets": ["1T46"]})
    assert [x for x in f if x.kind == "PRIMARY_SOURCE_UNDECLARED"], f


def test_prose_naming_a_paper_is_not_a_primary_source(tmp_path):
    """published_source held 'Gaudreault & Najmanovich 2015 JCIM Table 2' on six
    definitions. It names a paper; nothing can check it, and nothing did."""
    f = _prov(tmp_path, {"name": "p", "slug": "p", "targets": ["1T46"],
                         "published_source": "Hartshorn et al. 2007 JCIM Table 2"})
    assert [x for x in f if x.kind == "PRIMARY_SOURCE_UNDECLARED"], f


def test_an_explicit_unverified_admission_with_a_reason_passes(tmp_path):
    """A dataset with no verifiable source is not a failure -- but it must say so."""
    f = _prov(tmp_path, {"name": "p", "slug": "p", "targets": ["1T46"],
                         "primary_source": {"status": "UNVERIFIED", "reason": "none known"}})
    assert not f, f


def test_an_unverified_admission_without_a_reason_fails(tmp_path):
    f = _prov(tmp_path, {"name": "p", "slug": "p", "targets": ["1T46"],
                         "primary_source": {"status": "UNVERIFIED"}})
    assert [x for x in f if x.kind == "PRIMARY_SOURCE_UNDECLARED"], f


def test_a_registered_doi_passes_and_an_unregistered_one_fails(tmp_path):
    reg = _registry()
    good = next(d for d, v in reg.items() if v)
    bad = next(d for d, v in reg.items() if not v)
    base = {"name": "p", "slug": "p", "targets": ["1T46"]}
    assert not _prov(tmp_path, dict(base, primary_source={"doi": good}))
    f = _prov(tmp_path, dict(base, primary_source={"doi": bad}))
    assert [x for x in f if x.kind == "PRIMARY_SOURCE_UNRESOLVED"], f


def test_a_doi_absent_from_the_crossref_cache_fails(tmp_path):
    """Never doi.org -L: it returns HTTP 000 for valid DOIs and 404 for fakes."""
    absent = "10.9999" + "/" + "absent.control"   # assembled: see the gate's selftest
    ps = {"doi": absent}
    f = _prov(tmp_path, {"name": "p", "slug": "p", "targets": ["1T46"],
                         "primary_source": ps})
    assert [x for x in f if x.kind == "PRIMARY_SOURCE_UNRESOLVED"], f


def test_a_roster_csv_with_no_per_row_provenance_fails(tmp_path):
    f = _prov(tmp_path, csv_body="pdb_id,rmsd_threshold_A\n1T46,2.0\n2BSM,2.5\n")
    assert [x for x in f if x.kind == "ROW_PROVENANCE_UNDECLARED"], f


def test_a_roster_csv_with_per_row_provenance_passes(tmp_path):
    f = _prov(tmp_path, csv_body="pdb_id,source_locator\n1T46,table 2 p734\n2BSM,table 2 p735\n")
    assert not f, f


def test_a_provenance_column_empty_on_every_row_is_not_provenance(tmp_path):
    f = _prov(tmp_path, csv_body="pdb_id,source_locator\n1T46,\n2BSM,\n")
    assert [x for x in f if x.kind == "ROW_PROVENANCE_UNDECLARED"], f


def test_every_live_provenance_finding_is_quarantined_with_a_reason():
    import json
    raw = CDI.check_provenance(ROOT, _registry())[0]
    ledger = json.loads(CDI.DEFECTS_PATH.read_text())["quarantined"]
    unexplained = [
        "%s::%s::%s" % (f.path, f.kind, f.kind) for f in raw
        if not ledger.get("%s::%s::%s" % (f.path, f.kind, f.kind), {}).get("reason")]
    assert not unexplained, (
        "live provenance findings with no quarantine reason: %s" % unexplained)
