#!/usr/bin/env python3
"""check_dataset_identity.py - Fail-closed gate on dataset IDENTITY claims.

Motivation
----------
The 2026-09-10 dataset identity audit found that the repository's Astex Diverse
roster differs from the published set on 13 of 85 targets, that its ligand_id
column held the literal string ``LIG`` on all 85 rows, that its resolution_A
column held ``2.0`` on all 85 rows and disagreed with RCSB on 67, and that the
CASF-2016 C++ list contained two PDB codes that do not exist. None of it was
caught, because nothing in this repository ever asked the PDB whether a declared
code or a declared ligand was real.

This script is that gate. It is the identity counterpart of
``scripts/check_dois.py`` (the citation gate added by f005e77f) and deliberately
copies its shape: offline by default against a committed fixture, an explicit
``--refresh`` that will not write a cache unless controls behave, a ledger for
known defects that are not yet fixable, and a ``--selftest`` that proves the
gate is non-vacuous.

What it fails on
----------------
1. PDB_CODE_ABSENT   - a declared PDB code that the PDB does not serve.
                       Covers withdrawn depositions (3QGS, 3RP3), which return
                       404 and are absent from the obsoleted-entry list, so no
                       successor code exists.
2. CONSTANT_COLUMN   - a per-target column with ONE distinct value across all
                       rows. A constant column is a default, not a measurement.
                       Declare a legitimate constant in the dataset definition's
                       ``constant_columns:`` block with a justification and it is
                       accepted; leave it undeclared and this fires.
3. LIGAND_CODE_ABSENT- a declared ligand code that is not a non-polymer
                       component of its own entry. This is what ``LIG`` on all
                       85 Astex rows was, and what the astex_nonnative
                       ``ligand_id`` column was on 0/74 matching rows.
4. UNCACHED          - a declared code with no cached resolution. The cache is
                       stale: re-run with ``--refresh``.
5. DEFECT_LEDGER_DRIFT - a quarantined known defect that no longer reproduces,
                       or a quarantine entry for a file that no longer exists.
                       Stops the ledger from silently outliving the defect.

Known defects that are NOT yet fixed are quarantined in
``tests/fixtures/dataset_identity/known_defects.json`` with a reason, so the
gate fails on anything NEW while the outstanding defect stays enumerated and
visible instead of being silently tolerated.

Offline by default
------------------
``--refresh`` is the only mode that touches the network, and it refuses to write
a cache unless BOTH controls behave: a known-real entry must resolve WITH a
non-polymer component, and a known-withdrawn code must not resolve. A 404 from a
rate-limited or blocked host is not evidence that a code is fake, and this is how
that failure mode is excluded.

Usage
-----
    python3 scripts/check_dataset_identity.py            # offline gate (CI)
    python3 scripts/check_dataset_identity.py --refresh  # re-resolve, rewrite cache
    python3 scripts/check_dataset_identity.py --selftest # non-vacuity proof
    python3 scripts/check_dataset_identity.py --report   # print the full sweep

Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
SPDX-License-Identifier: Apache-2.0
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = ROOT / "tests" / "fixtures" / "dataset_identity"
CACHE_PATH = FIXTURE_DIR / "pdb_entry_cache.json"
DEFECTS_PATH = FIXTURE_DIR / "known_defects.json"

# Controls. The positive control is a released entry with a well-known ligand;
# the negative control is a WITHDRAWN deposition that has never been served.
POSITIVE_CONTROL = "1T46"          # must resolve, and must contain STI
POSITIVE_CONTROL_LIGAND = "STI"
NEGATIVE_CONTROL = "3QGS"          # WDRN 2011-01-24, never released

# A PDB id: digit followed by three alphanumerics. Anchored, because an
# unanchored version matches page numbers, years and four-digit potencies.
_PDB_CODE = re.compile(r"^[0-9][A-Za-z][A-Za-z0-9]{2}$")

# Columns that name a PDB entry, and columns that name a ligand component.
_PDB_COLUMNS = ("pdb_id", "target_pdb", "ligand_pdb", "receptor_id", "entry_id")
_LIGAND_COLUMNS = ("ligand_id", "ligand_code", "het_code", "ccd_code")

# Columns that are bookkeeping rather than per-target measurements, so a single
# distinct value in them is not a placeholder claim.
_NOT_A_MEASUREMENT = ("index", "notes", "basis", "status", "path", "sha256", "dir")

# The gate's own source and fixtures: a registry OF codes, not a declaration of
# a dataset. Scanning them makes the gate fail on itself. Same rationale as the
# citation gate's self-exclusion, and the same accepted limit: a genuine defect
# introduced into these files is not caught by this gate.
_SELF = {
    "scripts/check_dataset_identity.py",
    "tests/test_dataset_identity.py",
}

_SKIP_DIR_PARTS = {
    ".git", "build", "build-ci-fix", "build-audit", "node_modules", "__pycache__",
    ".venv", "venv", ".mypy_cache", ".pytest_cache", ".claude", ".cursor", ".gemini",
    ".trae", ".windsurf", ".qoder", ".roo", "docs",
}


class Finding(NamedTuple):
    kind: str
    path: str
    detail: str


class Declaration(NamedTuple):
    """One dataset declaration: where it is, the codes it names, its columns."""
    path: str
    kind: str                       # "csv" | "yaml" | "cpp"
    codes: List[str]                # declared PDB codes, uppercase
    columns: Dict[str, List[str]]   # column -> values (csv only)
    ligand_pairs: List[Tuple[str, str]]   # (entry_code, ligand_code)
    declared_constants: Dict[str, str]    # column -> justification


# ---------------------------------------------------------------------------
# Collecting declarations
# ---------------------------------------------------------------------------

def list_files(root: Path) -> List[str]:
    """Tracked files if git is usable, else a filtered walk. Never silently empty."""
    env = dict(os.environ, GIT_CONFIG_GLOBAL="/dev/null", GIT_CONFIG_SYSTEM="/dev/null")
    try:
        r = subprocess.run(["git", "ls-files", "-z"], cwd=str(root), env=env,
                           capture_output=True, text=True, timeout=120)
        if r.returncode == 0 and r.stdout.strip():
            return [f for f in r.stdout.split("\0") if f]
    except (OSError, subprocess.SubprocessError):
        pass
    # The fallback must select the SAME files git mode would, or the two modes
    # disagree on the scan surface and a quarantine entry drifts when the gate
    # runs on an exported tree. Suffix matching is not enough:
    # astex_diverse_set.csv.orig has suffix ".orig" and is still a roster that
    # code reads. Match on the full name.
    out: List[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in _SKIP_DIR_PARTS]
        for fn in filenames:
            low = fn.lower()
            if low.endswith((".csv", ".csv.orig", ".yaml", ".yml", ".cpp")):
                out.append(str((Path(dirpath) / fn).relative_to(root)))
    return out


def _yaml_constants(doc: dict) -> Dict[str, str]:
    """Read a definition's declared constant columns.

    Two shapes are accepted:
        constant_columns:
          - column: source_doi
            justification: "..."
    and the mapping form ``{column: justification}``. The list form is
    preferred because a mapping key ending in ``doi`` collides with the
    citation gate's field rule.
    """
    cc = doc.get("constant_columns")
    if isinstance(cc, dict):
        return {str(k): str(v) for k, v in cc.items()}
    if isinstance(cc, list):
        out = {}
        for item in cc:
            if isinstance(item, dict) and item.get("column"):
                out[str(item["column"])] = str(item.get("justification", ""))
        return out
    return {}


def collect(root: Path, files: Optional[Iterable[str]] = None
            ) -> Tuple[List[Declaration], Dict[str, int]]:
    if files is None:
        files = list_files(root)
    decls: List[Declaration] = []
    n_csv = n_yaml = n_cpp = 0

    for rel in files:
        rel_norm = rel.replace(os.sep, "/")
        if rel_norm in _SELF or rel_norm.startswith("tests/fixtures/dataset_identity/"):
            continue
        if any(part in _SKIP_DIR_PARTS for part in Path(rel_norm).parts[:-1]):
            continue
        p = root / rel
        if not p.is_file():
            continue
        suffix = p.suffix.lower()
        name = p.name.lower()

        # --- roster CSVs (including .csv.orig backups, which are still read) ---
        if suffix == ".csv" or name.endswith(".csv.orig"):
            try:
                with p.open(newline="") as fh:
                    rows = list(csv.DictReader(fh))
            except (OSError, UnicodeDecodeError, csv.Error):
                continue
            if not rows or not rows[0]:
                continue
            headers = [h for h in rows[0].keys() if h]
            pdb_cols = [h for h in headers if h.strip().lower() in _PDB_COLUMNS]
            if not pdb_cols:
                continue
            codes, pairs = [], []
            for r in rows:
                for c in pdb_cols:
                    v = (r.get(c) or "").strip()
                    if _PDB_CODE.match(v):
                        codes.append(v.upper())
                lig_cols = [h for h in headers if h.strip().lower() in _LIGAND_COLUMNS]
                anchor = (r.get("pdb_id") or r.get("ligand_pdb") or r.get("target_pdb") or "").strip()
                for lc in lig_cols:
                    lv = (r.get(lc) or "").strip()
                    # A ligand column may legitimately hold a PDB code (the
                    # entry the ligand comes from), not a component code.
                    if lv and not _PDB_CODE.match(lv) and _PDB_CODE.match(anchor):
                        pairs.append((anchor.upper(), lv.upper()))
            columns = {h: [(r.get(h) or "").strip() for r in rows] for h in headers}
            decls.append(Declaration(rel_norm, "csv", codes, columns, pairs, {}))
            n_csv += 1
            continue

        # --- dataset definition YAMLs ---
        if suffix in {".yaml", ".yml"}:
            try:
                import yaml  # lazy: the offline gate should not need it until here
                doc = yaml.safe_load(p.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001 - a malformed YAML is another gate's problem
                continue
            if not isinstance(doc, dict):
                continue
            tgts = doc.get("targets")
            if not isinstance(tgts, list):
                continue
            codes = [str(t).upper() for t in tgts if _PDB_CODE.match(str(t))]
            if not codes:
                continue
            nd = doc.get("non_dockable") or {}
            removed = doc.get("removed_nonexistent_codes") or []
            # Codes the definition itself flags as absent or unusable are
            # declarations ABOUT a defect, not declarations OF a target.
            excluded = {str(k).upper() for k in (nd.keys() if isinstance(nd, dict) else nd)}
            excluded |= {str(k).upper() for k in removed}
            codes = [c for c in codes if c not in excluded]
            decls.append(Declaration(rel_norm, "yaml", codes, {}, [], _yaml_constants(doc)))
            n_yaml += 1
            continue

        # --- hardcoded C++ lists (the 3QGS/3RP3 surface) ---
        if suffix == ".cpp":
            text = p.read_text(encoding="utf-8", errors="replace")
            hits: List[str] = []
            for m in re.finditer(r"\b(\w*_codes)\s*\(\s*\)\s*\{", text):
                seg = text[m.end():m.end() + 20000]
                depth = 1
                for i, ch in enumerate(seg):
                    if ch == "{":
                        depth += 1
                    elif ch == "}":
                        depth -= 1
                        if depth == 0:
                            seg = seg[:i]
                            break
                body = "\n".join(l for l in seg.splitlines() if not l.strip().startswith("//"))
                hits.extend(c.upper() for c in re.findall(r'"([0-9][A-Za-z][A-Za-z0-9]{2})"', body))
            if hits:
                decls.append(Declaration(rel_norm, "cpp", hits, {}, [], {}))
                n_cpp += 1

    # A roster CSV's legitimate constants are declared in the dataset DEFINITION
    # that owns it, not in the CSV (a CSV cannot carry a justification). Link
    # them: any yaml that names a .csv path lends that CSV its constant_columns.
    owned: Dict[str, Dict[str, str]] = {}
    for d in decls:
        if d.kind != "yaml" or not d.declared_constants:
            continue
        try:
            import yaml
            doc = yaml.safe_load((root / d.path).read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            continue

        def walk(node):
            if isinstance(node, dict):
                for v in node.values():
                    walk(v)
            elif isinstance(node, list):
                for v in node:
                    walk(v)
            elif isinstance(node, str) and node.strip().endswith(".csv"):
                owned.setdefault(node.strip().lstrip("./"), {}).update(d.declared_constants)

        walk(doc)
    if owned:
        decls = [d._replace(declared_constants=owned.get(d.path, d.declared_constants))
                 if d.kind == "csv" else d for d in decls]

    stats = {
        "declarations": len(decls),
        "csv_rosters": n_csv,
        "yaml_definitions": n_yaml,
        "cpp_lists": n_cpp,
        "codes_declared": sum(len(d.codes) for d in decls),
        "distinct_codes": len({c for d in decls for c in d.codes}),
        "columns_inspected": sum(len(d.columns) for d in decls),
        "ligand_fields_inspected": sum(len(d.ligand_pairs) for d in decls),
    }
    return decls, stats


# ---------------------------------------------------------------------------
# Online resolution (only under --refresh)
# ---------------------------------------------------------------------------

GRAPHQL = "https://data.rcsb.org/graphql"
_QUERY = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
  rcsb_entry_info{ resolution_combined experimental_method }
  nonpolymer_entities{ nonpolymer_comp{ chem_comp{ id } } } } }"""


def _fetch(ids: List[str]) -> Dict[str, dict]:
    import time
    import requests  # imported lazily: offline mode must not need it

    out: Dict[str, dict] = {}
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        for attempt in range(4):
            try:
                r = requests.post(GRAPHQL, json={"query": _QUERY, "variables": {"ids": chunk}},
                                  timeout=90)
                if r.status_code == 200:
                    data = (r.json().get("data") or {}).get("entries") or []
                    for e in data:
                        if not e:
                            continue
                        res = (e["rcsb_entry_info"] or {}).get("resolution_combined")
                        out[e["rcsb_id"].upper()] = {
                            "exists": True,
                            "resolution_A": (res[0] if res else None),
                            "experimental_method": (e["rcsb_entry_info"] or {}).get("experimental_method"),
                            "nonpolymer_comp_ids": sorted(
                                x["nonpolymer_comp"]["chem_comp"]["id"]
                                for x in (e["nonpolymer_entities"] or [])),
                        }
                    break
            except Exception:  # noqa: BLE001 - network layer
                pass
            time.sleep(1.5 * (attempt + 1))
        else:
            raise RuntimeError("RCSB unreachable for chunk starting %s" % chunk[0])
    return out


def refresh(root: Path) -> int:
    decls, stats = collect(root)
    codes = sorted({c for d in decls for c in d.codes} | {POSITIVE_CONTROL, NEGATIVE_CONTROL})
    print("[refresh] %d declarations (%d csv, %d yaml, %d cpp), %d distinct codes"
          % (stats["declarations"], stats["csv_rosters"], stats["yaml_definitions"],
             stats["cpp_lists"], len(codes)))
    fetched = _fetch(codes)

    pos = fetched.get(POSITIVE_CONTROL)
    neg = fetched.get(NEGATIVE_CONTROL)
    pos_ok = bool(pos) and POSITIVE_CONTROL_LIGAND in (pos.get("nonpolymer_comp_ids") or [])
    neg_ok = neg is None
    print("[refresh] positive control %s -> resolved=%s, contains %s=%s"
          % (POSITIVE_CONTROL, bool(pos), POSITIVE_CONTROL_LIGAND, pos_ok))
    print("[refresh] negative control %s -> resolved=%s (must be False)"
          % (NEGATIVE_CONTROL, not neg_ok))
    if not pos_ok or not neg_ok:
        print("VOID: controls did not behave. A 404 from a blocked or rate-limited host "
              "is not evidence that a code is fake. Cache NOT written.")
        return 2

    cache = {c: fetched.get(c, {"exists": False, "resolution_A": None,
                                "experimental_method": None, "nonpolymer_comp_ids": []})
             for c in codes}
    n_absent = sum(1 for v in cache.values() if not v["exists"])
    payload = {
        "_comment": "Committed PDB-entry cache for scripts/check_dataset_identity.py. "
                    "Refresh explicitly with --refresh; never hand-edit.",
        "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "oracle": "RCSB data API GraphQL entries{rcsb_entry_info, nonpolymer_entities}",
        "controls": {"positive": POSITIVE_CONTROL, "positive_ligand": POSITIVE_CONTROL_LIGAND,
                     "positive_ok": True, "negative": NEGATIVE_CONTROL, "negative_resolved": False},
        "counts": {"distinct": len(cache), "present": len(cache) - n_absent, "absent": n_absent},
        "entries": cache,
    }
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_PATH.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print("[refresh] wrote %s: %d codes, %d present, %d absent"
          % (CACHE_PATH.relative_to(root), len(cache), len(cache) - n_absent, n_absent))
    return 0


# ---------------------------------------------------------------------------
# Offline gate
# ---------------------------------------------------------------------------

def _is_measurement_column(col: str) -> bool:
    c = col.strip().lower()
    if not c:
        return False
    return not any(tok in c for tok in _NOT_A_MEASUREMENT)


def check(root: Path, cache: Dict[str, dict], defects: Dict[str, dict],
          files: Optional[Iterable[str]] = None) -> Tuple[List[Finding], dict]:
    decls, stats = collect(root, files)
    findings: List[Finding] = []
    seen_defects = set()
    n_code_checks = n_col_checks = n_lig_checks = 0

    def quarantined(path: str, kind: str, detail_key: str) -> bool:
        key = "%s::%s::%s" % (path, kind, detail_key)
        if key in defects:
            seen_defects.add(key)
            return True
        return False

    for d in decls:
        # Rule 1 / 4 - declared codes must exist
        for code in sorted(set(d.codes)):
            n_code_checks += 1
            entry = cache.get(code)
            if entry is None:
                if not quarantined(d.path, "UNCACHED", code):
                    findings.append(Finding("UNCACHED", d.path,
                                            "%s has no cached PDB resolution; run --refresh" % code))
            elif not entry.get("exists"):
                if not quarantined(d.path, "PDB_CODE_ABSENT", code):
                    findings.append(Finding("PDB_CODE_ABSENT", d.path,
                                            "%s is declared as a target but the PDB does not serve it" % code))

        # Rule 2 - no undeclared constant per-target column
        nrows = max((len(v) for v in d.columns.values()), default=0)
        if nrows > 1:
            for col, values in d.columns.items():
                if not _is_measurement_column(col):
                    continue
                n_col_checks += 1
                if len(set(values)) == 1 and col not in d.declared_constants:
                    if not quarantined(d.path, "CONSTANT_COLUMN", col):
                        findings.append(Finding(
                            "CONSTANT_COLUMN", d.path,
                            "column %r holds one value (%r) on all %d rows: that is a default, "
                            "not a measurement. Populate it, delete it, or declare it under "
                            "constant_columns with a justification."
                            % (col, values[0], nrows)))

        # Rule 3 - a declared ligand code must be a component of its own entry
        for entry_code, lig in d.ligand_pairs:
            n_lig_checks += 1
            entry = cache.get(entry_code)
            if entry is None or not entry.get("exists"):
                continue  # rule 1 already reported it
            if lig not in (entry.get("nonpolymer_comp_ids") or []):
                if not quarantined(d.path, "LIGAND_CODE_ABSENT", "%s/%s" % (entry_code, lig)):
                    findings.append(Finding(
                        "LIGAND_CODE_ABSENT", d.path,
                        "ligand %r is declared for %s, which has components %s"
                        % (lig, entry_code, entry.get("nonpolymer_comp_ids") or [])))

    # Rule 5 - the quarantine ledger must not outlive the defect
    for key, meta in defects.items():
        if key not in seen_defects:
            findings.append(Finding(
                "DEFECT_LEDGER_DRIFT", DEFECTS_PATH.name,
                "%s is quarantined as a known defect but no longer reproduces; drop the entry "
                "(reason on file: %s)" % (key, meta.get("reason", "?"))))

    stats.update({
        "code_checks": n_code_checks,
        "column_checks": n_col_checks,
        "ligand_checks": n_lig_checks,
        "cached": len(cache),
        "quarantined": len(defects),
        "quarantine_hits": len(seen_defects),
    })
    return findings, stats


# ---------------------------------------------------------------------------
# Non-vacuity self-test
# ---------------------------------------------------------------------------

def selftest(root: Path) -> int:
    """Prove the gate is non-vacuous: it must PASS a known-good tree and FAIL a
    known-bad tree ON EVERY RULE, and every rule must inspect a NON-ZERO number
    of fields on the real tree.

    The citation gate shipped with a rule that matched ZERO fields live while
    passing its own fixture - a green run that checked nothing. The denominator
    assertions below exist so that cannot happen here.
    """
    doc = json.loads(CACHE_PATH.read_text(encoding="utf-8")) if CACHE_PATH.exists() else {}
    cache = doc.get("entries", {})
    if not cache:
        print("VOID: no PDB entry cache at %s - cannot run the non-vacuity proof." % CACHE_PATH)
        return 2
    present = [c for c, v in cache.items() if v.get("exists") and v.get("nonpolymer_comp_ids")]
    absent = [c for c, v in cache.items() if not v.get("exists")]
    if not present:
        print("VOID: cache has no resolvable entry with a ligand; the pass arm is meaningless.")
        return 2
    if not absent:
        print("VOID: cache has no absent code, so the PDB_CODE_ABSENT arm cannot be exercised.")
        return 2

    good_code = POSITIVE_CONTROL if POSITIVE_CONTROL in present else present[0]
    good_lig = cache[good_code]["nonpolymer_comp_ids"][0]
    other = [c for c in present if c != good_code][:3]
    bad_code = absent[0]
    # A component id that is not in the entry. Assembled at runtime so a literal
    # cannot drift into the cache and silently disarm the arm.
    fake_lig = "Z" + "ZZ"

    rows_ok = ["pdb_id,ligand_id,resolution_A"]
    for i, c in enumerate([good_code] + other):
        lig = cache[c]["nonpolymer_comp_ids"][0]
        rows_ok.append("%s,%s,%s" % (c, lig, 1.5 + 0.1 * i))
    good_csv = "\n".join(rows_ok) + "\n"

    bad_fixtures = {
        "absent_code": "pdb_id,ligand_id,resolution_A\n%s,%s,1.5\n%s,%s,1.7\n"
                       % (bad_code, good_lig, good_code, good_lig),
        "constant_column": "pdb_id,ligand_id,resolution_A\n"
                           + "".join("%s,%s,2.0\n" % (c, cache[c]["nonpolymer_comp_ids"][0])
                                     for c in ([good_code] + other)),
        "absent_ligand": "pdb_id,ligand_id,resolution_A\n%s,%s,1.5\n%s,%s,1.7\n"
                         % (good_code, fake_lig, other[0] if other else good_code, good_lig),
    }

    results: Dict[str, List[Finding]] = {}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        (tmp / "good_set.csv").write_text(good_csv, encoding="utf-8")
        results["known_good"], _ = check(tmp, cache, {}, files=["good_set.csv"])
        for name, body in bad_fixtures.items():
            fn = "%s_set.csv" % name
            (tmp / fn).write_text(body, encoding="utf-8")
            results[name], _ = check(tmp, cache, {}, files=[fn])
        # Ledger drift: a quarantine entry that does not reproduce must fire.
        results["ledger_drift"], _ = check(
            tmp, cache, {"nonexistent.csv::CONSTANT_COLUMN::ghost": {"reason": "selftest"}},
            files=["good_set.csv"])
        # And a real defect that IS quarantined must NOT fire.
        qkey = "constant_column_set.csv::CONSTANT_COLUMN::resolution_A"
        quiet, _ = check(tmp, cache, {qkey: {"reason": "selftest"}},
                         files=["constant_column_set.csv"])
        results["_quarantine_suppresses"] = quiet

    print("[selftest] known-good fixture         -> %d findings (expected 0)"
          % len(results["known_good"]))
    # Each arm must fire its INTENDED rule. Counting findings is not enough: an
    # arm can be satisfied by an unrelated rule firing, which would let the rule
    # under test be silently dead.
    expected_kind = {
        "absent_code": "PDB_CODE_ABSENT",
        "constant_column": "CONSTANT_COLUMN",
        "absent_ligand": "LIGAND_CODE_ABSENT",
        "ledger_drift": "DEFECT_LEDGER_DRIFT",
    }
    arms = {k: v for k, v in results.items() if k not in ("known_good", "_quarantine_suppresses")}
    arms_ok = True
    for k, v in arms.items():
        kinds = {f.kind for f in v}
        want = expected_kind[k]
        hit = want in kinds
        arms_ok = arms_ok and hit
        print("[selftest] known-bad  %-18s -> %d findings, %s present: %s%s"
              % (k, len(v), want, hit, "" if hit else "   *** INTENDED RULE DID NOT FIRE ***"))
        for f in v:
            print("             %s: %s" % (f.kind, f.detail[:110]))
    print("[selftest] quarantined defect         -> %d findings (expected 0)"
          % len(results["_quarantine_suppresses"]))

    ok = (not results["known_good"]
          and all(len(v) > 0 for v in arms.values())
          and arms_ok
          and not results["_quarantine_suppresses"])

    # Denominators on the REAL tree: a rule that inspects nothing is vacuous
    # regardless of how its fixture behaves.
    _, live = collect(ROOT)
    live_findings, live_stats = check(ROOT, cache, {})
    print("[selftest] live denominators: %d code checks, %d column checks, %d ligand checks "
          "across %d declarations (%d csv / %d yaml / %d cpp)"
          % (live_stats["code_checks"], live_stats["column_checks"], live_stats["ligand_checks"],
             live_stats["declarations"], live_stats["csv_rosters"],
             live_stats["yaml_definitions"], live_stats["cpp_lists"]))
    for name, key in (("PDB_CODE_ABSENT/UNCACHED", "code_checks"),
                      ("CONSTANT_COLUMN", "column_checks"),
                      ("LIGAND_CODE_ABSENT", "ligand_checks")):
        if live_stats[key] == 0:
            print("VOID: rule %s inspected ZERO fields on the real tree. It is vacuous." % name)
            ok = False

    if not ok:
        print("VOID: the gate did not produce BOTH a pass on known-good and a failure on every "
              "known-bad arm with non-zero live denominators. It proves nothing in this state.")
        return 2
    print("[selftest] NON-VACUOUS: passes known-good, fails all %d known-bad arms, "
          "quarantine suppresses, and every rule inspects a non-zero number of live fields." % len(arms))
    return 0


# ---------------------------------------------------------------------------

def load_json(path: Path, default):
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--refresh", action="store_true", help="re-resolve online, rewrite cache")
    ap.add_argument("--selftest", action="store_true", help="prove the gate is non-vacuous")
    ap.add_argument("--report", action="store_true", help="print the full sweep table")
    args = ap.parse_args(argv)

    if args.refresh:
        return refresh(ROOT)
    if args.selftest:
        return selftest(ROOT)

    cache_doc = load_json(CACHE_PATH, {})
    cache = cache_doc.get("entries", {})
    defects = load_json(DEFECTS_PATH, {}).get("quarantined", {})
    if not cache:
        print("VOID: no PDB entry cache at %s. Run --refresh (needs network)." % CACHE_PATH)
        return 2

    findings, stats = check(ROOT, cache, defects)
    print("Dataset identity sweep - %d declarations (%d csv rosters, %d yaml definitions, "
          "%d cpp lists), %d distinct codes, %d code checks, %d column checks, %d ligand checks, "
          "%d cached, %d quarantined"
          % (stats["declarations"], stats["csv_rosters"], stats["yaml_definitions"],
             stats["cpp_lists"], stats["distinct_codes"], stats["code_checks"],
             stats["column_checks"], stats["ligand_checks"], stats["cached"], stats["quarantined"]))
    print("cache checked_utc=%s oracle=%s" % (cache_doc.get("checked_utc"), cache_doc.get("oracle")))

    if args.report:
        decls, _ = collect(ROOT)
        for d in sorted(decls, key=lambda x: x.path):
            print("  %-62s %-4s codes=%-4d cols=%-3d ligfields=%d"
                  % (d.path, d.kind, len(set(d.codes)), len(d.columns), len(d.ligand_pairs)))

    if findings:
        print("\nFAIL - %d finding(s):" % len(findings))
        for f in findings:
            print("  [%s] %s  %s" % (f.kind, f.path, f.detail))
        return 1
    print("PASS - every declared PDB code is served by the PDB, every declared ligand code is a "
          "component of its own entry, and no per-target column is an undeclared constant.")
    return 0


if __name__ == "__main__":
    sys.exit(main())