#!/usr/bin/env python3
"""check_run_receipt_codes.py - fail a run whose EXECUTED roster does not match
the canonical source it claims to implement.

Motivation
----------
Everything a RUN_RECEIPT.json recorded described HOW a run was configured --
matrix hash, binary hash, GA knobs, a full ProtocolConfig snapshot. None of it
recorded WHAT set of targets was enumerated. So a roster that was 72/85 of the
published Astex Diverse Set ran for months and not one receipt it wrote could
have contradicted the claim that it was the published set: the executed list was
never written down, and there was nothing to check it against.

Schema 2 of the receipt (LIB/RunReceipt.h) records the enumerated roster, the
sha256 of a rehashable sidecar holding it, and the canonical source path + sha256
that the engine's code list was GENERATED from
(LIB/generated/dataset_codes.inc, scripts/gen_dataset_codes.py). This script is
the check that makes those fields load-bearing.

What it fails on
----------------
  SIDECAR_SHA_MISMATCH  executed_codes_sha256 is not the sha256 of the sidecar,
                        or not the sha256 of the JSON array. Either the list or
                        the digest was edited after the run.
  SOURCE_SHA_MISMATCH   the canonical source in the repo today does not hash to
                        what the receipt recorded. The roster definition moved
                        under the run.
  ROSTER_MISMATCH       the executed codes are not a prefix-preserving subset of
                        the declared source's codes: the run executed something
                        the declared source does not contain.
  DECLARED_COUNT_MISMATCH the declared count in the receipt disagrees with the
                        source's actual row count.
  MISSING_FIELD         a schema-2 receipt with a roster field stripped or empty.

What it reports without failing
-------------------------------
  UNCHECKABLE_SCHEMA    schema_version < 2: the run predates roster recording.
                        This is NOT a pass -- it is the absence of evidence, and
                        --strict turns it into a failure.
  NO_DECLARED_SOURCE    the dataset's code list is not generated yet (e.g.
                        astex_nonnative). The executed list is still recorded and
                        its sidecar sha is still verified.
  INCOMPLETE_RUN        executed < declared. Legitimate (preparation dropped
                        targets) but always reported with both numbers, because a
                        silently short roster is how a 72/85 set passes for 85.

Usage
-----
    python3 scripts/check_run_receipt_codes.py <run_dir> [<run_dir> ...]
    python3 scripts/check_run_receipt_codes.py --strict <run_dir>   # schema<2 fails
    python3 scripts/check_run_receipt_codes.py --selftest           # non-vacuity proof
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional

ROOT = Path(__file__).resolve().parents[1]
MIN_SCHEMA = 2

# Fields a schema-2 receipt must carry non-empty. dataset_source_* are NOT here:
# a dataset with no generated code list legitimately leaves them blank, and that
# case is reported as NO_DECLARED_SOURCE rather than a missing field.
_REQUIRED = ("executed_codes", "executed_codes_sha256", "executed_codes_count")


class Finding(NamedTuple):
    kind: str
    run: str
    detail: str


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_source_codes(rel: str, repo: Path) -> Optional[List[str]]:
    """Read a canonical source the same way scripts/gen_dataset_codes.py does.

    Returns None when the source is absent from this checkout -- reported as a
    finding by the caller, never treated as an empty roster.
    """
    p = repo / rel
    if not p.is_file():
        return None
    if p.suffix.lower() == ".csv":
        with p.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        if not rows or "pdb_id" not in rows[0]:
            return None
        return [(r.get("pdb_id") or "").strip().upper() for r in rows]
    if p.suffix.lower() in {".yaml", ".yml"}:
        import yaml
        doc = yaml.safe_load(p.read_text(encoding="utf-8"))
        if not isinstance(doc, dict) or not isinstance(doc.get("targets"), list):
            return None
        return [str(x).strip().upper() for x in doc["targets"]]
    return None


def check_receipt(receipt_path: Path, repo: Path = ROOT,
                  strict: bool = False) -> tuple[List[Finding], Dict[str, object]]:
    run = str(receipt_path.parent)
    findings: List[Finding] = []
    info: Dict[str, object] = {"run": run, "checks": 0}

    try:
        doc = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as ex:
        return [Finding("UNREADABLE_RECEIPT", run, str(ex)[:160])], info

    schema = doc.get("schema_version")
    info["schema_version"] = schema
    info["dataset"] = doc.get("dataset")
    if not isinstance(schema, int) or schema < MIN_SCHEMA:
        kind = "UNCHECKABLE_SCHEMA" if not strict else "UNCHECKABLE_SCHEMA_STRICT"
        findings.append(Finding(
            kind, run,
            "schema_version=%r < %d: this run did not record its executed roster, so "
            "there is nothing to check it against. Absence of evidence, not a pass."
            % (schema, MIN_SCHEMA)))
        return findings, info

    for key in _REQUIRED:
        if key not in doc or doc.get(key) in (None, "", []):
            findings.append(Finding("MISSING_FIELD", run,
                                    "schema-2 receipt has no usable %r" % key))
    if findings:
        return findings, info

    executed = [str(c).upper() for c in doc["executed_codes"]]
    info["executed_count"] = len(executed)

    # 1. the count field must agree with the array
    info["checks"] = int(info["checks"]) + 1
    if int(doc.get("executed_codes_count", -1)) != len(executed):
        findings.append(Finding("SIDECAR_SHA_MISMATCH", run,
                                "executed_codes_count=%r but the array holds %d"
                                % (doc.get("executed_codes_count"), len(executed))))

    # 2. the digest must be the digest of the list, and of the sidecar on disk
    want = sha256_text("\n".join(executed) + "\n")
    got = str(doc.get("executed_codes_sha256", ""))
    info["checks"] = int(info["checks"]) + 1
    if got != want:
        findings.append(Finding(
            "SIDECAR_SHA_MISMATCH", run,
            "executed_codes_sha256=%s but sha256 of the recorded list is %s" % (got, want)))
    sidecar = doc.get("executed_codes_file")
    if sidecar:
        sp = Path(str(sidecar))
        if not sp.is_absolute():
            sp = receipt_path.parent / sp
        if sp.is_file():
            info["checks"] = int(info["checks"]) + 1
            have = sha256_file(sp)
            if have != got:
                findings.append(Finding(
                    "SIDECAR_SHA_MISMATCH", run,
                    "sidecar %s hashes to %s, receipt records %s" % (sp.name, have, got)))
        else:
            info["sidecar_absent"] = str(sp)

    # 3. what it claims to implement
    src_rel = str(doc.get("dataset_source_path") or "")
    if not src_rel:
        findings.append(Finding(
            "NO_DECLARED_SOURCE", run,
            "dataset=%r has no generated code list, so the executed roster names no "
            "canonical source. Its sidecar digest was still verified."
            % doc.get("dataset")))
        return findings, info

    info["dataset_source_path"] = src_rel
    info["dataset_source_provenance"] = doc.get("dataset_source_provenance")
    src = repo / src_rel
    if not src.is_file():
        findings.append(Finding("SOURCE_SHA_MISMATCH", run,
                                "declared source %s is absent from this checkout" % src_rel))
        return findings, info

    info["checks"] = int(info["checks"]) + 1
    live_sha = sha256_file(src)
    if live_sha != str(doc.get("dataset_source_sha256", "")):
        findings.append(Finding(
            "SOURCE_SHA_MISMATCH", run,
            "%s hashes to %s today; the run recorded %s. The roster definition changed "
            "after the run, so the run's claim about which set it executed cannot be "
            "confirmed from this checkout."
            % (src_rel, live_sha, doc.get("dataset_source_sha256"))))

    declared = read_source_codes(src_rel, repo)
    if declared is None:
        findings.append(Finding("SOURCE_SHA_MISMATCH", run,
                                "declared source %s is not readable as a roster" % src_rel))
        return findings, info
    info["declared_count_live"] = len(declared)

    info["checks"] = int(info["checks"]) + 1
    if int(doc.get("declared_codes_count", -1)) != len(declared):
        findings.append(Finding(
            "DECLARED_COUNT_MISMATCH", run,
            "receipt declares %r codes in %s; it holds %d today"
            % (doc.get("declared_codes_count"), src_rel, len(declared))))

    # 4. THE check: nothing executed that the declared source does not contain
    info["checks"] = int(info["checks"]) + 1
    extra = [c for c in executed if c not in set(declared)]
    if extra:
        findings.append(Finding(
            "ROSTER_MISMATCH", run,
            "%d executed code(s) are absent from the declared source %s: %s"
            % (len(extra), src_rel, extra[:10])))
    missing = [c for c in declared if c not in set(executed)]
    if missing:
        findings.append(Finding(
            "INCOMPLETE_RUN", run,
            "executed %d of %d declared codes from %s; %d not executed: %s"
            % (len(executed), len(declared), src_rel, len(missing), missing[:10])))

    return findings, info


_FAILING = {"UNREADABLE_RECEIPT", "MISSING_FIELD", "SIDECAR_SHA_MISMATCH",
            "SOURCE_SHA_MISMATCH", "ROSTER_MISMATCH", "DECLARED_COUNT_MISMATCH",
            "UNCHECKABLE_SCHEMA_STRICT"}


def is_failure(kind: str) -> bool:
    return kind in _FAILING


# ---------------------------------------------------------------------------
# Non-vacuity self-test
# ---------------------------------------------------------------------------

def _synth_receipt(codes: List[str], src_rel: str, repo: Path,
                   **over) -> dict:
    declared = read_source_codes(src_rel, repo) or []
    body = {
        "schema_version": 2,
        "run_id": "selftest", "dataset": "selftest",
        "executed_codes": codes,
        "executed_codes_count": len(codes),
        "executed_codes_sha256": sha256_text("\n".join(codes) + "\n"),
        "executed_codes_file": "executed_codes.txt",
        "declared_codes_count": len(declared),
        "dataset_source_path": src_rel,
        "dataset_source_sha256": sha256_file(repo / src_rel),
        "dataset_source_provenance": "ROSTER_CSV",
        "dataset_source_primary": "UNVERIFIED",
    }
    body.update(over)
    return body


def selftest(repo: Path = ROOT) -> int:
    """Prove the validator is non-vacuous: it must PASS a truthful receipt and
    FAIL each way of lying, with a non-zero number of checks actually performed."""
    src_rel = "benchmarks/astex_diverse/astex_diverse_set.csv"
    declared = read_source_codes(src_rel, repo)
    if not declared:
        print("VOID: cannot read %s; the proof has no roster to work with." % src_rel)
        return 2

    arms: Dict[str, tuple[dict, Optional[str]]] = {}
    arms["truthful"] = (_synth_receipt(declared, src_rel, repo), None)
    arms["tampered_digest"] = (
        _synth_receipt(declared, src_rel, repo,
                       executed_codes_sha256="0" * 64), "SIDECAR_SHA_MISMATCH")
    arms["tampered_count"] = (
        _synth_receipt(declared, src_rel, repo,
                       executed_codes_count=len(declared) + 1), "SIDECAR_SHA_MISMATCH")
    arms["stale_source_sha"] = (
        _synth_receipt(declared, src_rel, repo,
                       dataset_source_sha256="f" * 64), "SOURCE_SHA_MISMATCH")
    # A code the declared source does not contain. Built at runtime so a literal
    # cannot drift into the roster and silently disarm the arm.
    ghost = "9" + "ZZZ"
    arms["executed_a_code_not_declared"] = (
        _synth_receipt(declared + [ghost], src_rel, repo), "ROSTER_MISMATCH")
    arms["short_roster"] = (
        _synth_receipt(declared[:-3], src_rel, repo), "INCOMPLETE_RUN")
    arms["schema_1"] = ({"schema_version": 1, "dataset": "astex_diverse"},
                        "UNCHECKABLE_SCHEMA")
    arms["schema_2_stripped"] = ({"schema_version": 2, "dataset": "astex_diverse"},
                                 "MISSING_FIELD")

    ok = True
    total_checks = 0
    with tempfile.TemporaryDirectory() as td:
        for name, (body, want) in arms.items():
            d = Path(td) / name
            d.mkdir()
            (d / "RUN_RECEIPT.json").write_text(json.dumps(body, indent=2), encoding="utf-8")
            if "executed_codes" in body:
                (d / "executed_codes.txt").write_text(
                    "\n".join(body["executed_codes"]) + "\n", encoding="utf-8")
            findings, info = check_receipt(d / "RUN_RECEIPT.json", repo)
            kinds = {f.kind for f in findings}
            total_checks += int(info.get("checks", 0))
            if want is None:
                hit = not findings
            else:
                hit = want in kinds
            ok = ok and hit
            print("[selftest] %-30s checks=%-2s -> %-22s expected %-22s %s"
                  % (name, info.get("checks"), sorted(kinds) or "none",
                     want or "no findings", "OK" if hit else "*** ARM DID NOT FIRE ***"))

    # Denominator: an arm can "pass" while the validator compared nothing at all.
    if total_checks == 0:
        print("VOID: the validator performed ZERO comparisons across every arm.")
        ok = False
    if not ok:
        print("VOID: the validator did not pass the truthful receipt and fail every way of "
              "lying. It proves nothing in this state.")
        return 2
    print("[selftest] NON-VACUOUS: passes a truthful receipt, fails all %d lying arms, and "
          "performed %d field comparisons across them."
          % (len(arms) - 1, total_checks))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dirs", nargs="*", type=Path,
                    help="run output directories (or RUN_RECEIPT.json paths)")
    ap.add_argument("--strict", action="store_true",
                    help="treat a pre-schema-2 receipt as a failure, not a report")
    ap.add_argument("--selftest", action="store_true", help="prove non-vacuity")
    ap.add_argument("--repo", type=Path, default=ROOT,
                    help="repository root holding the canonical sources")
    args = ap.parse_args(argv)

    if args.selftest:
        return selftest(args.repo)
    if not args.run_dirs:
        ap.error("give at least one run directory, or --selftest")

    all_findings: List[Finding] = []
    n_checked = n_checks = 0
    for rd in args.run_dirs:
        rp = rd if rd.name == "RUN_RECEIPT.json" else rd / "RUN_RECEIPT.json"
        if not rp.is_file():
            all_findings.append(Finding("UNREADABLE_RECEIPT", str(rd), "no RUN_RECEIPT.json"))
            continue
        findings, info = check_receipt(rp, args.repo, strict=args.strict)
        n_checked += 1
        n_checks += int(info.get("checks", 0))
        print("%s  schema=%s dataset=%s executed=%s declared_live=%s checks=%s"
              % (rp.parent, info.get("schema_version"), info.get("dataset"),
                 info.get("executed_count"), info.get("declared_count_live"),
                 info.get("checks")))
        for f in findings:
            print("  [%s] %s" % (f.kind, f.detail))
        all_findings.extend(findings)

    fails = [f for f in all_findings if is_failure(f.kind)]
    reports = [f for f in all_findings if not is_failure(f.kind)]
    print("\n%d receipt(s), %d field comparisons, %d failure(s), %d report(s)"
          % (n_checked, n_checks, len(fails), len(reports)))
    if fails:
        print("FAIL - an executed roster does not match the source it claims to implement.")
        return 1
    if n_checks == 0:
        print("VOID - zero comparisons were performed. Nothing was verified.")
        return 2
    print("PASS - every checked run executed only codes its declared canonical source "
          "contains, and every recorded roster hashes to its recorded digest.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
