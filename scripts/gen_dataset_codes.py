#!/usr/bin/env python3
"""gen_dataset_codes.py - emit the engine's PDB code lists FROM their canonical source.

Motivation
----------
The same 85 Astex codes lived in six places: a hardcoded ``std::vector<std::string>``
in ``LIB/DatasetRunner.cpp``, two copies of a dataset YAML, a roster CSV, a target
manifest, and untracked build output. Nothing enforced agreement, and nothing
recorded which of them the engine actually read. The roster that ran was the C++
literal list, and it differed from the published Astex Diverse Set on 13 of 85
targets for months without a single check firing -- a gate comparing the YAMLs to
each other would have passed the whole time.

This script removes the hand-maintained copy. The C++ list is now GENERATED from a
single declared canonical source per dataset, and ``--check`` fails the build/test
if the committed generated file disagrees with what the canonical source produces.
The generated header also carries, for every list, the canonical source path and
its sha256, so a run receipt can record the source it claims to implement rather
than asserting it (see scripts/check_run_receipt_codes.py).

What is canonical, and what that word is worth here
---------------------------------------------------
"Canonical" means "the one file the generator reads" -- a mechanical guarantee that
the engine and the definition cannot diverge. It does NOT mean "verified against
the primary literature". Each entry carries an explicit ``provenance`` field:

  ROSTER_CSV          the code list is the roster CSV the benchmark harness reads.
  DEFINITION_YAML     the code list is the dataset definition's ``targets:`` block.
  EXTRACTED_FROM_CODE the code list was lifted out of the C++ literal because no
                      definition file held it. This is CIRCULAR provenance: the
                      source proves only what the code already said. It is labelled
                      so, and scripts/check_dataset_identity.py reports it UNVERIFIED.

The astex_diverse entry is deliberately pointed at the roster CSV that is LIVE
(72/85 of the published set), not at the corrected ``astex_diverse_hartshorn85.csv``.
Byte-identity with the previously hardcoded list is the acceptance criterion: a
codegen that silently changes the executed set is strictly worse than the
hardcoding it replaces. Switching rosters is a separate, deliberate commit --
change ``source`` here and the generated file, the receipt validator and the
identity gate all follow.

Usage
-----
    python3 scripts/gen_dataset_codes.py            # regenerate in place
    python3 scripts/gen_dataset_codes.py --check     # fail if committed != generated
    python3 scripts/gen_dataset_codes.py --print-json  # machine-readable registry
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple

ROOT = Path(__file__).resolve().parents[1]
OUT_REL = "LIB/generated/dataset_codes.inc"
BANNER_MARK = "// GENERATED -- DO NOT EDIT, run scripts/gen_dataset_codes.py"

# A PDB id: digit followed by three alphanumerics. Anchored -- an unanchored
# version matches years, page numbers and four-digit potencies.
_PDB_CODE = re.compile(r"^[0-9][A-Za-z][A-Za-z0-9]{2}$")


class CodeList(NamedTuple):
    """One generated code list and the single file it is generated from."""
    slug: str
    cpp_function: str          # the DatasetRunner method this list backs
    source: str                # repo-relative canonical source path
    reader: str                # "csv:<column>" | "yaml:<key>"
    expected_n: int            # tripwire: a reader that silently returns fewer
    provenance: str            # ROSTER_CSV | DEFINITION_YAML | EXTRACTED_FROM_CODE
    primary_source: str        # DOI of the publication, or "UNVERIFIED"
    note: str


REGISTRY: List[CodeList] = [
    CodeList(
        slug="astex_diverse",
        cpp_function="DatasetRunner::astex_diverse_codes",
        source="benchmarks/astex_diverse/astex_diverse_set.csv",
        reader="csv:pdb_id",
        expected_n=85,
        provenance="ROSTER_CSV",
        primary_source="10.1021/jm061277y",
        note=(
            "LIVE roster. It is 72/85 of the Hartshorn 2007 published set with 13 "
            "substitutions (dataset identity audit 2026-09-10). Kept as the source so "
            "this codegen is behaviour-preserving; the corrected roster is "
            "benchmarks/astex_diverse/astex_diverse_hartshorn85.csv and switching to it "
            "is a separate commit."
        ),
    ),
    CodeList(
        slug="hap2",
        cpp_function="DatasetRunner::hap2_codes",
        source="benchmarks/datasets/hap2.yaml",
        reader="yaml:targets",
        expected_n=59,
        provenance="DEFINITION_YAML",
        primary_source="10.1021/acs.jcim.5b00078",
        note="The definition YAML's targets: block already matched the C++ list exactly, in order.",
    ),
    CodeList(
        slug="casf2016",
        cpp_function="DatasetRunner::casf2016_codes",
        source="benchmarks/datasets/casf2016_core_codes.csv",
        reader="csv:pdb_id",
        expected_n=283,
        provenance="EXTRACTED_FROM_CODE",
        primary_source="UNVERIFIED",
        note=(
            "CIRCULAR: no definition file held these codes. Both copies of casf2016.yaml "
            "declare a 29-code tier-1 subset, not this list, so the source CSV was "
            "extracted from the C++ literal itself and proves only what the code already "
            "said. Separately: the audit found only 115 of these overlap the published "
            "CASF-2016 core set -- that defect is NOT addressed by this codegen."
        ),
    ),
]


# ---------------------------------------------------------------------------
# Reading canonical sources
# ---------------------------------------------------------------------------

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def read_codes(entry: CodeList) -> List[str]:
    """Read one canonical source. Raises on anything ambiguous -- never returns
    a short list quietly, because a reader that silently drops rows produces a
    generated file that looks authoritative and runs the wrong set."""
    path = ROOT / entry.source
    if not path.is_file():
        raise SystemExit("MISSING canonical source for %s: %s" % (entry.slug, entry.source))

    kind, _, key = entry.reader.partition(":")
    if kind == "csv":
        with path.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        if not rows:
            raise SystemExit("EMPTY canonical source for %s: %s" % (entry.slug, entry.source))
        if key not in rows[0]:
            raise SystemExit("column %r absent from %s (has %s)"
                             % (key, entry.source, sorted(k for k in rows[0] if k)))
        raw = [(r.get(key) or "").strip() for r in rows]
    elif kind == "yaml":
        import yaml
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        if not isinstance(doc, dict) or not isinstance(doc.get(key), list):
            raise SystemExit("%s: %r is not a list in %s" % (entry.slug, key, entry.source))
        raw = [str(x).strip() for x in doc[key]]
    else:
        raise SystemExit("unknown reader %r for %s" % (entry.reader, entry.slug))

    bad = [v for v in raw if not _PDB_CODE.match(v)]
    if bad:
        raise SystemExit("%s: %d value(s) in %s are not PDB ids: %s"
                         % (entry.slug, len(bad), entry.source, bad[:5]))
    codes = [v.upper() for v in raw]
    dupes = sorted({c for c in codes if codes.count(c) > 1})
    if dupes:
        raise SystemExit("%s: duplicate codes in %s: %s" % (entry.slug, entry.source, dupes))
    if len(codes) != entry.expected_n:
        raise SystemExit("%s: %s yielded %d codes, registry declares %d. If the roster "
                         "changed on purpose, update expected_n in the same commit."
                         % (entry.slug, entry.source, len(codes), entry.expected_n))
    return codes


# ---------------------------------------------------------------------------
# Emitting
# ---------------------------------------------------------------------------

def _c_str(s: str) -> str:
    return '"%s"' % s.replace("\\", "\\\\").replace('"', '\\"')


def render() -> str:
    lists: Dict[str, List[str]] = {}
    shas: Dict[str, str] = {}
    for e in REGISTRY:
        lists[e.slug] = read_codes(e)
        shas[e.slug] = sha256_file(ROOT / e.source)

    L: List[str] = []
    a = L.append
    a(BANNER_MARK)
    a("//")
    a("// Every PDB code list the engine executes, generated from ONE declared canonical")
    a("// source per dataset. Hand-editing this file is what the generator exists to stop:")
    a("// scripts/gen_dataset_codes.py --check regenerates and diffs, and fails with both")
    a("// sha256 values if the committed file disagrees with its source.")
    a("//")
    a("// Sources (path, sha256 of the source AT GENERATION TIME, provenance):")
    for e in REGISTRY:
        a("//   %-14s %s" % (e.slug, e.source))
        a("//   %-14s sha256=%s" % ("", shas[e.slug]))
        a("//   %-14s n=%d provenance=%s primary_source=%s"
          % ("", len(lists[e.slug]), e.provenance, e.primary_source))
    a("//")
    a("// provenance=EXTRACTED_FROM_CODE means the source was lifted out of the C++ literal")
    a("// because no definition file held the list. That is CIRCULAR: it proves only what")
    a("// the code already said. scripts/check_dataset_identity.py reports it UNVERIFIED.")
    a("// =============================================================================")
    a("")
    a("#ifndef FLEXAIDDS_GENERATED_DATASET_CODES_INC")
    a("#define FLEXAIDDS_GENERATED_DATASET_CODES_INC")
    a("")
    a("#include <cstddef>")
    a("#include <stdexcept>")
    a("#include <string>")
    a("#include <vector>")
    a("")
    a("namespace flexaids {")
    a("namespace generated {")
    a("")

    for e in REGISTRY:
        codes = lists[e.slug]
        a("// %s -- %s (%s)" % (e.slug, e.source, e.provenance))
        a("inline const std::vector<std::string>& %s_codes() {" % e.slug)
        a("    static const std::vector<std::string> v = {")
        for i in range(0, len(codes), 8):
            a("        " + ", ".join(_c_str(c) for c in codes[i:i + 8]) + ",")
        a("    };")
        a("    return v;")
        a("}")
        a("")

    a("/// What a generated list claims to implement: the canonical source it was read")
    a("/// from, that source's sha256 at generation time, and how trustworthy the")
    a("/// provenance is. Recorded into RUN_RECEIPT.json so a run's executed roster can")
    a("/// be checked against its declared source instead of taken on faith.")
    a("struct DatasetCodeSource {")
    a("    const char* slug;")
    a("    const char* cpp_function;")
    a("    const char* source_path;")
    a("    const char* source_sha256;")
    a("    const char* provenance;")
    a("    const char* primary_source;")
    a("    std::size_t count;")
    a("};")
    a("")
    a("inline const std::vector<DatasetCodeSource>& dataset_code_sources() {")
    a("    static const std::vector<DatasetCodeSource> v = {")
    for e in REGISTRY:
        a("        {%s, %s," % (_c_str(e.slug), _c_str(e.cpp_function)))
        a("         %s," % _c_str(e.source))
        a("         %s," % _c_str(shas[e.slug]))
        a("         %s, %s, %d}," % (_c_str(e.provenance), _c_str(e.primary_source),
                                     len(lists[e.slug])))
    a("    };")
    a("    return v;")
    a("}")
    a("")
    a("/// Codes for a dataset slug. THROWS on an unknown slug: an engine that cannot")
    a("/// name the declaration it is executing must not execute a list at all.")
    a("inline const std::vector<std::string>& dataset_codes(const std::string& slug) {")
    for e in REGISTRY:
        a("    if (slug == %s) return %s_codes();" % (_c_str(e.slug), e.slug))
    a("    throw std::runtime_error(\"dataset_codes: no generated code list declared for slug '\"")
    a("                             + slug + \"' (see scripts/gen_dataset_codes.py REGISTRY)\");")
    a("}")
    a("")
    a("/// Source record for a dataset slug. THROWS on an unknown slug, same rationale.")
    a("inline const DatasetCodeSource& dataset_code_source(const std::string& slug) {")
    a("    for (const auto& s : dataset_code_sources()) {")
    a("        if (slug == s.slug) return s;")
    a("    }")
    a("    throw std::runtime_error(\"dataset_code_source: no generated code list declared for \"")
    a("                             \"slug '\" + slug + \"'\");")
    a("}")
    a("")
    a("}  // namespace generated")
    a("}  // namespace flexaids")
    a("")
    a("#endif  // FLEXAIDDS_GENERATED_DATASET_CODES_INC")
    return "\n".join(L) + "\n"


def registry_json() -> dict:
    out = {"output": OUT_REL, "lists": []}
    for e in REGISTRY:
        codes = read_codes(e)
        out["lists"].append({
            "slug": e.slug,
            "cpp_function": e.cpp_function,
            "source": e.source,
            "source_sha256": sha256_file(ROOT / e.source),
            "reader": e.reader,
            "provenance": e.provenance,
            "primary_source": e.primary_source,
            "count": len(codes),
            "codes_sha256": hashlib.sha256(("\n".join(codes) + "\n").encode()).hexdigest(),
            "codes": codes,
            "note": e.note,
        })
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="do not write; fail if the committed file differs from generated")
    ap.add_argument("--print-json", action="store_true",
                    help="print the registry with codes and shas as JSON")
    args = ap.parse_args(argv)

    if args.print_json:
        json.dump(registry_json(), sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return 0

    want = render()
    out = ROOT / OUT_REL
    want_sha = hashlib.sha256(want.encode()).hexdigest()

    if args.check:
        if not out.is_file():
            print("FAIL - %s does not exist. Run scripts/gen_dataset_codes.py." % OUT_REL)
            return 1
        have = out.read_text(encoding="utf-8")
        have_sha = hashlib.sha256(have.encode()).hexdigest()
        if have == want:
            n = sum(s["count"] for s in registry_json()["lists"])
            print("PASS - %s matches its canonical sources (%d lists, %d codes, sha256=%s)"
                  % (OUT_REL, len(REGISTRY), n, want_sha))
            return 0
        print("FAIL - %s does not match what its canonical sources produce." % OUT_REL)
        print("  committed sha256: %s" % have_sha)
        print("  generated sha256: %s" % want_sha)
        import difflib
        diff = list(difflib.unified_diff(have.splitlines(), want.splitlines(),
                                         "committed", "generated", lineterm="", n=1))
        for line in diff[:40]:
            print("  %s" % line)
        if len(diff) > 40:
            print("  ... %d more diff lines" % (len(diff) - 40))
        print("  Run scripts/gen_dataset_codes.py to regenerate, and review the diff: a "
              "change here changes the set the engine executes.")
        return 1

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(want, encoding="utf-8")
    print("wrote %s (%d lists, %d codes, sha256=%s)"
          % (OUT_REL, len(REGISTRY), sum(len(read_codes(e)) for e in REGISTRY), want_sha))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
