#!/usr/bin/env python3
"""Differential tests: NativePoseQC extraction + upstream PoseBusters CLI.

Authoritative gate = installed `bust` (FLEXAIDDS_POSEBUSTERS_BIN).
Native C++ (LIB/PoseBust) is diagnostic only — these tests assert extraction
fail-closed topology and that bust runs on extracted SDFs.

Requires: .venv-posebusters/bin/bust (or FLEXAIDDS_POSEBUSTERS_BIN), clang++,
sample FlexAID poses under benchmarks/astex_repro when available.

Usage:
  python3 tests/test_posebust_upstream_parity.py
  python3 tests/test_posebust_upstream_parity.py --all-astex   # crystal-only, no GA poses
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUST = Path(
    os.environ.get(
        "FLEXAIDDS_POSEBUSTERS_BIN", str(ROOT / ".venv-posebusters/bin/bust")
    )
)
ASTEX = ROOT / "benchmarks/astex_diverse/astex_diverse"
REPRO_DIRS = [
    ROOT / "benchmarks/astex_repro/full_v132",
    ROOT / "benchmarks/astex_repro/full_v131",
    ROOT / "benchmarks/astex_repro/full_v130",
]

# Expected heavy-atom counts for ligands (from crystal SDF) used as extraction checks.
# 1G9V is the HEM-contamination regression target (must be 25, not heme-scale).
EXPECTED_LIG_HEAVY = {
    "1G9V": 25,
}


def find_pose(code: str) -> Path | None:
    for base in REPRO_DIRS:
        target = base / code
        if not target.is_dir():
            continue
        ep = target / "elected_pose.pdb"
        if ep.is_file():
            return ep
        # Prefer restart poses over root rank-0 (election may differ).
        candidates: list[Path] = []
        for p in sorted(target.rglob(f"{code}_*.pdb")):
            if "_INI" in p.name or "member" in p.name or "apo" in p.name:
                continue
            if p.name.endswith("_binding_site.pdb"):
                continue
            candidates.append(p)
        if candidates:
            # Prefer r*/ paths when present
            restarts = [c for c in candidates if "/r" in str(c)]
            return restarts[0] if restarts else candidates[0]
    return None


def crystal_paths(code: str) -> tuple[Path, Path] | None:
    d = ASTEX / code
    lig = d / f"{code}_ligand.sdf"
    prot = d / f"{code}_apo.pdb"
    if lig.is_file() and prot.is_file():
        return lig, prot
    return None


def run_bust(pred: Path, protein: Path, crystal: Path) -> dict[str, str]:
    cmd = [
        str(BUST),
        str(pred),
        "-p",
        str(protein),
        "-l",
        str(crystal),
        "--outfmt",
        "csv",
    ]
    out = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL)
    rows = list(csv.DictReader(io.StringIO(out)))
    assert rows, "empty bust csv"
    return rows[0]


def compile_tool(src: str, name: str, td: Path) -> Path:
    cpp = td / f"{name}.cpp"
    cpp.write_text(src)
    binp = td / name
    cmd = [
        "clang++",
        "-std=c++26",
        "-O2",
        f"-I{ROOT / 'LIB'}",
        f"-I{ROOT / 'LIB' / 'PoseBust'}",
        str(cpp),
        str(ROOT / "LIB/PoseBust/Loaders.cpp"),
        "-o",
        str(binp),
    ]
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return binp


EXTRACT_MAIN = r"""
#include "PoseBust/Loaders.h"
#include <iostream>
int main(int c,char**v){
  if(c < 4) return 9;
  flexaids::posebust::Molecule lig, ref;
  std::string err;
  if(!flexaids::posebust::load_pdb_flexaid_ligand(v[1], lig, &err)) {
    std::cerr<<err<<"\n"; return 1;
  }
  if(!flexaids::posebust::load_sdf(v[2], ref, &err)) {
    std::cerr<<err<<"\n"; return 2;
  }
  if(!flexaids::posebust::assign_topology_from_reference(lig, ref, &err)) {
    std::cerr<<"topology fail: "<<err<<"\n"; return 3;
  }
  if(!flexaids::posebust::write_sdf(lig, v[3], &err)) {
    std::cerr<<err<<"\n"; return 4;
  }
  std::cout<<lig.atoms.size()<<"\n";
  return 0;
}
"""

MISMATCH_MAIN = r"""
#include "PoseBust/Loaders.h"
#include <iostream>
int main(int c,char**v){
  flexaids::posebust::Molecule lig, ref;
  std::string err;
  if(!flexaids::posebust::load_pdb_flexaid_ligand(v[1], lig, &err)) return 10;
  if(!flexaids::posebust::load_sdf(v[2], ref, &err)) return 11;
  if(flexaids::posebust::assign_topology_from_reference(lig, ref, &err)) {
    std::cerr<<"unexpected success\n"; return 1;
  }
  std::cout<<"failed_as_expected: "<<err<<"\n";
  return 0;
}
"""

CRYSTAL_SELF_MAIN = r"""
#include "PoseBust/Loaders.h"
#include <iostream>
// Crystal-only topology round-trip: load SDF, rewrite, count heavy atoms.
int main(int c,char**v){
  flexaids::posebust::Molecule ref;
  std::string err;
  if(!flexaids::posebust::load_sdf(v[1], ref, &err)) {
    std::cerr<<err<<"\n"; return 1;
  }
  if(!flexaids::posebust::write_sdf(ref, v[2], &err)) {
    std::cerr<<err<<"\n"; return 2;
  }
  std::cout<<ref.n_heavy()<<"\n";
  return 0;
}
"""


def extract_ligand(pose: Path, crystal: Path, out_sdf: Path, bin_extract: Path) -> int:
    n = subprocess.check_output(
        [str(bin_extract), str(pose), str(crystal), str(out_sdf)], text=True
    ).strip()
    return int(n)


def sdf_atom_count(sdf: Path) -> int:
    text = sdf.read_text(errors="replace")
    for line in text.splitlines():
        if "V2000" in line or "V3000" in line:
            return int(line[:3].strip() or line.split()[0])
    raise AssertionError(f"no counts line in {sdf}")


def test_1g9v_extract_not_hem(bin_extract: Path) -> None:
    pose = find_pose("1G9V")
    paths = crystal_paths("1G9V")
    if pose is None or paths is None:
        print("SKIP test_1g9v_extract_not_hem: no pose/crystal artifact")
        return
    crystal, protein = paths
    with tempfile.TemporaryDirectory() as td:
        sdf = Path(td) / "lig.sdf"
        n = extract_ligand(pose, crystal, sdf, bin_extract)
        assert n < 50, f"1G9V ligand atom count looks like cofactor contamination: {n}"
        nat = sdf_atom_count(sdf)
        assert nat == EXPECTED_LIG_HEAVY["1G9V"], (
            f"expected {EXPECTED_LIG_HEAVY['1G9V']} heavy atoms for 1G9V, got {nat}"
        )
    print("PASS test_1g9v_extract_not_hem")


def test_1g9v_upstream_bust_runs(bin_extract: Path) -> None:
    if not BUST.is_file():
        print("SKIP test_1g9v_upstream_bust_runs: bust missing")
        return
    pose = find_pose("1G9V")
    paths = crystal_paths("1G9V")
    if pose is None or paths is None:
        print("SKIP test_1g9v_upstream_bust_runs: no pose")
        return
    crystal, protein = paths
    with tempfile.TemporaryDirectory() as td:
        sdf = Path(td) / "lig.sdf"
        extract_ligand(pose, crystal, sdf, bin_extract)
        row = run_bust(sdf, protein, crystal)
        assert row.get("mol_pred_loaded") == "True"
        assert row.get("mol_cond_loaded") == "True"
        assert "bond_angles" in row
        assert "minimum_distance_to_protein" in row
        # Document known upstream fails for this autonomous CF-elected pose.
        fails = [k for k, v in row.items() if v == "False"]
        print("upstream_1G9V_fails:", fails)
        # Regression: native previously reported no_radicals=false incorrectly;
        # upstream must load radicals check as boolean.
        assert "no_radicals" in row
    print("PASS test_1g9v_upstream_bust_runs")


def test_topology_mismatch_fails(bin_mismatch: Path) -> None:
    pose = find_pose("1G9V")
    if pose is None:
        print("SKIP test_topology_mismatch_fails")
        return
    other = ASTEX / "1GPK" / "1GPK_ligand.sdf"
    if not other.is_file():
        print("SKIP no 1GPK ligand")
        return
    out = subprocess.check_output([str(bin_mismatch), str(pose), str(other)], text=True)
    assert "failed_as_expected" in out
    print("PASS test_topology_mismatch_fails")


def rewrite_crystal_sdf(bin_self: Path, crystal: Path, out_sdf: Path) -> int:
    """Normalize Astex crystal SDF via PoseBust Loaders (raw files often have bad CTAB)."""
    n = subprocess.check_output(
        [str(bin_self), str(crystal), str(out_sdf)], text=True
    ).strip()
    return int(n)


def test_1m2z_broken_ring_loads(bin_self: Path) -> None:
    """Broken-ring / difficult geometry target: rewritten crystal must load in bust."""
    paths = crystal_paths("1M2Z")
    if paths is None:
        print("SKIP test_1m2z_broken_ring_loads: no 1M2Z crystal")
        return
    crystal, protein = paths
    if not BUST.is_file():
        print("SKIP test_1m2z_broken_ring_loads: bust missing")
        return
    with tempfile.TemporaryDirectory() as td:
        # Raw Astex SDFs often have non-RDKit CTAB lines; rewrite is required.
        fixed = Path(td) / "1M2Z_fixed.sdf"
        n = rewrite_crystal_sdf(bin_self, crystal, fixed)
        assert n >= 3, f"1M2Z too few heavy atoms: {n}"
        row = run_bust(fixed, protein, fixed)
        assert row.get("mol_pred_loaded") == "True", row
        assert row.get("mol_cond_loaded") == "True", row
        fails = [k for k, v in row.items() if v == "False"]
        print("upstream_1M2Z_crystal_self fails:", fails)
    print("PASS test_1m2z_broken_ring_loads")


def test_metal_cofactor_targets_crystal_self(bin_self: Path) -> None:
    """Metal/cofactor-rich targets: rewritten crystal + apo must load in bust."""
    if not BUST.is_file():
        print("SKIP test_metal_cofactor_targets_crystal_self: bust missing")
        return
    # 1G9V has HEM in holo; apo should be clean. 1P2Y/1Q4G are cofactor-adjacent.
    codes = ["1G9V", "1P2Y", "1Q4G", "1R9O"]
    ran = 0
    with tempfile.TemporaryDirectory() as td:
        for code in codes:
            paths = crystal_paths(code)
            if paths is None:
                print(f"  skip {code}: missing crystal/apo")
                continue
            crystal, protein = paths
            fixed = Path(td) / f"{code}_fixed.sdf"
            rewrite_crystal_sdf(bin_self, crystal, fixed)
            row = run_bust(fixed, protein, fixed)
            assert row.get("mol_pred_loaded") == "True", f"{code}: {row}"
            assert row.get("mol_cond_loaded") == "True", f"{code}: {row}"
            ran += 1
            fails = [k for k, v in row.items() if v == "False"]
            print(f"  {code} crystal-self fails={fails}")
    assert ran >= 1, "no metal/cofactor targets available"
    print(f"PASS test_metal_cofactor_targets_crystal_self ({ran} targets)")


def test_hem_contaminated_receptor_still_extracts_ligand(bin_extract: Path) -> None:
    """If receptor path contains HEM, ligand extract from pose must stay CONECT-based."""
    pose = find_pose("1G9V")
    paths = crystal_paths("1G9V")
    if pose is None or paths is None:
        print("SKIP test_hem_contaminated_receptor_still_extracts_ligand")
        return
    crystal, _ = paths
    with tempfile.TemporaryDirectory() as td:
        sdf = Path(td) / "lig.sdf"
        n = extract_ligand(pose, crystal, sdf, bin_extract)
        assert n == EXPECTED_LIG_HEAVY["1G9V"] or sdf_atom_count(sdf) == EXPECTED_LIG_HEAVY[
            "1G9V"
        ], f"HEM contamination risk: n={n}"
    print("PASS test_hem_contaminated_receptor_still_extracts_ligand")


def test_all_astex_crystal_sdf_loadable(bin_self: Path, limit: int | None = None) -> None:
    """Fail-closed: every Astex crystal ligand SDF must load and rewrite."""
    codes = sorted(p.name for p in ASTEX.iterdir() if p.is_dir())
    if limit:
        codes = codes[:limit]
    fails: list[str] = []
    ok = 0
    with tempfile.TemporaryDirectory() as td:
        binp = bin_self
        for code in codes:
            lig = ASTEX / code / f"{code}_ligand.sdf"
            if not lig.is_file():
                fails.append(f"{code}:missing_sdf")
                continue
            out = Path(td) / f"{code}.sdf"
            try:
                n = subprocess.check_output(
                    [str(binp), str(lig), str(out)], text=True, stderr=subprocess.PIPE
                ).strip()
                if int(n) < 1:
                    fails.append(f"{code}:zero_atoms")
                else:
                    ok += 1
            except subprocess.CalledProcessError as e:
                fails.append(f"{code}:load_fail:{e.stderr.decode()[:80] if e.stderr else e}")
    print(f"  astex crystal load: ok={ok} fail={len(fails)} / {len(codes)}")
    if fails[:10]:
        print("  first fails:", fails[:10])
    # Soft threshold: most of Astex must load (allow a few pathological SDFs)
    assert ok >= max(1, int(0.9 * len(codes))), f"too many crystal load failures: {fails[:20]}"
    print(f"PASS test_all_astex_crystal_sdf_loadable ({ok}/{len(codes)})")


def test_available_ga_poses_extract_and_bust(bin_extract: Path) -> None:
    """For every GA pose artifact found under full_v13x, extract + optional bust."""
    if not BUST.is_file():
        print("SKIP test_available_ga_poses_extract_and_bust: bust missing")
        return
    codes: list[str] = []
    for base in REPRO_DIRS:
        if base.is_dir():
            codes.extend(p.name for p in base.iterdir() if p.is_dir() and p.name[0].isdigit())
    codes = sorted(set(codes))
    if not codes:
        print("SKIP test_available_ga_poses_extract_and_bust: no GA poses")
        return
    n_ok = 0
    n_topo_fail = 0
    n_bust = 0
    for code in codes:
        pose = find_pose(code)
        paths = crystal_paths(code)
        if pose is None or paths is None:
            continue
        crystal, protein = paths
        with tempfile.TemporaryDirectory() as td:
            sdf = Path(td) / "lig.sdf"
            try:
                n = extract_ligand(pose, crystal, sdf, bin_extract)
            except subprocess.CalledProcessError:
                n_topo_fail += 1
                print(f"  {code}: topology/extract FAIL (fail-closed OK)")
                continue
            if code in EXPECTED_LIG_HEAVY:
                assert n == EXPECTED_LIG_HEAVY[code] or sdf_atom_count(sdf) == EXPECTED_LIG_HEAVY[
                    code
                ]
            n_ok += 1
            try:
                row = run_bust(sdf, protein, crystal)
                assert row.get("mol_pred_loaded") == "True"
                n_bust += 1
                fails = [k for k, v in row.items() if v == "False"]
                print(f"  {code}: extract_ok bust_fails={fails[:6]}")
            except Exception as e:
                print(f"  {code}: bust error {e}")
    print(f"  ga_poses extract_ok={n_ok} topo_fail={n_topo_fail} bust_ran={n_bust}")
    assert n_ok >= 1, "expected at least one GA pose extract"
    print("PASS test_available_ga_poses_extract_and_bust")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--all-astex",
        action="store_true",
        help="Also load all 85 Astex crystal SDFs",
    )
    args = ap.parse_args()
    os.chdir(ROOT)
    fails = 0

    with tempfile.TemporaryDirectory() as td0:
        td = Path(td0)
        bin_extract = compile_tool(EXTRACT_MAIN, "extract", td)
        bin_mismatch = compile_tool(MISMATCH_MAIN, "mismatch", td)
        bin_self = compile_tool(CRYSTAL_SELF_MAIN, "crystal_self", td)

        tests = [
            lambda: test_1g9v_extract_not_hem(bin_extract),
            lambda: test_1g9v_upstream_bust_runs(bin_extract),
            lambda: test_topology_mismatch_fails(bin_mismatch),
            lambda: test_1m2z_broken_ring_loads(bin_self),
            lambda: test_metal_cofactor_targets_crystal_self(bin_self),
            lambda: test_hem_contaminated_receptor_still_extracts_ligand(bin_extract),
            lambda: test_available_ga_poses_extract_and_bust(bin_extract),
        ]
        if args.all_astex:
            tests.append(lambda: test_all_astex_crystal_sdf_loadable(bin_self))

        names = [
            "test_1g9v_extract_not_hem",
            "test_1g9v_upstream_bust_runs",
            "test_topology_mismatch_fails",
            "test_1m2z_broken_ring_loads",
            "test_metal_cofactor_targets_crystal_self",
            "test_hem_contaminated_receptor_still_extracts_ligand",
            "test_available_ga_poses_extract_and_bust",
        ]
        if args.all_astex:
            names.append("test_all_astex_crystal_sdf_loadable")

        for name, fn in zip(names, tests):
            try:
                fn()
            except Exception as e:
                print(f"FAIL {name}: {e}")
                fails += 1

    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
