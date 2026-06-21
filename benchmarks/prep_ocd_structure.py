#!/usr/bin/env python3
"""
prep_ocd_structure.py — Download and prep a PDB complex for OCD mini-benchmark.

Creates per-structure directory with:
  {ID}.pdb               full RCSB PDB (complex with ligand)
  {ID}_apo.pdb           protein ATOM records only (no HETATM, no HOH)
  {ID}_ligand.sdf        main drug-like ligand as SDF (crystal geometry)
  {ID}_binding_site.pdb  Get_Cleft oracle cavity anchored on crystal ligand

Usage:
  python3 prep_ocd_structure.py <PDBID> <outdir> [--lig RESNAME] [--get-cleft /path/to/Get_Cleft]

Examples:
  python3 prep_ocd_structure.py 1DI8 /path/to/ocd_mini/cdk2
  python3 prep_ocd_structure.py 1AI8 /path/to/ocd_mini/thrombin --lig 4HB
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Residues to exclude when hunting for the main drug-like ligand
# ---------------------------------------------------------------------------
EXCLUDE_RESNAMES: frozenset[str] = frozenset({
    # Water / deuterium water
    "HOH", "WAT", "DOD", "H2O",
    # Common solvents / crystallization additives
    "DMS", "GOL", "EDO", "MPD", "ACT", "ACE", "ACY", "IPA", "MOH",
    "PEG", "1PE", "2PE", "P6G", "PG4", "PG6", "PGE", "15P",
    "TRS", "MES", "HEP", "BIS", "EPE", "BME", "BTB",
    "IMD", "TAM", "TAR", "CIT", "SUC", "FMT", "OXL", "MLI", "SCN",
    "PCA", "NH4", "CME", "CSO", "OXE", "FLC", "LI", "SAR",
    "SPM", "SPK", "SPD",   # polyamines
    # Ions / metals (single-atom HETATM)
    "MG", "ZN", "CA", "NA", "K", "CL", "MN", "FE", "CU", "NI", "CO",
    "BR", "IOD", "F", "SE", "PT", "AU", "HG", "CD", "BA", "RB", "CS", "SR",
    # Phosphate / sulfate
    "PO4", "SO4", "ACT",
    # Nucleotides (common cofactors; too large/flexible for typical drug-like OCD)
    # Comment out if you WANT ATP analogs as the ligand (e.g. for kinase APO-site test)
    # "ATP", "ADP", "AMP", "ANP", "ACP", "GNP", "GTP", "GDP",
    # Glycans (PTMs)
    "NAG", "FUC", "MAN", "GAL", "GLC", "BGC", "SIA",
    # Fatty acids / lipids
    "OLA", "PLM", "STE",
})

GET_CLEFT_DEFAULT = Path("/Users/lp.more/Projects/Get_Cleft/Get_Cleft")
RCSB_URL = "https://files.rcsb.org/download/{pdb}.pdb"


def download_pdb(pdb_id: str, dest: Path) -> Path:
    url = RCSB_URL.format(pdb=pdb_id.upper())
    out = dest / f"{pdb_id.upper()}.pdb"
    print(f"  Downloading {url} ...", end=" ", flush=True)
    urllib.request.urlretrieve(url, out)
    print(f"({out.stat().st_size // 1024} kB)")
    return out


def parse_hetatm_groups(pdb_path: Path) -> dict[tuple[str, str, int], list[str]]:
    """Return {(resname, chain, resnum): [lines]} for all HETATM residue groups."""
    groups: dict[tuple[str, str, int], list[str]] = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            chain   = line[21].strip() or " "
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            key = (resname, chain, resnum)
            groups.setdefault(key, []).append(line)
    return groups


def pick_ligand(groups: dict[tuple[str, str, int], list[str]],
                force_resname: str | None) -> tuple[tuple[str, str, int], list[str]] | None:
    """
    Select the main drug-like ligand group.
    - If force_resname given, use that.
    - Otherwise: exclude known solvents/ions, pick the group with the most heavy atoms.
    """
    if force_resname:
        hits = {k: v for k, v in groups.items() if k[0] == force_resname.upper()}
        if not hits:
            return None
        return max(hits.items(), key=lambda kv: len(kv[1]))

    candidates = {k: v for k, v in groups.items() if k[0] not in EXCLUDE_RESNAMES}
    if not candidates:
        return None
    return max(candidates.items(), key=lambda kv: len(kv[1]))


def write_apo(pdb_path: Path, out_path: Path) -> int:
    """Write protein ATOM records (no HETATM, no HOH) to out_path. Returns atom count."""
    count = 0
    with open(pdb_path) as src, open(out_path, "w") as dst:
        for line in src:
            if line.startswith("ATOM"):
                dst.write(line)
                count += 1
            elif line.startswith(("TER", "END")):
                dst.write(line)
    return count


def write_ligand_pdb(lines: list[str], out_path: Path) -> None:
    """Write ligand HETATM lines to a temp PDB fragment."""
    with open(out_path, "w") as f:
        for line in lines:
            f.write(line)
        f.write("END\n")


def convert_to_sdf(lig_pdb: Path, out_sdf: Path) -> bool:
    """Use obabel to convert a HETATM-fragment PDB → SDF."""
    cmd = ["obabel", str(lig_pdb), "-O", str(out_sdf), "--gen3D", "-h", "--perceive-charges"]
    # First try with --gen3D (may fail if already 3D); fall back without
    for attempt_cmd in [
        ["obabel", str(lig_pdb), "-O", str(out_sdf)],
    ]:
        result = subprocess.run(attempt_cmd, capture_output=True, text=True)
        if out_sdf.exists() and out_sdf.stat().st_size > 0:
            print(f"  SDF written: {out_sdf.name} ({out_sdf.stat().st_size} bytes)")
            return True
        print(f"  obabel attempt failed: {result.stderr.strip()[:200]}")
    return False


def build_anchor_string(resname: str, resnum: int, chain: str) -> str:
    """Format Get_Cleft -a argument."""
    if len(resname) < 3:
        resname = "-" * (3 - len(resname)) + resname
    chain_char = chain if (chain and chain.strip()) else "-"
    return f"{resname}{resnum}{chain_char}-"


def run_get_cleft(get_cleft: Path, pdb_path: Path, anchor: str, work_dir: Path) -> Path | None:
    pdb_base = pdb_path.name[:4]
    cmd = [str(get_cleft), "-p", str(pdb_path), "-o", pdb_base, "-a", anchor]
    print(f"  Get_Cleft: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR Get_Cleft rc={result.returncode}: {result.stderr[:300]}")
        return None
    matches = sorted(work_dir.glob("*_clf_*.pdb"))
    if not matches:
        print(f"  ERROR: no *_clf_*.pdb in {work_dir}; stdout={result.stdout[:200]}")
        return None
    return matches[0]


def prep_one(pdb_id: str, outdir: Path, get_cleft: Path,
             force_resname: str | None = None) -> dict | None:
    """
    Full prep pipeline for one PDB ID.
    Returns dict with file paths on success, None on failure.
    """
    pdb_id = pdb_id.upper()
    dest = outdir / pdb_id
    dest.mkdir(parents=True, exist_ok=True)

    pdb_path  = dest / f"{pdb_id}.pdb"
    apo_path  = dest / f"{pdb_id}_apo.pdb"
    lig_path  = dest / f"{pdb_id}_ligand.sdf"
    site_path = dest / f"{pdb_id}_binding_site.pdb"

    # Skip if fully prepped
    if pdb_path.exists() and apo_path.exists() and lig_path.exists() and site_path.exists():
        print(f"[{pdb_id}] Already prepped — skip")
        return {"pdb_id": pdb_id, "pdb": str(pdb_path), "apo": str(apo_path),
                "ligand_sdf": str(lig_path), "binding_site": str(site_path)}

    # 1. Download
    if not pdb_path.exists():
        try:
            download_pdb(pdb_id, dest)
        except Exception as e:
            print(f"[{pdb_id}] DOWNLOAD FAILED: {e}")
            return None

    # 2. Parse HETATM groups
    groups = parse_hetatm_groups(pdb_path)
    if not groups:
        print(f"[{pdb_id}] ERROR: no HETATM records found")
        return None

    # 3. Pick main ligand
    picked = pick_ligand(groups, force_resname)
    if picked is None:
        all_resnames = sorted({k[0] for k in groups})
        print(f"[{pdb_id}] ERROR: no drug-like ligand found. All HETATM: {all_resnames}")
        return None

    (resname, chain, resnum), lig_lines = picked
    print(f"[{pdb_id}] Ligand: {resname} chain={chain!r} resnum={resnum} ({len(lig_lines)} atoms)")

    # 4. Write apo PDB
    n_atoms = write_apo(pdb_path, apo_path)
    print(f"[{pdb_id}] Apo: {n_atoms} ATOM records → {apo_path.name}")

    # 5. Convert ligand → SDF
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        tmp_lig_pdb = Path(tmp.name)
    write_ligand_pdb(lig_lines, tmp_lig_pdb)
    if not convert_to_sdf(tmp_lig_pdb, lig_path):
        print(f"[{pdb_id}] ERROR: obabel SDF conversion failed")
        tmp_lig_pdb.unlink(missing_ok=True)
        return None
    tmp_lig_pdb.unlink(missing_ok=True)

    # 6. Run Get_Cleft for binding_site.pdb
    anchor = build_anchor_string(resname, resnum, chain)
    print(f"[{pdb_id}] Anchor string: {anchor}")
    with tempfile.TemporaryDirectory(prefix=f"gc_{pdb_id}_") as tmp:
        work = Path(tmp)
        import shutil
        shutil.copy2(pdb_path, work / pdb_path.name)
        clf_out = run_get_cleft(get_cleft, work / pdb_path.name, anchor, work)
        if clf_out is None:
            print(f"[{pdb_id}] ERROR: Get_Cleft failed")
            return None
        shutil.copy2(clf_out, site_path)
        print(f"[{pdb_id}] Binding site: {site_path.name} ({site_path.stat().st_size} bytes)")

    return {
        "pdb_id": pdb_id,
        "resname": resname,
        "chain": chain,
        "resnum": resnum,
        "pdb": str(pdb_path),
        "apo": str(apo_path),
        "ligand_sdf": str(lig_path),
        "binding_site": str(site_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pdb_id", help="4-letter PDB ID")
    parser.add_argument("outdir", type=Path, help="Output directory (e.g. benchmarks/ocd_mini/cdk2)")
    parser.add_argument("--lig", metavar="RESNAME", help="Force ligand residue name")
    parser.add_argument("--get-cleft", type=Path, default=GET_CLEFT_DEFAULT,
                        help=f"Path to Get_Cleft binary (default: {GET_CLEFT_DEFAULT})")
    args = parser.parse_args()

    if not args.get_cleft.is_file():
        print(f"ERROR: Get_Cleft not found: {args.get_cleft}", file=sys.stderr)
        return 1

    result = prep_one(args.pdb_id, args.outdir, args.get_cleft, force_resname=args.lig)
    if result is None:
        return 1

    import json
    print("\n=== RESULT ===")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
