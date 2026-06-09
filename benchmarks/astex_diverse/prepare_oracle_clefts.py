#!/usr/bin/env python3
"""
prepare_oracle_clefts.py — Oracle binding-site prep for all 85 Astex Diverse complexes.

For each complex:
  1. Download the original RCSB PDB (complex with crystal ligand).
  2. Find the crystal ligand's residue number + chain from HETATM records.
  3. Run Get_Cleft in oracle mode: -a <RESNAME><RESNUM><CHAIN>-
  4. Copy <base>_clf_1.pdb → astex_diverse/<PDB>/<PDB>_binding_site.pdb

Usage:
    python prepare_oracle_clefts.py --get-cleft /path/to/Get_Cleft [--dry-run] [--targets 1GM8 1GPK ...]
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Crystal-ligand residue names for all 85 Astex Diverse complexes.
# Source: DatasetRunner extraction from RCSB HETATM records (git history).
# These are the RCSB 3-letter residue codes for the docked ligand.
# ---------------------------------------------------------------------------
LIGAND_RESNAME: dict[str, str] = {
    "1G9V": "HEM",
    "1GM8": "SOX",
    "1GPK": "HUP",
    "1HNN": "SAH",
    "1HP0": "AD3",
    "1HQ2": "APC",
    "1IA1": "NDP",
    "1IGJ": "DGX",
    "1J3J": "NDP",
    "1JD0": "AZM",
    "1JJE": "BYS",
    "1K3U": "IAD",
    "1KE5": "LS1",
    "1KZK": "JE2",
    "1L2S": "STC",
    "1L7F": "BCZ",
    "1LPZ": "CMB",
    "1M2Z": "DEX",
    "1MEH": "IMP",
    "1MQ6": "XLD",
    "1N1M": "NAG",
    "1N2J": "PAF",
    "1N2V": "BDI",
    "1N46": "PFA",
    "1NAV": "IH5",
    "1OF1": "SCT",
    "1OF6": "DTY",
    "1OPK": "P16",
    "1OQ5": "CEL",
    "1OWE": "675",
    "1P2Y": "HEM",
    "1P62": "ADP",
    "1PMN": "984",
    "1Q1G": "MTI",
    "1Q41": "IXM",
    "1Q4G": "HEM",
    "1R1H": "BIR",
    "1R55": "097",
    "1R58": "AO5",
    "1R9O": "HEM",
    "1S19": "MC9",
    "1S3V": "TQD",
    "1SG0": "FAD",
    "1SJ0": "E4D",
    "1SQ5": "ADP",
    "1T40": "NAP",
    "1T46": "STI",
    "1T9B": "FAD",
    "1TT1": "KAI",
    # 1TW6 omitted from the HETATM map: its only HETATM are buffer/ion species
    # (BTB Bis-Tris, EDO, ZN, LI) — all blacklisted cofactors.  The cognate
    # ligand is the Smac-derived AVPI tetrapeptide bound in the IAP (IBM) groove,
    # stored as ATOM records (chain C).  See PEPTIDE_LIGAND below.  Anchoring the
    # cleft on BTB (the previous mapping) placed the oracle site ~39 Å away in a
    # crystallographic buffer pocket on the other BIR copy.
    "1TZ8": "DES",
    "1U1C": "BAU",
    "1U4D": "DBQ",
    "1UML": "FR4",
    "1UNL": "RRC",
    "1UOU": "CMU",
    "1V0P": "PVB",
    "1V48": "HA1",
    "1V4S": "MRK",
    "1VCJ": "IBA",
    "1W1P": "GIO",
    "1W2G": "THM",
    "1X8X": "TYR",
    "1XM6": "5RM",
    "1XOZ": "CIA",
    "1Y6B": "AAX",
    "1Y6R": "MTM",
    "1YGC": "905",
    "1YQY": "915",
    "1YV3": "ADP",
    "1YVF": "PH7",
    "1YWR": "LI9",
    "1Z95": "198",
    "2BM2": "PM2",
    "2BR1": "PFP",
    "2BSM": "BSM",
    "2BYS": "LOB",
    "2C3I": "IYZ",
    "2CET": "PGI",
    "2CGR": "GAS",
    "2D3U": "CCT",
    "2GBP": "BGC",
    "2HB1": "512",
    "2HR7": "P33",
    "2J62": "GSZ",
}

# Peptide ligands stored as ATOM records rather than HETATM small molecules.
# Maps PDB ID → (resname, resnum, chain) of an anchor residue inside the peptide
# (a residue buried in the binding groove works best as the Get_Cleft anchor).
# 1TW6: Smac AVPI tetrapeptide in chain C — anchor on Pro3, central to the motif.
PEPTIDE_LIGAND: dict[str, tuple[str, int, str]] = {
    "1TW6": ("PRO", 3, "C"),
}

RCSB_URL = "https://files.rcsb.org/download/{pdb}.pdb"

# Output directory: benchmarks/astex_diverse/astex_diverse/<PDB>/
SCRIPT_DIR = Path(__file__).parent.resolve()
ASTEX_DIR = SCRIPT_DIR / "astex_diverse"


def download_rcsb_pdb(pdb: str, dest: Path) -> Path:
    """Download complex PDB from RCSB into dest directory."""
    url = RCSB_URL.format(pdb=pdb)
    out_path = dest / f"{pdb}.pdb"
    print(f"  Downloading {url} ...", end=" ", flush=True)
    urllib.request.urlretrieve(url, out_path)
    print(f"({out_path.stat().st_size // 1024} kB)")
    return out_path


def find_ligand_anchor(pdb_path: Path, resname: str) -> tuple[int, str] | None:
    """
    Scan HETATM records for the first occurrence of resname.
    Returns (resnum, chain) or None if not found.
    PDB columns (0-indexed): resname cols 17-19, resnum cols 22-25, chain col 21.
    """
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            rnam = line[17:20].strip()
            if rnam == resname:
                chain = line[21].strip() or "-"
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    continue
                return resnum, chain
    return None


def build_anchor_string(resname: str, resnum: int, chain: str) -> str:
    """
    Format the Get_Cleft -a argument: RESNAMERESNUMCHAINALTLOC
    RESNAME = exactly 3 chars (pad with leading '-' for shorter names)
    CHAIN   = '-' for blank/space
    ALTLOC  = '-' (no alt loc for crystal ligand)
    """
    # Pad short resnames: Get_Cleft treats '-' as space when nSpaces > 0
    if len(resname) < 3:
        resname = "-" * (3 - len(resname)) + resname
    chain_char = chain if (chain and chain != " ") else "-"
    return f"{resname}{resnum}{chain_char}-"


def run_get_cleft(get_cleft: Path, pdb_path: Path, anchor: str, work_dir: Path) -> Path | None:
    """
    Run: Get_Cleft -p <pdb_path> -o <pdb_base> -a <anchor>
    Get_Cleft names output as <outbase>_<anchorinfo>_clf_1.pdb.
    We pass -o <pdb_base> so the file starts with the PDB ID, then glob.
    Returns the path to the output file, or None on failure.
    """
    pdb_base = pdb_path.name[:4]  # first 4 chars = PDB ID

    # Pass -o <pdb_base> so Get_Cleft prefixes output with the PDB ID.
    # Without -o, it defaults to "." which creates "._<anchor>_clf_1.pdb".
    cmd = [str(get_cleft), "-p", str(pdb_path), "-o", pdb_base, "-a", anchor]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: Get_Cleft exited {result.returncode}")
        print(f"  stderr: {result.stderr[:400]}")
        return None

    # Glob for the output cleft file.  Get_Cleft preserves the rank in the
    # filename (_clf_N.pdb where N is the cavity rank, not always 1).
    matches = sorted(work_dir.glob("*_clf_*.pdb"))
    if not matches:
        print(f"  ERROR: no *_clf_1.pdb found in {work_dir}")
        print(f"  stdout: {result.stdout[:400]}")
        return None

    clf_out = matches[0]
    print(f"  Get_Cleft output: {clf_out.name}")
    return clf_out


def find_peptide_anchor(pdb_path: Path, resname: str, resnum: int, chain: str
                        ) -> tuple[int, str] | None:
    """
    Confirm an ATOM-record peptide anchor residue exists.  Peptide ligands (e.g.
    1TW6 Smac AVPI) are stored as ATOM, not HETATM, so the HETATM scanner misses
    them.  Returns (resnum, chain) if the residue is present, else None.
    """
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if (line[17:20].strip() == resname
                    and line[21].strip() == chain.strip()
                    and line[22:26].strip() == str(resnum)):
                return resnum, chain
    return None


def process_one(pdb: str, get_cleft: Path, dry_run: bool) -> bool:
    """Process a single Astex complex. Returns True on success."""
    peptide = PEPTIDE_LIGAND.get(pdb)
    resname = peptide[0] if peptide else LIGAND_RESNAME.get(pdb)
    if resname is None:
        print(f"[{pdb}] SKIP — no ligand resname mapping")
        return False

    dest_dir = ASTEX_DIR / pdb
    dest_file = dest_dir / f"{pdb}_binding_site.pdb"

    print(f"\n[{pdb}] ligand={resname}")

    if dest_file.exists():
        print(f"  Already done: {dest_file} — skipping")
        return True

    with tempfile.TemporaryDirectory(prefix=f"gc_{pdb}_") as tmp:
        work = Path(tmp)

        # 1. Download fresh RCSB complex PDB
        try:
            pdb_path = download_rcsb_pdb(pdb, work)
        except Exception as e:
            print(f"  ERROR: download failed: {e}")
            return False

        # 2. Find ligand anchor (resnum + chain).  Peptide ligands are stored as
        #    ATOM records; small-molecule ligands as HETATM.
        if peptide:
            _, p_resnum, p_chain = peptide
            anchor_info = find_peptide_anchor(pdb_path, resname, p_resnum, p_chain)
            if anchor_info is None:
                print(f"  ERROR: no ATOM peptide anchor {resname}{p_resnum}{p_chain} in PDB")
                return False
        else:
            anchor_info = find_ligand_anchor(pdb_path, resname)
            if anchor_info is None:
                print(f"  ERROR: no HETATM {resname!r} found in downloaded PDB")
                return False
        resnum, chain = anchor_info
        anchor = build_anchor_string(resname, resnum, chain)
        print(f"  Anchor: {resname} resnum={resnum} chain={chain!r} → -a {anchor}")

        if dry_run:
            print(f"  DRY RUN: would run Get_Cleft -p {pdb}.pdb -a {anchor}")
            print(f"  DRY RUN: would copy _clf_1.pdb → {dest_file}")
            return True

        # 3. Run Get_Cleft
        clf_out = run_get_cleft(get_cleft, pdb_path, anchor, work)
        if clf_out is None:
            return False

        # 4. Copy _clf_1.pdb → <dest>/<PDB>_binding_site.pdb
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(clf_out, dest_file)
        print(f"  Written: {dest_file} ({dest_file.stat().st_size} bytes)")

    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--get-cleft", required=True, type=Path,
                        help="Path to compiled Get_Cleft binary")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without running them")
    parser.add_argument("--targets", nargs="+", metavar="PDB",
                        help="Subset of PDB IDs to process (default: all 85)")
    args = parser.parse_args()

    if not args.dry_run and not args.get_cleft.is_file():
        print(f"ERROR: Get_Cleft binary not found: {args.get_cleft}", file=sys.stderr)
        return 1

    known = set(LIGAND_RESNAME) | set(PEPTIDE_LIGAND)
    targets = [t.upper() for t in args.targets] if args.targets else sorted(known)
    unknown = [t for t in targets if t not in known]
    if unknown:
        print(f"ERROR: unknown PDB IDs: {unknown}", file=sys.stderr)
        return 1

    print(f"Processing {len(targets)} complexes (dry_run={args.dry_run})")
    print(f"Destination: {ASTEX_DIR}/")

    ok, fail = 0, 0
    for pdb in targets:
        if process_one(pdb, args.get_cleft, args.dry_run):
            ok += 1
        else:
            fail += 1

    print(f"\nDone: {ok} succeeded, {fail} failed")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
