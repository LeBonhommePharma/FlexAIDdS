#!/usr/bin/env python3
"""
fetch_itc_data.py — Build a unified ITC calibration table from many sources.

Casts the widest net for published protein–ligand (and host–guest) ITC
thermodynamics, normalizes everything to one schema, deduplicates, and reports
measurement variance for complexes measured more than once.

Unified schema (kcal/mol, one row per measurement):
    pdb_id, ligand_smiles, dH_kcal_mol, TdS_kcal_mol, dG_kcal_mol, T_K, source, doi

Sign convention:  dG = dH - TdS   (TdS is the entropy *contribution* T·ΔS).

Sources (see SOURCES registry below for access mode + URLs):
    scorpio     local curated CSV  (kcal/mol, has PDB IDs)            [FREE/local]
    bindingdb   BindingDB ITC TSV  (kJ/mol, has PDB IDs + SMILES)     [FREE/local+http]
    pdbbind     INDEX_general_PL_data  (affinity → dG only, no dH)    [MANUAL register]
    chembl      ChEMBL API (standard_type dH/Kd/Ka)                   [FREE API]
    csar        CSAR 2010/2012 NRC set                                [MANUAL register]
    freire      Velazquez-Campoy & Freire JACS/DDT compilations       [MANUAL paywalled SI]
    sampl       SAMPL4-7 host–guest ITC (CB/CD/octa-acid)             [FREE GitHub]
    biolip      BioLiP annotations                                    [FREE download]
    nist        IUPAC-NIST thermochemistry                            [MANUAL]
    csv         any file already in the unified schema                [local]

Usage
-----
    # Build unified table from everything available locally:
    python benchmarks/fetch_itc_data.py --sources scorpio bindingdb \
        --out benchmarks/itc_unified.csv

    # Add a manually-downloaded source file:
    python benchmarks/fetch_itc_data.py --sources pdbbind \
        --pdbbind-index ~/PDBbind/INDEX_general_PL_data.2020 \
        --out itc_pdbbind.csv

    # List every source and how to obtain it:
    python benchmarks/fetch_itc_data.py --list-sources
"""

import argparse
import csv
import os
import sys

KJ_PER_KCAL = 4.184
HOME = os.path.expanduser("~")

UNIFIED_COLS = ["pdb_id", "ligand_smiles", "dH_kcal_mol", "TdS_kcal_mol",
                "dG_kcal_mol", "T_K", "source", "doi"]


# ─────────────────────────────────────────────────────────────────────────────
# Source registry — access mode, URL, and where the local copy is expected.
# access: "local" | "http_free" | "api_free" | "manual_download"
# ─────────────────────────────────────────────────────────────────────────────
SOURCES = {
    "scorpio": {
        "access": "local",
        "url": "https://scorpio2.biophysics.ismb.lon.ac.uk/structure/itc_data/{id}",
        "local": os.path.join(HOME, "Documents/PhD/Programs/FlexAIDdS/data/SCORPIO/scorpio_itc_raw.csv"),
        "doi": "10.1016/j.jmb.2008.09.073",  # Olsson et al. JMB 384:1002 (2008)
        "notes": "Curated scrape (239 rows / 104 PDBs). Already kcal/mol with PDB IDs and T.",
    },
    "bindingdb": {
        "access": "http_free",
        "url": "https://www.bindingdb.org/bind/downloads.jsp",  # 'BindingDB_All' or ITC export
        "local": os.path.join(HOME, "Documents/PhD/Docking/BindingDB/BindingDB_ITC.tsv"),
        "doi": "10.1093/nar/gkv1072",
        "notes": "ITC subset TSV (kJ/mol). Filter: Measurement Type=ITC, Delta_H0 non-null, "
                 "PDB ID present, Temp=25C preferred. SMILES included.",
    },
    "pdbbind": {
        "access": "manual_download",
        "url": "http://www.pdbbind.org.cn (register; general+refined sets)",
        "local": None,  # pass --pdbbind-index
        "doi": "10.1021/jm048957q",
        "notes": "INDEX_general_PL_data.<year>: affinity only (-logKd/Ki) -> dG, NO dH/TdS. "
                 "Use for dG validation; pair with BindingDB for dH. Pass --pdbbind-index.",
    },
    "chembl": {
        "access": "api_free",
        "url": "https://www.ebi.ac.uk/chembl/api/data/activity?standard_type=dH&format=json",
        "local": None,
        "doi": "10.1093/nar/gkad1004",
        "notes": "REST API. assay_type=B, standard_type in (dH,Kd,Ka). Most entries lack a "
                 "matched crystal; join to PDBe SIFTS for pdb_id. Implemented as documented stub.",
    },
    "csar": {
        "access": "manual_download",
        "url": "http://www.csardock.org (CSAR-NRC HiQ 2010 / CSAR 2012)",
        "local": None,
        "doi": "10.1021/ci100366a",
        "notes": "Mostly dG; a subset has ITC dH. Download the affinity table, save as unified CSV.",
    },
    "freire": {
        "access": "manual_download",
        "url": "Velazquez-Campoy & Freire — JACS 2005/2006; Drug Discov. Today 2008 13:869",
        "local": None,
        "doi": "10.1016/j.drudis.2008.07.005",
        "notes": "Canonical drug-target thermodynamic signatures (HIV protease, statins, etc.). "
                 "Tables are in paywalled SI; transcribe to unified CSV and load via --source csv.",
    },
    "sampl": {
        "access": "http_free",
        "url": "https://github.com/samplchallenges (SAMPL4..SAMPL7 host_guest/ITC)",
        "local": None,
        "doi": "10.1007/s10822-018-0170-6",
        "notes": "Host–guest (CB7/8, β-CD, octa-acid) explicit ITC dH/TdS. No PDB id "
                 "(use host_guest label as pdb_id). Valid for dH/TdS scale calibration.",
    },
    "biolip": {
        "access": "http_free",
        "url": "https://zhanggroup.org/BioLiP/download.html",
        "local": None,
        "doi": "10.1093/nar/g1271",
        "notes": "Has PDB ids + ligand; binding-affinity field is sparse/mixed assay, few ITC dH. "
                 "Use to recover pdb_id/SMILES for affinity-only entries.",
    },
    "nist": {
        "access": "manual_download",
        "url": "https://www.nist.gov (IUPAC-NIST solubility/thermochem compilations)",
        "local": None,
        "doi": "",
        "notes": "Few protein–ligand entries; manual extraction to unified CSV.",
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# helpers
# ─────────────────────────────────────────────────────────────────────────────
def _f(v, scale=1.0):
    if v is None:
        return None
    s = str(v).strip()
    if s == "" or s.upper() in ("NA", "NAN", "NULL", "N/A"):
        return None
    # strip trailing unit tokens like "37.00 C"
    s = s.split()[0] if " " in s else s
    try:
        return float(s) * scale
    except ValueError:
        return None


def _first_pdb(field):
    s = (field or "").strip()
    if not s:
        return ""
    for tok in s.replace(";", ",").split(","):
        tok = tok.strip()
        if len(tok) == 4 and tok[0].isdigit():
            return tok.upper()
    return s.split(",")[0].strip().upper()


def _complete_dg(dH, TdS, dG):
    """Fill the missing one of (dH, TdS, dG) using dG = dH - TdS."""
    if dG is None and dH is not None and TdS is not None:
        dG = dH - TdS
    elif TdS is None and dH is not None and dG is not None:
        TdS = dH - dG
    elif dH is None and TdS is not None and dG is not None:
        dH = dG + TdS
    return dH, TdS, dG


# ─────────────────────────────────────────────────────────────────────────────
# parsers (yield unified-row dicts)
# ─────────────────────────────────────────────────────────────────────────────
def parse_scorpio(path):
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            dH = _f(row.get("dH_kcalmol"))
            TdS = _f(row.get("TdS_kcalmol"))
            dG = _f(row.get("dG_kcalmol"))
            dH, TdS, dG = _complete_dg(dH, TdS, dG)
            if dH is None and dG is None:
                continue
            yield {
                "pdb_id": (row.get("pdb_id") or "").strip().upper(),
                "ligand_smiles": "",  # SCORPIO has names, not SMILES
                "dH_kcal_mol": dH, "TdS_kcal_mol": TdS, "dG_kcal_mol": dG,
                "T_K": _f(row.get("temperature_K")) or 298.15,
                "source": "scorpio",
                "doi": row.get("source_url") or SOURCES["scorpio"]["doi"],
            }


def parse_bindingdb(path):
    sc = 1.0 / KJ_PER_KCAL
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            dH = _f(row.get("Delta_H0 (kJ/mol)"), sc)
            mTdS = _f(row.get("-T Delta_S0 (kJ/mol)"), sc)   # this is -T*dS
            TdS = (-mTdS) if mTdS is not None else None
            dG = _f(row.get("Delta_G0 (kJ/mol)"), sc)
            dH, TdS, dG = _complete_dg(dH, TdS, dG)
            if dH is None and dG is None:
                continue
            tc = _f(row.get("Temp (C)"))
            yield {
                "pdb_id": _first_pdb(row.get("PDB ID(s) for Ligand-Target Complex")),
                "ligand_smiles": (row.get("Ligand SMILES") or "").strip(),
                "dH_kcal_mol": dH, "TdS_kcal_mol": TdS, "dG_kcal_mol": dG,
                "T_K": (tc + 273.15) if tc is not None else 298.15,
                "source": "bindingdb",
                "doi": (row.get("Article DOI") or "").strip(),
            }


def parse_pdbbind_index(path):
    """INDEX_general_PL_data.<year>: '<code> <res> <year> <-logKd/Ki> <Kd=..> // ...'
    Affinity only -> dG = -RT ln(10) * pKd at 298 K. No dH/TdS."""
    R_kcal = 1.987204e-3
    RT_ln10 = R_kcal * 298.15 * 2.302585
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 4:
                continue
            try:
                pk = float(p[3])
            except ValueError:
                continue
            yield {
                "pdb_id": p[0].upper(), "ligand_smiles": "",
                "dH_kcal_mol": None, "TdS_kcal_mol": None,
                "dG_kcal_mol": -RT_ln10 * pk, "T_K": 298.15,
                "source": "pdbbind", "doi": SOURCES["pdbbind"]["doi"],
            }


def parse_unified_csv(path):
    """Read a file already in the unified schema (or a transcribed Freire/CSAR/SAMPL table)."""
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            low = {k.lower().strip(): v for k, v in row.items()}
            dH = _f(low.get("dh_kcal_mol"))
            TdS = _f(low.get("tds_kcal_mol"))
            dG = _f(low.get("dg_kcal_mol"))
            dH, TdS, dG = _complete_dg(dH, TdS, dG)
            yield {
                "pdb_id": (low.get("pdb_id") or "").strip().upper(),
                "ligand_smiles": (low.get("ligand_smiles") or "").strip(),
                "dH_kcal_mol": dH, "TdS_kcal_mol": TdS, "dG_kcal_mol": dG,
                "T_K": _f(low.get("t_k")) or 298.15,
                "source": (low.get("source") or "csv").strip(),
                "doi": (low.get("doi") or "").strip(),
            }


PARSERS = {
    "scorpio": parse_scorpio,
    "bindingdb": parse_bindingdb,
    "pdbbind": parse_pdbbind_index,
    "csv": parse_unified_csv,
}


# ─────────────────────────────────────────────────────────────────────────────
# dedup + variance
# ─────────────────────────────────────────────────────────────────────────────
def dedup_report(rows):
    """Group by (pdb_id, ligand_smiles); keep ALL rows, report multi-measurement variance."""
    groups = {}
    for r in rows:
        key = (r["pdb_id"], r["ligand_smiles"])
        groups.setdefault(key, []).append(r)

    def spread(vals):
        vals = [v for v in vals if v is not None]
        if len(vals) < 2:
            return 0.0
        m = sum(vals) / len(vals)
        return (sum((v - m) ** 2 for v in vals) / len(vals)) ** 0.5

    conflicts = []
    for key, g in groups.items():
        if len(g) > 1 and key[0]:
            conflicts.append({
                "pdb_id": key[0], "n": len(g),
                "sources": sorted({x["source"] for x in g}),
                "sd_dH": round(spread([x["dH_kcal_mol"] for x in g]), 2),
                "sd_dG": round(spread([x["dG_kcal_mol"] for x in g]), 2),
            })
    return conflicts


# ─────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sources", nargs="+", default=["scorpio", "bindingdb"],
                    help="Sources to integrate (default: scorpio bindingdb).")
    ap.add_argument("--out", default="benchmarks/itc_unified.csv",
                    help="Output unified CSV path.")
    ap.add_argument("--scorpio-csv", default=SOURCES["scorpio"]["local"])
    ap.add_argument("--bindingdb-tsv", default=SOURCES["bindingdb"]["local"])
    ap.add_argument("--pdbbind-index", default=None,
                    help="Path to INDEX_general_PL_data.<year> (manual download).")
    ap.add_argument("--csv-files", nargs="*", default=[],
                    help="Extra files already in unified schema (Freire/CSAR/SAMPL transcriptions).")
    ap.add_argument("--require-pdb", action="store_true",
                    help="Drop rows without a PDB id (needed for docking-based calibration).")
    ap.add_argument("--require-dh", action="store_true",
                    help="Drop rows without dH (needed for alpha calibration).")
    ap.add_argument("--list-sources", action="store_true")
    args = ap.parse_args()

    if args.list_sources:
        print(f"{'source':<11} {'access':<16} url / notes")
        print("-" * 90)
        for name, s in SOURCES.items():
            loc = "  local: " + (s["local"] if s["local"] and os.path.isfile(s["local"]) else "(not present)")
            print(f"{name:<11} {s['access']:<16} {s['url']}")
            print(f"{'':<28} {s['notes']}")
            print(f"{'':<28}{loc}")
        return

    rows = []
    for src in args.sources:
        try:
            if src == "scorpio":
                path = args.scorpio_csv
                if not os.path.isfile(path):
                    print(f"[skip] scorpio: {path} not found", file=sys.stderr); continue
                got = list(parse_scorpio(path))
            elif src == "bindingdb":
                path = args.bindingdb_tsv
                if not os.path.isfile(path):
                    print(f"[skip] bindingdb: {path} not found", file=sys.stderr); continue
                got = list(parse_bindingdb(path))
            elif src == "pdbbind":
                if not args.pdbbind_index or not os.path.isfile(args.pdbbind_index):
                    print("[skip] pdbbind: pass --pdbbind-index <INDEX_general_PL_data.YYYY> "
                          "(register at pdbbind.org.cn)", file=sys.stderr); continue
                got = list(parse_pdbbind_index(args.pdbbind_index))
            elif src in ("chembl", "csar", "freire", "sampl", "biolip", "nist"):
                s = SOURCES[src]
                print(f"[manual] {src}: {s['access']} — {s['url']}\n          {s['notes']}\n"
                      f"          Transcribe/export to the unified schema and pass via --csv-files.",
                      file=sys.stderr)
                continue
            else:
                print(f"[warn] unknown source '{src}'", file=sys.stderr); continue
            print(f"[ok] {src}: {len(got)} rows")
            rows.extend(got)
        except OSError as e:
            print(f"[err] {src}: {e}", file=sys.stderr)

    for cf in args.csv_files:
        if os.path.isfile(cf):
            got = list(parse_unified_csv(cf))
            print(f"[ok] csv {cf}: {len(got)} rows")
            rows.extend(got)
        else:
            print(f"[skip] csv {cf}: not found", file=sys.stderr)

    if args.require_pdb:
        rows = [r for r in rows if r["pdb_id"]]
    if args.require_dh:
        rows = [r for r in rows if r["dH_kcal_mol"] is not None]

    if not rows:
        sys.exit("No rows collected. Try --list-sources to see what's available.")

    conflicts = dedup_report(rows)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=UNIFIED_COLS)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in UNIFIED_COLS})

    n_pdb = sum(1 for r in rows if r["pdb_id"])
    n_dh = sum(1 for r in rows if r["dH_kcal_mol"] is not None)
    n_smiles = sum(1 for r in rows if r["ligand_smiles"])
    print(f"\nUnified table: {len(rows)} rows  "
          f"(with_pdb={n_pdb}, with_dH={n_dh}, with_SMILES={n_smiles})")
    print(f"Multi-measurement complexes (variance reported): {len(conflicts)}")
    for c in sorted(conflicts, key=lambda x: -x["sd_dH"])[:8]:
        print(f"  {c['pdb_id']}: n={c['n']} sources={'+'.join(c['sources'])} "
              f"sd_dH={c['sd_dH']} sd_dG={c['sd_dG']} kcal/mol")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
