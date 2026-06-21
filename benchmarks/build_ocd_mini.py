#!/usr/bin/env python3
"""
Build benchmark_crossdock_ocd_mini.json — same-protein OCD pairs from Astex Diverse.

Confirmed pairs (protein identity from CIF entity descriptions + organism source):
  - Factor Xa   : 1LPZ (human FXa) ← 1MQ6 (human FXa)         [Tier 1: same protein, same organism]
  - DHFR        : 1IA1 (C.albicans) ← 1S3V (H.sapiens)         [Tier 2: same enzyme, different organism]
  - Neuraminidase: 1L7F (Influenza A) ← 1VCJ (Influenza B)     [Tier 2: same enzyme, different strain]

Bidirectional pairs → 6 entries total.
rmsd_ref_sdf = receptor's own native ligand (in-frame, no superposition needed for same-protein).
"""
import json
from pathlib import Path

ASTEX = "/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse/astex_diverse"
OUT   = "/Users/lp.more/Projects/FlexAIDdS/benchmarks/datasets/benchmark_crossdock_ocd_mini.json"


def make_pair(index, receptor_id, donor_id, family, tier):
    rec = receptor_id.upper()
    don = donor_id.upper()
    return {
        "index":           index,
        "receptor_id":     rec,
        "ligand_id":       don,
        "family":          family,
        "tier":            tier,
        "same_protein":    True,
        # DatasetRunner core fields
        "receptor_pdb":    f"{ASTEX}/{rec}/{rec}_apo.pdb",
        "ligand_sdf":      f"{ASTEX}/{don}/{don}_ligand.sdf",
        "oracle_site_pdb": f"{ASTEX}/{rec}/{rec}_binding_site.pdb",
        # RMSD reference = receptor's own native ligand (same coordinate frame)
        "rmsd_ref_sdf":    f"{ASTEX}/{rec}/{rec}_ligand.sdf",
        # Metadata
        "receptor_name":   rec,
        "donor_name":      don,
        "note": (
            "Oracle cross-docking: receptor's binding_site.pdb used as oracle cleft. "
            "rmsd_ref_sdf = receptor native ligand (same frame as receptor_pdb). "
            "NO superposition needed — both are in the receptor coordinate system."
        ),
    }


pairs_spec = [
    # (receptor_id, donor_id, family, tier)
    ("1LPZ", "1MQ6", "Factor Xa",     1),   # same protein, same organism
    ("1MQ6", "1LPZ", "Factor Xa",     1),   # reverse
    ("1IA1", "1S3V", "DHFR",          2),   # C.albicans ← H.sapiens
    ("1S3V", "1IA1", "DHFR",          2),   # reverse
    ("1L7F", "1VCJ", "Neuraminidase", 2),   # Influenza A ← Influenza B
    ("1VCJ", "1L7F", "Neuraminidase", 2),   # reverse
]

pairs = [make_pair(i, r, d, fam, tier) for i, (r, d, fam, tier) in enumerate(pairs_spec)]

# Verify files exist
print("Checking file existence...")
ok = True
for p in pairs:
    for key in ["receptor_pdb", "ligand_sdf", "oracle_site_pdb", "rmsd_ref_sdf"]:
        path = Path(p[key])
        exists = path.exists()
        status = "OK" if exists else "MISSING"
        if not exists:
            ok = False
        print(f"  [{status}] {key}: {path.name}  ({p['receptor_id']}←{p['ligand_id']})")

manifest = {
    "schema_version":  1,
    "name":            "astex_ocd_mini",
    "description": (
        "Oracle cross-docking mini-benchmark: 6 bidirectional same-protein pairs from "
        "Astex Diverse (3 unique pairs × 2 directions). "
        "Factor Xa (Tier 1: same protein, same organism), "
        "DHFR (Tier 2: same enzyme, C.albicans/H.sapiens), "
        "Neuraminidase (Tier 2: same enzyme, Influenza A/B). "
        "Oracle mode: binding_site.pdb defines the search cleft. "
        "RMSD reference: receptor's own native ligand (no superposition needed)."
    ),
    "oracle_mode":       True,
    "n_pairs":           len(pairs),
    "astex_diverse_dir": ASTEX,
    "pairs":             pairs,
}

with open(OUT, "w") as f:
    json.dump(manifest, f, indent=2)

print(f"\n{'OK' if ok else 'ERRORS'}: {len(pairs)} pairs → {OUT}")
print(json.dumps(pairs[0], indent=2))
