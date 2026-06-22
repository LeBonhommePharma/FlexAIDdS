#!/usr/bin/env python3
"""
Deep-parse Astex Diverse CIFs to find same-protein pairs.
Extracts: entity descriptions, EC numbers, keywords, and checks file availability.
"""
import os, glob, re, json

astex_dir = "/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse/astex_diverse"

def parse_cif_value(s):
    s = s.strip()
    m = re.match(r"^'(.*)'$", s)
    if m: return m.group(1)
    m = re.match(r'^"(.*)"$', s)
    if m: return m.group(1)
    return s if s != "?" else ""

def extract_cif_info(cif_path):
    with open(cif_path, errors='replace') as f:
        content = f.read()

    # struct.title
    title = ""
    m = re.search(r"^_struct\.title\s+(.+)$", content, re.MULTILINE)
    if m: title = parse_cif_value(m.group(1))

    # EC numbers from _entity.pdbx_ec
    ec_numbers = re.findall(r"_entity\.pdbx_ec\s+(\S+)", content)
    ec_clean = [e.strip("'\"") for e in ec_numbers if e not in ("?", "")]

    # Keywords
    keywords = ""
    m = re.search(r"^_struct_keywords\.text\s+(.+)$", content, re.MULTILINE)
    if m: keywords = parse_cif_value(m.group(1))

    # Polymer entity descriptions via loop
    entities = []
    loop_pat = re.compile(
        r"loop_\n((?:_entity\.\S+\s*\n)+)([\s\S]*?)(?=\n#|\nloop_|\ndata_|\Z)",
        re.MULTILINE
    )
    for lm in loop_pat.finditer(content):
        header = lm.group(1)
        body   = lm.group(2)
        cols = re.findall(r"_entity\.(\S+)", header)
        if "pdbx_description" not in cols or "type" not in cols:
            continue
        desc_idx = cols.index("pdbx_description")
        type_idx = cols.index("type")
        ncols = len(cols)
        tokens = re.findall(r"'[^']*'|\"[^\"]*\"|\S+", body)
        for i in range(0, len(tokens) - ncols + 1, ncols):
            row = tokens[i:i + ncols]
            if len(row) < ncols: break
            etype = row[type_idx].strip("'\"")
            desc  = row[desc_idx].strip("'\"")
            if etype.lower() == "polymer":
                entities.append(desc)

    return title, ec_clean, keywords, entities

def check_files(pdb_dir):
    """Return dict of existing relevant files."""
    found = {}
    base = os.path.join(astex_dir, pdb_dir)
    pdb_id = pdb_dir.lower()
    for ext in ["_protein.pdb", "_ligand.sdf", "_ligand.mol2"]:
        f = os.path.join(base, f"{pdb_id}{ext}")
        if os.path.exists(f):
            found[ext] = f
    # Also check uppercase PDB id variants
    for ext in ["_protein.pdb", "_ligand.sdf", "_ligand.mol2"]:
        f = os.path.join(base, f"{pdb_dir.upper()}{ext}")
        if os.path.exists(f) and ext not in found:
            found[ext] = f
    # Check CIFs
    cifs = glob.glob(f"{base}/*.cif") + glob.glob(f"{base}/*.mmcif")
    if cifs:
        found["cif"] = cifs[0]
    return found

print("=== Astex Diverse Protein Identity Scan ===\n")

all_data = {}
for pdb_dir in sorted(os.listdir(astex_dir)):
    full_path = os.path.join(astex_dir, pdb_dir)
    if not os.path.isdir(full_path):
        continue
    pdb_id = pdb_dir.lower()
    cif_files = glob.glob(f"{full_path}/*.cif") + glob.glob(f"{full_path}/*.mmcif")
    title, ecs, kw, entities = "", [], "", []
    for cif in cif_files:
        try:
            title, ecs, kw, entities = extract_cif_info(cif)
            break
        except: pass

    files = check_files(pdb_dir)
    has_receptor = "_protein.pdb" in files
    has_ligand   = "_ligand.sdf" in files

    all_data[pdb_id] = {
        "title": title,
        "ec": ecs,
        "keywords": kw,
        "entities": entities,
        "has_receptor": has_receptor,
        "has_ligand": has_ligand,
        "files": files
    }

# Print summary for grouping
print(f"{'PDB':6} {'EC':15} {'Has_R':6} {'Has_L':6}  Entities")
print("-"*100)
for pid, d in sorted(all_data.items()):
    ec_str = ",".join(d["ec"]) if d["ec"] else "-"
    ent_str = " | ".join(d["entities"])[:80] if d["entities"] else "-"
    print(f"{pid:6} {ec_str:15} {str(d['has_receptor']):6} {str(d['has_ligand']):6}  {ent_str}")

# --- Find pairs by EC number ---
print("\n=== Pairs grouped by EC number ===")
from collections import defaultdict
ec_groups = defaultdict(list)
for pid, d in all_data.items():
    for ec in d["ec"]:
        if ec and ec != "?":
            ec_groups[ec].append(pid)

for ec, pids in sorted(ec_groups.items()):
    if len(pids) >= 2:
        print(f"EC {ec}: {pids}")

# --- Find pairs by entity name similarity ---
print("\n=== Pairs grouped by primary entity description ===")
ent_groups = defaultdict(list)
for pid, d in all_data.items():
    if d["entities"]:
        # Normalize: lowercase, strip spaces, first entity
        key = d["entities"][0].lower().strip()
        ent_groups[key].append(pid)

for ent, pids in sorted(ent_groups.items()):
    if len(pids) >= 2:
        print(f"'{ent}': {pids}")

# Save for next step
with open("/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_identities.json", "w") as f:
    json.dump(all_data, f, indent=2)
print("\nSaved -> astex_identities.json")
