#!/usr/bin/env python3
"""
Extract protein identities from Astex Diverse mmCIF files.
Outputs: pdb_id, title, polymer entities (comma-sep)
"""
import os, glob, re

astex_dir = "/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse/astex_diverse"
results = {}


def parse_cif_value(line):
    """Extract quoted or unquoted value from a CIF token string."""
    line = line.strip()
    m = re.match(r"'(.+?)'$", line)
    if m:
        return m.group(1)
    m = re.match(r'"(.+?)"$', line)
    if m:
        return m.group(1)
    return line


def extract_cif_info(cif_path):
    """Return (title, [polymer_descriptions]) from a CIF file."""
    with open(cif_path, errors='replace') as f:
        content = f.read()

    # --- struct.title (single-value field, value on same line) ---
    title = ""
    m = re.search(r"^_struct\.title\s+(.+)$", content, re.MULTILINE)
    if m:
        title = parse_cif_value(m.group(1))

    # --- _entity loop: find column indices, then parse data rows ---
    entities = []
    # Find all entity loops
    # Pattern: loop_ ... _entity.id ... _entity.pdbx_description ... <data rows>
    loop_pat = re.compile(
        r"loop_\n((?:_entity\.\S+\s*\n)+)((?:(?!loop_|^_|^#)[\s\S])*?)(?=\n#|\nloop_|\ndata_|\Z)",
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

        # Tokenize the body: respect quoted strings
        tokens = re.findall(r"'[^']*'|\"[^\"]*\"|\S+", body)

        # Group into rows of ncols
        for i in range(0, len(tokens) - ncols + 1, ncols):
            row = tokens[i:i + ncols]
            if len(row) < ncols:
                break
            etype = row[type_idx].strip("'\"")
            desc  = row[desc_idx].strip("'\"")
            if etype.lower() == "polymer":
                entities.append(desc)

    return title, entities


for pdb_dir in sorted(os.listdir(astex_dir)):
    full_path = os.path.join(astex_dir, pdb_dir)
    if not os.path.isdir(full_path):
        continue
    pdb_id = pdb_dir.lower()

    cif_files = glob.glob(f"{full_path}/*.cif") + glob.glob(f"{full_path}/*.mmcif")
    title, entities = "", []
    for cif in cif_files:
        try:
            title, entities = extract_cif_info(cif)
            if title or entities:
                break
        except Exception as e:
            print(f"  ERROR {pdb_id}: {e}")

    results[pdb_id] = {"title": title, "entities": entities}
    ent_str = " | ".join(entities) if entities else "(none)"
    print(f"{pdb_id}: {title[:70]}")
    print(f"        => {ent_str[:100]}")

print(f"\nTotal: {len(results)} entries")
