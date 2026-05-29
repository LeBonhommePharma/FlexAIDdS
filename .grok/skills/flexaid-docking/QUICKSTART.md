# Quickstart: Real FlexAIDδS Docking with Thermodynamics

This guide shows the fastest, most reliable way to run actual molecular docking + thermodynamic analysis using the `flexaid-docking` skill.

It assumes you want production-quality results (Voronoi CF/contact-function scoring + full ensemble thermodynamic ledger), not just smoke tests.

## How the Skill Interacts — Clarification Flow

The skill is designed to be conversational. When you trigger it with natural language (e.g. "dock biotin back into streptavidin", "redock the ligand from 1STP with thermodynamics", or "update the flexaid-docking skill"), it will ask focused clarifying questions **before** taking any action.

This prevents incorrect assumptions about organism, target chains, ligand source, self-docking vs cross-docking, thermodynamic depth, or special constraints.

**Example clarification flow:**

> **User:** Redock the co-crystallized ligand from PDB 1STP with full thermodynamic analysis.
>
> **Skill:** To do this correctly and reproducibly, I need to confirm a few details first:
>
> 1. Biological context — organism/species? (e.g. *Streptomyces avidinii* for streptavidin)
> 2. Is this self-docking (redocking the known biotin) or a different ligand?
> 3. Thermodynamic requirements — full ensemble ledger (Helmholtz F, entropy terms, Boltzmann weights) at a specific temperature? Include tENCoM vibrational entropy?
> 4. Any special constraints (covalent ligand, modified residues, specific chains only, NMR multi-model selection, bio-unit .pdb1 preference)?
> 5. Input preference — automatic download + intelligent splitting via the skill's tools, or do you already have prepared receptor/ligand files?
>
> Please reply with the details (or say "use defaults for 1STP biotin self-docking at 300 K" if that matches your intent). Once I have clarity I'll generate the exact validated workflow.

This behavior is enforced by the skill's activation rules.

## 1. Prerequisites

### A. Activate the skill's recommended Python environment
```bash
source ~/.venvs/flexaidds/bin/activate
```

**Tip:** For even shorter commands you can use the skill's `bin/` shortcuts (they are pure symlinks — identical behavior):
```bash
.grok/skills/flexaid-docking/bin/ensure-docking-data
```

### B. Ensure the FlexAIDδS binary is built
```bash
ls /path/to/your/build/FlexAIDδS
```

### C. Ensure the critical interaction matrices are available
```bash
python3 .grok/skills/flexaid-docking/scripts/ensure_docking_data.py
```

This command is **mandatory** before real docking. It will use the matrices bundled inside the skill or copy them from a known-good installation if you provide `--source`.

## 2. Prepare Clean Inputs

The skill strongly recommends high-quality, minimal inputs:

- **Receptor**: Clean PDB containing only the protein (or relevant chains), no waters or ligands.
- **Ligand**: MOL2 format with hydrogens and reasonable charges (Gasteiger is a good default).

Example preparation (using Open Babel):

```bash
# Extract ligand as MOL2 with hydrogens + Gasteiger charges
obabel ligand.pdb -O ligand.mol2 -h --partialcharge gasteiger
```

See `data/README.md` and the skill's guidance for more details on input hygiene.

## 3. Run a Docking Job

### Recommended: Use the high-level Python API

```python
import flexaidds as fd

results = fd.dock(
    receptor="receptor_clean.pdb",
    ligand="ligand.mol2",
    binary="/path/to/your/build/FlexAIDδS",
    compute_entropy=True,      # Enables full thermodynamic ledger
    temperature=300.0,
    # config="my_config.json", # optional
)

# Access thermodynamic results
for mode in results.rank_by_free_energy():
    thermo = mode.get_thermodynamics()
    print(f"Mode {mode.rank}: F = {thermo.free_energy:.2f} kcal/mol")
```

### Alternative: Direct binary invocation

```bash
/path/to/build/FlexAIDδS receptor_clean.pdb ligand.mol2 \
    -o my_docking_run \
    -c my_optional_config.json
```

Then load results afterward with `flexaidds.load_results("my_docking_run")`.

## 4. Analyze Results + Thermodynamics

After docking:

```python
import flexaidds as fd

# Load results
docking = fd.load_results("my_docking_run")

# Thermodynamic analysis is already computed if you used compute_entropy=True
for mode in docking.binding_modes:
    print(f"Mode {mode.rank}")
    print(f"  CF score (proxy): {mode.score:.2f}")
    if mode.thermodynamics:
        t = mode.thermodynamics
        print(f"  Free energy F:     {t.free_energy:.2f} kcal/mol")
        print(f"  Enthalpy H:        {t.mean_energy:.2f}")
        print(f"  Entropy term -TS:  {t.entropy_term:.2f}")
        print(f"  Boltzmann weight:  {mode.boltzmann_weight:.4g}")
```

## 5. Common Professional Workflows

- **Self-docking validation** (e.g. biotin into streptavidin): Use crystal pose as reference for RMSD.
- **Ensemble thermodynamics**: Always set `compute_entropy=True` when you care about ranking by free energy rather than raw CF score.
- **Reproducibility**: Pin the exact FlexAIDδS binary + matrix files used. The skill's `ensure_docking_data.py --source` mechanism helps here.
- **CI / automated runs**: Use `--check` mode of the ensure script and run the validator before launching jobs.

## 6. Troubleshooting

| Problem                        | Likely Cause                              | Solution |
|--------------------------------|-------------------------------------------|----------|
| Binary complains about missing `MC_*.dat` | Matrices not next to binary | Run `ensure_docking_data.py` |
| No thermodynamic data in results | `compute_entropy` was not enabled | Re-run with the flag or Python kwarg |
| Poor ranking vs experiment     | Using raw CF score instead of free energy | Always use the thermodynamic ledger for final ranking |
| Binary not found               | Wrong path or not built with Metal/CPU support | Check your build and pass correct `--binary` |

## 7. Redocking from a Single PDB ID (New)

For the common case of redocking a cocrystallized ligand back into its target:

```bash
python3 .grok/skills/flexaid-docking/scripts/redock_from_pdb.py 1STP
```

This experimental helper will:
- Download the PDB
- Attempt to extract protein vs ligand
- Prepare inputs
- Run docking + basic thermodynamics

**Always inspect the generated receptor and ligand files.** The splitting is heuristic.

## 8. Next Steps

- Read the full skill documentation: `SKILL.md`
- See advanced guidance: `references/flexaid-docking-guidance.md`
- Explore the thermodynamic models: `python/flexaidds/thermodynamics.py`

The `flexaid-docking` skill is designed to make the above workflows safe, reproducible, and correctly scoped between scoring proxies and real statistical mechanics.

Run with confidence.