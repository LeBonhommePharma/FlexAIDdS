# ITC Thermodynamic Calibration — Methodology & Runbook

Make FlexAIDdS's enthalpy and entropy terms **quantitatively comparable to ITC**
(Isothermal Titration Calorimetry) by fitting two physically-interpretable
scaling factors, α (enthalpy) and β (entropy).

```
ΔH_pred  = α · (T_eff · H_vct_raw) = α · H_vct     → fit α against ΔH_ITC
TΔS_pred = β · TdS_vib                              → fit β against TΔS_ITC
ΔG_pred  = ΔH_pred − TΔS_pred                       → validated against ΔG_ITC (no extra fit)
```

Result: `T_eff_calibrated = α · T_eff_current`, `tds_vib_scale = β`.
Sign convention throughout: **ΔG = ΔH − TΔS** (TΔS is the entropy *contribution*).

---

## Components

| File | Role |
|------|------|
| `benchmarks/fetch_itc_data.py` | Unify ITC data from many sources → one schema; dedup; variance report. |
| `benchmarks/calibrate_itc.py` | Fit α/β from FlexAIDdS results vs ITC; per-source r/RMSE; optional plots. |
| `benchmarks/tests/test_itc_calibration.py` | Offline test suite (synthetic fixtures; no network, no docking). |
| `LIB/DatasetRunner.cpp` | Emits `g_bind,h_vct,h_vct_raw,n_heavy,tds_shannon,tds_vib` into the aggregate `*_results.csv` when `FLEXAIDDS_THERMO_CSV=1`. |

Unified schema (kcal/mol):
`pdb_id,ligand_smiles,dH_kcal_mol,TdS_kcal_mol,dG_kcal_mol,T_K,source,doi`

---

## Step 1 — Build the ITC reference table

```bash
# Local sources (SCORPIO curated CSV + BindingDB ITC TSV):
python benchmarks/fetch_itc_data.py --sources scorpio bindingdb \
    --out benchmarks/itc_unified.csv

# Calibration-ready subset (docking needs a structure; α needs ΔH):
python benchmarks/fetch_itc_data.py --sources scorpio bindingdb \
    --require-pdb --require-dh --out benchmarks/itc_calib.csv

# Pull ChEMBL ΔH activities over the network (free API):
python benchmarks/fetch_itc_data.py --sources chembl --download --chembl-limit 500 \
    --out /tmp/itc_chembl.csv

# See every source and how to obtain it (free/API vs manual paywalled SI):
python benchmarks/fetch_itc_data.py --list-sources
```

Sources: **scorpio**, **bindingdb** (local/http, working parsers), **pdbbind**
(`--pdbbind-index`, ΔG-only), **chembl** (`--download`), and
**csar/freire/sampl/biolip/nist** (transcribe to unified schema, load via
`--csv-files`). Duplicate measurements of a complex are all kept; per-complex
standard deviation is reported as experimental uncertainty.

Current local yield (SCORPIO + BindingDB): **944 rows** — 496 with PDB, 914 with
ΔH, 705 with SMILES; 69 multi-measurement complexes.

## Step 2 — Dock the calibration set

Dock each complex with the thermo hooks enabled so the decomposition is written:

```bash
export FLEXAIDDS_THERMO=1        # emit [THERMO] stdout line
export FLEXAIDDS_THERMO_CSV=1    # add thermo columns to aggregate *_results.csv
# ... run the usual DatasetRunner / benchmark harness (oracle-site mode is fine —
#     calibration is about the energy decomposition at the native pose) ...
```

`calibrate_itc.py` reads, per PDB id, whichever of these exists (in priority):
aggregate `*_results.csv`, per-target `result.csv`, per-target `stdout.log`
`[THERMO]` line, or `result.json`.

## Step 3 — Fit α and β

```bash
python benchmarks/calibrate_itc.py \
    --itc-csv benchmarks/itc_calib.csv --source unified \
    --results <results_dir> \
    --t-eff-current <T_eff used in the runs> \
    --plot-dir figures/itc \
    --out itc_calibration_report.json
```

`--source` also accepts native formats directly: `scorpio` (scorpio_itc_raw.csv),
`bindingdb` (ITC TSV), or `csv` (simple `pdb_id,dH_kcal_mol,TdS_kcal_mol,dG_kcal_mol`).

Output: fitted α, β, `T_eff_calibrated`, `tds_vib_scale`; Pearson r + RMSE for
ΔH / TΔS / ΔG **separately** and **per source** (so you can check that α, β fit on
one source generalize to another); ready-to-paste env lines; and per-source
3-panel PNGs (ΔH/TΔS/ΔG predicted-vs-ITC) when matplotlib is present.

### Recommended split
Fit on **BindingDB** (has ΔH/TΔS + SMILES + PDB), validate ΔG generalization on
**SCORPIO** (structure-prepped). ΔG uses the fitted α, β with no further tuning —
it is the falsifiable check that the decomposition is physical rather than two
independent rescalings. Watch for enthalpy–entropy compensation: if ΔG correlates
but ΔH and TΔS don't, the terms are compensating, not individually calibrated.

## Step 4 — Apply the calibration (future engine wiring)

```bash
FLEXAIDDS_THERMO_CALIBRATE=1
FLEXAIDDS_ALPHA_ENTHALPY=<α from step 3>
FLEXAIDDS_BETA_ENTROPY=<β from step 3>
```
(Reading these in the scorer is the follow-up step, after the G_bind fix lands.)

---

## Testing

Fully offline — no network, no docking, no numpy/matplotlib required:

```bash
pytest benchmarks/tests/test_itc_calibration.py -q      # 14 tests
python benchmarks/tests/test_itc_calibration.py         # same, no pytest needed
```

Covered: unit conversion (kJ↔kcal), BindingDB −TΔS sign handling, ΔG derivation,
ChEMBL parsing/filtering, PDBbind index, dedup variance, Pearson/RMSE,
through-origin scale fit, duplicate averaging, `[THERMO]` log + aggregate-CSV
extraction, an end-to-end fit that recovers a known α/β, and graceful plot
skipping when matplotlib is absent. Wired into CI's `pure_python_results` job.

## Notes / caveats

- **Units:** unified/scorpio are kcal/mol; bindingdb kJ→kcal auto-converted;
  `--source csv` honours `--itc-units {kcal,kJ}`.
- **SCORPIO has no SMILES** (ligand names only) → its dedup key is `(pdb_id, "")`.
  Fine for the pdb_id-keyed calibration join; map names→SMILES if you need
  SMILES-level dedup across sources.
- **Astex Diverse / PDBbind index** carry Kd/ΔG only — insufficient for α/β
  (ΔH/TΔS); usable for ΔG validation.
- **Temperature:** prefer 298 K; correct or exclude 37 °C entries.
