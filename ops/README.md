# ops/ — Merge Gate and Scoring Guardrails

## CF Scoring-Regression Gate

**Script**: `ops/gates/cf_gate_probe_cf.sh`  
**Manifest**: `ops/gates/panel_manifest.tsv`  
**Instrument**: `build/probe_cf` (built from `tools/probe_cf.cpp`, commit 81acfed6+)

The gate scores the native crystal pose and best decoy pose for each panel target
via the real engine CF path (`vcfunction()` / `score_native_pose()`), computes
`ΔCF = cf_total(native) − cf_total(decoy)`, and fails if the inverted fraction
(ΔCF > tol) exceeds MAX_INV_FRAC (default 1/8).

### Running the gate

```bash
# Build probe_cf if stale
cmake --build build --target probe_cf

# Run the gate
bash ops/gates/cf_gate_probe_cf.sh ops/gates/panel_manifest.tsv
# Exit 0 = PASS, 1 = FAIL (inverted fraction too high), 2 = no data
```

### CI trigger paths

The gate runs on any PR that touches:
```
LIB/vcfunction.cpp
LIB/top.cpp
LIB/read_input.cpp
LIB/soft_wall.h
LIB/ic2cf.cpp
LIB/ProcessLigand/SybylTyper.cpp
data/MC_st0r5.2_6.dat
*.dat
```

### Panel manifest

`ops/gates/panel_manifest.tsv` is tab-separated: `pdb <TAB> receptor.pdb <TAB> native_pose.sdf <TAB> decoy_pose.pdb`

Commented rows (`#`) have no archived decoy yet — populate from the next full
benchmark run output (`best-CF non-native pose per target`).

---

## SCORING_PROVENANCE.json — Required on all scorer PRs

**Any PR that modifies scoring sources** (matrix files, `sas_weight`, hbond flags,
`r0`, soft-wall cutoff, WAL cap logic, or `vcfunction.cpp`) **must include** a
populated `SCORING_PROVENANCE.json` file in the PR root.

Use `ops/SCORING_PROVENANCE_template.json` as the starting point:

```bash
cp ops/SCORING_PROVENANCE_template.json SCORING_PROVENANCE.json
# Fill in all fields, then git add SCORING_PROVENANCE.json
```

### Required fields

| Field | Description |
|-------|-------------|
| `matrix_md5` | MD5 of the `.dat` scoring matrix in use |
| `matrix_file` | Path to the matrix file relative to repo root |
| `sas_weight` | Current `sas_weight` value in `read_input.cpp` / config |
| `hbond_flags` | Object with all active hbond flag names → values |
| `r0` | Voronoi contact `r0` parameter |
| `soft_wall_cutoff` | Soft-wall energy cutoff value |
| `wal_cap_state` | `"on"` / `"off"` / `"conditional:<condition>"` |
| `git_commit` | Full SHA of the commit being reviewed |
| `binary_sha256` | SHA-256 of `build/FlexAIDdS` or `build/probe_cf` used for validation |
| `timestamp` | ISO-8601 timestamp of when provenance was recorded |

### Rationale

The df2f36c58 regression landed because a `.dat` matrix content-swap was
invisible to git review (the file had `skip-worktree` set). SCORING_PROVENANCE
ensures the exact matrix MD5 and all scoring parameters are on record for every
scorer PR, making silent parameter changes impossible to miss in review.

---

## Pre-commit hook

`.git/hooks/pre-commit` rejects any commit where a tracked file under `LIB/`
or matching `*.dat` has `assume-unchanged` or `skip-worktree` set via
`git update-index`. This is enforced locally; CI enforces the gate independently.
