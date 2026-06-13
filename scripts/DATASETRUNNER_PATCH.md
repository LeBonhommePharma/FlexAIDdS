# DatasetRunner Integration Patch — Rotamer Pre-Relaxation

## Files to modify

1. `LIB/DatasetRunner.h`  — add field to `DockingConfig`
2. `LIB/DatasetRunner.cpp` — add `#include`, insert hook before cmd builder
3. `CMakeLists.txt`        — add `LIB/receptor_prep.cpp` to source list

---

## Patch 1 — DatasetRunner.h : add field to DockingConfig (line ~226)

Find the end of the `DockingConfig` struct (after `bool force_rigid{false};`),
insert ONE new field before the closing `};`:

```cpp
    /// Binding-site rotamer pre-relaxation (Option 3 — apo-strain fix).
    /// Greedy Dunbrack rotamer search on pocket sidechains before docking.
    /// Default false: behaviour unchanged from current baseline.
    /// Enable in v43+ DatasetRunner injection for CF_false_minimum targets.
    bool   receptor_rotamer_prep{false};
};
```

**Exact insertion point** (after line 224 `bool   force_rigid{false};`):
```
224:    bool   force_rigid{false};
225:    // ↓↓↓ INSERT HERE ↓↓↓
226:    bool   receptor_rotamer_prep{false};  // Option 3 apo-strain fix
227:};
```

---

## Patch 2 — DatasetRunner.cpp : add #include near top of file

Near the other LIB/ includes at the top of DatasetRunner.cpp, add:

```cpp
#include "receptor_prep.h"
```

---

## Patch 3 — DatasetRunner.cpp : insert prep hook in run() (before cmd builder)

**Context**: In `DatasetRunner::run()`, find the block that builds the docking
command (around line 4888):

```cpp
            std::ostringstream cmd;
            // Oracle LOCCLF: pass binding site PDB via env var so top.cpp
            // can skip SURFNET auto-detection and load the oracle spheres.
            if (!entry.binding_site_path.empty() && fs::exists(entry.binding_site_path)) {
                cmd << "FLEXAIDDS_ORACLE_SITE='" << entry.binding_site_path << "' ";
```

**Insert immediately before** `std::ostringstream cmd;`:

```cpp
            // ── Rotamer pre-relaxation (Option 3: apo-strain fix) ──────────
            // Pre-relax binding-site sidechains before docking to reduce
            // CF.wal false penalties on near-native poses in apo structures.
            // Gated on config.receptor_rotamer_prep (default false).
            std::string effective_receptor = entry.receptor_path;
            if (config.receptor_rotamer_prep &&
                !entry.binding_site_path.empty() &&
                fs::exists(entry.binding_site_path))
            {
                std::string prepped = out_dir + "/" + entry.pdb_id + "_prepped.pdb";
                bool need_prep =
                    !fs::exists(prepped) ||
                    fs::last_write_time(prepped) <
                        fs::last_write_time(entry.receptor_path) ||
                    fs::last_write_time(prepped) <
                        fs::last_write_time(entry.binding_site_path);

                if (need_prep) {
                    int n_mod = receptor_prep::prep_receptor_rotamers(
                        entry.receptor_path,
                        entry.binding_site_path,
                        prepped);
                    if (n_mod >= 0) {
                        effective_receptor = prepped;
                        std::cerr << "  [PREP] " << entry.pdb_id
                                  << ": rotamer-prepped " << n_mod
                                  << " pocket residues → " << prepped << "\n";
                    } else {
                        std::cerr << "  [PREP-WARN] " << entry.pdb_id
                                  << ": prep_receptor_rotamers() failed, "
                                     "using unmodified apo receptor\n";
                    }
                } else {
                    effective_receptor = prepped;
                    std::cerr << "  [PREP] " << entry.pdb_id
                              << ": using cached rotamer-prepped receptor\n";
                }
            }
            // ── end rotamer pre-relaxation ──────────────────────────────────
```

Then, ONE line further down (currently line 4906):

```cpp
                << "'" << entry.receptor_path << "' "
```

Change to:

```cpp
                << "'" << effective_receptor << "' "
```

That is the **only** change to the existing command-builder logic.

---

## Patch 4 — CMakeLists.txt : add receptor_prep.cpp to FLEXAID_SOURCES

Find the `FLEXAID_SOURCES` list (the long set(...) block).  It already contains
entries like:

```cmake
    LIB/build_rotamers.cpp
    LIB/read_rotlib.cpp
    LIB/read_rotobs.cpp
```

Add adjacent to those (anywhere in the list is fine):

```cmake
    LIB/receptor_prep.cpp
```

Also add `receptor_prep.cpp` to the `flexaid_core` library target if it exists
separately in the CMakeLists (check for `add_library(flexaid_core ...)`).

---

## Patch 5 — Provenance logging (optional but recommended)

In the existing dock_config.json writer (around line 4656), add a line to record
whether prep was active:

```cpp
               << "    \"receptor_rotamer_prep\": "
               << (config.receptor_rotamer_prep ? "true" : "false") << ",\n"
```

Insert adjacent to the `"force_rigid"` line.

---

## Enabling in v43 DatasetRunner injection

In the benchmark runner call site (wherever `DockingConfig config;` is constructed
before calling `runner.run(entries, config)`), set:

```cpp
config.receptor_rotamer_prep = true;   // enable Option 3
```

Or in whatever YAML/JSON deserialisation reads the benchmark config:

```yaml
receptor_rotamer_prep: true
```

---

## Summary of changed/added files

| File | Type | Change |
|------|------|--------|
| `LIB/receptor_prep.h`      | **NEW** | Public API header |
| `LIB/receptor_prep.cpp`    | **NEW** | Full implementation (~320 lines) |
| `LIB/DatasetRunner.h`      | **MODIFY** | +1 field to DockingConfig |
| `LIB/DatasetRunner.cpp`    | **MODIFY** | +1 `#include`, +30 lines in run(), 1 line changed |
| `CMakeLists.txt`           | **MODIFY** | +1 source file to FLEXAID_SOURCES |
| `scripts/prep_receptor.py` | Reference  | Python equivalent (kept for standalone testing) |
