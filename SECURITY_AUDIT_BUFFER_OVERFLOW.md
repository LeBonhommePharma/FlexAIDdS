# Buffer Overflow Security Audit Report

**Date:** 2026-03-14
**Scope:** Full codebase — C++, CUDA, Metal, Objective-C++, Python bindings
**Auditor:** Automated static analysis

> Status note, May 2026: this document is a historical audit input, not a live
> closure report. Some findings have since been fixed or partially bounded.
> Current status belongs in `docs/KNOWN_LIMITATIONS.md` and source-level tests.

---

## Executive Summary

The FlexAIDdS codebase contains **widespread use of unsafe C string functions** (`sprintf`, `strcpy`, `strcat`, `sscanf`) across the core C++ library. Most buffer overflow risks are **MEDIUM** severity because input is typically from trusted configuration files rather than untrusted user input. However, several patterns present **HIGH** risk if input validation assumptions change.

**Total findings: 28 distinct vulnerability patterns across 25+ files.**

---

## CRITICAL Findings (0)

No critical findings. The codebase does not expose network services or parse untrusted binary formats.

---

## HIGH Severity Findings (7)

### H-1: `modify_pdb.cpp:29,86,102,112` — Stack buffer overflow via line accumulation

```cpp
char lines[50][100]; // store residue lines  (line 29)
strcpy(lines[nlines++], buffer);             (lines 86, 102, 112)
```

**Risk:** `nlines` is incremented without bounds checking against the hard limit of 50. A PDB residue with >50 ATOM lines (e.g., a modified amino acid or unusual ligand labeled as ATOM) overflows the `lines` array.

**Fix:** Add `if (nlines >= 50) { fprintf(stderr, "..."); continue; }` before each `strcpy(lines[nlines++], ...)`.

---

### H-2: `modify_pdb.cpp:190-192` — `strncpy` + `sprintf` + `strcat` on uninitialized `newline`

```cpp
char newline[100];
strncpy(newline, lines[i], 6);        // copies 6 bytes, does NOT null-terminate
sprintf(&newline[6], "%5d", ++(*wrote)); // writes at offset 6
strcat(&newline[11], &lines[i][11]);  // strcat needs null-terminated destination
```

**Risk:** `strncpy(dst, src, 6)` does not null-terminate if `src >= 6` chars. `sprintf` writes at offset 6, producing `"ATOM  12345"` (11 bytes). Then `strcat(&newline[11], ...)` appends the rest of the PDB line. If `lines[i]` is 99+ chars, the result exceeds the 100-byte `newline` buffer.

**Fix:** Use `snprintf` with explicit length. Initialize `newline` to zeros.

---

### H-3: `read_input.cpp:85-132` — Unbounded `strcpy` from parsed config lines

```cpp
if(strcmp(field,"PDBNAM") == 0){strcpy(pdb_name, &buffer[7]);}
if(strcmp(field,"INPLIG") == 0){strcpy(lig_file, &buffer[7]);}
if(strcmp(field,"DEPSPA") == 0){strcpy(FA->dependencies_path, &buffer[7]);}
if(strcmp(field,"STATEP") == 0){strcpy(FA->state_path, &buffer[7]);}
if(strcmp(field,"TEMPOP") == 0){strcpy(FA->temp_path, &buffer[7]);}
```

**Risk:** `buffer` is `char[MAX_PATH__*2]` (510 bytes). The destination buffers (`pdb_name`, `FA->state_path`, etc.) are `char[MAX_PATH__]` (255 bytes). A crafted config file with a path >254 chars overflows the destination.

**Fix:** Use `strncpy(pdb_name, &buffer[7], MAX_PATH__ - 1); pdb_name[MAX_PATH__ - 1] = '\0';`

---

### H-4: `read_input.cpp:97-98` — Unbounded array indexing for `optline` / `flexscline`

```cpp
char optline[MAX_PAR][MAX_PATH__];    // MAX_PAR = 100
char flexscline[MAX_PAR][MAX_PATH__];
if(strcmp(field,"OPTIMZ") == 0){strcpy(optline[nopt++], buffer);}
if(strcmp(field,"FLEXSC") == 0){strcpy(flexscline[nflexsc++], buffer);}
```

**Risk:** If the config file contains >100 OPTIMZ or FLEXSC lines, `nopt` or `nflexsc` exceeds `MAX_PAR`, causing a stack buffer overflow. The bounds check at line 529 (`if(i==MAX_PAR)`) is in the wrong loop (it checks during processing, not during reading).

**Fix:** Add `if (nopt >= MAX_PAR) continue;` and `if (nflexsc >= MAX_PAR) continue;` before the `strcpy` calls.

---

### H-5: `read_input.cpp:87-91` — `sscanf %s` into tiny fixed buffers

```cpp
sscanf(buffer, "%s %s", a, FA->metopt);    // FA->metopt is char[3]
sscanf(buffer, "%s %s", a, FA->bpkenm);    // FA->bpkenm is char[3]
sscanf(buffer, "%s %s", a, FA->complf);    // FA->complf is char[4]
sscanf(buffer, "%s %s", a, FA->vcontacts_self_consistency); // char[6]
```

**Risk:** `sscanf %s` writes until whitespace, with no length limit. A value like `"METOPT GENETIC_ALGORITHM"` writes 18 bytes into a 3-byte buffer.

**Fix:** Use width-limited format specifiers: `sscanf(buffer, "%s %2s", a, FA->metopt);`

---

### H-6: `gaboom.cpp:1804,1809,1812,1814` — `sscanf %s` into small GA config buffers

```cpp
sscanf(buffer, "%s %s", field, GB->pop_init_method);  // char[9]
strcpy(GB->pop_init_file, &buffer[16]);                // char[MAX_PATH__]
sscanf(buffer, "%s %s", field, GB->fitness_model);     // char[9]
sscanf(buffer, "%s %s", field, GB->rep_model);         // char[9]
```

**Risk:** `pop_init_method`, `fitness_model`, `rep_model` are `char[9]`. Any value >=9 chars overflows. `pop_init_file` uses `strcpy` from a fixed offset which can overflow if the line format deviates.

**Fix:** Use `sscanf(buffer, "%s %8s", field, GB->pop_init_method);` with width limits.

---

### H-7: `gaboom.cpp:76-95` — Path concatenation without overflow check

```cpp
char PAUSEFILE[MAX_PATH__];   // 255 bytes
strcpy(PAUSEFILE, FA->state_path);
strcat(PAUSEFILE, "/.pause");
```

**Risk:** If `FA->state_path` is close to 255 chars, appending `"/.pause"` (7 bytes) overflows `PAUSEFILE`. Same pattern for `ABORTFILE`, `STOPFILE`, `UPDATEFILE`.

**Fix:** Use `snprintf(PAUSEFILE, MAX_PATH__, "%s/.pause", FA->state_path);`

---

## MEDIUM Severity Findings (14)

### M-1: `DensityPeak_Cluster.cpp:43-44,404-526` — Repeated `strcat` into `remark[MAX_REMARK]`

`remark` is `char[5000]`. Multiple `sprintf` + `strcat` calls append REMARK lines. Each tmpremark can be ~120 bytes. With many optimizable residues, the total can approach 5000 bytes. No running-length check.

**Fix:** Track `remark` length and use `snprintf` / bounds check before each `strcat`.

---

### M-2: `BindingMode.cpp:261-332,344-404` — Same `remark` accumulation pattern

Same as M-1. Two separate functions both build `remark[MAX_REMARK]` via repeated `strcat`.

---

### M-3: `cluster.cpp:23-260` — Same `remark` accumulation pattern

Same pattern as M-1/M-2.

---

### M-4: `FOPTICS.cpp:244-306` — Same `remark` accumulation pattern

Same pattern.

---

### M-5: `top.cpp:31-510` — Same `remark` accumulation pattern

Same pattern.

---

### M-6: `read_input.cpp:182-307` — Multiple `strcpy`/`strcat` path constructions

Path concatenation throughout the function (e.g., `tmpprotname`, `emat`, `deftyp`, `rotlib_file`, `rotobs_file`) uses `strcpy`/`strcat` without checking combined length against `MAX_PATH__`.

---

### M-7: `gaboom.cpp:300-321` — `sprintf` into `gridfilename` then `strcat`

```cpp
sprintf(gridfilename, "/grid.%d.prt.pdb", i+1);
strcpy(gridfile, FA->temp_path);
strcat(gridfile, gridfilename);
```

If `temp_path` + filename > 255, overflow occurs.

---

### M-8: `read_coor.cpp:64` — `strcpy` into `atom.name[5]`

```cpp
strcpy((*atoms)[FA->atm_cnt].name, atm_typ);  // name is char[5]
```

`atm_typ` is `char[5]`, extracted via `strncpy` from PDB columns. Safe only if the PDB format is well-formed.

---

### M-9: `read_lig.cpp:140` — `strcpy` into `residue.name[4]`

```cpp
strcpy((*residue)[FA->res_cnt].name, rnam);  // name is char[4]
```

`rnam` is `char[4]` from PDB parsing. Off-by-one risk if null terminator not set.

---

### M-10: `read_flexscfile.cpp:64` — `strcpy` into `flex_res.name[4]`

```cpp
strcpy(FA->flex_res[FA->nflxsc].name, resname);
```

`resname` comes from `sscanf %s` which is unbounded.

---

### M-11: `shortest_path.cpp:50,82` — `strcpy` from `std::string` into `shortpath`

```cpp
strcpy(residue->shortpath[i][j], ss.str().c_str());
```

There is a length check at line 78 (`MAX_SHORTEST_PATH*6`), but `strcpy` at line 50 has no check for the simpler `i==j` case.

---

### M-12: `metal_eval.mm:278,299-302` — Metal gene buffer hardcoded to 256 genes

```cpp
ctx->buf_genes_f = [ctx->device newBufferWithLength:(size_t)max_pop * 256 * sizeof(float) ...];
// Later:
for (int g = 0; g < n_genes; ++g)
    genes_f[c * n_genes + g] = (float)h_genes[c * n_genes + g];
```

Buffer allocated for 256 genes per chromosome. If `n_genes > 256`, the buffer is too small, but the write loop uses `n_genes` (unbounded). Also, the write pattern `c * n_genes` differs from the allocation pattern `c * 256`, which could cause out-of-bounds writes when `n_genes != 256`.

**Fix:** Allocate `max_pop * max_genes * sizeof(float)` and use consistent indexing.

---

### M-13: `cuda_eval.cu:276-279` — `pop_size` not validated against `max_pop`

```cpp
const size_t gene_bytes = (size_t)pop_size * n_genes * sizeof(double);
CUDA_CHECK(cudaMemcpy(ctx->d_genes, h_genes, gene_bytes, cudaMemcpyHostToDevice));
```

If `pop_size > ctx->max_pop`, the `cudaMemcpy` writes past the allocated device buffer.

**Fix:** Add `assert(pop_size <= ctx->max_pop);` or validation.

---

### M-14: `maps.cpp:12` — `sprintf` into `char_coor[9]`

```cpp
char char_coor[9];
sprintf(char_coor, "%8.3f", coor[i]);
```

`%8.3f` can produce >8 chars for large values (e.g., `-12345.678` is 10 chars). Overflow.

**Fix:** Use `snprintf(char_coor, sizeof(char_coor), "%8.3f", coor[i]);`

---

## LOW Severity Findings (7)

### L-1: `read_rotlib.cpp:81` — `sscanf %s` into `rotamer.name[9]`

Rotamer name from library file. Low risk as library files are trusted.

### L-2: `read_conect.cpp:13` — `char number[6]` with `sscanf %d`

PDB CONECT records. Standard format limits this to 5 digits.

### L-3: `SdfReader.cpp:75` — `sscanf %63s` into `mol_name[64]`

Width-limited. Safe.

### L-4: `Mol2Reader.cpp:146` — `sscanf %63s` into `mol_name[64]`

Width-limited. Safe.

### L-5: `tencom_entropy_diff.cpp:80` — `strncpy` into `ca.resname[3]`

Uses `substr(17, 3)`. Length bounded by substr.

### L-6: `python/bindings/core_bindings.cpp:67,197,298` — `snprintf` into `char buf[256]`

Uses `snprintf(buf, sizeof(buf), ...)`. Safe.

### L-7: `benchmark_tencom.cpp:45-54` — `strncpy` with `sizeof(ca.name) - 1`

Uses safe pattern with explicit null terminator.

---

## Recommendations

### Immediate Actions (HIGH findings)

1. **Replace all `sprintf` with `snprintf`** across the codebase
2. **Add bounds checks** in `modify_pdb.cpp` for `nlines` counter
3. **Add width specifiers to `sscanf %s`** formats matching destination buffer sizes
4. **Validate `nopt`/`nflexsc` counters** before writing to `optline`/`flexscline` arrays
5. **Use `snprintf` for all path concatenation** instead of `strcpy`/`strcat` chains

### Systematic Improvements

6. **Replace `strcpy`/`strcat` chains** with `snprintf` throughout path construction code
7. **Add running length tracking** to `remark` buffer accumulation in cluster/binding mode output
8. **Validate `pop_size <= max_pop`** in CUDA/Metal batch evaluation entry points
9. **Fix Metal gene buffer sizing** to use `max_genes` instead of hardcoded 256
10. **Consider `std::string`** for path and remark buffer management in new code

### Files Requiring Most Attention

| File | # of Issues | Priority |
|------|-------------|----------|
| `LIB/read_input.cpp` | 5 | HIGH |
| `LIB/modify_pdb.cpp` | 2 | HIGH |
| `LIB/gaboom.cpp` | 3 | HIGH |
| `LIB/DensityPeak_Cluster.cpp` | 1 | MEDIUM |
| `LIB/BindingMode.cpp` | 1 | MEDIUM |
| `LIB/cluster.cpp` | 1 | MEDIUM |
| `LIB/FOPTICS.cpp` | 1 | MEDIUM |
| `LIB/top.cpp` | 1 | MEDIUM |
| `LIB/metal_eval.mm` | 1 | MEDIUM |
| `LIB/cuda_eval.cu` | 1 | MEDIUM |
| `LIB/maps.cpp` | 1 | MEDIUM |

---

## Scope Notes

- **Python code** (`python/flexaidds/`): No buffer overflow risks. Python handles memory management automatically.
- **pybind11 bindings** (`core_bindings.cpp`): Uses `snprintf` with `sizeof(buf)` consistently. Safe.
- **`_core.cpp`**: Uses `strncpy` with explicit sizes. Safe.
- **ShannonThermoStack**: Uses `std::vector` throughout. No raw buffer issues.
- **Test files**: Some use `strncpy` for test fixtures; acceptable in test context.
- **GPU kernels (CUDA/Metal)**: Kernel-internal memory is bounded by shared memory declarations. Host-side buffer management has validation gaps (M-12, M-13).
