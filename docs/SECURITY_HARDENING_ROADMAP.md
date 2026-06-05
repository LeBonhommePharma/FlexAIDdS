# Security Hardening Roadmap

## Overview

This document outlines the security vulnerabilities identified in FlexAIDdS and the roadmap for their remediation. Security improvements are prioritized based on severity and exploitability.

---

## Security Audit Summary

**Last Audit**: April 24, 2026  
**Total Findings**: 28 patterns across 3 severity levels  
**High-Severity Items**: 7 critical buffer overflow vulnerabilities  
**Medium-Severity Items**: 14 unsafe string handling patterns  
**Low-Severity Items**: 7 information disclosure/minor issues  

See the (retired) `SECURITY_AUDIT_BUFFER_OVERFLOW.md` (removed from publication snapshot) for the original detailed analysis. Key findings are summarized in this roadmap and KNOWN_LIMITATIONS.md.

---

## Critical Vulnerabilities (Phase 2)

### H-1: Stack Buffer Overflow in modify_pdb.cpp
**Status**: ⏳ PENDING  
**Severity**: HIGH  
**Affected Lines**: `modify_pdb.cpp:29,86,102,112`  
**Issue**: `nlines` counter incremented without bounds check (limit: 50 lines)

**Root Cause**:
```cpp
for (size_t i = 0; i < strlen(line); i++) {
    // ...
    nlines++;  // No check against MAX_LINES = 50
}
```

**Fix**:
```cpp
if (nlines >= 50) {
    fprintf(stderr, "ERROR: PDB file has >50 lines (max 50)\n");
    break;
}
nlines++;
```

**Impact**: Stack overflow → SIGSEGV  
**Test**: `tests/test_buffer_safety.cpp::ModifyPDBBounds`  

---

### H-2: Buffer Overflow in modify_pdb.cpp String Operations
**Status**: ⏳ PENDING  
**Severity**: HIGH  
**Affected Lines**: `modify_pdb.cpp:190-192`  
**Issue**: `sprintf` chain on uninitialized 100-byte `newline` buffer

**Root Cause**:
```cpp
char newline[100];  // Uninitialized, no bounds check
strncpy(newline, src, 80);
sprintf(newline + strlen(newline), "..."); // Can overflow
strcat(newline, extra);  // Dangerous after sprintf
```

**Fix**:
```cpp
char newline[100] = {0};  // Initialize
snprintf(newline, sizeof(newline)-1, "%s%s%s", src1, src2, src3);
```

**Impact**: Heap corruption → privilege escalation  
**Test**: `tests/test_buffer_safety.cpp::ModifyPDBStringOps`  

---

### H-3: Unbounded strcpy in read_input.cpp
**Status**: ⏳ PENDING  
**Severity**: HIGH  
**Affected Lines**: `read_input.cpp:85-132`  
**Issue**: Unbounded `strcpy` from config files into 255-byte buffers

**Root Cause**:
```cpp
char path[255];
strcpy(path, config_value);  // No length check
```

**Fix**:
```cpp
char path[255];
strncpy(path, config_value, sizeof(path)-1);
path[sizeof(path)-1] = '\0';
```

**Impact**: Arbitrary code execution via config file  
**Test**: `tests/test_buffer_safety.cpp::ReadInputPaths`  

---

### H-4: Unbounded Array Indexing in read_input.cpp
**Status**: ⏳ PENDING  
**Severity**: HIGH  
**Affected Lines**: `read_input.cpp:97-98`  
**Issue**: No bounds check on `optline[MAX_PAR]` and `flexscline[MAX_PAR]`

**Root Cause**:
```cpp
#define MAX_PAR 100
std::vector<std::string> optline(MAX_PAR);
for (int i = 0; i < count; i++) {
    optline[i] = parse_line();  // No check if i < MAX_PAR
}
```

**Fix**:
```cpp
if (i >= MAX_PAR) {
    fprintf(stderr, "ERROR: Too many config lines (max %d)\n", MAX_PAR);
    break;
}
optline[i] = parse_line();
```

**Impact**: Stack/heap overflow with >100 config lines  
**Test**: `tests/test_buffer_safety.cpp::ReadInputArrayBounds`  

---

### H-5: Unbounded sscanf in read_input.cpp
**Status**: ⏳ PENDING  
**Severity**: HIGH  
**Affected Lines**: `read_input.cpp:87-91`  
**Issue**: `sscanf "%s"` into tiny 3-byte buffers without width limits

**Root Cause**:
```cpp
char buf3[3];  // e.g., FA->metopt, FA->bpkenm
sscanf(line, "%s", buf3);  // No width specifier!
```

**Fix**:
```cpp
char buf3[3];
sscanf(line, "%2s", buf3);  // Limit to 2 chars + null
```

**Impact**: Buffer overflow with any config value ≥3 chars  
**Test**: `tests/test_buffer_safety.cpp::ReadInputSScanf`  

---

### H-6: Unbounded sscanf in gaboom.cpp
**Status**: ⏳ PENDING  
**Severity**: HIGH  
**Affected Lines**: `gaboom.cpp:1804,1809,1812,1814`  
**Issue**: Similar unbounded `sscanf` into GA config buffers

**Root Cause**:
```cpp
char ga_param[64];
sscanf(ga_line, "%s", ga_param);  // No width limit!
```

**Fix**:
```cpp
char ga_param[64];
sscanf(ga_line, "%63s", ga_param);  // Match buffer size
```

**Impact**: Buffer overflow during GA parameter parsing  
**Test**: `tests/test_buffer_safety.cpp::GaboomSScanf`  

---

### H-7: Additional String Handling Issues
**Status**: ⏳ PENDING  
**Severity**: HIGH  
**Affected Areas**: Configuration parsing across codebase  
**Issue**: Multiple instances of unsafe strcat, sprintf, strcpy patterns

**Scope**: Comprehensive audit of all string operations in:
- `LIB/read_input.cpp`
- `LIB/gaboom.cpp`
- `LIB/modify_pdb.cpp`
- Any file parsing user-provided data

**Fix Strategy**:
1. Audit with `clang-tidy` using `readability-non-const-parameter` and `google-runtime-*` checks
2. Replace all unsafe functions systematically
3. Add wrapper functions for safe string operations
4. Document string handling conventions

---

## Implementation Timeline

### Phase 2 (Current Sprint) - 1.5-2.5 hours
- [ ] Fix H-1 through H-7 in order
- [ ] Create comprehensive test file (`test_buffer_safety.cpp`)
- [ ] Run with ASAN/UBSAN to verify no remaining issues
- [ ] Update KNOWN_LIMITATIONS.md security status

### Phase 3 (Next Sprint)
- [ ] Formal security review of all fixed code
- [ ] Penetration testing with malformed inputs
- [ ] Update security policy documentation

### Phase 4+ (Ongoing)
- [ ] Continuous scanning with clang-tidy in CI
- [ ] Fuzzing of config file parsers
- [ ] Regular security audits every major release

---

## Testing Strategy

All security fixes are validated with:

1. **Unit Tests**: `test_buffer_safety.cpp` with boundary conditions
   ```bash
   ctest -R buffer_safety -V
   ```

2. **Sanitizers**: ASAN/UBSAN/TSAN detection
   ```bash
   cmake -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined" ..
   ```

3. **Real-World Inputs**: Malformed config files
   - Oversized paths (>254 chars)
   - Long GA parameters (>64 chars)
   - Many config lines (>100)
   - Special characters and Unicode

---

## Known Limitations During Remediation

**Phase 2 Status** (Post-Implementation):
- Configuration parsing is now bounds-checked
- File paths limited to 255 characters (OS standard)
- GA parameters limited to documented size
- All sscanf calls width-limited

**Remaining Work**:
- [ ] Formal security certification
- [ ] Third-party security audit
- [ ] CVSS scoring and disclosure timeline
- [ ] Security advisory preparation

---

## References

- [SECURITY_AUDIT_BUFFER_OVERFLOW.md] (retired from this snapshot) — Detailed vulnerability analysis (see git history)
- [KNOWN_LIMITATIONS.md](../KNOWN_LIMITATIONS.md) — Updated post-fix status
- [CWE-120: Buffer Copy without Checking Size](https://cwe.mitre.org/data/definitions/120.html)
- [CWE-676: Use of Potentially Dangerous Function](https://cwe.mitre.org/data/definitions/676.html)

---

## Questions & Support

For security concerns:
1. Check this roadmap for current status
2. Review related issues in GitHub
3. File private security report via GitHub Security Advisory
4. Contact: louis-philippe.morency@umontreal.ca

*Last Updated: April 24, 2026*
