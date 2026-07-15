# FastOPTICS MinPts — literature validation

**Status:** implementation in `LIB/FastOPTICS_cluster.cpp`, constants in `LIB/ga_constants.h`.  
**Log tag:** `[FO-MINPTS]`  
**Production rule:** **one** FastOPTICS + BindingPopulation pass. The old triple MinPts ladder (×1.5 scales) was testing-only in the legacy repo and is **not** used.

## Peer-reviewed sources

| Source | Venue | MinPts guidance used |
|--------|--------|----------------------|
| **Ankerst, Breunig, Kriegel, Sander** — *OPTICS* | ACM SIGMOD Record **28**(2):49–60 (1999), [doi:10.1145/304181.304187](https://doi.org/10.1145/304181.304187) | MinPts is the primary density parameter. Experiments: good results for **MinPts ∈ [10, 20]**. Larger MinPts smooths the reachability plot and weakens single-link chaining. |
| **Ester et al.** — DBSCAN | KDD 1996 | Default **MinPts = 4** for 2-D → floor when N allows. |
| **Sander et al.** — GDBSCAN | Data Min. Knowl. Disc. **2**:169–194 (1998) | For **dim > 2**, **MinPts ≈ 2 · dim**. |

## Mapping into FlexAIDdS (single pass)

1. **Effective dimension** (Sander): `dim_eff = clamp(6 + fdih, 2, 20)` with `fdih` from `FA->resligand->fdih`; else IC gene count.
2. **One MinPts** (Sander + Ankerst + Ester + feasibility): `sander = 2·dim_eff`, prefer Ankerst **[10, 20]** when `nChrom ≥ 20`, else Ester floor / `nChrom/3` cap.
3. **Diversity softener only** if CF diversity &lt; 5% (toward Ester 4).
4. **No multi-scale re-run.** Super-cluster energy pre-filter and minibatch sampling (if enabled) are *not* additional FO pose clusterings.

## Runtime check

```text
[FO-MINPTS] literature=Ankerst1999[10-20]+Sander1998(2*dim)+Ester1996(floor4) ... minPts=N (single FO pass; Ankerst band [10,20])
Size of Population is K Binding Modes (minPts=N).
-- end of FastOPTICS_cluster --
```

Exactly **one** “Size of Population” line per clustering call.
