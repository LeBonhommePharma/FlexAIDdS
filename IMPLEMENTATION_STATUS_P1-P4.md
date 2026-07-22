# FlexAIDdS priority work orders — implementation status (P1–P4)

**Branch:** `feat/cluster-rep-medoid` off HEAD `3e674479c` (no merge to main).
**Session:** Claude (Cowork) executing the consolidated handoff.
**Landed here:** P3 (`d5c33b349`), P2 (`89d8dcd3e`).
**Delivered as plans (not committed):** P1, P4 — see §3/§4 for why and the exact recipe.

---

## 0. What actually happened vs. the handoff

The handoff referenced four detailed work orders as Cowork artifacts. Only **P2**
(`WORKORDER_clustering_medoid.md`) survived to disk; the P3/P1/P4 artifact links
were dead on this device (no artifacts present). So:

- **P2** was implemented to its full written spec, with every anchor re-verified
  against live source at `3e674479c`.
- **P3** was reconstructed from the handoff summary + live source (small, ~15 LOC,
  low ambiguity) and is flagged as a reconstruction in its commit + code comment.
- **P1 and P4** are delivered as **anchor-grounded implementation plans**, not code.
  Two independent reasons: (a) their detailed work orders are gone — writing 250–550
  LOC of invented integration surface would be exactly the fabrication the swarm's
  rules forbid; (b) both are **gated on measurement this environment cannot run** —
  the P2 work order itself states a single 2000-gen dock OOMs the host and the
  committer must *not* run the benchmark, and P1's whole payoff must be *measured*
  on the P2-enriched pool. Unbuilt, unmeasured code committed as "done" would be
  worse than an honest plan.

**Environment limits (stated plainly):** no `cmake` in PATH and the native build is
macOS; this bridge is a Linux VM. I ran `g++ -std=c++20 -fsyntax-only` against the
real include tree (all three modified TUs pass) but did **not** produce a working
binary or run any dock/benchmark. Bit-identity and RMSD gates are OPS-owned exactly
as the P2 work order prescribes.

---

## 1. P3 — `FLEXAIDDS_COM_FLOOR` soft lower clamp  ·  commit `d5c33b349`

**File:** `LIB/vcfunction.cpp`, immediately after the `FLEXAIDDS_VCT_NORM`
intensive-com block (Lever A, already shipped).

**What it does.** The favorable (negative) `CF.com` channel is unbounded below, so
an overpacked non-native pose can drive `com → −∞` and swamp every attractive term,
leaving no headroom for a downstream orientation-aware rescorer (P1) to out-vote it.
This installs a soft floor at `−F`:

```
softfloor(x) = −F + F·softplus((x+F)/F),   softplus(z)=max(z,0)+log1p(e^−|z|)
  x ≫ −F  → x     (near-identity; pose order by com preserved)
  x → −∞  → −F    (bounded; no single term swamps the CF sum)
  softfloor′ ∈ (0,1]   (monotone ⇒ rank-preserving)
```

**Gate.** `FLEXAIDDS_COM_FLOOR=F` (F>0). **Default-OFF**: unset or F≤0 ⇒ block
skipped ⇒ bit-identical. It is an **enabler, not an accuracy fix** — the commit
message says so.

**⚠ Reconstruction flag.** The exact soft-floor functional form was in the missing
P3 work order. The softplus form here is the standard monotone-bounded realization
of the handoff spec ("soft floor at −F, rank-preserving + bounding"). **Confirm F
and the squashing form against the original work order before the OPS canary run.**

**OPS acceptance:** (1) `FLEXAIDDS_COM_FLOOR` unset ⇒ byte-identical poses vs HEAD.
(2) Set F (e.g. the value from the original WO) and confirm the com channel is
bounded on the 9,500-pose set the handoff cites, with rank preserved above the floor.

---

## 2. P2 — `FLEXAIDDS_CLUSTER_REP` election + `.pop.tsv` dump  ·  commit `89d8dcd3e`

Implements `WORKORDER_clustering_medoid.md` and **supersedes the clustering portion
of `3e674479c`**. New files/edits: `LIB/ClusterRepMode.h` (gate), `LIB/cluster.cpp`
(IP-1 + IP-5 + REMARKs), `LIB/DensityPeak_Cluster.cpp` (IP-2).

**The gate (single source of truth, `ClusterRepMode.h`):**

| `FLEXAIDDS_CLUSTER_REP` | behavior | CF/leader | DensityPeak |
|---|---|---|---|
| unset / `lowcf` | **DEFAULT, bit-identical** | lowest-CF head | `Representative` (lowest-CF) |
| `medoid` | pure **unweighted** geometric medoid (≥3 members) | ✅ | falls back to lowcf |
| `bmedoid` | Boltzmann-CF-weighted medoid (HEAD variant, **ablation only**, ≥2) | ✅ | falls back to lowcf |
| `center` | density-peak center | n/a (lowcf) | `Center` |

**Three HEAD defects corrected:** (1) default-ON → default-OFF; (2) Boltzmann-CF
weighting → pure geometric `medoid` (bmedoid retained for ablation only); (3) the
DP `#define OUTPUT_CLUSTER_CENTER true→false` hardcode → runtime
`output_cluster_center = (mode==center)`, re-exposing the center as an *option* with
lowest-CF as the default. Legacy `FLEXAIDDS_MEDOID_REFINE` is honored only when
explicitly non-zero (aliases to `bmedoid` + deprecation notice); its default-ON
behavior is removed.

**Invariants held:** `Clus_ACF` / between-cluster ranking untouched (representative-
independent — the election only changes *which member* is emitted, not cluster
order). Provenance REMARKs (`cluster_rep_mode`, `cluster_rep_shifted`) emit **only
for non-default modes**, so the `lowcf` PDB stays byte-identical (acceptance gate #1).

**IP-5 dump (`<prefix>_rN.pop.tsv`):** gated on `refstructure==1` **AND**
`FLEXAIDDS_DUMP_POP=1`; default `.rrd/.cad/.mcf/.pdb` byte-unchanged. Columns:
`idx  cluster  rmsd_to_head  rmsd_raw  rmsd_sym  cf_total  cf_com  cf_wal  pose_id  is_elected`.
`is_elected`/`pose_id` are the join that answers "was the near-native population
pose the one elected?" — the instrument that unblocks P1. Per-chrom `com/wal` come
from an audit-only `ic2cf` re-score (never on the benchmark hot path).

**Deferred (per work order, low value for v1):** IP-3 (FO/BindingMode medoid mode —
FO already elects a consensus center), the DP-path `.pop.tsv` (only the CF/leader
default path dumps), and IP-4 (InStream — not the benchmark clustering path).

**OPS acceptance gates (unchanged from the work order §7):**
1. **Bit-identity (MANDATORY):** `FLEXAIDDS_CLUSTER_REP` + `FLEXAIDDS_MEDOID_REFINE`
   both unset ⇒ byte-identical `_j.pdb`/`.cad`/`.mcf` vs a HEAD build with
   `FLEXAIDDS_MEDOID_REFINE=0`, canary 1G9V/1SJ0/1OPK/1M2Z/2HB1, `FLEXAID_SEED=12345`.
   (It will NOT match HEAD's *default*, because HEAD's default is ON — that's the bug.)
2. `ctest` green with the gate unset.
3. Mode smoke on 1G9V: `=medoid` emits `REMARK cluster_rep_mode=medoid` + a
   `[MEDOID_REFINE]` line for ≥1 cluster; `=lowcf` emits neither and equals gate #1;
   `=center` under `clustering_algorithm=DP` emits the density center.
4. (informational) elected-pose RMSD delta — expected ≈ 0 on today's pool (§8 of WO).
5. `FLEXAIDDS_DUMP_POP=1` on 1OF6 ⇒ `.pop.tsv` has ≥1 `rmsd_sym<3` row; confirm from
   `is_elected` whether any sub-3 Å pose was elected.

---

## 3. P1 — KORP-PL orientation-dependent KB rescorer (PLAN, not committed)

**Why plan-only:** detailed WO gone; and it is **gated on P2 measurement** — a
rank-0 rescorer can only re-rank *already-emitted* representatives, so it is
provably neutral until the P2 `.pop.tsv` audit confirms near-native poses now reach
the election pool. Landing it before that measurement is untestable.

**Licensing verdict (from handoff, decisive):** KORP-PL is binary-only / not
redistributable ⇒ it must be an **optional external runtime executable FlexAIDdS
shells out to**, exactly like the PoseBusters `bust` integration. No vendored
copyleft/ARR, no new build-time dependency.

**Verified anchors to build on:**
- Shell-out + stdout capture template: `DatasetRunner::exec_cmd_output` (popen,
  `LIB/DatasetRunner.cpp:2078-2099`) and the tracked/timeout variant
  `exec_dock`/`fork_exec` (`:2057-2075`). Copy this pattern; do **not** invent a new
  process layer.
- Emitted-pose REMARK schema to append to: `cluster.cpp` emit loop (rank-0 `_0.pdb`),
  `REMARK CF.*` block ~`:447-475`.

**Recipe:**
1. New env gate `FLEXAIDDS_KORP_EXE=/path/to/korp-pl` (+ optional
   `FLEXAIDDS_KORP_WEIGHT`). Unset ⇒ **default-OFF, bit-identical** — no rescore.
2. After the representatives are elected and written (post-`write_pdb` in the emit
   loop, or as a post-pass over the `_j.pdb` set), for each emitted pose shell out:
   `"$FLEXAIDDS_KORP_EXE" --receptor <rec> --ligand <pose_j.pdb>` → capture score via
   the `exec_cmd_output` template.
3. Parse the scalar KORP-PL score; store as `REMARK CF.korp=<..>` on each pose.
4. Rank-0 re-election: reorder the emitted representatives by a blended objective
   `rank_key = CF.app + w·KORP` (or pure KORP if `w` sentinel), re-emit the `_j.pdb`
   ordering. Keep the original CF ordering when the exe is absent/fails (fail-open).
5. **Pre-registered G3 kill criterion** (from handoff): if, on the P2-enriched pool,
   pooled `Spearman(rank_key, RMSD)` does not exceed the CF baseline (the handoff's
   offline proof: CF ρ≈+0.02 vs distance-proxy +0.65 on 40 poses), **revert** — the
   term is not recovering orientation signal and must not ship.

**Effort:** ~250–330 LOC, external runtime only. **Do not start until P2 gate #5
shows near-native poses in the elected/near-elected pool.**

---

## 4. P4 — Solis-Wets local refinement of top-K (PLAN, not committed)

**Why plan-only:** detailed WO gone; and the handoff itself flags this as "the least
offline-validatable — it genuinely needs the engine to run," which this environment
cannot do. It also carries a hard merge gate that can only be checked by running.

**Verified anchors to build on:**
- `CF.wal`-only refinement objective + the nearest existing "promote within a mode by
  a single CF channel" precedent: `BindingMode.cpp:575-643`
  (`FLEXAIDDS_PB_AWARE_PROMOTION` picks min `CF.wal` within a mode). P4 generalizes
  this from *selection* to *local optimization* against `cf.wal`.
- Per-pose scoring entry point for the trust-region evals: `ic2cf(...)`
  (`LIB/ic2cf.cpp`), which fills `FA->optres[].cf.wal`; `get_cf_evalue` for the scalar.
- Emitted representative set + gene→IC mapping: `cluster.cpp` emit loop
  (`FA->opt_par[k] = chrom[Clus_TOP[j]].genes[k].to_ic`).

**Recipe (from handoff summary):**
1. Env gate `FLEXAIDDS_REFINE_LOCAL` — **default-OFF, bit-identical**.
2. For each emitted rank-j representative, run Solis-Wets on the rotational/torsional
   DOFs **against `cf.wal` only**, **translation frozen**, under a **15°/25° per-DOF
   trust region + 2 Å centroid cap**, budget ~200 evals/pose (handoff: median 118 used).
3. **Hard merge gate (agent-flagged, blocking):** elected rank-0 symmetry-corrected
   RMSD must **NOT regress** vs the un-refined pick. The ρ≈0 trap applies locally —
   refining toward the wrong minimum can worsen top-1. Gate is checkable only by
   running the canary; **do not merge on a regression.**

**Effort:** ~220 LOC, no new deps, but **requires a working engine build + a canary
run to clear the non-regression gate.** Independent code path from P1/P2/P3.

---

## 5. Recommended sequence for the committer/OPS (unchanged intent)

1. **Build** this branch natively (macOS) and run **P2 gate #1** (bit-identity) and
   **P3** default-off byte check. These prove the two landed changes are safe defaults.
2. Run **P2 gate #5** (`FLEXAIDDS_DUMP_POP=1`) to measure whether near-native poses
   reach the election pool. This is the decision point for P1.
3. Confirm the **P3 `F`** and **soft-floor form** against the original P3 work order.
4. Only then implement **P1** (§3) and measure the ρ lift on the enriched pool
   against the G3 kill criterion.
5. Implement **P4** (§4) in parallel with the elected-RMSD non-regression gate enforced.

Nothing here is merged to `main`. Delete the branch to discard.
