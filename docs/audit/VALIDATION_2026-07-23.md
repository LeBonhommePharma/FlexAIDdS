# Validation of the Claude Science 24 h FlexAIDdS Audit

**Date:** 2026-07-23
**Branch:** `feat/smart-water-retention` (HEAD `cf428ab28`, 7 commits ahead of `main`)
**Mode:** read-only — no source files were modified. `ctest` was attempted but no test
binaries are currently built in any `build*/` directory (see Claim 6).

Every finding below cites the file and line it was derived from.

---

## Summary table

| # | Claim | Verdict |
|---|---|---|
| 1 | COM_FLOOR stable-softplus is numerically correct | **Confirmed** (code correct; the doc-comment formula above it is wrong) |
| 2 | Blow-up is protein-oxygen driven, not water-driven | **Refuted — inverted.** Retained waters are typed as **carbon** (`assign_types.cpp:112`). The "C×O.3" contact is almost certainly water. |
| 3 | VCT_NORM=1 does not tame the blow-up | **Confirmed, and stronger than stated** — VCT_NORM *amplifies* it whenever contact count < 100 |
| 4 | Thermo gate is inert for selection | **Confirmed exactly** (and `dG_eff` is a population scalar, so it could never rank poses) |
| 5 | N.2→N.ar remap is unconditional | **Was true at HEAD; already fixed in the uncommitted working tree** |
| 6 | PoseBusters schema-pin is the lone ctest failure | **Mechanism confirmed; fix already applied uncommitted and is correct, not masking** |
| 7 | Smart-water retention is default-OFF and bit-identical | **Confirmed and provable from code.** But it leaks crystal-ligand coordinates into receptor prep. |
| 8 | The 3-arm campaign is the decisive next gate | **Refuted as designed** — the script cannot execute on this machine, and no arm isolates the one lever that works |

**Two blockers that matter more than any individual claim:**

1. **`assign_types.cpp:112` types every crystallographic water oxygen as VCT row 1 (`C.1`).**
   Row 1's two most attractive partners are the two most negative entries in the whole
   scoring matrix. This is the mechanism behind the com blow-up, and it re-frames Claims
   2, 3, 7 and 8.
2. **The fixes for Claims 5 and 6 exist only as uncommitted working-tree edits.** They are
   unbuilt, untested, and one `git checkout` away from being lost.

---

## Claim 1 — COM_FLOOR stable softplus

**Verdict: the code is correct. The comment documenting it is not.**

`LIB/vcfunction.cpp:884–899`:

```cpp
if(const char* com_floor_env = std::getenv("FLEXAIDDS_COM_FLOOR")){
    const double F = std::atof(com_floor_env);
    if(F > 0.0){
        for(int j=0; j<FA->num_optres; ++j){
            double& com = FA->optres[j].cf.com;
            const double z = com + F;
            com = -F + (z > 0.0 ? z + std::log1p(std::exp(-z))
                                : std::log1p(std::exp(z)));
        }
    }
}
```

### (a) Is this the correct stable soft-floor at −F?

Yes. It computes `softfloor(x) = −F + softplus(x + F)` using the standard branch-stable
softplus. Both branches evaluate `exp` only on a non-positive argument, so the argument is
in `(0, 1]` and can never overflow. The two properties that matter hold:

- `x → −∞` ⇒ `z → −∞` ⇒ `log1p(exp(z)) → 0` ⇒ `com → −F`. Bounded below.
- `x ≫ −F` ⇒ `z ≫ 0` ⇒ `log1p(exp(−z)) → 0` ⇒ `com → x`. Near-identity.
- `softfloor′(x) = σ(z) ∈ (0,1)`. Strictly monotone, so ordering *by com* is preserved.

The near-identity is far tighter than the comment suggests: the residual is
`log1p(exp(−z))`, which is below `2×10⁻⁹` once `z > 20`, i.e. for any `com > −F + 20`. At
`F = 500` that means every healthy pose (`com ≈ −130`) is passed through untouched to
nine decimal places. That is exactly the behaviour you want from an enabler.

### The documentation defect

Line 871 states the implemented form is:

> `softfloor(x) = −F + F·softplus((x + F)/F)`

The code implements `−F + softplus(x + F)` — unit scale, not `F`-scaled. **The code is the
better of the two**, and the comment should be corrected rather than the code. The
documented `F`-scaled version would be badly behaved: at `x = 0, F = 500` it returns
`−500 + 500·softplus(1) = +156`, i.e. it would corrupt every healthy `com` value in the
population. The unit-scale form is the right choice; only the comment is stale.

The commit message for `4805de1d8` records a smoke test (`−3000→−500, −130→−130, 0→0,
209→209, 300→300, 1000→1000`) whose outputs are consistent with the unit-scale form and
inconsistent with the `F`-scaled form — further evidence the code, not the comment, is the
intended artifact.

### (b) Remaining edge cases

| Input | Result | Assessment |
|---|---|---|
| `z = 0` (`com = −F`) | `−F + log1p(1) = −F + 0.693` | Continuous, no branch discontinuity. Both branches agree in the limit. |
| `F = 0` or negative | Block skipped by `if(F > 0.0)` | Safe. |
| `com = −1e300` | `exp(−1e300) = 0`, `log1p(0) = 0` → `−F` | Correct, no NaN. |
| `com = +inf` | Stays `+inf` | Pre-existing pathology, not introduced here. |
| `com = NaN` | Propagates | Same. |
| `FLEXAIDDS_COM_FLOOR="abc"` | `atof` → 0.0 → **silently disabled** | Real usability hazard: a typo in the launch script produces a silently un-floored run that looks like a valid arm. |

**Recommendations (all low-risk):**

1. Fix the comment at line 871 to state `−F + softplus(x + F)`.
2. Reject unparseable `FLEXAIDDS_COM_FLOOR` loudly instead of silently disabling — use
   `strtod` with endptr checking and `fprintf(stderr, ...)` on failure. Given that the
   value of the whole campaign depends on this env var being live, silent-off is the
   single worst failure mode available here.
3. Emit the effective `F` once per dock in the provenance JSON so a completed run can be
   audited after the fact.

### (c) One structural caveat on rank preservation

The comment's "monotone ⇒ rank-preserving" claim is true *of the com channel alone*. Total
CF is `com + sas + wal + elec + hbond + …`, so squashing `com` changes the relative weight
of every other channel and **does** reorder poses by total CF. That is the intended effect
— but the code comment should say so, because "rank-preserving" reads as a much stronger
guarantee than what is delivered.

---

## Claim 2 — "Water blow-up is protein-oxygen driven, not water-driven"

**Verdict: refuted, and the causality is inverted. This is the most consequential finding
in this review.**

### The matrix is bounded; the area multiplier is not

`LIB/vcfunction.cpp:576–584`:

```cpp
double contribution = 0.0;
if(energy_matrix->weight){
    if(FA->normalize_area){ contribution = yval*area/surfA; }
    else                  { contribution = yval*area; }
}
```

With `energy_matrix->weight` set, `get_yval` (`vcfunction.cpp:1311`) returns the raw matrix
scalar. I scanned `MC_st0r5.2_6.dat` directly: **the most negative entry in the entire
matrix is `2-2 = −198.9`.** No pair value exceeds ~199 in magnitude. `FA->normalize_area`
defaults to `false` (`config_defaults.h:26`, `top.cpp:591`), so the live path is
`contribution = yval * area` — extensive in the Voronoi facet area.

So a single contact at `−3254` requires an area of about `3254 / 180.8 ≈ 18 Å²`, which is
a perfectly ordinary Voronoi facet. **Nothing is numerically broken.** The unboundedness
is not in the potential's *values* — it is in the `value × area × Σcontacts` construction,
which has no reference state and no per-contact ceiling. Science's phrase "unbounded in
the NRG matrix" is not right; the matrix is tightly bounded. The correct statement is that
the *functional form* is unbounded and the matrix units are arbitrary ROC-trained scores,
not kcal/mol.

### The "C" in "C×O.3" is a water molecule

`LIB/assign_types.cpp:108–114`:

```cpp
// Set a Hydrophilic type to water molecules
for(k=1;k<=FA->res_cnt;k++){
    rot=residue[k].rot;
    for(i=residue[k].fatm[rot];i<=residue[k].latm[rot];i++){
        if(strcmp(residue[atoms[i].ofres].name,"HOH") == 0){atoms[i].type = 1;}
    }
}
```

Every retained crystallographic water oxygen is assigned `type = 1`. Under the canonical
VCT numbering documented at `top.cpp:60–67`, **row 1 is `C.1`** — sp carbon. Row 40
(`SOLVENT`) exists and is *entirely zero* in `MC_st0r5.2_6.dat`, so waters are not merely
mistyped, they are typed into a row that was never trained on solvent.

Now look at what row 1 is worth against oxygen:

```
1-13 (C.1 × O.2)  = −198.3   ← the 2nd most attractive pair in the whole matrix
1-14 (C.1 × O.3)  = −180.8
1- 3 (C.1 × C.3)  = −162.9
```

Science's reported signature is *one dominant `C × O.3` contact at ~−3254*. Under this
typing, `C.1 × O.3` is precisely the pair you get from **a retained water oxygen touching
a ligand or protein sp3 oxygen** — and it draws from one of the strongest cells in the
matrix. The mechanism logger in `ops/run_astex85_twoarm.sh:110` reads the contact *type
pair* out of the pose REMARKs; if the REMARK reports canonical type indices, a water shows
up as "C" and is indistinguishable from a genuine carbon. Science read the label and
concluded "protein carbon"; the label is wrong.

The repository already knows this at the prose level. `LIB/modify_pdb.cpp:49–52`:

> "Without this filter, every low-B-factor water in the structure becomes a receptor atom
> and the GA can bury the ligand inside the solvent shell, harvesting unbounded
> complementarity from sub-Ångström ligand-O ⋯ HOH-O contacts."

and `DatasetRunner.cpp:5843–5846` records `1JD0` retaining 156 waters with `CF = −4269`
against an expected `~−50`. The water hypothesis was the original hypothesis; the typing
bug is *why* it looks like a carbon problem in the instrumentation.

This is the same class of defect already recorded in project memory as the "ligand
type-numbering bug" — a reader writing a non-canonical type index that the matrix then
reads as a different element.

### What this means for the audit

- Science's Type I / Type II split in the campaign header
  (`ops/run_astex85_twoarm.sh:8–17`) is likely **one mechanism, not two**. Arm C strips
  all waters; if the blow-up survives Arm C, that is evidence for a genuine protein-oxygen
  Type II. If it does not, Type II never existed.
- `1SG0` should be re-examined specifically: dump the residue name of the atom on the
  dominant contact. If it is `HOH`, Claim 2 is settled.

**Recommendations:**

1. **Highest priority.** Decide what waters should be typed as and make it explicit. Row 40
   (`SOLVENT`) being all-zero means it cannot be used as-is; the honest options are
   (a) type water O as `O.3` (row 14) and accept that it is scored as an ordinary sp3
   oxygen, or (b) keep row 1 but document it as a deliberate hydrophilic proxy and retrain
   that row. What is not defensible is the current state, where a comment says
   "hydrophilic" and the matrix says "sp carbon".
2. Extend the per-pose instrumentation from `cf428ab28` to emit the **residue name and
   element** of the dominant contact's two atoms, not just the VCT type indices. The whole
   Type I/Type II confusion is downstream of the logger printing an index whose meaning is
   wrong.
3. Add a regression test asserting `assign_types` produces the intended type for an `HOH`
   oxygen, so this cannot silently drift again.

---

## Claim 3 — VCT_NORM=1 does not tame the blow-up

**Verdict: confirmed, and the canary's 0/5 is fully explained by the code. VCT_NORM does
not merely fail to bound the sum — it makes the pathological pose worse.**

`LIB/vcfunction.cpp:854–861`:

```cpp
if(FA->vct_normalize_contacts){
    constexpr double VCT_NREF = 100.0;
    for(int j=0; j<FA->num_optres; ++j){
        if(vct_ncon[j] > 0){
            FA->optres[j].cf.com *= VCT_NREF / (double)vct_ncon[j];
        }
    }
}
```

Three independent reasons this cannot work:

**1. The multiplier is greater than 1 for any contact count below 100.** `vct_ncon[j]`
counts com-contributing Voronoi contacts for one optimizable residue — for a ligand this
is typically tens. At `N = 20` the multiplier is `×5`: a `com` of `−3254` becomes
`−16 270`. The rescale that was added to keep the term's magnitude comparable
(`VCT_NREF`) is precisely what converts a normalizer into an amplifier in the regime that
matters.

**2. It normalizes by the wrong quantity.** The blow-up is concentrated in *one* contact
carrying ~95 % of the total (per Science's 1SG0 figure). Dividing a sum by its cardinality
does nothing to a sum dominated by a single term. Worse, a ligand buried in a solvent
shell has *fewer, larger* Voronoi facets than a ligand in a normal pocket — so the
pathological pose gets a *smaller* `N`, hence a *larger* multiplier. The normalization is
anti-correlated with the pathology.

**3. `com` is extensive in area, not in count.** `contribution = yval * area`. The natural
extensive variable is buried surface area, so `Σ area` — not `N` — is the only
denominator that makes the quantity intensive. This mirrors the already-recorded finding
that ACF is extensive and regressed 1P2Y.

### Normalization schemes that would work

Ranked by risk-adjusted expected value:

**(a) Per-contact energy cap — recommended first move.** Directly mirrors the per-contact
wall ceiling already in the codebase (`WAL_CONTACT_CAP = 50`, `vcfunction.cpp:543`), which
is documented as having restored the GA gradient. The com channel has no equivalent:

```cpp
// at vcfunction.cpp:632, replacing `cfs->com += contribution;`
cfs->com += (contribution < -COM_CONTACT_CAP) ? -COM_CONTACT_CAP : contribution;
```

This is strictly better than the global soft floor for the observed failure mode, because
it removes the *single* runaway contact while leaving the other 40 healthy contacts fully
expressed. The global `COM_FLOOR` at `−500` cannot distinguish "one absurd contact" from
"forty good ones", and once a pose saturates the floor its com gradient vanishes entirely.
The cap keeps a gradient everywhere. Precedent, symmetry with the wall channel, and a
one-line diff make this the highest-value experiment available.

**(b) Buried-area normalization.** `com / Σ area` gives mean complementarity per Å² of
contact — genuinely intensive, and immune to the "fewer larger facets" inversion that
breaks count normalization. Rescale by a reference area (~300 Å² for a drug-like ligand)
rather than a reference count.

**(c) Per-heavy-atom normalization.** `com / n_heavy`. Crude but ligand-size-invariant, and
the machinery already exists — `ThermodynamicEngine` computes `H_vct` per heavy atom
(`gaboom.cpp:1276–1277`, `thermo_result.n_heavy_atoms`).

**(d) Sigmoid saturation per contact.** `−E_max · tanh(|c| / E_max)`. Smooth and
differentiable, but it distorts every contact rather than only the outliers, and it adds a
parameter. Prefer (a) unless the GA turns out to be sensitive to the cap's kink.

Whichever is chosen, **VCT_NORM should be retired, not merely defaulted off.** Leaving a
live env var that silently multiplies com by `100/N` is a trap for the next campaign.

---

## Claim 4 — Thermo gate is inert for selection

**Verdict: confirmed exactly, with an additional reason Science did not give that makes
the claim even stronger.**

### (a) Ordering of QuickSort vs thermo compute

`LIB/gaboom.cpp`:

- **line 1238** — `QuickSort((*chrom),0,GB->num_chrom-1,true);` — ranking fixed.
- **line 1241** — `if (FA->thermo_engine_enabled && FA->thermo_engine != nullptr) {` — thermo block *begins*.
- **line 1308** — `FA->thermo_result = FA->thermo_engine->compute(...)`.
- **lines 1340–1352** — `printf("[THERMO3] dG_eff=...")`.

The sort strictly precedes the compute. Confirmed.

Note the sort is at line **1238**, not 705 as the campaign script's header comment claims
(`ops/run_astex85_twoarm.sh:32`). Line 705 is a different, earlier sort inside the
generation loop. The conclusion is unaffected — every `QuickSort` call site precedes line
1241 — but the citation in the script header is wrong.

### (b) No selection consumer

`grep -rn "dG_eff" LIB/` returns exactly four files: `ThermodynamicEngine.cpp` (computes
it), `ThermodynamicEngine.h` (declares it), and `gaboom.cpp` (prints it). There is no read
of `thermo_result.dG_eff` in any clustering, ranking, or election path. The header states
this explicitly at `ThermodynamicEngine.h:62`:

> "DIAGNOSTIC ONLY — dG_eff does not affect pose selection, at any flag setting."

`FLEXAIDDS_THERMO_SCORE=1` gates only `apply_gate` (`ThermodynamicEngine.h:35`), which
overwrites the *reported* `dG_eff` with a `+1000` sentinel. That value is then printed and
discarded.

### The stronger argument

`ThermoResult::dG_eff` is **a single `float` for the whole population** — `<CF> − T·H`
where `<CF>` is `mean_CF` and `H` is the Shannon entropy of the population's energy
histogram (`ThermodynamicEngine.cpp:155`). It is not a per-pose quantity. Even if it were
wired into `QuickSort` tomorrow, it would add the same constant to every chromosome and
change nothing. Claim 4 is not "currently unwired" — it is **structurally incapable of
ranking poses in its present formulation.** This matters because it means "wire dG_eff
into selection" is not a small change; it requires defining a per-pose thermodynamic
quantity first.

### (c) Is wiring it in scientifically sound, and where?

Not as it stands, for the reason above. A defensible path, in order:

1. **Define a per-pose ΔG.** The population entropy `H` must be decomposed into a per-pose
   term. The natural object is a per-*cluster* entropy: for binding mode `m` with member
   poses `{i}`, `S_m = −Σ p_i ln p_i` over Boltzmann weights within the cluster. This is
   already close to what `BindingMode.cpp` computes.
2. **Apply it at cluster election, not at GA selection.** The GA's job is to find minima;
   entropy is a property of a basin, not of a point, and a single pose has no entropy. The
   right integration point is where `cluster.cpp` elects the representative of each mode
   and where modes are ranked against each other — i.e. `G_m = <CF>_m − T·S_m`, ranking
   *modes*, with the representative still chosen by CF within the mode.
3. **Never let it touch the GA fitness.** Injecting a population-level entropy into
   per-chromosome fitness creates a feedback loop: fitness depends on the population's
   diversity, which selection then changes, which changes fitness. Given the recorded
   gen-10 energy-Shannon collapse, this codebase is already fragile in exactly that way.

**Regression risk.** Low if confined to mode ranking (the current mode ranking is
CF-only, so an A/B is clean and cheap). High if placed in the GA loop. Project memory
already records that near-native poses are consistently *worse* in CF than the elected
false minimum by −2 to −86 units, and that selection is CF-bound — an entropy term that
favours broad basins is one of the few principled discriminators that could close that
gap, which is a real argument for doing (2) properly rather than dismissing the thermo
stack.

---

## Claim 5 — N.2 → N.ar remap is unconditional

**Verdict: true of committed HEAD; already fixed in the working tree, but uncommitted.**

`git diff LIB/top.cpp` shows an *uncommitted* change:

```diff
-	if (!strcmp(s, "N.2"))   return 10;  // N.ar — sp2 imine is an acceptor; ...
+	// N.2 keeps its own canonical row 7. It is NOT remapped to N.ar (10):
+	// aromaticity is already perceived upstream — SybylTyper.cpp tests
+	// in_aromatic_ring() first and emits N.ar for any N in an aromatic ring
+	if (!strcmp(s, "N.2"))   return 7;
```

with a matching uncommitted change in `LIB/Mol2Reader.cpp:43` keeping the two readers in
sync. So Science's claim was accurate against `HEAD`, and someone has since applied the
fix locally without committing it.

### (c) Are rows 7 and 10 degenerate? No — the remap materially changed scoring.

I counted non-zero entries directly in `MC_st0r5.2_6.dat`:

| Row | Type | Non-zero entries |
|---|---|---|
| 7 | N.2 | **9** |
| 8 | N.3 | **0** (dead) |
| 10 | N.ar | **14** |

They are emphatically not degenerate. Comparing the shared partners:

| Partner | Row 7 (N.2) | Row 10 (N.ar) |
|---|---|---|
| 11 (N.am) | −195.7 | −157.3 |
| 12 (N.pl3) | −122.3 | −26.35 |
| 13 (O.2) | −87.38 | −15.0 |
| 14 (O.3) | −144.8 | −189.6 |
| 15 (O.co2) | −197.6 | −125.3 |
| 18 (S.3) | +28.13 | −190.1 |

The `N × O.2` interaction differs by a factor of **5.8** and the `N × S.3` interaction
**flips sign** (+28 repulsive vs −190 strongly attractive). Remapping N.2→N.ar was not a
cosmetic substitution — it rewrote the entire nitrogen interaction profile for every sp2
imine nitrogen in the benchmark. Amidines, guanidines, imines, azomethines, and
non-aromatic C=N in fused systems are common in the Astex set; these are exactly the
groups affected.

The fix's own justification is also sound and I verified the premise: aromaticity is
perceived upstream, so an atom reaching the `N.2` branch has already been judged
non-aromatic and an "only when aromatic" guard would be dead code. Row 7 being live (9
entries) means the dead-row argument that legitimately justifies `N.3 → N.am` does not
transfer.

### The remaining contested remap

`N.3 → 11 (N.am)` is *still* unconditional, and this one is genuinely contestable. Every
sp3 amine — including every basic aliphatic amine, the most common protonatable group in
CNS drugs — is scored with the amide-nitrogen row. The justification (row 8 is all-zero) is
factually correct and I confirmed it, but the consequence is that FlexAIDdS has **no
parameters at all for sp3 amine nitrogen**, and silently substitutes a chemically distinct
type. For a psychopharmacology-targeted engine this is a first-order gap: protonated
tertiary amines are the defining pharmacophore of most aminergic ligands.

**Recommendations:**

1. **Commit the N.2 fix immediately** and re-run the Astex-85 baseline. It is a scoring
   change to a live matrix row and invalidates any prior comparison across it.
2. Treat dead row 8 as a **training-data gap, not a typing problem.** The right fix is to
   parameterize N.3, not to keep aliasing it. Until then, log a one-time warning per dock
   listing how many ligand atoms were type-substituted, so downstream results carry the
   caveat.
3. Add a unit test asserting `sybyl_name_to_canonical_vct` and
   `Mol2Reader::sybyl_to_flexaid_type` agree on every SYBYL name. These two tables have now
   drifted apart at least once, and drift means SDF and MOL2 inputs score the same molecule
   differently.

---

## Claim 6 — PoseBusters schema-pin is the lone ctest failure

**Verdict: mechanism confirmed by inspection; the fix is already applied uncommitted and is
correct. It masks nothing.** I could not run `ctest` — see caveat below.

### (a) The code

`git diff LIB/PoseBust/BustCli.cpp` shows an uncommitted removal:

```diff
         "minimum_distance_to_protein",
-        "no_protein_clashes",
         "volume_overlap_with_protein",
```

### (b) The failure mechanism

`tests/test_posebust.cpp:156–163` defines `synthetic_full_pb_header()`, whose column list
does **not** contain `no_protein_clashes`. With that name in
`mandatory_pb_check_columns()`, `TEST(BustCliSchema, PassesWithFullMandatorySet)` at line
172 must fail with `mandatory_checks_missing` — the schema pin fails closed on a column
the fixture never provides. That is exactly the failure Science reports, and the removal
resolves it.

**Caveat on verification:** I ran `ctest` in both `build/` and `build_test/` and every test
reported `Not Run` — no test binaries are compiled in any `build*/` directory. So I confirm
the failure *by inspection of the fixture and the pin*, not by execution. The claim that it
is the **lone** failure is therefore unverified; it needs
`cmake -DBUILD_TESTING=ON && cmake --build . && ctest --output-on-failure` to stand.

### (c) Is dropping the column correct, or is it masking a real check?

**Correct, and it masks nothing.** I checked the installed PoseBusters package in
`.venv-posebusters/`. `posebusters/config/redock.yml` defines the protein-clash checks as:

```
line 160:  no_clashes: "Minimum distance to protein"
line 231:  no_volume_clash: "Volume overlap with protein"
```

There is **no `no_protein_clashes` column** anywhere in the redock suite. The name was
fabricated — a plausible-sounding column that upstream never emits. Both genuine
protein-clash checks (`minimum_distance_to_protein` and `volume_overlap_with_protein`)
remain in the mandatory list after the fix, so protein-clash gating is fully intact. The
pin is *strengthened* by the removal: previously it could never pass against real
PoseBusters output either, meaning the pin was failing closed on every genuine run.

**Recommendations:**

1. Commit the fix with the `redock.yml` line references as the justification.
2. Generate the fixture header **from** `mandatory_pb_check_columns()` rather than
   hand-writing it in the test, so a pin bump cannot desynchronize from its own fixture
   again.
3. Add a check that verifies the pin against a real `bust --outfmt csv` header captured
   from the pinned PoseBusters version, so an upstream rename is caught at test time rather
   than at campaign time.

---

## Claim 7 — Smart-water retention is default-OFF and bit-identical

**Verdict: confirmed, and provable from code rather than empirically. But the feature has a
methodological problem Science did not raise.**

### Provable default-off

Three independent gates, all defaulting to disabled:

1. `LIB/config_defaults.h:162` — `{"binding_site_water_radius", V(0.0f)}`.
2. `LIB/config_parser.cpp:309` — reads with default `0.0f`.
3. `LIB/DatasetRunner.cpp:5849–5852`:
   ```cpp
   const bool smart_water = (smart_water_env != nullptr && std::string(smart_water_env) == "1");
   double bs_water_radius = smart_water ? 4.5 : 0.0;
   ```

And the consuming guard, `LIB/modify_pdb.cpp:278–279`:

```cpp
if(keep_structural_waters && remove_water && !exclude_het &&
   binding_site_water_radius > 0.0f && !oracle_lig.empty())
```

With `binding_site_water_radius == 0.0f` the entire block is skipped, `bs_water_keep`
stays empty and `bs_water_filter_active` stays `false`. No branch downstream of it can
observe a difference. Bit-identity is a property of the control flow, not an empirical
observation — Claim 7 is provable, and stronger than Science stated it.

### The problem Science missed: crystal-pose leakage into receptor preparation

`LIB/modify_pdb.cpp:271–276`:

```cpp
std::string oracle_lig = (oracle_ligand_path != nullptr) ? oracle_ligand_path : "";
if(oracle_lig.empty()) {
    const char* env = getenv("FLEXAIDDS_RMSDST");
    if(env != nullptr) oracle_lig = env;
}
```

and the retention criterion at line 54: `d(HOH-O, nearest crystal-ligand heavy atom) ≤ radius`.

**Smart water retention selects which waters to keep using the crystal ligand's
coordinates** — falling back to `FLEXAIDDS_RMSDST`, the RMSD reference pose. The comment at
lines 269–271 is candid that the runner "deliberately leaves [reference_ligand] empty to
keep crystal coordinates out of the GA", then reaches for the RMSD reference instead.

Keeping crystal coordinates out of the GA but using them to prepare the receptor is still
information leakage: the receptor handed to the search encodes where the answer is. Arm A
and Arm B of the campaign both enable this. Any success-rate improvement they show is not
attributable to "better water handling" in a way that would transfer to a prospective
target, where no crystal ligand exists. This is the same class of issue as the recorded
seed-echo inflation finding.

**Recommendations:**

1. **Report smart-water results as oracle-assisted**, clearly separated from blind numbers.
   The engine already distinguishes oracle from blind modes elsewhere; this path should
   inherit that discipline.
2. Provide a **prospective variant** that selects bridging waters using the detected cavity
   (`CavityDetect/`) rather than the crystal ligand. Same two criteria — inside the cavity,
   H-bonded to protein — with a blind cavity center. That version is publishable; the
   current one is a diagnostic.
3. The fallback at `modify_pdb.cpp:283` (no readable ligand → plain B-factor retention)
   silently changes the protocol mid-campaign per target. It logs to stderr, which is
   good; it should also be recorded in the per-target provenance JSON so the analysis can
   segregate those targets.

---

## Claim 8 — The N=85 three-arm campaign is the decisive next gate

**Verdict: refuted as currently designed. The script cannot run on this machine, and even
if it could, no arm isolates the one lever the canary says works.**

### Blocker 1: the script aborts before doing anything

`ops/run_astex85_twoarm.sh:63–64`:

```bash
exec 9>"${LOCKFILE}"
flock -n 9 || { echo "[ABORT] another campaign holds the lock: ${LOCKFILE}" >&2; exit 1; }
```

`command -v flock` on this machine returns **nothing**. `flock` is a util-linux tool and is
not present in macOS by default. The `||` catches the `command not found` (exit 127) and
the script exits 1 with a misleading "another campaign holds the lock" message. **The
three-arm campaign has never run and cannot run as written.**

Fix: use a portable lock. `mkdir "${LOCKFILE}.d"` is atomic on every POSIX filesystem:

```bash
LOCKDIR="/tmp/flexaidds_benchmark.lock.d"
mkdir "${LOCKDIR}" 2>/dev/null || { echo "[ABORT] another campaign holds ${LOCKDIR}" >&2; exit 1; }
trap 'rmdir "${LOCKDIR}"' EXIT
```

### Blocker 2: arm failures are silently swallowed

Each arm is wrapped as `( ... ) || { echo "[WARN] ARM X failed rc=$?" ; }` (lines 148, 160,
172). But `run_arm` ends with `com_summary`, which ends with `echo` — so `run_arm` always
returns 0, the subshell always exits 0, and the `||` handler is unreachable. A crashed arm
produces an empty output directory and a clean-looking log. Capture and propagate the
runner's `rc` explicitly: `return "${rc}"` at the end of `run_arm`.

### Blocker 3 (scientific): no arm isolates COM_FLOOR

This is the fatal design problem given the canary result. The arms are:

| Arm | Water | VCT_NORM | COM_FLOOR | Thermo |
|---|---|---|---|---|
| A | smart | off | off | off |
| B | smart | **on** | 500 | on (inert) |
| C | stripped | **on** | 500 | off |

The canary found **VCT_NORM = 0/5 and COM_FLOOR = 5/5**. Both arms that carry the com fix
carry it as `VCT_NORM=1 + COM_FLOOR=500` bundled together — and per Claim 3 these two are
not independent: `VCT_NORM` multiplies `com` by `100/N` **before** `COM_FLOOR` clamps it
(`vcfunction.cpp:854` runs before `:884`). With `N < 100` that pre-multiplication inflates
`com`, driving more poses into the floor's saturated region where the com gradient is
zero. The composition is actively worse than `COM_FLOOR` alone, and neither B nor C can
measure `COM_FLOOR` on its own. **The campaign as written cannot deliver the verdict it was
designed to deliver.**

Two further design notes:

- The header's "thermo gate is INERT" reasoning (lines 30–36) is **correct** — validated in
  Claim 4 — but it means Arm B spends a full 85-target run on three env vars that provably
  do nothing. The B−C delta is confounded anyway: B and C differ in *both* water treatment
  and thermo flags. Since thermo is inert, the delta is interpretable as water — but only
  because of an argument the script itself has to make in a comment. A clean design would
  not need the disclaimer.
- One thing the script gets right that is worth preserving: water retention is applied in
  `modify_pdb` at engine init from the per-job config, not at cache-prep time, so the
  shared `--cache` directory does **not** leak water settings between arms. I checked this
  specifically because it would have invalidated everything; it is fine.

### Redesigned campaign

Drop VCT_NORM entirely (per Claim 3 it should be retired, not tested). Make water and
com-taming a clean 2×2 with the smallest arm count that answers the question:

| Arm | Water | com taming | Question answered |
|---|---|---|---|
| **A0** | all waters kept (status quo) | none | Baseline. Reproduces the blow-up; anchors everything. |
| **A1** | all waters kept | `COM_FLOOR=500` only | Does the floor alone recover success? (canary says yes, 5/5) |
| **A2** | **all waters stripped** | none | Is the blow-up water-driven? Given Claim 2, this is now *the* decisive arm. |
| **A3** | all waters stripped | `COM_FLOOR=500` | Do the two levers compose, or is one redundant? |

`A2` is the arm that adjudicates Claim 2 and it does not currently exist — Arm C confounds
water-stripping with com-taming, so a null result in Arm C cannot distinguish "water wasn't
the problem" from "the com fix already handled it". Smart-water retention should be dropped
from the decisive campaign entirely and run separately as an oracle-assisted study, per
Claim 7.

If the water typing fix from Claim 2 lands first, add `A2'` (waters kept, typed as `O.3`),
which is the *scientifically correct* configuration and may make the whole com-taming
question moot.

---

# Enhancement recommendations

## 1. What pushes FlexAIDdS past the 91.8 % (v127) ceiling

Ordered by expected value per unit of effort, grounded in what this review and the recorded
project history actually show.

**(a) Fix water typing — `assign_types.cpp:112`.** One line. Waters currently score against
the two most attractive cells in the matrix while masquerading as carbon. Everything
downstream — the com blow-up, the Type I/Type II taxonomy, the smart-water feature, the
three-arm campaign — is built on top of this defect. Fix it before running anything else,
because it may dissolve the problem the campaign was designed to study.

**(b) Per-contact com cap.** Section "Claim 3 → (a)". One line, exact precedent in the
wall channel, and it targets the observed single-dominant-contact failure mode more
precisely than a global floor.

**(c) Close the selection gap, not the search gap.** Project memory is unambiguous that
near-native poses are *found* and then *not elected* — 22 found vs 8 selected in v102, with
the near-native always worse in CF by −2 to −86. At a 91.8 % ceiling the remaining ~7
targets are almost certainly selection failures, not search failures. The levers are
(i) mode-level ranking by `G_m = <CF>_m − T·S_m` (Claim 4 recommendation 2), and
(ii) a consensus or orientation-aware rescorer applied *only* to mode representatives,
where it is cheap and cannot perturb the GA.

**(d) Parameterize the dead matrix rows.** Row 8 (N.3) is entirely zero, and row 40
(SOLVENT) is entirely zero. Two of the most common chemical entities in a binding site —
sp3 amine nitrogen and water — have no parameters at all. This is the deepest limitation
found in this review, and unlike everything else it is a data problem, not a code problem.

**(e) Symmetry-aware RMSD.** Already flagged in the recorded methodology audit. At >90 %
success, several of the remaining failures may be scoring artifacts of RMSD rather than
docking failures. Cheap to check and it can only improve the honest number.

A caution: at 91.8 % on Astex-85, seven targets separate you from the ceiling and
run-to-run GA noise is recorded at ~2 Å with roughly six targets sitting within 0.5 Å of
the 2 Å cutoff. **Any change below about ±4 targets is not measurable in a single run.**
Every experiment above needs `OMP_NUM_THREADS=1` plus fixed seeds (the recorded determinism
protocol) and, ideally, repeats — otherwise you will be attributing noise.

## 2. A principled bound on CF.com, better than a hard soft-floor

The soft floor treats a symptom. The disease is that `com = Σ (matrix_value × area)` has
**no reference state** — there is no configuration whose com is defined to be zero, so the
only thing constraining the magnitude is how much surface the pose can bury. A pose that
buries more always scores better, without limit. This is the same structural flaw that
made ACF extensive and regressed 1P2Y.

The principled fix is a **reference-state normalization**, in increasing order of ambition:

**Tier 1 — per-contact cap (ship this).** `min(contribution, −E_max)` per contact, `E_max`
around 200–400. Bounds the outlier, preserves the gradient everywhere else, mirrors the
wall channel. One line.

**Tier 2 — buried-area normalization.** `com_normalized = com × (A_ref / Σ area)` with
`A_ref ≈ 300 Å²`. Makes the term intensive in the *correct* extensive variable, and unlike
count normalization it cannot invert on the pathological pose.

**Tier 3 — an actual reference state.** For each contact, subtract the expectation for that
atom-type pair over a random-contact ensemble:
`ΔE_ij = E_ij − <E>_ij`. This is what every knowledge-based potential (DrugScore, PMF,
ITScore) does, and it is what makes those functions size-transferable. It requires
computing per-type-pair contact-frequency expectations over a reference database — a real
project, but it is the difference between a heuristic and a potential of mean force.

Given the recorded null results on matrix magnitude rescaling (`FA_matrix_v1`: 27/85 both
arms, r ≈ 0), do **not** expect Tier 3 to move accuracy by itself. Its value is that it
makes the score comparable across ligand sizes, which is what a docking engine needs for
prospective virtual screening — a different and arguably more important goal than one more
Astex target.

## 3. The right thermodynamic integration for dG_eff

**Not FEP, not canonical MC. A per-mode entropy correction at cluster election.**

FEP and canonical MC both require sampling that the GA does not produce — the GA generates
a biased, non-Boltzmann population by construction, so any free energy computed from it is
not a free energy. The current `dG_eff = <CF> − T·H` already commits this error: `H` is the
Shannon entropy of a *GA population's* energy histogram, which reflects the GA's
convergence state (and, per the recorded gen-10 collapse, its failure modes) far more than
it reflects the physical density of states.

The simpler correction that is actually defensible:

1. **Cluster poses into modes** (already done in `BindingMode.cpp`).
2. **Within each mode, compute** `S_m = −k Σ p_i ln p_i` over Boltzmann weights
   `p_i ∝ exp(−β·CF_i)` for member poses. This *is* a configurational entropy estimate over
   a basin, which is a physically meaningful object even from biased sampling, because it
   is a property of the local density of poses rather than the global population.
3. **Rank modes by** `G_m = <CF>_m − T·S_m`. Elect the representative within each mode by
   CF, exactly as today.
4. **Validate against ITC data before believing it.** The ITC calibration framework already
   exists in this repo (recorded in project memory, with ~944 unified BindingDB rows). A
   thermodynamic term that does not improve correlation against measured ΔG/ΔH/ΔS is not
   worth its regression risk. This is the single highest-value use of that dataset.

The key discipline: **keep it out of the GA fitness.** Entropy is a basin property. Putting
a population-derived quantity into per-individual fitness creates the exact
diversity-fitness feedback loop this codebase has already been bitten by.

## 4. Scoring-architecture changes within reach

**(a) Desolvation — the largest missing physical term, and the one most relevant to the
water problem.** The current `SAS` channel is a raw solvent-accessible-surface penalty
(`vcfunction.cpp:777–791`) with no atom-type dependence: burying a charged oxygen and
burying an aliphatic carbon cost the same. An atom-type-weighted desolvation penalty
`ΔG_desolv = Σ σ_type · ΔASA_i` (Eisenberg–McLachlan style, five to seven atom classes) is
a contained change to an existing channel, and it is the term that would *physically*
prevent burying a ligand in a solvent shell — the failure this whole audit is about.
Recommended as the highest-value architectural addition.

**(b) Distance-dependent dielectric — already implemented, and worth an audit.**
`vcfunction.cpp:655`:

```cpp
double E_elec = KCOULOMB * qA * qB / (FA->dielectric * dist * dist);
```

The comment above says `E_elec = (332.0637·qA·qB)/(eps·r)` with `eps = dielectric·r`, which
gives `1/r²` — and the code implements `1/r²`. Code and intent agree, so this is not a bug.
But note that the *comment on line 654* labels it "distance-dependent dielectric" while the
form used is the sigmoidal-free simple DDD; if `FA->dielectric` is being set to a bulk value
like 78 rather than a DDD coefficient (typically 4), the term is being scaled ~20× too
small. Worth one grep of where `FA->dielectric` is set before the next campaign.

**(c) Directional hydrogen bonding — present but disabled during search.**
`vcfunction.cpp:664–666` gates the Gaussian H-bond term on
`use_hbond_search && !hbond_rank_rescore`, i.e. by default it contributes at rank time but
not during the GA. Since H-bond geometry is the primary thing distinguishing a native pose
from a shape-complementary false minimum, and Claim 5 shows the nitrogen typing that feeds
it was wrong until the uncommitted fix, this deserves a clean A/B **after** the N.2 fix
lands. The recorded v51/v52 history shows a previous attempt over-corrected and collapsed
`cf_native`; the lesson there was about the *penalty* formulation, not about whether
directionality belongs in search.

**(d) What not to do.** Do not revisit matrix magnitude rescaling or `r0`/VCT geometry
tuning. Two independent, well-documented experiments (`FA_matrix_v1`, `lever2_r0`) returned
clean nulls — deepening the CF wells 18× moved zero poses. The lever is selection and
physics coverage, not the shape of the existing potential.

---

# Immediate actions

Ordered. The first three are cheap and unblock everything else.

1. **Commit the working-tree fixes.** `LIB/top.cpp`, `LIB/Mol2Reader.cpp` (N.2 → row 7) and
   `LIB/PoseBust/BustCli.cpp` (schema pin) are correct, validated above, and currently
   exist only as uncommitted edits. They are one stray `git checkout` from being lost.
2. **Build with `-DBUILD_TESTING=ON` and run `ctest --output-on-failure`.** No test binaries
   exist in any `build*/` directory right now, so the "lone failure" claim is unverified and
   the project is flying without its test suite.
3. **Fix `assign_types.cpp:112`** (water typing) and correct the stale comment at
   `vcfunction.cpp:871` (soft-floor formula).
4. **Fix the campaign script**: replace `flock` with `mkdir`-based locking, propagate
   `run_arm`'s exit code, correct the `gaboom.cpp:705` citation to `:1238`, and remove the
   stale `setsid` comment (`setsid` is absent from macOS and is not used in the script).
5. **Redesign the campaign as A0/A1/A2/A3** (Claim 8), dropping VCT_NORM and moving
   smart-water to a separate oracle-assisted study.
6. **Re-baseline Astex-85 after the N.2 and water-typing fixes.** Both change live matrix
   rows; every number predating them, including the 91.8 % v127 figure, is measured against
   a different scoring function.
