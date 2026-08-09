# entropy.help

## Concept Note for a Planned Public Audit Layer

**Status: draft and unverified.** No completed entropy.help audit, deposited
ensemble, or provenance-backed quantitative result is present in this repository
as of this revision. The cases below are proposed evaluation targets, not
published findings.

Molecular docking mixes several distinct objects that are too often described
with the same thermodynamic vocabulary. A search engine may rank poses with an
empirical or contact-function score. A finite sampled ensemble can support a
reproducible weighting calculation. Experimental binding free energy additionally
depends on physical energy calibration, state definitions, concentration,
solvent, and bound/unbound reference terms. Those layers are related, but they
are not interchangeable.

entropy.help is a proposal for making those distinctions inspectable.

### Proposed Scope

For a deposited pose ensemble, an audit would attempt to reconstruct the exact
finite-sample calculation declared by the producer:

- the sampled states and multiplicities;
- the energy or score supplied for each state;
- the weighting convention and temperature-like parameter;
- the source revision, binary/input digests, command, seeds, and timestamps; and
- any separate configurational, vibrational, solvation, concentration, or
  reference-state terms.

For physically calibrated energies, a canonical calculation may use

**F_config = -kT ln Z_sampled**

and

**S_config = -k_B sum_i p_i ln p_i.**

For FlexAIDdS GA output, however, the CF/contact-function score is a ranking
proxy unless a separate calibration establishes physical energy units. Applying
the same algebra to CF values can produce a useful finite-sample *score-space
diagnostic*, but it does not by itself establish physical F, S, or binding ΔG.

### What a Public Record Would Require

A future audit may be promoted from `PLANNED_UNVERIFIED` only when the repository
contains, at minimum:

1. a machine-readable JSON report;
2. a human-readable Markdown summary;
3. a separate provenance record with source and binary/input identity; and
4. the ensemble or durable receipt needed to verify every reported digest and
   recompute the declared quantities.

The repository validator rejects missing artifact links, placeholder signatures
or digests, and completion/publication/reproducibility language without that
evidence. Passing the validator establishes artifact presence and claim hygiene;
it is not, by itself, scientific validation of the method.

### Planned, Unverified Candidate Cases

The following candidates are retained as a prospective test matrix. They have no
result values or outcome claims in the current repository:

1. μ-Opioid receptor + fentanyl — proposed psychopharmacology case study.
2. HIV-1 protease + darunavir — proposed flexible-ligand case study.
3. CDK2 + dinaciclib — proposed kinase case study.
4. BACE1 + verubecestat — proposed alternative-pocket case study.
5. ITC-187 — proposed comparison set, conditional on a licensed, traceable
   experimental-data source and matched thermodynamic definitions.
6. CASF-2016 — proposed docking/scoring benchmark, with pose generation and
   thermodynamic claims reported separately.
7. Thrombin + dabigatran — proposed negative-control case.

These entries are hypotheses about useful validation coverage. They are not
evidence that entropy-aware ranking rescues poses, improves affinity correlation,
or reproduces experiment.

### Request a Candidate Audit

The coordination issue accepts candidate datasets and methodology discussion:

https://github.com/LeBonhommePharma/FlexAIDdS/issues/219

Submitting a request does not create a completed audit or guarantee publication.
Any future result must remain scoped to its deposited artifacts and must preserve
the boundary between CF/contact-function scoring proxies, ensemble-derived
estimates, and experimentally validated thermodynamic quantities.

---

*entropy.help* is a draft open-science initiative seeded by the FlexAIDdS
project. The point is not to manufacture authority; it is to make claims
falsifiable by attaching them to inspectable inputs, calculations, and receipts.
