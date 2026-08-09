# entropy.help Documentation

This directory contains draft documentation, templates, and a planning registry for the proposed entropy.help thermodynamic-claims audit workflow.

**Evidence status:** no completed, published, validated, or reproducible entropy.help audit is present in this checkout. `PLANNED_UNVERIFIED` and `EXAMPLE_UNVERIFIED` records are non-claiming scaffolds.

## Contents

- `MANIFESTO.md` — Draft statement of purpose and scientific boundaries
- `THERMODYNAMIC_OUTPUT_SCHEMA.md` — Proposed, not-yet-certified output contract
- `audit-report-template.md` — Non-claiming human-readable report template
- `audit-report-example.json` — Synthetic null-valued schema example; it is unsigned and is not an audit result
- `audits/audits.json` — Planning registry; entries stay `PLANNED_UNVERIFIED` until artifacts exist

## Status

These artifacts are pre-publication planning material. They do not establish that a TotalSampledPartitionFunction path is integrated, that CF/contact-function scores have physical energy units, or that any entropy correction improves docking or affinity prediction.

A status of complete, published, validated, or reproducible requires a provenance record: the registry entry has to link an on-disk JSON report, a Markdown summary, and a separate provenance JSON. The underlying ensemble or durable receipt has to be available wherever a digest or quantitative result is claimed.

Run the fail-closed claims check with:

```bash
python3 scripts/validate_thermo_claims.py
```

The validator checks artifact presence and obvious placeholder/fake cryptographic fields. A passing result is a documentation/provenance gate, not proof of physical thermodynamic validity.

## Coordination

All work tracked in the public GitHub issue:  
https://github.com/LeBonhommePharma/FlexAIDdS/issues/219

Contributions, corrections, and audit requests are welcome.
