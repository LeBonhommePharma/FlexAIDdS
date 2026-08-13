## Summary

What changed and why.

## Change class (pick **one** primary)

- [ ] **Fix** — bug; default path unchanged, or fail-closed
- [ ] **Science enhancement** — opt-in; env-gated **OFF** (`flexaids::env_bool`)
- [ ] **Engine/code** — LIB/src behavior that is not a science gate
- [ ] **Docs / audit / swarm pack** — `docs/swarm/`, `docs/audit/` (own PR; do not mix with LIB)
- [ ] **Benchmark harness / frozen referee**
- [ ] **CI / hygiene**

Run `python3 scripts/classify_diff.py` on this branch. If it prints `VERDICT: MIXED`, split the PR.

## Science-Impact

- [ ] none (cannot move docked coordinates or CF ranking)
- [ ] gated OFF — flag name: `FLEXAIDDS_…`
- [ ] default-path (requires METHODOLOGY.md §1 parity)

## Scope

- [ ] Core 1.0 supported surface
- [ ] Experimental surface only
- [ ] Documentation only
- [ ] CI / release engineering
- [ ] Security hardening
- [ ] Benchmark / reproducibility

## Release impact

- [ ] affects CLI behavior
- [ ] affects Python package behavior
- [ ] affects configuration schema or defaults
- [ ] affects benchmark or metric reporting
- [ ] affects installation instructions
- [ ] no user-facing release impact

## Validation

Describe how this was validated.

- [ ] unit tests
- [ ] integration tests
- [ ] smoke benchmark
- [ ] documentation review
- [ ] manual validation

## Security and safety

- [ ] no new third-party dependency
- [ ] third-party dependency added and documented
- [ ] no known security-sensitive parsing change
- [ ] security-sensitive parsing change reviewed

## Reproducibility

- [ ] no benchmark claim involved
- [ ] benchmark claim updated with reproducibility artifact
- [ ] benchmark language remains preliminary

## Notes

Anything reviewers should pay attention to.
