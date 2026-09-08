# Retraction: the 1P2Y 4.16 A -> 0.89 A improvement is not attributable to the
# coordinate guard (or to the local-copy fix)

**Affects commit `2f5365ac` "fix(vcontacts): pristine-coordinate guard ..."**,
whose message states: *"IMPACT IS NOT NOISE. With the guard on, 1P2Y goes from
rmsd 4.1597 A (FAILURE) ..."* and cites pose hash `e19e6e5562661333`.

## What was actually measured, later, with a control

A three-arm test on the post-rebase engine ran gates-off, coordinate-guard and
local-copy on 1P2Y with identical settings:

| arm | rc | poses | pose hash | rmsd |
|---|---|---|---|---|
| gates off | 0 | 87 | `83e4146955d66105` | 0.8939 |
| coordinate guard | 0 | 87 | `83e4146955d66105` | 0.8939 |
| local copy | 0 | 87 | `83e4146955d66105` | 0.8939 |

**All three are byte-identical, including the control.** Gates-off already
produces 0.8939 A, so on this engine the guard changes nothing on 1P2Y and the
improvement cannot be attributed to it.

## Why the original claim reads as wrong today -- and what is NOT wrong about it

**CORRECTION (second pass, external audit).** An earlier version of this note
said the 4.1597 A baseline "was never accompanied by a gates-off control on the
same engine." **That was false, and it unfairly impugned `2f5365ac`'s
methodology.** The original experiment was properly controlled: a single build,
engine `2bbd841ca4f192af`, ran an env-gated three-arm test on 1P2Y --

| arm | poses | pose hash | bestCF | rmsd |
|---|---|---|---|---|
| guard off | 128 | `866b89c0d2ea86ed` | -254.86286 | 4.1597 |
| guard on | 128 | `e19e6e5562661333` | -231.46644 | 0.8939 |
| guard on + diag | 128 | `e19e6e5562661333` | -231.46644 | 0.8939 |

Same build, same target, same settings, gate toggled. On **that** engine the
guard genuinely produced 4.1597 A -> 0.8939 A, and `2f5365ac` was entitled to
say so.

**What is actually retracted is the generalisation, not the measurement.** The
effect is ENGINE-DEPENDENT. On the post-rebase engine the defect no longer
manifests -- gates-off already succeeds -- so the guard is a no-op there and
nothing on the current engine can be credited to it. The cause of the change
came in from upstream, not from either fix.

Two operational consequences:

* Do not cite 4.1597 -> 0.8939 as evidence for the guard **on the current
  engine**. It is evidence about `2bbd841ca4f192af` only.
* `e19e6e5562661333` is not reproducible today and must not be used as a
  regression target; the current engine's value for that configuration is
  `83e4146955d66105`.

## What still stands

Nothing about the *code* argument is retracted. Perturbing shared `atoms[]`
coordinates under OpenMP is a data race, and the restore is skipped on the
`return -1` path, so a perturbation can outlive the call. Those are defects
whether or not they move any pose today. The local-copy commit removes the cause
and is behaviour-preserving on 1P2Y, which is all the measurement supports.

## Scope of the correction

- `2f5365ac` -- message left intact (it is history); this note is the correction.
- The local-copy commit was rewritten before push to remove a claim that the
  three-arm test reproduced `e19e6e5562661333` "bit for bit". It did not; every
  arm produced `83e4146955d66105`, and because the control matched too, the test
  could not have failed and is not evidence either fix works.

Found by external audit against the session archive.
