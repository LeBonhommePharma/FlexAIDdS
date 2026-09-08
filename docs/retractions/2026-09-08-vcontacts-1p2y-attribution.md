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

## Why the original claim was wrong

The 4.1597 A "failure" baseline was measured on a **pre-rebase engine** and was
never accompanied by a gates-off control on the same engine. A before/after
across two different engines with no control does not isolate the guard. The
`e19e6e5562661333` hash belongs to that older engine and is not reproducible on
the current one -- it is not a target any present-day run should be checked
against.

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
