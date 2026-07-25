# W2.2 — COM_BURIAL_CAP=-130 is **not** product default

**Decision:** Do **not** merge or ship `FLEXAIDDS_COM_BURIAL_CAP=-130` as default.

**Reasons (on-disk):**
- Run `v_comcap_fixed_20260724_133920` is **UNCITABLE** (OOM workers=6, 10/85, empty git_commit).
- Cap value ≈ 1G9V native com (−129) — single-target tuning.
- Per-optres/total CF still reached −1144 class — does not globally bound com.
- Live autonomous CF magnitudes already moderated without CAP=-130.

**Allowed later:** redesigned **global** com bound behind env, multi-target oracle, full SCORING_PROVENANCE, after wall (W2.1) clarity.

**COM_FLOOR / VCT_NORM:** optional serial canary only after wall design; not packaged with Softβ+72d7+autonomous.
