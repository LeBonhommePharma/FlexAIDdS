# ELECTION_V135 a posteriori

```
# ELECTION_V135 near-miss a posteriori (FINAL)
OUT: /Users/lp.more/flexaidds_results/election_v135_near_miss_20260726_225823
ALL_ARMS_DONE: True
experiment: ELECTION_V135_1N1M_near_miss
one_variable: FLEXAIDDS_ELECTION_V135=1 (+ defaults SCORE_TAU=25, INCLUDE_SINGLETONS) vs control unset
magnitude_floor: 1N1M elect_rmsd <= 2.5 OR (actual_elect - pool_near) reduced by >= 1.0 A; no wipeout (elect not worse by >0.5 on both)
codes: ['1N1M', '1L7F'] restarts: 5 matrix: 9dc93717dfed0698006d88dd6a9627bc

## arm_control receipt
  restarts=5 matrix_md5=9dc93717dfed0698006d88dd6a9627bc
  election_v135=False election_score_tau=0.0
  consensus_scorer=False no_sec=True
  binary_sha256=a3fa78c1ad8e7778715cd56c21f9604d51f469213d656804c23920214acc4255

## arm_v135 receipt
  restarts=5 matrix_md5=9dc93717dfed0698006d88dd6a9627bc
  election_v135=True election_score_tau=25.0
  consensus_scorer=False no_sec=True
  binary_sha256=a3fa78c1ad8e7778715cd56c21f9604d51f469213d656804c23920214acc4255

## Metrics
arm        code      elect      BCR   Δelect           cf  restart
control    1L7F     3.9907   3.9233             -157.7286        1
control    1N1M     6.3999   4.0427              -99.3141        2
v135       1L7F     3.9907   3.9233  +0.0000    -157.7286        1
v135       1N1M     6.3999   4.0427  +0.0000     -99.3141        4

## Judgment vs a priori floors
1N1M control_elect=6.3999 v135_elect=6.3999 gap_shrink=+0.0000
1L7F control_elect=3.9907 v135_elect=3.9907 d=+0.0000
accept_1N1M_elect_floor=False  (elect<=2.5 or gap_shrink>=1.0)
wipeout=False
ACCEPT_ELECTION_V135=False

## L4-ish log markers
  control: V135/election mentions=116
  v135: V135/election mentions=222

## Next flip
rule: ELECTION_V135_null → prioritize G4.3 mutation or new selection architecture
priority_order: [G4.3_mutation, new_search_arch, scoring_non_burial_optional]
NOTE: full-85 still blocked (Phase-4 sampling ACCEPT not met)

WROTE accept.txt

```
