#!/bin/bash
set -euo pipefail
SCRATCH="${1:-/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-86a3d1efec00/implementer}"
mkdir -p "$SCRATCH"
BASELINE="f37b02fc"

echo "=== run audit ==="
python3 scripts/v88_gate_audit.py "$SCRATCH" > /dev/null
AUDIT="$SCRATCH/gate_audit.json"
echo "audit: $(cat $AUDIT | head -c 200)"

echo "=== generate v88 ref ==="
python3 - <<PY > "$SCRATCH/v88_91pct_reference.txt"
import json
a = json.load(open("$AUDIT"))
print("v88 Astex 91.4% reference (historical claim from v88 run)")
print(f"Claim: {a['historical_claim']['rate']}% ({a['historical_claim']['n']}/85)")
print(f"Dict loose <2: {a['v88_dict']['loose_lt2']}/85 ({100*a['v88_dict']['loose_lt2']/a['v88_dict']['n']:.1f}%)")
print(f"Dict strict 0<rmsd<2: {a['v88_dict']['strict_0_lt2']}/85 ({100*a['v88_dict']['strict_0_lt2']/a['v88_dict']['n']:.1f}%)")
print(f"Seed echoes (0.00): {a['v88_dict']['seed_echo']}")
print(a['strict_note'])
print("Traceable to Published v88 numbers dict in REPRODUCIBILITY.md")
PY

echo "=== generate observed rate ==="
python3 - <<PY > "$SCRATCH/current_observed_rate.json"
import json, datetime
a = json.load(open("$AUDIT"))
c = a['current_strict']
print(json.dumps({"n":c['n'], "succ":c['succ'], "rate":c['rate'], "gate":a['gate'], "note":"NATIVE=0.0 + strict gate + raw C++ (fresh context)", "timestamp": datetime.datetime.now().isoformat()}, indent=2))
PY

echo "=== generate success verify ==="
python3 - <<PY > "$SCRATCH/success_rate_verify.txt"
import json
a = json.load(open("$AUDIT"))
print("=== VERIF (analysis, NATIVE=0, strict gate)")
print(f"v88 historical claim: {a['historical_claim']['rate']}% ({a['historical_claim']['n']}/85)")
print(f"Dict strict 0<rmsd<2: {a['v88_dict']['strict_0_lt2']}/85")
print(f"Current strict: {a['current_strict']['succ']}/85 ({a['current_strict']['rate']}%)")
print(a['strict_note'])
PY

echo "=== fresh grep for ban (live) ==="
(
grep -n -i -E 'NATIVE_SEED_FRAC|seeding|seed-echo|forbidden|0\.0' benchmarks/BENCHMARK_STANDARD.md scripts/rate_slice_verify.sh scripts/reproduce_astex85.sh | head -20
) > "$SCRATCH/seeding_ban_evidence.txt"

echo "=== ctest full output ==="
(
for t in build_verify/test_statmech build_verify/test_binding_mode_statmech; do
  [ -x $t ] && echo "=== $t ===" && $t 2>&1 | cat
done
) > "$SCRATCH/ctest_pytest_results.txt"

echo "=== scope check ==="
git diff --name-only $BASELINE..HEAD > "$SCRATCH/analysis_scope_check.txt" || true
python3 - <<PY
import sys
allowed = {"REPRODUCIBILITY.md", "BENCHMARK_STANDARD.md"}
with open("$SCRATCH/analysis_scope_check.txt") as f:
  bad = [l.strip() for l in f if l.strip() and not any(a in l for a in ["scripts/", "REPRO", "BENCHMARK", "plan.md", "docs/"])]
print("bad files (if any):", bad)
print("PASS" if not bad else "FAIL")
PY > "$SCRATCH/scope_pass.txt" 2>&1 || true

echo "=== clean stale ==="
rm -f "$SCRATCH/v88_reproduced_results.csv" 2>/dev/null || true

echo "=== done, ls ==="
ls -l "$SCRATCH"/*.txt "$SCRATCH"/*.json 2>/dev/null | cat
