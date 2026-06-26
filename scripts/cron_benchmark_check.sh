#!/usr/bin/env bash
# cron_benchmark_check.sh — lightweight cron watchdog for benchmark campaigns.
#
# - Logs a one-line status snapshot every invocation
# - Restarts monitor_benchmark_campaign.py if it died
# - Triggers post-run audits if a campaign hit 85/85 but audits missing
#
# Safe: never signals benchmark_datasets or FlexAIDdS processes.
#
# Crontab (every 5 min):
#   */5 * * * * /path/to/scripts/cron_benchmark_check.sh >> ~/flexaidds_results/cron_benchmark_check.log 2>&1
#
set -euo pipefail

REPO="/Users/lp.more/.grok/worktrees/projects-flexaidds/opus-48-science-fixes"
V111="${HOME}/flexaidds_results/v111_science_20260626_0613"
BASE="${HOME}/flexaidds_results/baseline_8196829_audit"
MONITOR_LOG="${HOME}/flexaidds_results/benchmark_monitor.log"
MONITOR_PY="${REPO}/scripts/monitor_benchmark_campaign.py"
TOTAL=85

utc() { date -u +"%Y-%m-%dT%H:%M:%SZ"; }

count_done() {
  local root="$1"
  find "$root" -maxdepth 2 -name result.csv 2>/dev/null | wc -l | tr -d ' '
}

count_sub2() {
  local root="$1"
  python3 - "$root" <<'PY'
import sys, glob, csv, os
root = sys.argv[1]
n = 0
for rp in glob.glob(os.path.join(root, "*/result.csv")):
    try:
        with open(rp) as f:
            row = next(csv.DictReader(f))
        r = float(row.get("rmsd_hungarian") or row.get("rmsd_to_crystal") or -1)
        if 0 <= r < 2:
            n += 1
    except Exception:
        pass
print(n)
PY
}

snapshot() {
  local label="$1" root="$2"
  if [[ ! -d "$root" ]]; then
    echo "$(utc) [cron] ${label}: MISSING"
    return
  fi
  local done sub2
  done=$(count_done "$root")
  sub2=$(count_sub2 "$root")
  local fail=$(( done - sub2 ))
  [[ $fail -lt 0 ]] && fail=0
  local pct
  pct=$(python3 -c "print(f'{100*${done}/${TOTAL}:.1f}')")
  echo "$(utc) [cron] ${label}: ${done}/${TOTAL} (${pct}%) sub2=${sub2} fail=${fail}"
}

# ── 1. Status snapshot ───────────────────────────────────────────────────────
snapshot "v111_science" "$V111"
snapshot "baseline_8196829" "$BASE"

# ── 2. Ensure long-running monitor is alive ──────────────────────────────────
if ! pgrep -f "monitor_benchmark_campaign.py" >/dev/null 2>&1; then
  echo "$(utc) [cron] monitor DEAD — restarting"
  nohup python3 "$MONITOR_PY" --interval 90 >> "$MONITOR_LOG" 2>&1 &
  echo "$(utc) [cron] monitor restarted pid=$!"
fi

# ── 3. Safety-net audits if monitor missed completion ───────────────────────
for pair in "v111_science:${V111}" "baseline_8196829:${BASE}"; do
  label="${pair%%:*}"
  root="${pair##*:}"
  [[ -d "$root" ]] || continue
  done=$(count_done "$root")
  marker="${root}/.monitor_audit_done"
  if [[ "$done" -ge "$TOTAL" && ! -f "$marker" ]]; then
    echo "$(utc) [cron] ${label} complete — triggering audits"
    python3 "${REPO}/scripts/failure_classify.py" "$root" || true
    python3 "${REPO}/scripts/cf_ground_truth_audit.py" "$root" || true
    date -u +"%Y-%m-%dT%H:%M:%SZ" > "$marker"
  fi
done

# ── 4. Milestone alerts (25/50/75/85) ───────────────────────────────────────
MILESTONE_FILE="${HOME}/flexaidds_results/.cron_milestones.log"
for pair in "v111:${V111}" "baseline:${BASE}"; do
  tag="${pair%%:*}"
  root="${pair##*:}"
  [[ -d "$root" ]] || continue
  done=$(count_done "$root")
  for m in 25 50 75 85; do
    key="${tag}_${m}"
    if [[ "$done" -ge "$m" ]] && ! grep -q "^${key}:" "$MILESTONE_FILE" 2>/dev/null; then
      echo "${key}:$(utc)" >> "$MILESTONE_FILE"
      echo "$(utc) [cron] MILESTONE ${tag} reached ${m}/85"
    fi
  done
done