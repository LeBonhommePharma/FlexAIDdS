#!/usr/bin/env bash
# ============================================================================
# OVERNIGHT CHAIN v2 — stability-gated. Runs, in order:
#   1. contention probe (CON-first, then ISO)      ~1.2 h
#   2. water ablation, 3-wide, resumes from 11     ~20.5 h
#
# WHY v2 REPLACES v1's ABSOLUTE-QUIET GATE.
# v1 required load1 < 4.0 on 3 consecutive polls. MEASURED over 10+ polls the
# minimum load1 on this box was 7.70, with NOTHING writing to disk anywhere on
# the volume, no fresh .o files and no docking output. LP reports the other
# agent runs 3 threads, which accounts for ~3 of ~9.5; the remaining ~6.5 is a
# machine floor I could not attribute or remove. So 4.0 was UNREACHABLE, and
# v1 would have burned its 8 h deadline launching nothing. A gate that cannot
# pass is as useless as one that cannot fail.
#
# THE CORRECT PRECONDITION IS A STABLE BASELINE, NOT A SILENT ONE.
# What breaks attribution in a contention experiment is a baseline that VARIES
# between waves -- then a pose divergence cannot be pinned on my own
# concurrency. A steady baseline still leaves a real manipulation: the ISO arm
# adds 1 engine process x 3 OMP threads, the CON arm adds 3 x 3. Every cell
# receipt records load1, so the contrast is reported from measurement.
#
# GATE = over STABLE_SAMPLES consecutive polls: coefficient of variation of
#        load1 < STABLE_CV, AND mean load1 < LOAD_CEILING (a saturated box
#        would make BOTH arms contention-bound and flatten the contrast),
#        AND no file written in the probe batch (the orphaned engine from the
#        interrupted launch must be finished).
#
# THE BUG THIS FIXES BEFORE ANYTHING RUNS. The interrupted launch left
# ISO_1P2Y_A mid-flight with 30 PARTIAL pose files and no DONE receipt (its
# parent was killed, so cell() never wrote one -- correctly recorded as not-run
# rather than partial). But cell() skips a cell only when DONE exists, so the
# ISO arm would re-run INTO THE SAME DIRECTORY, and posehash() globs every
# *.pdb there -- hashing orphan poses together with fresh ones, in the arm
# meant to be the clean control. Any incomplete cell dir is RENAMED ASIDE with
# a timestamp first. Renamed, never deleted: these dirs carry provenance.
#
# STOP LEVER: touch $R/state/CHAIN2_STOP   (its own name, so v1's CHAIN_STOP
# retires v1 without also stopping v2.)
# ============================================================================
set -u
R=/Users/lp.more/flexaidds_results
LOG=$R/state/overnight_chain_v2.log
LOAD_CEILING=15.0
STABLE_CV=0.35
STABLE_SAMPLES=5
POLL=60
DEADLINE=$(( $(date +%s) + 8*3600 ))

say(){ echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }

say "=== OVERNIGHT CHAIN v2 ARMED ==="
say "  stage 1: run_contention_probe_confirst.sh   (CON first, then ISO)"
say "  stage 2: run_water_ablation_w3.sh           (168 cells, 3-wide, resumes from 11)"
say "  gate: load1 cv < $STABLE_CV over $STABLE_SAMPLES polls AND mean < $LOAD_CEILING AND batch not writing"
say "  stop lever: touch $R/state/CHAIN2_STOP"

# ------------------------------------------------- wait for a STABLE baseline
SAMPLES=""
passed=0
while [ "$(date +%s)" -lt "$DEADLINE" ]; do
  if [ -f "$R/state/CHAIN2_STOP" ]; then say "STOP sentinel before stage 1 — exiting"; exit 0; fi
  P=$(ls -d "$R"/contention_probe_* 2>/dev/null | sort | tail -1)
  a=$(find "$P" -type f -exec stat -f %m {} + 2>/dev/null | sort -n | tail -1)
  l1=$(uptime | sed 's/.*load average[s]*: //' | awk '{print $1}' | tr -d ,)
  mem=$(vm_stat | awk '/page size of/{ps=$8} /Pages free/{f=$3} /Pages inactive/{i=$3} /Pages purgeable/{p=$3}
        END{gsub(/\./,"",f);gsub(/\./,"",i);gsub(/\./,"",p); printf "%.2f",(f+i+p)*ps/1073741824}')
  sleep "$POLL"
  b=$(find "$P" -type f -exec stat -f %m {} + 2>/dev/null | sort -n | tail -1)
  batch_quiet=0
  [ "${b:-0}" -le "${a:-0}" ] && batch_quiet=1

  SAMPLES=$(echo $SAMPLES $l1 | awk -v k="$STABLE_SAMPLES" \
            '{s="";for(i=(NF>k?NF-k+1:1);i<=NF;i++)s=s" "$i; print s}')
  n=$(echo $SAMPLES | wc -w | tr -d ' ')
  stats=$(echo $SAMPLES | awk '{mn=$1;mx=$1;t=0
          for(i=1;i<=NF;i++){if($i<mn)mn=$i; if($i>mx)mx=$i; t+=$i}
          m=t/NF; printf "%.2f %.2f %.2f %.3f", mn, mx, m, (m>0?(mx-mn)/m:9)}')
  set -- $stats
  smin=$1; smax=$2; smean=$3; scv=$4
  stable=0
  if [ "$n" -ge "$STABLE_SAMPLES" ]; then
    stable=$(awk -v cv="$scv" -v c="$STABLE_CV" -v m="$smean" -v ceil="$LOAD_CEILING" \
             'BEGIN{print (cv<c && m<ceil)?1:0}')
  fi
  say "  poll load1=$l1 mem=${mem}GB batch_quiet=$batch_quiet | window n=$n min=$smin max=$smax mean=$smean cv=$scv -> stable=$stable"
  if [ "$stable" -eq 1 ] && [ "$batch_quiet" -eq 1 ]; then
    say "GATE PASSED — baseline stable: mean load1 $smean, cv $scv over $n polls; batch quiet"
    say "  NOTE FOR THE WRITE-UP: this baseline is NOT idle. Per-cell load1 is recorded in"
    say "  every receipt, and the ISO-vs-CON contrast must be reported against those figures."
    echo "$smean" > "$R/state/contention_baseline_load"
    passed=1
    break
  fi
done

if [ "$passed" -ne 1 ]; then
  say "DEADLINE: 8 h elapsed and the baseline never stabilised (cv < $STABLE_CV, mean < $LOAD_CEILING)."
  say "  Nothing was launched. Re-arm by hand, or raise STABLE_CV if the box is genuinely this noisy."
  exit 0
fi

# --------------------------------------------- rename aside any incomplete cell
P=$(ls -d "$R"/contention_probe_* 2>/dev/null | sort | tail -1)
if [ -n "${P:-}" ] && [ -d "$P/run" ]; then
  TS=$(date -u +%Y%m%d_%H%M%S)
  for d in "$P/run"/*; do
    [ -d "$d" ] || continue
    case "$(basename "$d")" in *.orphan_*) continue ;; esac
    if [ ! -f "$d/DONE" ]; then
      np=$(find "$d" -name '*.pdb' 2>/dev/null | grep -vc _INI || echo 0)
      mv "$d" "${d}.orphan_$TS"
      say "  RENAMED ASIDE (no DONE, $np partial poses; would have contaminated posehash): $(basename "$d")"
    fi
  done
fi

# ------------------------------------------------------------------- stage 1
if [ -f "$R/state/CHAIN2_STOP" ]; then say "STOP sentinel — not starting stage 1"; exit 0; fi
say "=== STAGE 1: CONTENTION PROBE ==="
bash "$R/state/run_contention_probe_confirst.sh"
rc1=$?
say "=== STAGE 1 EXIT rc=$rc1 ==="
P=$(ls -d "$R"/contention_probe_* 2>/dev/null | sort | tail -1)
say "  cells: $(find "$P/run" -name DONE 2>/dev/null | wc -l | tr -d ' ')/12"

# ------------------------------------------------------------------- stage 2
if [ -f "$R/state/CHAIN2_STOP" ]; then say "STOP sentinel — not starting stage 2"; exit 0; fi
if [ "$rc1" -ne 0 ]; then
  say "STAGE 1 FAILED rc=$rc1 — NOT starting the 20.5 h water ablation on a failed precondition."
  say "  (stage 1's gates refuse on a wrong binary or a non-overlapping CON wave; either is a real stop.)"
  exit "$rc1"
fi
say "=== STAGE 2: WATER ABLATION 3-WIDE ==="
bash "$R/state/run_water_ablation_w3.sh"
rc2=$?
say "=== STAGE 2 EXIT rc=$rc2 ==="
say "=== OVERNIGHT CHAIN v2 COMPLETE ==="
