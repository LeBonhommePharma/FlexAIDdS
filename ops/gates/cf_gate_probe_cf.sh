#!/bin/bash
# CF scoring-regression merge gate — wraps the in-repo build/probe_cf instrument.
# Supersedes the standalone python re-implementation: scores via the ENGINE's real
# vcfunction()/score_native_pose(), not a proxy. For each panel target it scores the
# NATIVE pose and the best DECOY pose, computes ΔCF = cf_total(native) - cf_total(decoy),
# and FAILS if the inverted fraction (ΔCF>tol) exceeds MAX_INV_FRAC.
#
# GATEFIX v2 (ligand + config):
#   DEFECT 1 — PDB decoys need --ligand (crystal SDF topology) or probe_cf returns
#              empty cf_total (rc=256) and the gate silently skips the target.
#   DEFECT 2 — without --config, LOCCLF/optres pocket pruning never applies and CF
#              is ~200× inflated (whole-receptor optres). Production-equivalent CF
#              requires per-target dock_config.json on BOTH native and decoy calls.
#
# Usage: cf_gate_probe_cf.sh <panel_manifest.tsv> [tol] [max_inv_frac]
#   manifest rows (tab-sep, 5 fields):
#     pdb <TAB> receptor.pdb <TAB> native_pose.sdf <TAB> decoy_pose.pdb <TAB> dock_config.json
# Requires: build/probe_cf (built from tools/probe_cf.cpp, commit 81acfed6+).
set -u
# Resolve repo root from this script so the gate is portable.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"
PROBE="${FLEXAIDDS_PROBE_CF:-$REPO/build/probe_cf}"
MANIFEST="${1:?manifest tsv required}"
TOL="${2:-0.0}"
MAX_INV_FRAC="${3:-0.125}"

[ -x "$PROBE" ] || { echo "GATE ERROR: $PROBE not built (cmake --build build --target probe_cf)"; exit 2; }

# $1=receptor $2=pose $3=ligand_topology_sdf $4=dock_config.json  -> echoes cf_total
cf() {
  local rec="$1" pose="$2" lig="$3" cfg="$4"
  local -a args=(--receptor "$rec" --pose "$pose")
  # PDB poses lack Sybyl typing; always supply crystal/topology SDF as --ligand.
  case "${pose##*.}" in
    pdb|PDB|ent|ENT)
      if [ -z "$lig" ] || [ ! -f "$lig" ]; then
        echo "GATE ERROR: PDB pose requires ligand topology SDF (got: '${lig:-empty}')" >&2
        return 2
      fi
      args+=(--ligand "$lig")
      ;;
  esac
  if [ -z "$cfg" ] || [ ! -f "$cfg" ]; then
    echo "GATE ERROR: dock_config.json required for production CF (got: '${cfg:-empty}')" >&2
    return 2
  fi
  args+=(--config "$cfg")
  "$PROBE" "${args[@]}" 2>/dev/null \
    | sed -n 's/.*"cf_total": *\([-0-9.][-0-9.]*\).*/\1/p' | head -1
}

n=0; inv=0; rows=""
while IFS=$'\t' read -r pdb rec nat dec cfg || [ -n "${pdb:-}" ]; do
  [ -z "${pdb:-}" ] && continue
  case "$pdb" in \#*) continue;; esac
  # 5th field required for active rows (GATEFIX v2 defect 2).
  if [ -z "${cfg:-}" ]; then
    echo "GATE ERROR: $pdb missing dock_config.json (5th manifest field)" >&2
    exit 2
  fi
  # Resolve relative paths against REPO (manifest may use absolute or repo-relative).
  for v in rec nat dec cfg; do
    eval "p=\$$v"
    case "$p" in
      /*) ;;
      *) eval "$v=\"\$REPO/\$p\"" ;;
    esac
  done
  for p in "$rec" "$nat" "$dec" "$cfg"; do
    [ -f "$p" ] || { echo "GATE ERROR: $pdb missing file: $p" >&2; exit 2; }
  done

  cn=$(cf "$rec" "$nat" "$nat" "$cfg") || true
  cd_=$(cf "$rec" "$dec" "$nat" "$cfg") || true
  if [ -z "$cn" ] || [ -z "$cd_" ]; then
    echo "  skip $pdb (probe_cf gave no cf_total)"
    continue
  fi
  d=$(echo "$cn - $cd_" | bc -l)
  bad=$(echo "$d > $TOL" | bc -l)
  n=$((n+1)); [ "$bad" = 1 ] && inv=$((inv+1))
  rows="$rows$(printf '  %-8s dCF=%+8.2f (native %s vs decoy %s)%s\n' "$pdb" "$d" "$cn" "$cd_" "$([ "$bad" = 1 ] && echo '  INVERTED')")\n"
done < "$MANIFEST"

[ "$n" -eq 0 ] && { echo "GATE ERROR: no targets scored"; exit 2; }
frac=$(echo "scale=3; $inv / $n" | bc -l)
printf '%b' "$rows"
echo "CF-REGRESSION GATE | n=$n inverted=$inv frac=$frac (threshold $MAX_INV_FRAC)"
pass=$(echo "$frac <= $MAX_INV_FRAC" | bc -l)
if [ "$pass" = 1 ]; then echo "  VERDICT: PASS"; exit 0; else echo "  VERDICT: FAIL"; exit 1; fi
