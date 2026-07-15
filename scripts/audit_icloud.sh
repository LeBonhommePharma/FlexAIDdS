#!/usr/bin/env bash
# =============================================================================
# audit_icloud.sh — Robust iCloud Drive audit for FlexAIDdS
#
# Scans primary iCloud locations for common sync problems:
#   - Conflicted copies (Finder "conflicted copy")
#   - Placeholder files (*.icloud — evicted / not downloaded)
#   - Zero-byte files
#   - Suspiciously small files inside archive/ directories
#   - Recent (last 5 min) file activity without a .verified marker
#   - Quota / usage summaries (du + df)
#
# Usage:
#   ./scripts/audit_icloud.sh
#   ./scripts/audit_icloud.sh --json
#   ./scripts/audit_icloud.sh --fix     # reports only (read-only for now)
#
# Environment:
#   FLEXAIDDS_ICLOUD   Primary base (if unset: .../FlexAIDdS_benchmarks under iCloud)
#
# Always also scans sibling FlexAIDdS / flexaidds variants under the iCloud container.
#
# This script is read-only / idempotent. --fix currently only reports.
# Apache-2.0 (c) 2026 NRGlab, Université de Montréal
# =============================================================================

set -euo pipefail

# ─── Colors & helpers ─────────────────────────────────────────────────────────

if [[ -t 1 ]]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    CYAN='\033[0;36m'
    BOLD='\033[1m'
    NC='\033[0m'
else
    RED='' GREEN='' YELLOW='' CYAN='' BOLD='' NC=''
fi

info()  { $QUIET || printf "${CYAN}[INFO]${NC}  %s\n" "$*"; }
ok()    { $QUIET || printf "${GREEN}[OK]${NC}    %s\n" "$*"; }
warn()  { $QUIET || printf "${YELLOW}[WARN]${NC}  %s\n" "$*"; }
err()   { printf "${RED}[ERR]${NC}   %s\n" "$*" >&2; }

usage() {
    cat <<'EOF'
audit_icloud.sh — FlexAIDdS iCloud sync auditor

Options:
  --json     Emit machine-readable JSON summary (no human report)
  --fix      Request remediation hints (currently read-only: only reports)
  -h, --help Show this help

Environment variables respected:
  FLEXAIDDS_ICLOUD   Override primary iCloud base directory
                     (falls back to standard iCloud/FlexAIDdS_benchmarks)

Examples:
  FLEXAIDDS_ICLOUD=/path/to/mybench ./scripts/audit_icloud.sh --json
  ./scripts/audit_icloud.sh
EOF
}

# ─── Arg parsing ──────────────────────────────────────────────────────────────

JSON_MODE=false
FIX_MODE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --json) JSON_MODE=true; shift ;;
        --fix)  FIX_MODE=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) err "Unknown argument: $1"; usage; exit 1 ;;
    esac
done

QUIET=false
if $JSON_MODE; then QUIET=true; fi

# ─── Resolve scan locations (env var first, no hardcoded usernames) ───────────

ICLOUD_CONTAINER="$HOME/Library/Mobile Documents/com~apple~CloudDocs"

# Default per task spec + safe_archive convention
DEFAULT_BASE="${ICLOUD_CONTAINER}/FlexAIDdS_benchmarks"
PRIMARY="${FLEXAIDDS_ICLOUD:-$DEFAULT_BASE}"

declare -a SCAN_LOCS=()

add_loc() {
    local p="$1"
    [[ -z "$p" ]] && return
    # Only add if it exists or is the explicitly requested primary (so we can report missing)
    if [[ -d "$p" || "$p" == "$PRIMARY" ]]; then
        SCAN_LOCS+=("$p")
    fi
}

add_loc "$PRIMARY"

# Always consider the explicitly set FLEXAIDDS_ICLOUD (even if different)
if [[ -n "${FLEXAIDDS_ICLOUD:-}" ]]; then
    add_loc "$FLEXAIDDS_ICLOUD"
fi

# Standard sibling variants (FlexAIDdS, flexaidds, etc.)
add_loc "${ICLOUD_CONTAINER}/FlexAIDdS"
add_loc "${ICLOUD_CONTAINER}/flexaidds"
add_loc "${ICLOUD_CONTAINER}/FlexAIDdS_benchmarks"

# Discover other flex* / benchmark dirs at container root (robust to renames)
if [[ -d "$ICLOUD_CONTAINER" ]]; then
    while IFS= read -r -d '' d; do
        add_loc "$d"
    done < <(find "$ICLOUD_CONTAINER" -maxdepth 1 -type d \
        \( -iname '*flexaidds*' -o -iname '*flexaid*' -o -iname '*benchmarks*' \) \
        -not -path '*/.*' -print0 2>/dev/null || true)
fi

# Deduplicate while preserving order (space-safe, no word-split)
if [[ ${#SCAN_LOCS[@]} -gt 0 ]]; then
    tmp_dedup="$(mktemp)"
    printf "%s\n" "${SCAN_LOCS[@]}" > "$tmp_dedup"
    SCAN_LOCS=()
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -n "$line" ]] && SCAN_LOCS+=("$line")
    done < <(awk '!seen[$0]++' "$tmp_dedup" 2>/dev/null || cat "$tmp_dedup")
    rm -f "$tmp_dedup"
fi

if [[ ${#SCAN_LOCS[@]} -eq 0 ]]; then
    warn "No iCloud locations found or resolvable."
    SCAN_LOCS=("$PRIMARY")   # still report on the target primary
fi

# ─── Data collection (read-only) ──────────────────────────────────────────────

TOTAL_CONFLICTS=0
TOTAL_PLACEHOLDERS=0
TOTAL_ZEROS=0
TOTAL_SMALL_ARCHIVE=0
TOTAL_RECENT_UNVER=0
TOTAL_RECENT_UNVER_DIRS=0

declare -a CONFLICT_EXAMPLES=()
declare -a PLACEHOLDER_EXAMPLES=()
declare -a ZERO_EXAMPLES=()
declare -a SMALL_ARCHIVE_EXAMPLES=()
declare -a RECENT_UNVER_EXAMPLES=()

# Collect usage info
declare -a USAGE_LINES=()

scan_location() {
    local root="$1"
    local label="${2:-$root}"

    if [[ ! -d "$root" ]]; then
        $QUIET || warn "Location not present: $root"
        return 0
    fi

    info "Scanning: $root"

    # Conflicted copies (files and dirs)
    local conf_count=0
    local conf_list
    conf_list=$(find "$root" \( -iname '*conflicted*' -o -iname '*(*conflicted*' \) -print 2>/dev/null | head -30 || true)
    if [[ -n "$conf_list" ]]; then
        conf_count=$(printf "%s\n" "$conf_list" | wc -l | tr -d ' ')
        while IFS= read -r line; do
            CONFLICT_EXAMPLES+=("[$label] $line")
        done <<< "$conf_list"
    fi
    TOTAL_CONFLICTS=$((TOTAL_CONFLICTS + conf_count))

    # Placeholders (*.icloud)
    local ph_count=0
    local ph_list
    ph_list=$(find "$root" -name '*.icloud' -print 2>/dev/null | head -30 || true)
    if [[ -n "$ph_list" ]]; then
        ph_count=$(printf "%s\n" "$ph_list" | wc -l | tr -d ' ')
        while IFS= read -r line; do
            PLACEHOLDER_EXAMPLES+=("[$label] $line")
        done <<< "$ph_list"
    fi
    TOTAL_PLACEHOLDERS=$((TOTAL_PLACEHOLDERS + ph_count))

    # Zero-byte files
    local z_count=0
    local z_list
    z_list=$(find "$root" -type f -size 0 -print 2>/dev/null | head -20 || true)
    if [[ -n "$z_list" ]]; then
        z_count=$(printf "%s\n" "$z_list" | wc -l | tr -d ' ')
        while IFS= read -r line; do
            ZERO_EXAMPLES+=("[$label] $line")
        done <<< "$z_list"
    fi
    TOTAL_ZEROS=$((TOTAL_ZEROS + z_count))

    # Empty (zero-byte) files discovered inside archive trees (truly suspicious for incomplete archives).
    # Separate from global zero count for emphasis in "recent archives".
    local sa_count=0
    local sa_list=""
    if command -v find >/dev/null 2>&1; then
        while IFS= read -r adir; do
            [[ -d "$adir" ]] || continue
            local empties
            empties=$(find "$adir" -type f -size 0 -print 2>/dev/null | head -10 || true)
            if [[ -n "$empties" ]]; then
                sa_count=$((sa_count + $(printf "%s\n" "$empties" | wc -l | tr -d ' ')))
                sa_list+="${empties}"$'\n'
            fi
        done < <(find "$root" -type d \( -iname '*archive*' -o -iname '*archiv*' \) -print 2>/dev/null | head -15 || true)
    fi
    TOTAL_SMALL_ARCHIVE=$((TOTAL_SMALL_ARCHIVE + sa_count))
    if [[ -n "$sa_list" ]]; then
        while IFS= read -r line; do
            [[ -n "$line" ]] && SMALL_ARCHIVE_EXAMPLES+=("[$label] $line")
        done <<< "$sa_list"
    fi

    # Recent (last 5 minutes) modified files without .verified marker
    local ru_count=0
    local ru_list
    ru_list=$(find "$root" -type f -mmin -5 ! -iname '*.verified' -print 2>/dev/null | head -20 || true)
    if [[ -n "$ru_list" ]]; then
        ru_count=$(printf "%s\n" "$ru_list" | wc -l | tr -d ' ')
        while IFS= read -r line; do
            RECENT_UNVER_EXAMPLES+=("[$label] $line")
        done <<< "$ru_list"
    fi
    TOTAL_RECENT_UNVER=$((TOTAL_RECENT_UNVER + ru_count))

    # Recently touched directories that lack a .verified marker at top level
    local rud_count=0
    if command -v find >/dev/null 2>&1; then
        while IFS= read -r d; do
            if [[ -d "$d" ]]; then
                if [[ ! -f "$d/.verified" ]]; then
                    rud_count=$((rud_count + 1))
                    if [[ ${#RECENT_UNVER_EXAMPLES[@]} -lt 25 ]]; then
                        RECENT_UNVER_EXAMPLES+=("[$label DIR] $d (no .verified)")
                    fi
                fi
            fi
        done < <(find "$root" -type d -mmin -5 -print 2>/dev/null | head -30 || true)
    fi
    TOTAL_RECENT_UNVER_DIRS=$((TOTAL_RECENT_UNVER_DIRS + rud_count))

    # Usage for this location
    local du_line
    du_line=$(du -sh "$root" 2>/dev/null || echo "N/A $root")
    USAGE_LINES+=("$du_line")
}

# Perform scans
for loc in "${SCAN_LOCS[@]}"; do
    scan_location "$loc"
done

# Volume / quota summaries (best effort)
df_lines=$(df -h "$HOME" 2>/dev/null | cat || echo "df unavailable")
# Try to get iCloud-specific if visible (often on the Data volume)
df_ic=$(df -h "$ICLOUD_CONTAINER" 2>/dev/null | tail -1 || true)
if [[ -n "$df_ic" ]]; then
    df_lines="$df_lines"$'\n'"iCloud container: $df_ic"
fi

# ─── Reporting ────────────────────────────────────────────────────────────────

timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

if $JSON_MODE; then
    # Generate clean, properly escaped JSON via Python (handles paths with unicode/spaces)
    export AUDIT_TS="$timestamp"
    export C_CONFLICTS="$TOTAL_CONFLICTS"
    export C_PLACE="$TOTAL_PLACEHOLDERS"
    export C_ZERO="$TOTAL_ZEROS"
    export C_SMALL="$TOTAL_SMALL_ARCHIVE"
    export C_RECENT_F="$TOTAL_RECENT_UNVER"
    export C_RECENT_D="$TOTAL_RECENT_UNVER_DIRS"
    export FIX_MODE="$FIX_MODE"

    # Write arrays to temp files for python consumption (space/unicode safe)
    AUDIT_LOCS_FILE=$(mktemp)
    AUDIT_USAGE_FILE=$(mktemp)
    AUDIT_DF_FILE=$(mktemp)
    EX_CONF_FILE=$(mktemp)
    EX_PH_FILE=$(mktemp)
    EX_Z_FILE=$(mktemp)
    EX_SA_FILE=$(mktemp)
    EX_RU_FILE=$(mktemp)

    export AUDIT_LOCS_FILE AUDIT_USAGE_FILE AUDIT_DF_FILE EX_CONF_FILE EX_PH_FILE EX_Z_FILE EX_SA_FILE EX_RU_FILE

    printf "%s\n" "${SCAN_LOCS[@]:-}" > "$AUDIT_LOCS_FILE"
    printf "%s\n" "${USAGE_LINES[@]:-}" > "$AUDIT_USAGE_FILE"
    printf "%s\n" "$df_lines" > "$AUDIT_DF_FILE"
    printf "%s\n" "${CONFLICT_EXAMPLES[@]:0:8}" > "$EX_CONF_FILE"
    printf "%s\n" "${PLACEHOLDER_EXAMPLES[@]:0:8}" > "$EX_PH_FILE"
    printf "%s\n" "${ZERO_EXAMPLES[@]:0:5}" > "$EX_Z_FILE"
    printf "%s\n" "${SMALL_ARCHIVE_EXAMPLES[@]:0:5}" > "$EX_SA_FILE"
    printf "%s\n" "${RECENT_UNVER_EXAMPLES[@]:0:8}" > "$EX_RU_FILE"

    python3 - <<'PYEOF'
import json, os, sys
data = {
  "timestamp": os.environ.get("AUDIT_TS", ""),
  "scanned_locations": [],
  "counts": {
    "conflicted_copies": int(os.environ.get("C_CONFLICTS", 0)),
    "placeholders_icloud": int(os.environ.get("C_PLACE", 0)),
    "zero_byte_files": int(os.environ.get("C_ZERO", 0)),
    "small_files_in_archives": int(os.environ.get("C_SMALL", 0)),
    "recent_unverified_files": int(os.environ.get("C_RECENT_F", 0)),
    "recent_dirs_without_verified": int(os.environ.get("C_RECENT_D", 0)),
  },
  "usage": {"du_summaries": [], "df": "unavailable"},
  "examples": {"conflicts": [], "placeholders": [], "zeros": [], "small_archive": [], "recent_unverified": []},
  "fix_mode": os.environ.get("FIX_MODE", "false").lower() == "true",
  "notes": "Read-only audit. --fix is currently a no-op (reporting only)."
}

for k, env in [("scanned_locations", "AUDIT_LOCS_FILE"), ("du_summaries", "AUDIT_USAGE_FILE")]:
    p = os.environ.get(env)
    if p and os.path.exists(p):
        with open(p) as f:
            vals = [ln.rstrip("\n") for ln in f if ln.strip()]
            if k == "du_summaries":
                data["usage"]["du_summaries"] = vals
            else:
                data[k] = vals

p = os.environ.get("AUDIT_DF_FILE")
if p and os.path.exists(p):
    with open(p) as f:
        data["usage"]["df"] = f.read().strip()

exmap = {"conflicts":"EX_CONF_FILE", "placeholders":"EX_PH_FILE", "zeros":"EX_Z_FILE", "small_archive":"EX_SA_FILE", "recent_unverified":"EX_RU_FILE"}
for k, env in exmap.items():
    p = os.environ.get(env)
    if p and os.path.exists(p):
        with open(p) as f:
            data["examples"][k] = [ln.rstrip("\n") for ln in f if ln.strip()]

print(json.dumps(data, indent=2, ensure_ascii=False))
PYEOF

    # cleanup temps
    rm -f "$AUDIT_LOCS_FILE" "$AUDIT_USAGE_FILE" "$AUDIT_DF_FILE" "$EX_CONF_FILE" "$EX_PH_FILE" "$EX_Z_FILE" "$EX_SA_FILE" "$EX_RU_FILE" 2>/dev/null || true
    exit 0
fi

# Human-readable report
echo "======================================================================"
echo "  FlexAIDdS iCloud Audit Report"
echo "  Generated: $timestamp (UTC)"
echo "======================================================================"
echo ""
echo "${BOLD}Scanned locations:${NC}"
for l in "${SCAN_LOCS[@]}"; do
    if [[ -d "$l" ]]; then
        printf "  • %s\n" "$l"
    else
        printf "  • %s  ${YELLOW}(not present)${NC}\n" "$l"
    fi
done
echo ""

echo "${BOLD}Detection summary:${NC}"
printf "  Conflicted copies (\"*conflicted*\" etc.):   %s\n" "$TOTAL_CONFLICTS"
printf "  Placeholder files (*.icloud):                 %s\n" "$TOTAL_PLACEHOLDERS"
printf "  Zero-byte files:                              %s\n" "$TOTAL_ZEROS"
printf "  Empty files inside archive trees:                   %s\n" "$TOTAL_SMALL_ARCHIVE"
printf "  Files modified <5min ago (no .verified):      %s\n" "$TOTAL_RECENT_UNVER"
printf "  Recent dirs lacking .verified marker:         %s\n" "$TOTAL_RECENT_UNVER_DIRS"
echo ""

# Examples sections (only if >0)
print_examples() {
    local title="$1"; shift
    local -a arr=()
    if [[ $# -gt 0 ]]; then
        arr=("$@")
    fi
    if [[ ${#arr[@]} -gt 0 ]]; then
        echo "${BOLD}${title}:${NC}"
        local shown=0
        for e in "${arr[@]}"; do
            printf "  - %s\n" "$e"
            shown=$((shown+1))
            [[ $shown -ge 8 ]] && { echo "    ... (more omitted)"; break; }
        done
        echo ""
    fi
}

print_examples "Conflicted copies (examples)" "${CONFLICT_EXAMPLES[@]:-}"
print_examples "Placeholder (*.icloud) examples" "${PLACEHOLDER_EXAMPLES[@]:-}"
print_examples "Zero-byte file examples" "${ZERO_EXAMPLES[@]:-}"
print_examples "Small files inside archive trees" "${SMALL_ARCHIVE_EXAMPLES[@]:-}"
print_examples "Recent activity without .verified marker" "${RECENT_UNVER_EXAMPLES[@]:-}"

echo "${BOLD}Storage usage:${NC}"
for u in "${USAGE_LINES[@]}"; do
    printf "  %s\n" "$u"
done
echo ""
echo "${BOLD}Volume info (df):${NC}"
printf "%s\n" "$df_lines" | sed 's/^/  /'
echo ""

echo "======================================================================"
echo "${BOLD}Recommendations:${NC}"
if [[ $TOTAL_CONFLICTS -gt 0 ]]; then
    echo "  • Conflicted copies detected: Open in Finder, choose the version to keep,"
    echo "    then delete the conflicting copy. Re-run audit after."
fi
if [[ $TOTAL_PLACEHOLDERS -gt 0 ]]; then
    echo "  • *.icloud placeholders: Files are not fully downloaded (Optimize Storage)."
    echo "    Right-click → Download Now, or in Finder iCloud Drive settings disable"
    echo "    optimization for large benchmark data. Consider 'brctl download'."
fi
if [[ $TOTAL_ZEROS -gt 0 || $TOTAL_SMALL_ARCHIVE -gt 0 ]]; then
    echo "  • Zero or tiny files in archives: Possible incomplete writes or sync cuts."
    echo "    Re-archive using scripts/safe_archive_to_icoud.py --verify and ensure"
    echo "    .verified marker + manifest are present."
fi
if [[ $TOTAL_RECENT_UNVER -gt 0 || $TOTAL_RECENT_UNVER_DIRS -gt 0 ]]; then
    echo "  • Recent modifications without .verified: Sync may still be in flight."
    echo "    Wait a few minutes, ensure writer process wrote the marker, then re-audit."
    echo "    Never treat a directory as complete until .verified exists."
fi
if [[ $TOTAL_CONFLICTS -eq 0 && $TOTAL_PLACEHOLDERS -eq 0 && $TOTAL_ZEROS -eq 0 && $TOTAL_RECENT_UNVER -eq 0 ]]; then
    ok "No obvious iCloud sync red flags detected in scanned locations."
fi
echo ""
echo "  Run with --json for machine-readable output (CI / monitoring)."
echo "  Use FLEXAIDDS_ICLOUD=/path to target a specific tree."
echo "  This script is read-only. --fix currently prints this report only."
echo "======================================================================"

if $FIX_MODE; then
    warn "--fix was supplied. Current implementation performs no mutations (reporting only)."
    info "Future versions may suggest or (with confirmation) perform safe cleanups."
fi

exit 0
