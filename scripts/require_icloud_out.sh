#!/usr/bin/env bash
# require_icloud_out.sh — fail if $1 (or $OUT) is not under iCloud Drive.
#
# Usage:
#   source scripts/require_icloud_out.sh
#   require_icloud_out "$OUT"
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0

require_icloud_out() {
  local path="${1:-${OUT:-}}"
  if [[ -z "$path" ]]; then
    echo "REFUSE: empty output path" >&2
    return 91
  fi
  # Local-first anti-hang: allow APFS OUT when explicitly enabled
  if [[ "${FLEXAIDDS_ALLOW_LOCAL_OUT:-0}" == "1" ]]; then
    case "$path" in
      *"/Mobile Documents/com~apple~CloudDocs/"*)
        echo "WARN: local-first mode but OUT is still on CloudDocs: $path" >&2
        echo "  Prefer \$FLEXAIDDS_LOCAL_ROOT (see use_local_first_benchmark_storage.sh)" >&2
        ;;
    esac
    return 0
  fi
  case "$path" in
    *"/Mobile Documents/com~apple~CloudDocs/"*)
      return 0
      ;;
    *)
      echo "REFUSE: benchmark OUT must be under iCloud Drive (or set FLEXAIDDS_ALLOW_LOCAL_OUT=1)" >&2
      echo "  got: $path" >&2
      echo "  expected prefix: \$HOME/Library/Mobile Documents/com~apple~CloudDocs/" >&2
      echo "  local-first: source scripts/use_local_first_benchmark_storage.sh" >&2
      echo "  iCloud-only: source scripts/use_icloud_benchmark_storage.sh" >&2
      return 91
      ;;
  esac
}

# If executed (not sourced) with an argument, check and exit
if [[ "${BASH_SOURCE[0]:-}" == "${0:-}" ]]; then
  require_icloud_out "${1:-}"
  exit $?
fi
