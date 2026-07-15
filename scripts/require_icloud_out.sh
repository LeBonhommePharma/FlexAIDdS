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
  case "$path" in
    *"/Mobile Documents/com~apple~CloudDocs/"*)
      return 0
      ;;
    *)
      echo "REFUSE: benchmark OUT must be under iCloud Drive" >&2
      echo "  got: $path" >&2
      echo "  expected prefix: \$HOME/Library/Mobile Documents/com~apple~CloudDocs/" >&2
      echo "  set FLEXAIDDS_RESULTS via: source scripts/use_icloud_benchmark_storage.sh" >&2
      return 91
      ;;
  esac
}

# If executed (not sourced) with an argument, check and exit
if [[ "${BASH_SOURCE[0]:-}" == "${0:-}" ]]; then
  require_icloud_out "${1:-}"
  exit $?
fi
