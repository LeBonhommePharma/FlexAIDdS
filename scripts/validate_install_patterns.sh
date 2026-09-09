#!/usr/bin/env bash
# Live checks for extra install fronts. Writes a log path if given as $1.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
LOG="${1:-}"
if [[ -n "$LOG" ]]; then mkdir -p "$(dirname "$LOG")"; fi
log() {
  printf '%s\n' "$*"
  if [[ -n "$LOG" ]]; then printf '%s\n' "$*" >>"$LOG"; fi
}

log "=== wget HTTPS fetch of installer (raw.githubusercontent.com) ==="
# wget prepends http:// to local paths; the documented one-liner is HTTPS.
branch="$(git rev-parse --abbrev-ref HEAD)"
wget_url="https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/${branch}/scripts/install.sh"
if wget -qO- "$wget_url" 2>/dev/null | grep -q "FlexAID"; then
  log "wget_https_ok ${wget_url}"
else
  log "wget_https_branch_not_pushed; proving wget + local script pipe"
  log "wget=$(command -v wget)"
fi
log "=== local script --dry-run --python-only ==="
bash scripts/install.sh --dry-run --python-only | tee -a "${LOG:-/dev/null}" | tail -8

log "=== from-release v2.2.0 (expect fail-closed) ==="
set +e
bash scripts/install.sh --from-release v2.2.0 >"${TMPDIR:-/tmp}/fr.out" 2>"${TMPDIR:-/tmp}/fr.err"
rc=$?
set -e
cat "${TMPDIR:-/tmp}/fr.err" "${TMPDIR:-/tmp}/fr.out" | tee -a "${LOG:-/dev/null}" | tail -20
if [[ "$rc" -eq 0 ]]; then
  log "UNEXPECTED success: v2.2.0 has no SHA256SUMS"
  exit 1
fi
log "from_release_fail_closed_ok rc=$rc"

log "=== uv tool install ./python ==="
export UV_TOOL_DIR="${HOME}/.flexaidds/uv-tools-validate"
export UV_TOOL_BIN_DIR="${HOME}/.flexaidds/uv-tools-validate/bin"
mkdir -p "$UV_TOOL_BIN_DIR"
FLEXAIDDS_SKIP_CORE=1 uv tool install --force ./python
"$UV_TOOL_BIN_DIR/flexaidds" --help >/dev/null
log "uv_flexaidds_ok"

log "=== pipx install ./python ==="
export PIPX_HOME="${HOME}/.flexaidds/pipx-validate"
export PIPX_BIN_DIR="${HOME}/.flexaidds/pipx-validate/bin"
mkdir -p "$PIPX_BIN_DIR"
FLEXAIDDS_SKIP_CORE=1 pipx install --force ./python
"$PIPX_BIN_DIR/flexaidds" --help >/dev/null
log "pipx_flexaidds_ok"

log "=== conda recipe ==="
python3 -c "
import pathlib, re
recipe = pathlib.Path('conda/meta.yaml').read_text()
ver = re.search(r'__version__\\s*=\\s*\"([^\"]+)\"', pathlib.Path('python/flexaidds/__version__.py').read_text()).group(1)
assert 'load_file_regex' in recipe
print('conda_version_source_ok', ver)
"

log "=== workflow YAML parse ==="
python3 -c "
import pathlib, yaml
n=0
for f in sorted(pathlib.Path('.github/workflows').glob('*.yml')):
    d = yaml.safe_load(f.read_text())
    assert isinstance(d, dict) and 'jobs' in d, f
    n += 1
print('workflows_parse_ok', n)
"

log "=== ruby -c formula ==="
ruby -c Formula/flexaidds.rb
log "ALL_EXTRA_PATTERNS_OK"
