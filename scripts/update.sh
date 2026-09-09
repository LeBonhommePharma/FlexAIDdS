#!/usr/bin/env bash
# Refresh every FlexAID∆S front that is already on this machine.
#
#   curl -fsSL https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/main/scripts/update.sh | bash
#   wget -qO-  https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/main/scripts/update.sh | bash
#   bash scripts/update.sh --dry-run
#   python -m flexaidds --self-update
#
# Detects Homebrew HEAD engine, ~/.flexaidds/venv, uv tool, pipx, and a local
# GHCR image. Updates all of them in one invocation. Does not invent a PyPI
# upgrade while the package is unpublished.
set -euo pipefail

REPO_HTTPS="https://github.com/LeBonhommePharma/FlexAIDdS.git"
TAP_NAME="lebonhommepharma/flexaidds"
FORMULA="${TAP_NAME}/flexaidds"
PY_GIT="git+${REPO_HTTPS}#subdirectory=python"
VENV="${FLEXAIDDS_VENV:-${HOME}/.flexaidds/venv}"
GHCR="ghcr.io/lebonhommepharma/flexaidds"

DRY_RUN=0
DO_ENGINE=1
DO_PYTHON=1
DO_DOCKER=1

usage() {
  cat <<'EOF'
FlexAID∆S unified updater

Usage: scripts/update.sh [--dry-run] [--engine-only] [--python-only] [--docker-only]

Refreshes every detected install front:
  - macOS Homebrew engine: brew upgrade --fetch-HEAD (fallback: uninstall + --HEAD)
  - Python: venv pip, uv tool, and/or pipx (git+subdirectory, SKIP_CORE=1)
  - Docker: docker pull ghcr.io/lebonhommepharma/flexaidds:latest if that image exists locally

Does not call unpublished `pip install --upgrade flexaidds` as the only path.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=1; shift ;;
    --engine-only) DO_PYTHON=0; DO_DOCKER=0; shift ;;
    --python-only) DO_ENGINE=0; DO_DOCKER=0; shift ;;
    --docker-only) DO_ENGINE=0; DO_PYTHON=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

run() {
  printf '+ %s\n' "$*"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    return 0
  fi
  "$@"
}

have() { command -v "$1" >/dev/null 2>&1; }

echo "FlexAID∆S updater  dry_run=${DRY_RUN}"

ec=0

update_engine() {
  [[ "$DO_ENGINE" -eq 1 ]] || return 0
  if ! have brew; then
    echo "engine: brew not found — skip"
    return 0
  fi
  if ! brew list --formula lebonhommepharma/flexaidds/flexaidds >/dev/null 2>&1 \
     && ! brew list --formula flexaidds >/dev/null 2>&1; then
    echo "engine: flexaidds formula not installed — skip (install: scripts/install.sh --engine-only)"
    return 0
  fi
  run brew tap "$TAP_NAME" "$REPO_HTTPS" || true
  run brew update || true
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "+ brew upgrade --fetch-HEAD ${FORMULA}"
    echo "+ FlexAIDdS --help"
    return 0
  fi
  if ! brew upgrade --fetch-HEAD "$FORMULA"; then
    echo "engine: upgrade --fetch-HEAD failed; reinstall --HEAD"
    run brew uninstall --force "$FORMULA" || true
    run brew install --HEAD "$FORMULA"
  fi
  if have FlexAIDdS; then
    FlexAIDdS --help >/dev/null
  elif [[ -x /opt/homebrew/bin/FlexAIDdS ]]; then
    /opt/homebrew/bin/FlexAIDdS --help >/dev/null
  fi
  echo "engine: updated"
}

update_python() {
  [[ "$DO_PYTHON" -eq 1 ]] || return 0
  local did=0
  if [[ -x "${VENV}/bin/pip" ]]; then
    run env FLEXAIDDS_SKIP_CORE=1 "${VENV}/bin/pip" install --upgrade --force-reinstall "$PY_GIT"
    did=1
  fi
  if have uv && uv tool list 2>/dev/null | grep -q '^flexaidds '; then
    run env FLEXAIDDS_SKIP_CORE=1 uv tool upgrade flexaidds || \
      run env FLEXAIDDS_SKIP_CORE=1 uv tool install --force "$PY_GIT"
    did=1
  fi
  if have pipx && pipx list 2>/dev/null | grep -qi flexaidds; then
    run env FLEXAIDDS_SKIP_CORE=1 pipx upgrade flexaidds || \
      run env FLEXAIDDS_SKIP_CORE=1 pipx install --force "$PY_GIT"
    did=1
  fi
  if [[ "$did" -eq 0 ]]; then
    echo "python: no venv/uv/pipx flexaidds detected — skip"
    echo "        install: scripts/install.sh --python-only"
    return 0
  fi
  echo "python: updated"
}

update_docker() {
  [[ "$DO_DOCKER" -eq 1 ]] || return 0
  if ! have docker; then
    echo "docker: docker not on PATH — skip"
    return 0
  fi
  if ! docker info >/dev/null 2>&1; then
    echo "docker: daemon not running — skip"
    return 0
  fi
  if ! docker image inspect "${GHCR}:latest" >/dev/null 2>&1; then
    echo "docker: ${GHCR}:latest not present locally — skip (no unpublished pull)"
    return 0
  fi
  run docker pull "${GHCR}:latest"
  echo "docker: pulled ${GHCR}:latest"
}

update_engine || ec=1
update_python || ec=1
update_docker || ec=1

if [[ "$ec" -ne 0 ]]; then
  echo "update finished with errors." >&2
  exit "$ec"
fi
echo "done."
exit 0
