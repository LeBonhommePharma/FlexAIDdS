#!/usr/bin/env bash
# FlexAID∆S first-shot installer (engine and Python package are separate).
#
# Review first:
#   curl -fsSL https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/main/scripts/install.sh | less
# Run:
#   curl -fsSL https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/main/scripts/install.sh | bash
# Flags (when piping, pass them after -s --):
#   curl -fsSL …/scripts/install.sh | bash -s -- --python-only
#   curl -fsSL …/scripts/install.sh | bash -s -- --engine-only
#   curl -fsSL …/scripts/install.sh | bash -s -- --dry-run
#
# Default: Python package always; native engine on macOS when Homebrew is present
# (brew install --HEAD). Does not call unpublished `pip install flexaidds` and
# does not call stable `brew install` without --HEAD.
set -euo pipefail

REPO_HTTPS="https://github.com/LeBonhommePharma/FlexAIDdS.git"
TAP_NAME="lebonhommepharma/flexaidds"
FORMULA="${TAP_NAME}/flexaidds"
PY_GIT="git+${REPO_HTTPS}#subdirectory=python"
VENV="${FLEXAIDDS_VENV:-${HOME}/.flexaidds/venv}"

WANT_ENGINE=1
WANT_PYTHON=1
DRY_RUN=0

usage() {
  cat <<'EOF'
FlexAID∆S first-shot installer

Usage:
  scripts/install.sh [--all | --engine-only | --python-only] [--dry-run]

  --all            Native engine (macOS Homebrew --HEAD) + Python package (default)
  --engine-only    Homebrew FlexAIDdS / tENCoM only
  --python-only    flexaidds Python package only (git+subdirectory, SKIP_CORE=1)
  --dry-run        Print commands; do not install

The native engine and the flexaidds Python package are separate. Installing
one does not give you the other. flexaidds is not on public PyPI. Homebrew
stable v2.0.3 is refused (no provenance JSON); first-shot is --HEAD.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --all) WANT_ENGINE=1; WANT_PYTHON=1; shift ;;
    --engine-only) WANT_ENGINE=1; WANT_PYTHON=0; shift ;;
    --python-only) WANT_ENGINE=0; WANT_PYTHON=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

run() {
  printf '+ %s\n' "$*"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    return 0
  fi
  "$@"
}

have() {
  command -v "$1" >/dev/null 2>&1
}

os="$(uname -s 2>/dev/null || echo unknown)"

echo "FlexAID∆S installer"
echo "  os=${os}"
echo "  engine=${WANT_ENGINE} python=${WANT_PYTHON} dry_run=${DRY_RUN}"
echo "  engine = native FlexAIDdS binary (Homebrew --HEAD on macOS)"
echo "  python = flexaidds analysis package (not the docking engine)"

install_engine() {
  if [[ "$os" != Darwin ]]; then
    echo "Native engine curl path is macOS Homebrew today."
    echo "On ${os}: use CMake or Docker — see docs/INSTALL.md"
    echo "  git clone ${REPO_HTTPS}"
    echo "  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel"
    return 0
  fi
  if ! have brew; then
    echo "Homebrew is required for the macOS engine first-shot." >&2
    echo "Install Homebrew from https://brew.sh then re-run this script." >&2
    echo "Or build from source (docs/INSTALL.md)." >&2
    return 1
  fi
  run brew tap "$TAP_NAME" "$REPO_HTTPS"
  # tap trust is a no-op when already trusted; required on Homebrew 6+ when
  # HOMEBREW_REQUIRE_TAP_TRUST is set.
  run brew trust --formula "$FORMULA" || true
  run brew install --HEAD "$FORMULA"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "+ FlexAIDdS --help"
    return 0
  fi
  if have FlexAIDdS; then
    FlexAIDdS --help >/dev/null
    echo "engine: FlexAIDdS --help ok ($(command -v FlexAIDdS))"
  elif [[ -x /opt/homebrew/bin/FlexAIDdS ]]; then
    /opt/homebrew/bin/FlexAIDdS --help >/dev/null
    echo "engine: /opt/homebrew/bin/FlexAIDdS --help ok (not first on PATH)"
  elif [[ -x /usr/local/bin/FlexAIDdS ]]; then
    /usr/local/bin/FlexAIDdS --help >/dev/null
    echo "engine: /usr/local/bin/FlexAIDdS --help ok (not first on PATH)"
  else
    echo "FlexAIDdS not on PATH after brew install --HEAD" >&2
    return 1
  fi
}

install_python() {
  if ! have python3; then
    echo "python3 >= 3.9 is required for the flexaidds package." >&2
    return 1
  fi
  if ! have git; then
    echo "git is required (pip installs from GitHub, not PyPI)." >&2
    return 1
  fi
  run mkdir -p "$(dirname "$VENV")"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "+ python3 -m venv ${VENV}"
    echo "+ ${VENV}/bin/python -m pip install -U pip"
    echo "+ FLEXAIDDS_SKIP_CORE=1 ${VENV}/bin/pip install ${PY_GIT}"
    echo "+ ${VENV}/bin/python -c 'import flexaidds as fd; print(fd.__version__)'"
    echo "+ ${VENV}/bin/flexaidds --help"
    return 0
  fi
  if [[ ! -x "${VENV}/bin/python" ]]; then
    run python3 -m venv "$VENV"
  fi
  run "${VENV}/bin/python" -m pip install -U pip
  run env FLEXAIDDS_SKIP_CORE=1 "${VENV}/bin/pip" install "$PY_GIT"
  ver="$("${VENV}/bin/python" -c "import flexaidds as fd; print(fd.__version__)")"
  "${VENV}/bin/flexaidds" --help >/dev/null
  echo "python: flexaidds ${ver} ok (${VENV}/bin/flexaidds)"
  echo "Activate with:  source ${VENV}/bin/activate"
}

ec=0
if [[ "$WANT_ENGINE" -eq 1 ]]; then
  install_engine || ec=1
fi
if [[ "$WANT_PYTHON" -eq 1 ]]; then
  install_python || ec=1
fi

if [[ "$ec" -ne 0 ]]; then
  echo "install finished with errors (see above)." >&2
  exit "$ec"
fi
echo "done."
exit 0
