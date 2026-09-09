#!/usr/bin/env bash
# FlexAID∆S first-shot installer (engine and Python package are separate).
#
# Review:
#   curl -fsSL https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/main/scripts/install.sh | less
#   wget -qO- https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/main/scripts/install.sh | less
# Run:
#   curl -fsSL …/scripts/install.sh | bash
#   wget -qO-  …/scripts/install.sh | bash
# Flags (after -s -- when piping):
#   --python-only | --engine-only | --uv | --pipx | --from-release [tag]
#   --from-archive FILE --sha256 HEX --prefix DIR | --dry-run
set -euo pipefail

REPO_HTTPS="https://github.com/LeBonhommePharma/FlexAIDdS.git"
REPO_SLUG="LeBonhommePharma/FlexAIDdS"
TAP_NAME="lebonhommepharma/flexaidds"
FORMULA="${TAP_NAME}/flexaidds"
PY_GIT="git+${REPO_HTTPS}#subdirectory=python"
VENV="${FLEXAIDDS_VENV:-${HOME}/.flexaidds/venv}"
PREFIX="${FLEXAIDDS_PREFIX:-${HOME}/.local}"

WANT_ENGINE=1
WANT_PYTHON=1
DRY_RUN=0
PY_BACKEND="pip" # pip | uv | pipx
FROM_RELEASE=""
FROM_ARCHIVE=""
EXPECT_SHA256=""

usage() {
  cat <<'EOF'
FlexAID∆S first-shot installer

Usage:
  scripts/install.sh [options]

  --all              Homebrew --HEAD engine (macOS) + Python (default)
  --engine-only      Native engine only
  --python-only      flexaidds Python package only
  --uv               Python via `uv tool install` (git+subdirectory)
  --pipx             Python via `pipx install` (git+subdirectory)
  --from-release [TAG]
                     Download GitHub Release tarball + SHA256SUMS.txt (fail-closed)
  --from-archive FILE
                     Install a local tarball (pair with --sha256)
  --sha256 HEX       Expected SHA-256 of --from-archive
  --prefix DIR       Unpack prefix for release/archive (default: ~/.local)
  --dry-run          Print commands; do not install

Does not call unpublished `pip install flexaidds`.
Does not call stable Homebrew without --HEAD.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --all) WANT_ENGINE=1; WANT_PYTHON=1; shift ;;
    --engine-only) WANT_ENGINE=1; WANT_PYTHON=0; shift ;;
    --python-only) WANT_ENGINE=0; WANT_PYTHON=1; shift ;;
    --uv) PY_BACKEND="uv"; shift ;;
    --pipx) PY_BACKEND="pipx"; shift ;;
    --from-release)
      WANT_ENGINE=1
      WANT_PYTHON=0
      if [[ $# -ge 2 && "$2" != --* ]]; then FROM_RELEASE="$2"; shift 2; else FROM_RELEASE="latest"; shift; fi
      ;;
    --from-archive)
      WANT_ENGINE=1
      WANT_PYTHON=0
      FROM_ARCHIVE="${2:?--from-archive needs a file}"
      shift 2
      ;;
    --sha256) EXPECT_SHA256="${2:?}"; shift 2 ;;
    --prefix) PREFIX="${2:?}"; shift 2 ;;
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

have() { command -v "$1" >/dev/null 2>&1; }

file_sha256() {
  if have shasum; then
    shasum -a 256 "$1" | awk '{print $1}'
  else
    sha256sum "$1" | awk '{print $1}'
  fi
}

os="$(uname -s 2>/dev/null || echo unknown)"
arch="$(uname -m 2>/dev/null || echo unknown)"

echo "FlexAID∆S installer"
echo "  os=${os} arch=${arch} backend=${PY_BACKEND} dry_run=${DRY_RUN}"
echo "  engine=${WANT_ENGINE} python=${WANT_PYTHON}"

install_engine_brew() {
  if [[ "$os" != Darwin ]]; then
    echo "Native Homebrew engine is macOS-only. On ${os} use CMake, Docker, or --from-release."
    echo "  git clone ${REPO_HTTPS}"
    echo "  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel"
    return 0
  fi
  if ! have brew; then
    echo "Homebrew is required for the macOS engine first-shot." >&2
    echo "Install Homebrew from https://brew.sh then re-run." >&2
    return 1
  fi
  run brew tap "$TAP_NAME" "$REPO_HTTPS"
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
    echo "engine: /opt/homebrew/bin/FlexAIDdS --help ok"
  elif [[ -x /usr/local/bin/FlexAIDdS ]]; then
    /usr/local/bin/FlexAIDdS --help >/dev/null
    echo "engine: /usr/local/bin/FlexAIDdS --help ok"
  else
    echo "FlexAIDdS not on PATH after brew install --HEAD" >&2
    return 1
  fi
}

release_asset_name() {
  case "${os}-${arch}" in
    Darwin-arm64|Darwin-aarch64) echo "flexaidds-macos.tar.gz" ;;
    Darwin-x86_64) echo "flexaidds-macos.tar.gz" ;;
    Linux-x86_64|Linux-amd64) echo "flexaidds-linux.tar.gz" ;;
    MINGW*|MSYS*|CYGWIN*|Windows*) echo "flexaidds-windows.zip" ;;
    *) echo "flexaidds-macos.tar.gz" ;;
  esac
}

install_from_archive() {
  local archive="$1"
  local expect="$2"
  if [[ ! -f "$archive" && "$DRY_RUN" -eq 0 ]]; then
    echo "archive not found: $archive" >&2
    return 1
  fi
  if [[ -z "$expect" ]]; then
    echo "--from-archive/--from-release requires a SHA-256 (SHA256SUMS.txt or --sha256)." >&2
    return 1
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "+ sha256 ${archive} == ${expect}"
    echo "+ tar/unzip ${archive} -> ${PREFIX}"
    return 0
  fi
  local got
  got="$(file_sha256 "$archive")"
  if [[ "$got" != "$expect" ]]; then
    echo "SHA-256 mismatch for ${archive}" >&2
    echo "  expected ${expect}" >&2
    echo "  got      ${got}" >&2
    return 1
  fi
  mkdir -p "${PREFIX}/bin"
  case "$archive" in
    *.zip)
      unzip -o "$archive" -d "${PREFIX}/.flexaidds-unpack"
      ;;
    *)
      tar -xzf "$archive" -C "${PREFIX}" --strip-components=0
      ;;
  esac
  # Release layout is dist/FlexAIDdS (see .github/workflows/release.yml).
  local candidate
  for candidate in \
      "${PREFIX}/dist/FlexAIDdS" \
      "${PREFIX}/FlexAIDdS" \
      "${PREFIX}/bin/FlexAIDdS" \
      "${PREFIX}/.flexaidds-unpack/dist/FlexAIDdS" \
      "${PREFIX}/.flexaidds-unpack/FlexAIDdS"
  do
    if [[ -f "$candidate" ]]; then
      chmod +x "$candidate"
      ln -sf "$candidate" "${PREFIX}/bin/FlexAIDdS"
      echo "engine: installed ${candidate} -> ${PREFIX}/bin/FlexAIDdS"
      "${PREFIX}/bin/FlexAIDdS" --help >/dev/null || true
      return 0
    fi
  done
  echo "tarball hashed OK but FlexAIDdS binary not found under ${PREFIX}" >&2
  return 1
}

install_from_release() {
  local tag="$1"
  local asset
  asset="$(release_asset_name)"
  local base="https://github.com/${REPO_SLUG}/releases"
  local sums_url asset_url
  if [[ "$tag" == "latest" ]]; then
    sums_url="${base}/latest/download/SHA256SUMS.txt"
    asset_url="${base}/latest/download/${asset}"
  else
    sums_url="${base}/download/${tag}/SHA256SUMS.txt"
    asset_url="${base}/download/${tag}/${asset}"
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "+ curl -fsSL ${sums_url}"
    echo "+ curl -fsSL ${asset_url}"
    echo "+ verify sha256 of ${asset}"
    return 0
  fi
  local tmp
  tmp="$(mktemp -d "${TMPDIR:-/tmp}/flexaidds-rel.XXXXXX")"
  if ! curl -fsSL "$sums_url" -o "${tmp}/SHA256SUMS.txt"; then
    echo "No SHA256SUMS.txt on GitHub Release ${tag}." >&2
    echo "Current tags may have zero assets; this path is live after the next tagged release that runs .github/workflows/release.yml." >&2
    rm -rf "$tmp"
    return 1
  fi
  curl -fsSL "$asset_url" -o "${tmp}/${asset}"
  local expect
  expect="$(awk -v a="${asset}" '$2==a || $2=="*"a {print $1; exit}' "${tmp}/SHA256SUMS.txt")"
  if [[ -z "$expect" ]]; then
    echo "${asset} not listed in SHA256SUMS.txt" >&2
    rm -rf "$tmp"
    return 1
  fi
  install_from_archive "${tmp}/${asset}" "$expect"
  local rc=$?
  rm -rf "$tmp"
  return "$rc"
}

install_python_pip() {
  if ! have python3; then echo "python3 >= 3.9 required" >&2; return 1; fi
  if ! have git; then echo "git required (not on PyPI)" >&2; return 1; fi
  run mkdir -p "$(dirname "$VENV")"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "+ python3 -m venv ${VENV}"
    echo "+ ${VENV}/bin/python -m pip install -U pip"
    echo "+ FLEXAIDDS_SKIP_CORE=1 ${VENV}/bin/pip install ${PY_GIT}"
    echo "+ ${VENV}/bin/flexaidds --help"
    return 0
  fi
  [[ -x "${VENV}/bin/python" ]] || run python3 -m venv "$VENV"
  run "${VENV}/bin/python" -m pip install -U pip
  run env FLEXAIDDS_SKIP_CORE=1 "${VENV}/bin/pip" install "$PY_GIT"
  ver="$("${VENV}/bin/python" -c "import flexaidds as fd; print(fd.__version__)")"
  "${VENV}/bin/flexaidds" --help >/dev/null
  echo "python: flexaidds ${ver} via pip (${VENV}/bin/flexaidds)"
  echo "Activate:  source ${VENV}/bin/activate"
}

install_python_uv() {
  if ! have uv; then
    echo "uv not on PATH. Install from https://docs.astral.sh/uv/ then re-run --uv" >&2
    return 1
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "+ FLEXAIDDS_SKIP_CORE=1 uv tool install ${PY_GIT}"
    echo "+ uvx --from ${PY_GIT} flexaidds --help"
    return 0
  fi
  run env FLEXAIDDS_SKIP_CORE=1 uv tool install --force "$PY_GIT"
  uv tool run --from "$PY_GIT" flexaidds --help >/dev/null || \
    flexaidds --help >/dev/null
  echo "python: flexaidds via uv tool"
}

install_python_pipx() {
  if ! have pipx; then
    echo "pipx not on PATH. pip install pipx  OR  brew install pipx" >&2
    return 1
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "+ FLEXAIDDS_SKIP_CORE=1 pipx install ${PY_GIT}"
    echo "+ flexaidds --help"
    return 0
  fi
  run env FLEXAIDDS_SKIP_CORE=1 pipx install --force "$PY_GIT"
  flexaidds --help >/dev/null
  echo "python: flexaidds via pipx"
}

install_python() {
  case "$PY_BACKEND" in
    uv) install_python_uv ;;
    pipx) install_python_pipx ;;
    *) install_python_pip ;;
  esac
}

ec=0
if [[ -n "$FROM_ARCHIVE" ]]; then
  install_from_archive "$FROM_ARCHIVE" "$EXPECT_SHA256" || ec=1
elif [[ -n "$FROM_RELEASE" ]]; then
  install_from_release "$FROM_RELEASE" || ec=1
elif [[ "$WANT_ENGINE" -eq 1 ]]; then
  install_engine_brew || ec=1
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
