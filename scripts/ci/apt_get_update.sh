#!/usr/bin/env bash
# GitHub-hosted Ubuntu runners occasionally ship broken Microsoft apt repos
# (NOSPLIT / invalid clearsign). Disable them before apt-get update.
set -euo pipefail

if [[ "$(uname -s)" != "Linux" ]]; then
  exit 0
fi

shopt -s nullglob
for repo in /etc/apt/sources.list.d/*microsoft* /etc/apt/sources.list.d/*azure* /etc/apt/sources.list.d/*msedge*; do
  sudo mv "$repo" "${repo}.disabled" 2>/dev/null || sudo rm -f "$repo"
done

sudo apt-get update -qq -o Acquire::Retries=3