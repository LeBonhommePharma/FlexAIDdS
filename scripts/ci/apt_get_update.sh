#!/usr/bin/env bash
# GitHub-hosted Ubuntu runners occasionally ship broken Microsoft apt repos
# (NOSPLIT / invalid clearsign). Disable them before apt-get update.
set -euo pipefail

if [[ "$(uname -s)" != "Linux" ]]; then
  exit 0
fi

for repo in \
  /etc/apt/sources.list.d/microsoft-prod.list \
  /etc/apt/sources.list.d/azure-cli.list \
  /etc/apt/sources.list.d/msedge.list; do
  if [[ -f "$repo" ]]; then
    sudo mv "$repo" "${repo}.disabled" || sudo rm -f "$repo"
  fi
done

sudo apt-get update -qq