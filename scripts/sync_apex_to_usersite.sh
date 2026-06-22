#!/usr/bin/env bash
# Sync full site/ to LeBonhommePharma/lebonhommepharma.github.io (thebonhomme.com CNAME).
# User-site Pages serves the apex domain; gh-pages alone does not update /.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SITE="$ROOT/site"
USER_REPO="${USER_SITE_REPO:-LeBonhommePharma/lebonhommepharma.github.io}"
WORKDIR="/tmp/usersite-sync"
TOKEN="${GITHUB_TOKEN:-}"

if [ -z "$TOKEN" ]; then
  echo "GITHUB_TOKEN is required to sync $USER_REPO" >&2
  exit 1
fi

rm -rf "$WORKDIR"
git clone "https://x-access-token:${TOKEN}@github.com/${USER_REPO}.git" "$WORKDIR"

rsync -a --delete \
  --exclude '.git' \
  "$SITE/" "$WORKDIR/"

cd "$WORKDIR"
git add -A

if git diff --staged --quiet; then
  echo "user-site sync: no changes"
else
  git -c user.name="github-actions[bot]" \
      -c user.email="github-actions[bot]@users.noreply.github.com" \
      commit -m "Sync full site/ from FlexAIDdS (corporate homepage + product paths)"
  git push origin HEAD:main
  echo "user-site sync: pushed"
fi

rm -rf "$WORKDIR"