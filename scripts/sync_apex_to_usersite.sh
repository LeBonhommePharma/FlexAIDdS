#!/usr/bin/env bash
# Sync apex homepage assets to LeBonhommePharma/lebonhommepharma.github.io.
# thebonhomme.com/ is served from the user-site repo (custom domain CNAME),
# not from FlexAIDdS/gh-pages — keep index.html + assets + CNAME in sync after deploys.

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

APEX_FILES=(index.html app.js style.css theme.css theme.js CNAME)

rm -rf "$WORKDIR"
git clone "https://x-access-token:${TOKEN}@github.com/${USER_REPO}.git" "$WORKDIR"

for file in "${APEX_FILES[@]}"; do
  cp "$SITE/$file" "$WORKDIR/$file"
done
mkdir -p "$WORKDIR/assets"
rsync -a "$SITE/assets/" "$WORKDIR/assets/"
mkdir -p "$WORKDIR/FlexAIDdS"
rsync -a --delete "$SITE/FlexAIDdS/" "$WORKDIR/FlexAIDdS/"

cd "$WORKDIR"
git add index.html app.js style.css theme.css theme.js assets/ FlexAIDdS/

if git diff --staged --quiet; then
  echo "user-site apex sync: no changes"
else
  git -c user.name="github-actions[bot]" \
      -c user.email="github-actions[bot]@users.noreply.github.com" \
      commit -m "Sync apex homepage from FlexAIDdS/site (Mol* viewer)"
  git push origin HEAD:main
  echo "user-site apex sync: pushed"
fi

rm -rf "$WORKDIR"