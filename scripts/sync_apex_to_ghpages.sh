#!/usr/bin/env bash
# Publish FlexAID∆S product files to gh-pages ROOT.
# GitHub mounts this repo's gh-pages at thebonhomme.com/FlexAIDdS/ — root index.html
# must be the product page, NOT the corporate homepage from site/index.html.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SITE="$ROOT/site"
WORKTREE="/tmp/gh-pages-sync"

git -C "$ROOT" fetch origin gh-pages

if [ -d "$WORKTREE" ]; then
  git -C "$ROOT" worktree remove --force "$WORKTREE" 2>/dev/null || rm -rf "$WORKTREE"
fi

git -C "$ROOT" worktree add "$WORKTREE" origin/gh-pages

# Product site only (served at /FlexAIDdS/ on the apex domain).
rsync -a --delete \
  --exclude '.git' \
  "$SITE/FlexAIDdS/" "$WORKTREE/"

# Never claim the apex custom domain from this repo.
rm -f "$WORKTREE/CNAME"

cd "$WORKTREE"
git add -A

if git diff --staged --quiet; then
  echo "gh-pages publish: no changes"
else
  git -c user.name="github-actions[bot]" \
      -c user.email="github-actions[bot]@users.noreply.github.com" \
      commit -m "Publish FlexAID∆S product to gh-pages root (/FlexAIDdS/ mount)"
  git push origin HEAD:gh-pages
  echo "gh-pages publish: pushed"
fi

git -C "$ROOT" worktree remove --force "$WORKTREE"