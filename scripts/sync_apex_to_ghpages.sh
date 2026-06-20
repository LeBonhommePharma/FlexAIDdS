#!/usr/bin/env bash
# Publish the full site/ tree to the gh-pages branch.
# GitHub Pages project workflow artifacts only update /FlexAIDdS/ on the custom
# domain, leaving stale orphan files at thebonhomme.com/. Serving from gh-pages
# keeps apex (/) and /FlexAIDdS/ paths in sync.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SITE="$ROOT/site"
WORKTREE="/tmp/gh-pages-sync"

git -C "$ROOT" fetch origin gh-pages

if [ -d "$WORKTREE" ]; then
  git -C "$ROOT" worktree remove --force "$WORKTREE" 2>/dev/null || rm -rf "$WORKTREE"
fi

git -C "$ROOT" worktree add "$WORKTREE" origin/gh-pages

# Replace published tree with the current site/ contents.
find "$WORKTREE" -mindepth 1 -maxdepth 1 ! -name '.git' -exec rm -rf {} +
rsync -a --delete "$SITE/" "$WORKTREE/"

cd "$WORKTREE"
git add -A

if git diff --staged --quiet; then
  echo "gh-pages publish: no changes"
else
  git -c user.name="github-actions[bot]" \
      -c user.email="github-actions[bot]@users.noreply.github.com" \
      commit -m "Publish site/ to gh-pages (Mol* viewer + apex sync)"
  git push origin HEAD:gh-pages
  echo "gh-pages publish: pushed"
fi

git -C "$ROOT" worktree remove --force "$WORKTREE"