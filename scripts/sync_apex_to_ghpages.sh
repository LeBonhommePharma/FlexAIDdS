#!/usr/bin/env bash
# Publish site/ to the gh-pages branch (backup / workflow artifact).
# Do NOT publish CNAME here — only lebonhommepharma.github.io may claim
# thebonhomme.com. A duplicate CNAME on gh-pages hijacks /FlexAIDdS/ paths.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SITE="$ROOT/site"
WORKTREE="/tmp/gh-pages-sync"

git -C "$ROOT" fetch origin gh-pages

if [ -d "$WORKTREE" ]; then
  git -C "$ROOT" worktree remove --force "$WORKTREE" 2>/dev/null || rm -rf "$WORKTREE"
fi

git -C "$ROOT" worktree add "$WORKTREE" origin/gh-pages

# Replace published tree with the current site/ contents (keep worktree git metadata).
rsync -a --delete --exclude '.git' --exclude 'CNAME' "$SITE/" "$WORKTREE/"
# Ensure gh-pages never steals the custom domain from the user-site repo.
rm -f "$WORKTREE/CNAME"

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