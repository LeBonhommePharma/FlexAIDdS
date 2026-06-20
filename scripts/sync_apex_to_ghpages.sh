#!/usr/bin/env bash
# Sync apex homepage assets from site/ to the gh-pages branch root.
# GitHub Pages project sites serve workflow artifacts under /FlexAIDdS/ while
# thebonhomme.com/ is still backed by the legacy gh-pages branch — keep both in sync.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SITE="$ROOT/site"
WORKTREE="/tmp/gh-pages-sync"

APEX_FILES=(
  index.html
  app.js
  style.css
  theme.css
  theme.js
)

STALE_ROOT_FILES=(
  app.jsx
  components.jsx
  sections.jsx
  styles.css
  colors_and_type.css
)

git -C "$ROOT" fetch origin gh-pages

if [ -d "$WORKTREE" ]; then
  git -C "$ROOT" worktree remove --force "$WORKTREE" 2>/dev/null || rm -rf "$WORKTREE"
fi

git -C "$ROOT" worktree add "$WORKTREE" origin/gh-pages

for file in "${STALE_ROOT_FILES[@]}"; do
  rm -f "$WORKTREE/$file"
done

for file in "${APEX_FILES[@]}"; do
  cp "$SITE/$file" "$WORKTREE/$file"
done

mkdir -p "$WORKTREE/assets"
rsync -a --delete "$SITE/assets/" "$WORKTREE/assets/"

cd "$WORKTREE"
git add -A

if git diff --staged --quiet; then
  echo "gh-pages apex sync: no changes"
else
  git -c user.name="github-actions[bot]" \
      -c user.email="github-actions[bot]@users.noreply.github.com" \
      commit -m "Sync apex homepage from site/ (Mol* viewer)"
  git push origin HEAD:gh-pages
  echo "gh-pages apex sync: pushed"
fi

git -C "$ROOT" worktree remove --force "$WORKTREE"