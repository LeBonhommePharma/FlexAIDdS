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
# --checksum: git worktree add stamps dest mtimes as "now", which is newer
# than the just-patched source. rsync -a then skips equal-size files, so a
# 4-digit commit count change (2553→2560) never reached Pages.
rsync -a --delete --checksum \
  --exclude '.git' \
  "$SITE/FlexAIDdS/" "$WORKTREE/"

if ! cmp -s "$SITE/FlexAIDdS/index.html" "$WORKTREE/index.html"; then
  echo "gh-pages publish: site/FlexAIDdS/index.html was not copied" >&2
  exit 1
fi

# Never claim the apex custom domain from this repo.
rm -f "$WORKTREE/CNAME"
# Branch-deploy fallback: GitHub's Jekyll builder must not mangle the product tree.
touch "$WORKTREE/.nojekyll"
# Install the Pages deploy workflow onto this branch so it can be dispatched
# with --ref gh-pages (environment policy allows gh-pages + master, not main).
mkdir -p "$WORKTREE/.github/workflows"
cp "$ROOT/.github/workflows/pages.yml" "$WORKTREE/.github/workflows/pages.yml"

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