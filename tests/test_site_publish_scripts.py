"""Guards for the FlexAIDdS → Pages / apex stats publish path.

These are pure file/regex checks. They exist because a green Deploy Site
workflow previously left https://thebonhomme.com/FlexAIDdS/ frozen: it
pushed gh-pages (branch snapshot) but Pages build_type is "workflow",
deploy-pages on main is rejected by the github-pages environment (allowed
branches: gh-pages + master), and the usersite step swallowed a 403 after
an rsync --delete of the whole apex tree.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import update_site_stats as stats  # noqa: E402
import verify_live_site_stats as verify  # noqa: E402

WORKFLOW = (REPO_ROOT / ".github/workflows/update-site.yml").read_text(encoding="utf-8")
PAGES = (REPO_ROOT / ".github/workflows/pages.yml").read_text(encoding="utf-8")
USERSITE = (REPO_ROOT / "scripts/sync_apex_to_usersite.sh").read_text(encoding="utf-8")
GHPAGES = (REPO_ROOT / "scripts/sync_apex_to_ghpages.sh").read_text(encoding="utf-8")


class TestWorkflowDeploysPages(unittest.TestCase):
    def test_pages_workflow_calls_deploy_pages(self) -> None:
        self.assertIn("actions/deploy-pages@", PAGES)
        self.assertIn("actions/upload-pages-artifact@", PAGES)
        self.assertIn("github-pages-deploy", PAGES)

    def test_pages_workflow_runs_on_allowed_refs(self) -> None:
        self.assertIn("branches: [gh-pages]", PAGES)
        self.assertIn("workflow_dispatch:", PAGES)
        self.assertIn("name: github-pages", PAGES)

    def test_update_site_does_not_attach_github_pages_environment(self) -> None:
        # environment: github-pages rejects branch main (policy is gh-pages + master).
        self.assertNotIn("name: github-pages", WORKFLOW)
        self.assertNotIn("actions/deploy-pages@", WORKFLOW)

    def test_update_site_dispatches_pages_on_gh_pages_ref(self) -> None:
        self.assertIn("gh workflow run pages.yml --ref gh-pages", WORKFLOW)
        self.assertIn("actions: write", WORKFLOW)

    def test_does_not_swallow_usersite_failure_after_destructive_rsync(self) -> None:
        self.assertNotIn("user-site sync failed", WORKFLOW)

    def test_verifies_live_markers(self) -> None:
        self.assertIn("verify_live_site_stats.py", WORKFLOW)
        self.assertIn("thebonhomme.com/FlexAIDdS/", WORKFLOW)
        self.assertIn("lebonhommepharma.github.io/FlexAIDdS/", WORKFLOW)


class TestUsersiteScriptIsSurgical(unittest.TestCase):
    def test_does_not_rsync_delete_site_onto_usersite(self) -> None:
        self.assertNotIn('"$SITE/" "$WORKDIR/"', USERSITE)
        self.assertNotIn("rsync -a --delete", USERSITE)

    def test_requires_explicit_usersite_token(self) -> None:
        self.assertIn("USER_SITE_TOKEN", USERSITE)
        self.assertNotIn('TOKEN="${GITHUB_TOKEN:-}"', USERSITE)

    def test_patches_via_update_site_stats(self) -> None:
        self.assertIn("--patch-tree", USERSITE)
        self.assertIn("--from-json", USERSITE)


class TestGhpagesPublish(unittest.TestCase):
    def test_still_publishes_product_root(self) -> None:
        self.assertIn('"$SITE/FlexAIDdS/" "$WORKTREE/"', GHPAGES)
        self.assertIn('rm -f "$WORKTREE/CNAME"', GHPAGES)

    def test_writes_nojekyll(self) -> None:
        self.assertIn(".nojekyll", GHPAGES)

    def test_installs_pages_workflow_on_gh_pages(self) -> None:
        self.assertIn(".github/workflows/pages.yml", GHPAGES)
        self.assertIn("pages.yml", GHPAGES)


class TestPatchHtmlMarkers(unittest.TestCase):
    SAMPLE = (
        '<span id="stat-commits" hidden>1744</span>\n'
        '<span id="stat-langs" hidden>15</span>\n'
        '<span id="stat-stars" hidden>9</span>\n'
        '<span id="last-updated" hidden>2026-07-14</span>\n'
        '<span id="latest-release" hidden>v2.0.2</span>\n'
    )

    def test_replaces_every_marker(self) -> None:
        out = stats.patch_html_markers(
            self.SAMPLE,
            2553,
            15,
            stars=11,
            release="v2.2.0",
            last_updated="2026-09-08",
        )
        self.assertIn(">2553<", out)
        self.assertIn(">15<", out)
        self.assertIn(">11<", out)
        self.assertIn(">2026-09-08<", out)
        self.assertIn(">v2.2.0<", out)
        self.assertNotIn("1744", out)
        self.assertNotIn("v2.0.2", out)

    def test_apply_snapshot_to_tree_is_surgical(self) -> None:
        snapshot = {
            "commits": 2553,
            "languageCount": 15,
            "lastUpdated": "2026-09-08",
            "stars": 11,
            "latestRelease": "v2.2.0",
            "languages": [
                {"id": "cpp", "name": "C++", "percent": 50.7, "color": "#f34b7d"},
            ],
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "FlexAIDdS").mkdir()
            (root / "flexaid-ds").mkdir()
            (root / "assets").mkdir()
            (root / "rive").mkdir()
            (root / "FlexAIDdS" / "index.html").write_text(self.SAMPLE, encoding="utf-8")
            (root / "flexaid-ds" / "index.html").write_text(self.SAMPLE, encoding="utf-8")
            (root / "rive" / "keep-me.txt").write_text("atlas\n", encoding="utf-8")
            changed = stats.apply_snapshot_to_tree(str(root), snapshot)
            self.assertTrue(any("FlexAIDdS/index.html" in path for path in changed))
            self.assertTrue((root / "rive" / "keep-me.txt").is_file())
            payload = json.loads((root / "assets" / "repo-stats.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["commits"], 2553)
            self.assertEqual(payload["stars"], 11)
            self.assertEqual(payload["latestRelease"], "v2.2.0")
            self.assertEqual(payload["lastUpdated"], "2026-09-08")
            html = (root / "FlexAIDdS" / "index.html").read_text(encoding="utf-8")
            self.assertIn(">2553<", html)
            self.assertIn(">v2.2.0<", html)

    def test_verify_parse_markers(self) -> None:
        found = verify.parse_markers(self.SAMPLE)
        self.assertEqual(found["stat-commits"], "1744")
        self.assertEqual(found["stat-stars"], "9")
        self.assertEqual(found["latest-release"], "v2.0.2")
        self.assertEqual(found["last-updated"], "2026-07-14")


if __name__ == "__main__":
    unittest.main()
