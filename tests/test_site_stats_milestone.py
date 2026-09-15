"""Guard the FlexAIDdS Repository Stats 2,500th-commit milestone chrome.

Live https://thebonhomme.com/FlexAIDdS/ is published from site/FlexAIDdS/
(gh-pages snapshot + pages.yml). The badge threshold and the CSS that keeps
the gold card on the same baseline as C++26 / Apache 2.0 / Source Languages
must not drift back to 1500/2000 or in-flow badge layout (card sitting high).
"""

from __future__ import annotations

import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
JSX = (REPO_ROOT / "site" / "FlexAIDdS" / "sections.jsx").read_text(encoding="utf-8")
CSS = (REPO_ROOT / "site" / "FlexAIDdS" / "styles.css").read_text(encoding="utf-8")
APEX_CSS = (REPO_ROOT / "site" / "style.css").read_text(encoding="utf-8")


class TestProductMilestone(unittest.TestCase):
    def test_threshold_is_2500(self) -> None:
        self.assertIn("commits >= 2500", JSX)
        self.assertIn("2,500th commit", JSX)
        self.assertNotIn("commits >= 1500", JSX)
        self.assertNotIn("commits >= 2000", JSX)
        self.assertNotIn("1,500th", JSX)
        self.assertNotIn("2,000th", JSX)

    def test_stats_row_stretches_items_to_one_baseline(self) -> None:
        row = CSS.split(".stats-row", 1)[1].split(".stat-item {", 1)[0]
        self.assertIn("align-items: stretch", row)
        milestone = CSS.split(".stat-item--milestone", 1)[1].split(
            ".stat-milestone-badge", 1
        )[0]
        self.assertIn("position: relative", milestone)
        self.assertIn("translate(-50%, -50%)", CSS)

    def test_badge_is_out_of_flow(self) -> None:
        badge_block = CSS.split(".stat-milestone-badge", 1)[1].split(".stat-item .stat-value", 1)[0]
        self.assertIn("position: absolute", badge_block)
        self.assertNotIn("margin-bottom", badge_block)
        self.assertNotIn("display: inline-block", badge_block)


class TestApexStatsCss(unittest.TestCase):
    def test_apex_badge_matches_product_overlay(self) -> None:
        self.assertIn("2500-commit milestone", APEX_CSS)
        badge_block = APEX_CSS.split(".stat-milestone-badge", 1)[1].split(".stat-value", 1)[0]
        self.assertIn("position: absolute", badge_block)
        self.assertIn("translate(-50%, -50%)", badge_block)


if __name__ == "__main__":
    unittest.main()
