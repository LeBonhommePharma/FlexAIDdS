"""
Lightweight pytest wrapper for the FlexAIDdS grok-build skill packaging validator.

Run:
    python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
    python3 -m pytest tests/test_flexaid_skill.py -q
"""
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
VALIDATOR = REPO_ROOT / ".grok" / "skills" / "flexaid-docking" / "scripts" / "validate_skill.py"


def test_skill_validator_exists_and_executable():
    """The validator script must exist and be runnable."""
    assert VALIDATOR.exists(), f"Missing validator at {VALIDATOR}"
    # On Unix the +x bit is nice but not required for python3 invocation
    assert VALIDATOR.suffix == ".py"


def test_skill_validator_passes():
    """Running the validator must exit 0 and report VALIDATION PASSED."""
    proc = subprocess.run(
        [sys.executable, str(VALIDATOR)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=30,
    )
    combined = (proc.stdout or "") + (proc.stderr or "")
    assert proc.returncode == 0, (
        f"Validator failed with code {proc.returncode}\n"
        f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    assert "VALIDATION PASSED" in combined, (
        "Validator did not emit success message.\n" + combined
    )


def test_skill_metadata_and_aliases_present():
    """Quick smoke that SKILL.md still carries the required frontmatter + aliases."""
    skill_md = REPO_ROOT / ".grok" / "skills" / "flexaid-docking" / "SKILL.md"
    assert skill_md.exists()
    text = skill_md.read_text(encoding="utf-8")
    assert "name: flexaid-docking" in text
    assert "description:" in text
    for alias in ("/FlexAid docking", "/FlexAidDS", "FlexAIDdS", "FlexAID∆S"):
        assert alias in text, f"Missing documented alias {alias}"


if __name__ == "__main__":
    # Allow direct execution: python tests/test_flexaid_skill.py
    pytest.main([__file__, "-q", "--tb=line"])
