"""Tests for issue-16: /simulate skill 3D support.

The /simulate skill is interpreted by the agent at runtime, so these are
documentation-presence guards over `.claude/skills/simulate/SKILL.md`: they ensure
the instructions for 3D config generation, the confirmation table, and the 2D-only
warning remain present, and that the 2D path is not removed (scenario 4).
"""

from __future__ import annotations

from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SKILL = _REPO_ROOT / ".claude" / "skills" / "simulate" / "SKILL.md"


@pytest.fixture(scope="module")
def skill_text() -> str:
    return _SKILL.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Scenario 1 + 2: 3D config generation + confirmation table
# ---------------------------------------------------------------------------


class TestThreeDInstructions:
    def test_references_3d_template(self, skill_text):
        assert "fouling_supply_3d.json" in skill_text

    def test_mentions_dimensions_and_nz(self, skill_text):
        assert "dimensions" in skill_text
        assert "nz" in skill_text

    def test_has_3d_example_block(self, skill_text):
        # A generation example that sets dimensions: 3.
        assert '"dimensions": 3' in skill_text

    def test_references_3d_geometry(self, skill_text):
        assert "four_cylinder_3d.json" in skill_text


# ---------------------------------------------------------------------------
# Scenario 3: warning for 2D-only templates/geometries under 3D intent
# ---------------------------------------------------------------------------


class TestTwoDOnlyWarning:
    def test_warns_about_2d_only_assets(self, skill_text):
        # Some warning tying the 2D-only geometry/template to 3D misuse.
        assert "2D" in skill_text or "2 次元" in skill_text or "二次元" in skill_text
        # The 2D staggered geometry is explicitly flagged as 2D-only somewhere
        # in the 3D guidance.
        assert "four_cylinder_staggered.json" in skill_text


# ---------------------------------------------------------------------------
# Scenario 4: 2D path unchanged (still documented)
# ---------------------------------------------------------------------------


class TestTwoDStillPresent:
    def test_2d_base_template_still_referenced(self, skill_text):
        assert "../templates/fouling_supply.json" in skill_text

    def test_2d_default_behaviour_documented(self, skill_text):
        # The skill still describes the default (2D) flow.
        assert "fouling_supply" in skill_text
