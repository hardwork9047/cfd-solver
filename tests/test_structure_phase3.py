"""Tests for Issue #3 — src/bin abolition, src/tools creation, configs reorganisation.

Acceptance scenarios:
1. src/bin/ does not exist
2. src/tools/ scripts are importable (no syntax / import errors)
3. configs/lbm_dem/cases/ contains the smoke config
4. configs/lbm_dem/{geometries,materials,templates,sweeps}/ exist
5. All lbm_dem case configs load correctly via SimulationConfig.from_json()
6. run_lbm_dem_sweep.py exists in src/runners/
"""

from __future__ import annotations

import ast
import py_compile
import tempfile
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
_CONFIGS = _REPO_ROOT / "configs"


# ---------------------------------------------------------------------------
# Scenario 1: src/bin/ must not exist
# ---------------------------------------------------------------------------

class TestBinAbolished:
    def test_src_bin_does_not_exist(self):
        assert not (_SRC / "bin").exists(), (
            "src/bin/ still exists — it should have been removed in phase 3"
        )


# ---------------------------------------------------------------------------
# Scenario 2: src/tools/ scripts compile cleanly
# ---------------------------------------------------------------------------

class TestToolsCompile:
    @pytest.fixture
    def tool_scripts(self):
        tools_dir = _SRC / "tools"
        scripts = list(tools_dir.glob("*.py"))
        assert scripts, "src/tools/ is empty or missing"
        return scripts

    def test_tools_dir_exists(self):
        assert (_SRC / "tools").is_dir(), "src/tools/ does not exist"

    @pytest.mark.parametrize("script", list((_SRC / "tools").glob("*.py")) if (_SRC / "tools").is_dir() else [])
    def test_script_compiles(self, script: Path):
        """Each script in src/tools/ must compile without syntax errors."""
        py_compile.compile(str(script), doraise=True)

    def test_tools_has_scripts(self):
        tools_dir = _SRC / "tools"
        scripts = list(tools_dir.glob("*.py"))
        assert len(scripts) >= 5, (
            f"Expected at least 5 scripts in src/tools/, found {len(scripts)}: {[s.name for s in scripts]}"
        )


# ---------------------------------------------------------------------------
# Scenario 3: smoke config is in configs/lbm_dem/cases/
# ---------------------------------------------------------------------------

class TestConfigsCases:
    def test_smoke_config_in_cases(self):
        smoke = _CONFIGS / "lbm_dem" / "cases" / "membrane_pressure_periodic_smoke.json"
        assert smoke.exists(), (
            "membrane_pressure_periodic_smoke.json should be in configs/lbm_dem/cases/"
        )

    def test_cases_dir_exists(self):
        assert (_CONFIGS / "lbm_dem" / "cases").is_dir()

    def test_no_stray_json_at_lbm_dem_root(self):
        """No loose .json files should remain at configs/lbm_dem/ root."""
        stray = [f for f in (_CONFIGS / "lbm_dem").glob("*.json")]
        assert stray == [], f"Stray JSON at lbm_dem root: {[f.name for f in stray]}"


# ---------------------------------------------------------------------------
# Scenario 4: configs/lbm_dem sub-directories exist
# ---------------------------------------------------------------------------

class TestConfigsDirStructure:
    @pytest.mark.parametrize("subdir", ["geometries", "materials", "templates", "sweeps", "cases"])
    def test_subdir_exists(self, subdir: str):
        path = _CONFIGS / "lbm_dem" / subdir
        assert path.is_dir(), f"configs/lbm_dem/{subdir}/ does not exist"


# ---------------------------------------------------------------------------
# Scenario 5: all lbm_dem case configs load via SimulationConfig.from_json()
# ---------------------------------------------------------------------------

class TestConfigLoading:
    @pytest.mark.parametrize(
        "config_path",
        list((_CONFIGS / "lbm_dem" / "cases").glob("*.json"))
        if (_CONFIGS / "lbm_dem" / "cases").is_dir() else [],
    )
    def test_case_loads(self, config_path: Path):
        """Each case config must load without errors (includes extends resolution)."""
        from particulate_flow.io.config import SimulationConfig
        cfg = SimulationConfig.from_json(config_path)
        assert cfg.values, f"Empty values after loading {config_path.name}"


# ---------------------------------------------------------------------------
# Scenario 6: sweep runner exists in src/runners/
# ---------------------------------------------------------------------------

class TestSweepRunner:
    def test_sweep_runner_exists(self):
        runner = _SRC / "runners" / "run_lbm_dem_sweep.py"
        assert runner.exists(), "src/runners/run_lbm_dem_sweep.py does not exist"

    def test_sweep_runner_compiles(self):
        runner = _SRC / "runners" / "run_lbm_dem_sweep.py"
        if runner.exists():
            py_compile.compile(str(runner), doraise=True)

    def test_sweep_runner_points_to_existing_lbm_dem_runner(self):
        """The RUNNER constant in run_lbm_dem_sweep.py must point to an existing file."""
        import importlib.util
        sweep_path = _SRC / "runners" / "run_lbm_dem_sweep.py"
        if not sweep_path.exists():
            pytest.skip("sweep runner not found")
        spec = importlib.util.spec_from_file_location("run_lbm_dem_sweep", sweep_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        assert mod.RUNNER.exists(), f"RUNNER path {mod.RUNNER} does not exist"
