"""Tests for Issue #2 — runner entry-point refactor.

Acceptance scenarios:
1. build_lbm_dem_solver() exists and returns a FastLBMDEM
2. solver section in SimulationConfig maps fluid_method etc.
3. accelerator section maps fluid→fluid_accelerator, compute→compute_accelerator
4. build_dem_packing_solver() exists and returns a DEMPackingSimulation
5. run_lbm_dem.py contains no inline FastLBMDEM construction
6. build_lbm_dem_solver uses subpackage imports
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Acceptance Scenario 1 & 6: builder module and subpackage imports
# ---------------------------------------------------------------------------

class TestBuilderImport:
    def test_builder_module_importable(self):
        from particulate_flow import builder
        assert builder is not None

    def test_build_lbm_dem_solver_exists(self):
        from particulate_flow.builder import build_lbm_dem_solver
        assert callable(build_lbm_dem_solver)

    def test_build_dem_packing_solver_exists(self):
        from particulate_flow.builder import build_dem_packing_solver
        assert callable(build_dem_packing_solver)

    def test_builder_uses_subpackage_lbm(self):
        """builder.py must import from particulate_flow.lbm or subpackages."""
        builder_path = Path("src/particulate_flow/builder.py")
        source = builder_path.read_text()
        assert "particulate_flow.lbm" in source or "from .lbm" in source or "from ..lbm" in source

    def test_builder_uses_subpackage_dem(self):
        builder_path = Path("src/particulate_flow/builder.py")
        source = builder_path.read_text()
        assert "particulate_flow.dem" in source or "from .dem" in source or "from ..dem" in source


# ---------------------------------------------------------------------------
# Acceptance Scenario 2: solver section in SimulationConfig
# ---------------------------------------------------------------------------

class TestSolverSection:
    def test_solver_section_accepted(self):
        from particulate_flow.io.config import SimulationConfig
        cfg = SimulationConfig.from_mapping({"solver": {"fluid_method": "lbm-bgk-guo"}})
        assert cfg.values["fluid_method"] == "lbm-bgk-guo"

    def test_solver_particle_method(self):
        from particulate_flow.io.config import SimulationConfig
        cfg = SimulationConfig.from_mapping({"solver": {"particle_method": "dem-hertz"}})
        assert cfg.values["particle_method"] == "dem-hertz"

    def test_solver_particle_fluid_coupling(self):
        from particulate_flow.io.config import SimulationConfig
        cfg = SimulationConfig.from_mapping({
            "solver": {"particle_fluid_coupling": "immersed_boundary"}
        })
        assert cfg.values["particle_fluid_coupling"] == "immersed_boundary"

    def test_solver_section_in_section_keys(self):
        from particulate_flow.io.config import SECTION_KEYS
        assert "solver" in SECTION_KEYS


# ---------------------------------------------------------------------------
# Acceptance Scenario 3: accelerator section mapping
# ---------------------------------------------------------------------------

class TestAcceleratorSection:
    def test_accelerator_fluid_maps_to_fluid_accelerator(self):
        from particulate_flow.io.config import SimulationConfig
        cfg = SimulationConfig.from_mapping({"accelerator": {"fluid": "numpy"}})
        assert cfg.values["fluid_accelerator"] == "numpy"

    def test_accelerator_compute_maps_to_compute_accelerator(self):
        from particulate_flow.io.config import SimulationConfig
        cfg = SimulationConfig.from_mapping({"accelerator": {"compute": "numpy"}})
        assert cfg.values["compute_accelerator"] == "numpy"

    def test_accelerator_auto(self):
        from particulate_flow.io.config import SimulationConfig
        cfg = SimulationConfig.from_mapping({"accelerator": {"fluid": "auto", "compute": "auto"}})
        assert cfg.values["fluid_accelerator"] == "auto"
        assert cfg.values["compute_accelerator"] == "auto"

    def test_accelerator_section_in_section_keys(self):
        from particulate_flow.io.config import SECTION_KEYS
        assert "accelerator" in SECTION_KEYS

    def test_numerics_backward_compat_fluid_method(self):
        """fluid_method in numerics section still works."""
        from particulate_flow.io.config import SimulationConfig
        cfg = SimulationConfig.from_mapping({"numerics": {"fluid_method": "lbm-trt-guo"}})
        assert cfg.values["fluid_method"] == "lbm-trt-guo"

    def test_numerics_backward_compat_fluid_accelerator(self):
        from particulate_flow.io.config import SimulationConfig
        cfg = SimulationConfig.from_mapping({"numerics": {"fluid_accelerator": "numba"}})
        assert cfg.values["fluid_accelerator"] == "numba"


# ---------------------------------------------------------------------------
# Acceptance Scenario 4: build_dem_packing_solver returns DEMPackingSimulation
# ---------------------------------------------------------------------------

class TestBuildDEMPackingSolver:
    def test_build_dem_packing_solver_returns_simulation(self):
        from particulate_flow.builder import build_dem_packing_solver
        from particulate_flow.dem.packing import DEMPackingSimulation

        args = argparse.Namespace(
            nx=50, ny=50,
            n_particles=5,
            particle_radius=2.0,
            particle_density=1.0,
            gravity=9.81,
            e_n=0.5, e_t=0.3, mu_friction=0.3,
            k_n=1000.0, k_t=500.0,
            packing_steps=10,
            packing_dt=1e-4,
            cylinders=None,
        )
        sim = build_dem_packing_solver(args)
        assert isinstance(sim, DEMPackingSimulation)


# ---------------------------------------------------------------------------
# Acceptance Scenario 5: run_lbm_dem.py has no inline solver construction
# ---------------------------------------------------------------------------

class TestRunnerEntryPointOnly:
    def test_run_lbm_dem_no_inline_fastlbmdem_construction(self):
        """run_lbm_dem.py must not contain FastLBMDEM(nx= ... ) inline construction."""
        runner_path = Path("src/runners/run_lbm_dem.py")
        source = runner_path.read_text()
        # The inline construction pattern is FastLBMDEM(nx=
        assert "FastLBMDEM(nx=" not in source, (
            "run_lbm_dem.py still contains inline FastLBMDEM construction; "
            "use build_lbm_dem_solver(args) instead"
        )

    def test_run_lbm_dem_calls_builder(self):
        runner_path = Path("src/runners/run_lbm_dem.py")
        source = runner_path.read_text()
        assert "build_lbm_dem_solver" in source

    def test_run_dem_packing_no_inline_construction(self):
        runner_path = Path("src/runners/run_dem_packing.py")
        source = runner_path.read_text()
        assert "DEMPackingSimulation(" not in source, (
            "run_dem_packing.py still contains inline DEMPackingSimulation construction"
        )

    def test_run_dem_packing_calls_builder(self):
        runner_path = Path("src/runners/run_dem_packing.py")
        source = runner_path.read_text()
        assert "build_dem_packing_solver" in source
