"""Tests for subpackage structure (Acceptance Scenarios 1, 2, 3)."""

from __future__ import annotations

import pytest


# ---------------------------------------------------------------------------
# Acceptance Scenario 1: サブパッケージ import が動作する
# ---------------------------------------------------------------------------

class TestLBMSubpackageImports:
    def test_lbm_constants(self):
        from particulate_flow.lbm.constants import C, W, OPPOSITE, Q, CS2
        assert Q == 9
        assert len(W) == 9
        assert len(C) == 9

    def test_lbm_operators(self):
        from particulate_flow.lbm.operators import equilibrium, guo_forcing
        assert callable(equilibrium)
        assert callable(guo_forcing)

    def test_lbm_collision(self):
        from particulate_flow.lbm import collision
        assert hasattr(collision, "FLUID_METHODS")

    def test_lbm_boundary(self):
        from particulate_flow.lbm import boundary
        assert boundary is not None

    def test_lbm_macroscopic(self):
        from particulate_flow.lbm import macroscopic
        assert macroscopic is not None

    def test_lbm_kernels_importable(self):
        from particulate_flow.lbm import kernels  # Numba optional
        assert kernels is not None


class TestDEMSubpackageImports:
    def test_dem_solver(self):
        from particulate_flow.dem.solver import DEMSolver
        assert DEMSolver is not None

    def test_dem_packing(self):
        from particulate_flow.dem.packing import DEMPackingSimulation, PackingMetrics
        assert DEMPackingSimulation is not None

    def test_dem_contact(self):
        from particulate_flow.dem import contact
        assert contact is not None

    def test_dem_particle_manager(self):
        from particulate_flow.dem import particle_manager
        assert particle_manager is not None

    def test_dem_kernels_importable(self):
        from particulate_flow.dem import kernels
        assert kernels is not None


class TestIBMSubpackageImports:
    def test_ibm_coupling(self):
        from particulate_flow.ibm import coupling
        assert coupling is not None

    def test_ibm_kernels_importable(self):
        from particulate_flow.ibm import kernels
        assert kernels is not None


class TestGeometrySubpackageImports:
    def test_geometry_pore(self):
        from particulate_flow.geometry.pore import PoreGeometry, Cylinder
        assert PoreGeometry is not None
        assert Cylinder is not None


class TestIOSubpackageImports:
    def test_io_config(self):
        from particulate_flow.io.config import SimulationConfig
        assert SimulationConfig is not None

    def test_io_paths(self):
        from particulate_flow.io.paths import program_results_dir
        assert callable(program_results_dir)

    def test_io_verification(self):
        from particulate_flow.io.verification import VerificationResult, run_fluid_verification
        assert VerificationResult is not None

    def test_io_visualization(self):
        from particulate_flow.io.visualization import plot_fields, plot_particles
        assert callable(plot_fields)
        assert callable(plot_particles)


# ---------------------------------------------------------------------------
# Acceptance Scenario 2: 既存の公開 API が後方互換で動作する
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    def test_top_level_lbm_dem_solver(self):
        from particulate_flow import LBMDEMSolver
        assert LBMDEMSolver is not None

    def test_top_level_fast_solver(self):
        from particulate_flow import FastLBMDEM
        assert FastLBMDEM is not None

    def test_top_level_dem_solver(self):
        from particulate_flow import DEMSolver
        assert DEMSolver is not None

    def test_top_level_geometry(self):
        from particulate_flow import PoreGeometry, Cylinder
        assert PoreGeometry is not None
        assert Cylinder is not None

    def test_top_level_config(self):
        from particulate_flow import SimulationConfig
        assert SimulationConfig is not None

    def test_top_level_packing(self):
        from particulate_flow import DEMPackingSimulation, PackingMetrics
        assert DEMPackingSimulation is not None

    def test_top_level_plot_functions(self):
        from particulate_flow import plot_fields, plot_particles
        assert callable(plot_fields)
        assert callable(plot_particles)

    def test_top_level_verification(self):
        from particulate_flow import VerificationResult, run_fluid_verification
        assert VerificationResult is not None


# ---------------------------------------------------------------------------
# Acceptance Scenario 3: NumPy/Numba バックエンド切り替えが機能する
# ---------------------------------------------------------------------------

class TestAcceleratorConstants:
    def test_fluid_accelerators_defined(self):
        from particulate_flow.lbm.collision import FLUID_ACCELERATORS
        assert "numpy" in FLUID_ACCELERATORS
        assert "numba" in FLUID_ACCELERATORS
        assert "auto" in FLUID_ACCELERATORS

    def test_compute_accelerators_defined(self):
        from particulate_flow.lbm.constants import COMPUTE_ACCELERATORS
        assert "numpy" in COMPUTE_ACCELERATORS
        assert "numba" in COMPUTE_ACCELERATORS
