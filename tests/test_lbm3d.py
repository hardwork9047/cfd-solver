"""Tests for issue-9: 3D D3Q15 LBM solver with Lees-Edwards BC."""

from __future__ import annotations

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_3d_solver(
    *,
    nx: int = 16,
    ny: int = 32,
    nz: int = 8,
    le_shear_rate: float = 0.001,
    **kwargs,
):
    from particulate_flow.lbm3d import LBMDEMSolver3D

    return LBMDEMSolver3D(
        nx=nx,
        ny=ny,
        nz=nz,
        le_shear_rate=le_shear_rate,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# D3Q15 lattice constants
# ---------------------------------------------------------------------------

class TestD3Q15Constants:
    def test_15_directions(self):
        from particulate_flow.lbm3d import C3, W3, Q3

        assert Q3 == 15
        assert C3.shape == (15, 3)
        assert W3.shape == (15,)

    def test_weights_sum_to_one(self):
        from particulate_flow.lbm3d import W3

        assert pytest.approx(W3.sum(), rel=1e-10) == 1.0

    def test_opposite_directions(self):
        from particulate_flow.lbm3d import C3, OPPOSITE3

        for i, j in enumerate(OPPOSITE3):
            np.testing.assert_array_equal(C3[i] + C3[j], [0, 0, 0])

    def test_cs2_is_one_third(self):
        from particulate_flow.lbm3d import CS2_3

        assert CS2_3 == pytest.approx(1.0 / 3.0)


# ---------------------------------------------------------------------------
# LBMDEMSolver3D construction
# ---------------------------------------------------------------------------

class TestLBMDEMSolver3DInit:
    def test_f_shape(self):
        sim = _make_3d_solver(nx=8, ny=16, nz=4)
        assert sim.f.shape == (15, 8, 16, 4)

    def test_step_count_zero(self):
        sim = _make_3d_solver()
        assert sim.step_count == 0

    def test_get_fields_returns_four_arrays(self):
        sim = _make_3d_solver(nx=8, ny=16, nz=4)
        result = sim.get_fields()
        assert len(result) == 4
        rho, ux, uy, uz = result
        assert rho.shape == (8, 16, 4)
        assert ux.shape == (8, 16, 4)
        assert uy.shape == (8, 16, 4)
        assert uz.shape == (8, 16, 4)


# ---------------------------------------------------------------------------
# Advance / step
# ---------------------------------------------------------------------------

class TestAdvance:
    def test_step_count_increments(self):
        sim = _make_3d_solver()
        sim.advance(5)
        assert sim.step_count == 5

    def test_advance_does_not_crash(self):
        sim = _make_3d_solver(nx=8, ny=16, nz=4)
        sim.advance(10)

    def test_rho_stays_close_to_one(self):
        sim = _make_3d_solver(nx=8, ny=16, nz=4)
        sim.advance(20)
        rho, *_ = sim.get_fields()
        np.testing.assert_allclose(rho, 1.0, atol=1e-4)


# ---------------------------------------------------------------------------
# Linear shear profile (L2 < 1%)
# ---------------------------------------------------------------------------

class TestLinearShearProfile:
    def test_ux_linear_in_y_after_spinup(self):
        """Starting from analytical profile, ux should stay linear (L2 < 1%)."""
        ny = 32
        nx = 16
        nz = 8
        shear_rate = 0.001

        from particulate_flow.lbm3d import LBMDEMSolver3D

        sim = LBMDEMSolver3D(
            nx=nx,
            ny=ny,
            nz=nz,
            le_shear_rate=shear_rate,
            init_analytical=True,  # start from ux = gamma_dot * (y - ny/2)
        )
        sim.advance(200)

        _, ux, _, _ = sim.get_fields()
        # Expected: ux(y) = shear_rate * (y - (ny-1)/2), averaged over x and z
        ux_mean = ux.mean(axis=(0, 2))  # shape (ny,)
        y = np.arange(ny, dtype=float)
        ux_analytical = shear_rate * (y - (ny - 1) / 2.0)
        # Avoid zero-denominator in relative error
        norm = np.sqrt(np.mean(ux_analytical**2))
        l2_rel = np.sqrt(np.mean((ux_mean - ux_analytical) ** 2)) / norm
        assert l2_rel < 0.01, f"L2 relative error {l2_rel:.4f} >= 1%"

    def test_uy_uz_near_zero_at_steady_state(self):
        """Transverse velocity components should be negligible in pure shear."""
        from particulate_flow.lbm3d import LBMDEMSolver3D

        sim = LBMDEMSolver3D(
            nx=8, ny=16, nz=4,
            le_shear_rate=0.001,
            init_analytical=True,
        )
        sim.advance(100)
        _, _, uy, uz = sim.get_fields()
        assert np.abs(uy).max() < 1e-6
        assert np.abs(uz).max() < 1e-6


# ---------------------------------------------------------------------------
# Config plumbing
# ---------------------------------------------------------------------------

class TestConfigDimensions:
    def test_dimensions_2_in_config(self):
        """dimensions=2 passes through SimulationConfig."""
        from particulate_flow.io.config import SimulationConfig

        cfg = SimulationConfig.from_mapping(
            {"domain": {"nx": 40, "ny": 30, "dimensions": 2}}
        )
        assert cfg.values.get("dimensions") == 2

    def test_dimensions_3_in_config(self):
        """dimensions=3 passes through SimulationConfig."""
        from particulate_flow.io.config import SimulationConfig

        cfg = SimulationConfig.from_mapping(
            {"domain": {"nx": 20, "ny": 32, "nz": 8, "dimensions": 3}}
        )
        assert cfg.values.get("dimensions") == 3
        assert cfg.values.get("nz") == 8

    def test_builder_returns_3d_solver_for_dimensions_3(self):
        """build_lbm_dem_solver returns LBMDEMSolver3D when dimensions=3."""
        import argparse
        from particulate_flow.builder import build_lbm_dem_solver
        from particulate_flow.lbm3d import LBMDEMSolver3D

        args = argparse.Namespace(
            nx=16, ny=32, nz=8,
            dimensions=3,
            reynolds_number=10.0,
            u_max=0.0,
            le_shear_rate=0.001,
            le_shear_axis=0,
            le_boundary_axis=1,
            le_interpolation_order=3,
        )
        solver = build_lbm_dem_solver(args)
        assert isinstance(solver, LBMDEMSolver3D)

    def test_builder_returns_2d_solver_for_dimensions_2(self):
        """build_lbm_dem_solver returns FastLBMDEM when dimensions=2 (regression)."""
        import argparse
        from particulate_flow.builder import build_lbm_dem_solver
        from particulate_flow.fast_solver import FastLBMDEM

        args = argparse.Namespace(
            nx=40, ny=30,
            dimensions=2,
            reynolds_number=100.0,
            u_max=0.05,
            reynolds_length=None,
            flow_condition=None,
            flow_control=False,
            flow_control_gain=0.2,
            y_boundary="wall",
            streamwise_boundary="pressure",
            pressure_drop=1e-4,
            rho_out=1.0,
            fluid_method="lbm-bgk-guo",
            fluid_accelerator="numpy",
            compute_accelerator="numpy",
            particle_method="dem-hertz",
            particle_search="cell_list",
            particle_fluid_coupling="point_force",
            ibm_stiffness=1.0,
            ibm_marker_spacing=1.0,
            n_particles=0,
            particle_radius=3.0,
            radius_variation=0.0,
            density_ratio=2.0,
            gravity=0.0,
            dem_substeps=4,
            rolling_friction=False,
            sliding_friction=0.5,
            tangential_damping=0.4,
            rolling_friction_coeff=0.05,
            rolling_damping=0.2,
            particle_attraction=False,
            particle_repulsion=False,
            attraction_strength=1e-3,
            repulsion_strength=1e-3,
            attraction_cutoff=3.0,
            repulsion_cutoff=3.0,
            attraction_min_gap=0.05,
            repulsion_min_gap=0.05,
            porous_resistance=False,
            porous_resistance_coeff=0.0,
            particle_volume_fraction=None,
            particle_source="initial",
            total_steps=100,
            warmup_steps=0,
            cylinder=False,
            cyl_x=None,
            cyl_y=None,
            cyl_r=6.0,
            cylinder_spec=None,
            le_shear_rate=0.0,
            le_shear_axis=0,
            le_boundary_axis=1,
            le_interpolation_order=3,
            surface_roughness=0.0,
        )
        solver = build_lbm_dem_solver(args)
        assert isinstance(solver, FastLBMDEM)
