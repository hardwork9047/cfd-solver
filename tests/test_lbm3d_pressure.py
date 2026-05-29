"""Tests for issue-14: 3D pressure-driven flow vertical slice.

Covers acceptance scenarios 1, 5, 6 from issue #14 (the fluid-side path and
regression), plus guard tests that the not-yet-implemented particle stack
(IBM/DEM/obstacles/inlet — follow-ups #17-#20) raises informatively rather than
silently no-opping.
"""

from __future__ import annotations

import argparse

import numpy as np
import pytest

from particulate_flow.lbm3d import CS2_3, LBMDEMSolver3D

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pressure_solver(*, nx=24, ny=12, nz=8, pressure_drop=2e-3, **kwargs):
    """Build a 3D pressure-driven solver (no shear, x=pressure)."""
    return LBMDEMSolver3D(
        nx=nx,
        ny=ny,
        nz=nz,
        Re=10.0,
        u_max=0.05,
        le_shear_rate=0.0,
        streamwise_boundary="pressure",
        pressure_drop=pressure_drop,
        rho_out=1.0,
        **kwargs,
    )


def _builder_args_3d_pressure(**overrides):
    args = argparse.Namespace(
        nx=24,
        ny=12,
        nz=8,
        dimensions=3,
        reynolds_number=10.0,
        u_max=0.05,
        le_shear_rate=0.0,
        le_shear_axis=0,
        le_boundary_axis=1,
        le_interpolation_order=3,
        streamwise_boundary="pressure",
        pressure_drop=2e-3,
        rho_out=1.0,
    )
    for k, v in overrides.items():
        setattr(args, k, v)
    return args


# ---------------------------------------------------------------------------
# Scenario 1: 3D pressure-driven flow works
# ---------------------------------------------------------------------------


class TestPressureDrivenFlow:
    def test_rho_in_exceeds_rho_out(self):
        """Inlet density is set above outlet by pressure_drop / cs^2."""
        sim = _pressure_solver(pressure_drop=3e-3)
        assert sim.rho_in == pytest.approx(1.0 + 3e-3 / CS2_3)
        assert sim.rho_in > sim.rho_out

    def test_positive_bulk_velocity_develops(self):
        """A positive pressure drop drives net +x flow."""
        sim = _pressure_solver(pressure_drop=2e-3)
        sim.advance(300)
        _, ux, _, _ = sim.get_fields()
        # Bulk streamwise velocity in the interior should be positive.
        bulk_ux = ux[1:-1].mean()
        assert bulk_ux > 1e-5, f"expected positive bulk ux, got {bulk_ux:.2e}"

    def test_density_decreases_inlet_to_outlet(self):
        """Pressure (density) falls monotonically from inlet to outlet."""
        sim = _pressure_solver(pressure_drop=2e-3)
        sim.advance(300)
        rho, *_ = sim.get_fields()
        rho_x = rho.mean(axis=(1, 2))  # profile along x
        # Inlet plane denser than outlet plane.
        assert rho_x[0] > rho_x[-1]
        # Broadly monotone decreasing (allow tiny non-monotonic numerical wiggle).
        diffs = np.diff(rho_x)
        assert (diffs <= 1e-6).mean() > 0.8

    def test_stays_finite_and_stable(self):
        sim = _pressure_solver(pressure_drop=2e-3)
        sim.advance(300)
        rho, ux, uy, uz = sim.get_fields()
        for field in (rho, ux, uy, uz):
            assert np.all(np.isfinite(field))
        assert np.abs(ux).max() < 1.0  # subsonic / stable


# ---------------------------------------------------------------------------
# Scenario 5: builder dispatch forwards pressure params
# ---------------------------------------------------------------------------


class TestBuilderPressureForwarding:
    def test_builder_returns_3d_pressure_solver(self):
        from particulate_flow.builder import build_lbm_dem_solver

        solver = build_lbm_dem_solver(_builder_args_3d_pressure())
        assert isinstance(solver, LBMDEMSolver3D)
        assert solver.rho_in > solver.rho_out
        assert solver.streamwise_boundary == "pressure"


# ---------------------------------------------------------------------------
# Scenario 6: 2D + existing 3D LE behaviour unchanged
# ---------------------------------------------------------------------------


class TestRegression:
    def test_default_periodic_unchanged(self):
        """Default streamwise_boundary stays periodic (LE path untouched)."""
        sim = LBMDEMSolver3D(nx=8, ny=16, nz=4, le_shear_rate=0.001)
        assert sim.streamwise_boundary == "periodic"
        # LE shear run still works.
        sim.advance(20)
        rho, *_ = sim.get_fields()
        np.testing.assert_allclose(rho, 1.0, atol=1e-3)


# ---------------------------------------------------------------------------
# Scenarios 2-4: particle stack must fail loudly (not silently absent)
# ---------------------------------------------------------------------------


class TestParticleStackStubs:
    def test_particles_raise_not_implemented(self):
        with pytest.raises(NotImplementedError, match=r"#1[78]"):
            LBMDEMSolver3D(nx=16, ny=12, nz=8, n_particles=5)

    def test_cylinders_now_build_solid_mask(self):
        # Cylinders are implemented as of issue #19: they build a solid mask
        # rather than raising. (Kept here to document the stub's retirement.)
        sim = LBMDEMSolver3D(
            nx=16,
            ny=12,
            nz=8,
            cylinders=[(8.0, 6.0, 3.0)],
        )
        assert sim.solid.any()

    def test_inlet_source_now_supported_with_pressure(self):
        # Inlet injection is implemented as of issue #20 (requires pressure flow).
        # (Kept here to document the stub's retirement.)
        sim = LBMDEMSolver3D(
            nx=16,
            ny=12,
            nz=8,
            streamwise_boundary="pressure",
            pressure_drop=1e-3,
            particle_source="left_inlet",
            particle_fluid_coupling="immersed_boundary",
            source_volume_fraction=0.05,
        )
        assert sim.particle_source == "left_inlet"
        assert sim.dem is not None and sim.dem.n_p == 0

    def test_coupling_raises_not_implemented(self):
        with pytest.raises(NotImplementedError, match=r"#17"):
            LBMDEMSolver3D(
                nx=16,
                ny=12,
                nz=8,
                particle_fluid_coupling="ibm",
            )
