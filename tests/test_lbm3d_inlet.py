"""Tests for issue-20: 3D left-inlet particle injection/removal.

Acceptance scenarios:
1. Particles appear at the inlet plane over time and advect downstream.
2. Particles past the outlet are removed.
3. Injected volume fraction tracks the configured source fraction.

Plus DEM3D add/remove unit tests and a periodic-mode guard.
"""

from __future__ import annotations

import numpy as np
import pytest

from particulate_flow.dem.contact3d import DEM3D
from particulate_flow.lbm3d import LBMDEMSolver3D


def _inlet_solver(*, source_volume_fraction=0.05, pressure_drop=4e-3, **kwargs):
    nx, ny, nz = 48, 16, 16
    return LBMDEMSolver3D(
        nx=nx,
        ny=ny,
        nz=nz,
        Re=10.0,
        u_max=0.05,
        streamwise_boundary="pressure",
        pressure_drop=pressure_drop,
        rho_out=1.0,
        particle_source="left_inlet",
        particle_fluid_coupling="immersed_boundary",
        particle_radius=2.0,
        density_ratio=2.0,
        gravity=0.0,
        source_volume_fraction=source_volume_fraction,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# DEM3D dynamic array unit tests
# ---------------------------------------------------------------------------


class TestDEM3DAddRemove:
    def _solver(self, n=2):
        pos = np.array([[float(i), 5.0, 5.0] for i in range(n)])
        return DEM3D(
            pos=pos,
            vel=np.zeros((n, 3)),
            radii=np.full(n, 2.0),
            nx=20,
            ny=20,
            nz=20,
        )

    def test_add_grows_all_arrays(self):
        sim = self._solver(2)
        sim.add_particles(
            np.array([[10.0, 10.0, 10.0]]),
            np.array([[0.1, 0.0, 0.0]]),
            np.array([3.0]),
        )
        assert sim.n_p == 3
        for arr in (sim.pos, sim.vel, sim.omega):
            assert arr.shape == (3, 3)
        assert sim.radii.shape == (3,)
        assert sim.masses.shape == (3,)
        assert sim.inertias.shape == (3,)
        # New particle's omega is zero and mass matches the sphere formula.
        np.testing.assert_allclose(sim.omega[2], 0.0)
        expected_mass = sim.density_ratio * (4.0 / 3.0) * np.pi * 3.0**3
        assert sim.masses[2] == pytest.approx(expected_mass)

    def test_remove_keeps_complement(self):
        sim = self._solver(3)
        sim.remove_particles(np.array([False, True, False]))
        assert sim.n_p == 2
        assert sim.pos.shape == (2, 3)
        # The surviving particles are indices 0 and 2 (x = 0 and 2).
        np.testing.assert_allclose(sim.pos[:, 0], [0.0, 2.0])

    def test_empty_ops_are_noops(self):
        sim = self._solver(2)
        sim.add_particles(np.empty((0, 3)), np.empty((0, 3)), np.empty(0))
        assert sim.n_p == 2
        sim.remove_particles(np.zeros(2, dtype=bool))
        assert sim.n_p == 2

    def test_remove_rejects_wrong_size_mask(self):
        sim = self._solver(3)
        with pytest.raises(ValueError, match="mask size"):
            sim.remove_particles(np.array([True, False]))  # size 2 != n_p 3

    def test_add_then_remove_round_trip(self):
        sim = self._solver(2)
        sim.add_particles(
            np.array([[10.0, 10.0, 10.0]]),
            np.zeros((1, 3)),
            np.array([2.0]),
        )
        assert sim.n_p == 3
        # Remove the just-added particle; original two remain consistent.
        sim.remove_particles(np.array([False, False, True]))
        assert sim.n_p == 2
        np.testing.assert_allclose(sim.pos[:, 0], [0.0, 1.0])
        assert sim.masses.shape == (2,) and sim.inertias.shape == (2,)


# ---------------------------------------------------------------------------
# Scenario 1: particles appear and advect downstream
# ---------------------------------------------------------------------------


class TestInjectionAndAdvection:
    def test_particles_appear_over_time(self):
        sim = _inlet_solver()
        assert sim.dem.n_p == 0
        sim.advance(300)
        assert sim.dem.n_p > 0, "inlet should inject particles once flow develops"
        assert np.all(np.isfinite(sim.dem.pos))

    def test_injected_particles_advect_downstream(self):
        sim = _inlet_solver()
        sim.advance(200)
        x_mean_early = sim.dem.pos[:, 0].mean() if sim.dem.n_p else 0.0
        sim.advance(200)
        # Mean x should increase as particles move toward the outlet.
        assert sim.dem.n_p > 0
        assert sim.dem.pos[:, 0].mean() >= x_mean_early

    def test_particles_enter_near_inlet(self):
        sim = _inlet_solver()
        sim.advance(120)
        if sim.dem.n_p:
            # Newly injected particles start near x=0 (within a few radii).
            assert sim.dem.pos[:, 0].min() < 8.0


# ---------------------------------------------------------------------------
# Scenario 2: outflow removal
# ---------------------------------------------------------------------------


class TestOutflowRemoval:
    def test_particles_past_outlet_removed(self):
        sim = _inlet_solver()
        sim.advance(1200)
        # No surviving particle should be fully past the outlet.
        if sim.dem.n_p:
            assert np.all(sim.dem.pos[:, 0] - sim.dem.radii <= sim.nx)
        # Removal actually happened over a long run.
        assert sim.removed_particles > 0

    def test_population_bounded(self):
        sim = _inlet_solver()
        sim.advance(800)
        n_mid = sim.dem.n_p
        sim.advance(800)
        # Population reaches a steady-ish state, not unbounded growth.
        assert sim.dem.n_p <= n_mid + 40


# ---------------------------------------------------------------------------
# Scenario 3: injected volume tracks phi
# ---------------------------------------------------------------------------


class TestVolumeFractionTracking:
    def test_higher_phi_injects_more(self):
        lo = _inlet_solver(source_volume_fraction=0.02)
        hi = _inlet_solver(source_volume_fraction=0.10)
        lo.advance(400)
        hi.advance(400)
        assert hi.injected_particle_volume > lo.injected_particle_volume

    def test_injected_volume_tracks_budget(self):
        sim = _inlet_solver(source_volume_fraction=0.05)
        sim.advance(400)
        # Injected solid volume should be close to phi * cumulative inlet flux,
        # within one particle volume of slack (budget is consumed in quanta).
        one_particle = (4.0 / 3.0) * np.pi * 2.0**3
        expected = 0.05 * sim.cumulative_inlet_flow_volume
        assert sim.injected_particle_volume <= expected + 1e-9
        assert expected - sim.injected_particle_volume < one_particle + 1e-9


# ---------------------------------------------------------------------------
# Guard + regression
# ---------------------------------------------------------------------------


class TestGuards:
    def test_left_inlet_requires_pressure_boundary(self):
        with pytest.raises(ValueError, match="pressure"):
            LBMDEMSolver3D(
                nx=32,
                ny=16,
                nz=16,
                le_shear_rate=0.0,
                streamwise_boundary="periodic",
                particle_source="left_inlet",
                particle_fluid_coupling="immersed_boundary",
                source_volume_fraction=0.05,
            )

    def test_left_inlet_requires_immersed_boundary_coupling(self):
        with pytest.raises(ValueError, match="immersed_boundary"):
            LBMDEMSolver3D(
                nx=32,
                ny=16,
                nz=16,
                streamwise_boundary="pressure",
                pressure_drop=1e-3,
                particle_source="left_inlet",
                particle_fluid_coupling="none",
                source_volume_fraction=0.05,
            )

    def test_left_inlet_requires_source_fraction(self):
        with pytest.raises(ValueError, match="source_volume_fraction"):
            LBMDEMSolver3D(
                nx=32,
                ny=16,
                nz=16,
                streamwise_boundary="pressure",
                pressure_drop=1e-3,
                particle_source="left_inlet",
                particle_fluid_coupling="immersed_boundary",
                source_volume_fraction=None,
            )
