"""Tests for issue-17: 3D IBM particle-fluid coupling.

Acceptance scenarios:
1. A particle in a 3D pressure-driven flow experiences fluid drag.
2. Newton's 3rd-law back-reaction is applied to the fluid body force.
3. 2D IBM behaviour unchanged.

Plus a DEM3D external-force hook test and coupled-run stability.
"""

from __future__ import annotations

import numpy as np

from particulate_flow.dem.contact3d import DEM3D
from particulate_flow.lbm3d import LBMDEMSolver3D


def _coupled_solver(*, pressure_drop=2e-3, n_particles=1, **kwargs):
    """Build a 3D pressure-driven solver with IBM-coupled particles."""
    nx, ny, nz = 32, 16, 16
    # One particle near the inlet, centred in y-z.
    pos = np.array([[8.0, ny / 2.0, nz / 2.0]])
    return LBMDEMSolver3D(
        nx=nx,
        ny=ny,
        nz=nz,
        Re=10.0,
        u_max=0.05,
        streamwise_boundary="pressure",
        pressure_drop=pressure_drop,
        rho_out=1.0,
        n_particles=n_particles,
        particle_fluid_coupling="immersed_boundary",
        particle_positions=pos,
        particle_radius=2.0,
        density_ratio=2.0,
        gravity=0.0,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Scenario 1: particle experiences fluid drag
# ---------------------------------------------------------------------------


class TestParticleDrag:
    def test_particle_advects_downstream(self):
        sim = _coupled_solver(pressure_drop=3e-3)
        x0 = sim.dem.pos[0, 0]
        vx0 = sim.dem.vel[0, 0]
        sim.advance(400)
        assert sim.dem.vel[0, 0] > vx0, "drag should accelerate particle in +x"
        assert sim.dem.pos[0, 0] > x0, "particle should move downstream"
        assert np.all(np.isfinite(sim.dem.pos))

    def test_quiescent_fluid_leaves_particle_at_rest(self):
        # No pressure drop, no gravity → no drag → particle stays put.
        sim = _coupled_solver(pressure_drop=0.0)
        p0 = sim.dem.pos.copy()
        sim.advance(100)
        np.testing.assert_allclose(sim.dem.pos, p0, atol=1e-6)


# ---------------------------------------------------------------------------
# Scenario 2: Newton's 3rd-law back-reaction on the fluid
# ---------------------------------------------------------------------------


class TestBackReaction:
    def test_body_force_field_populated_near_particle(self):
        sim = _coupled_solver(pressure_drop=3e-3)
        sim.advance(20)
        # The IBM body-force field must be nonzero somewhere near the particle.
        assert np.any(np.abs(sim.Fx) > 0.0)
        # Force concentrated near the particle, not spread across the whole domain.
        cx = int(round(sim.dem.pos[0, 0]))
        near = np.abs(sim.Fx[max(cx - 4, 0) : cx + 4]).sum()
        total = np.abs(sim.Fx).sum()
        assert near > 0.5 * total

    def test_reaction_balances_spread_force(self):
        # The force spread to the fluid must equal minus the reaction on particles
        # (Newton's 3rd law), to tolerance, for the IBM exchange in one step.
        sim = _coupled_solver(pressure_drop=3e-3)
        spread, reaction = sim._ibm_force_audit()
        np.testing.assert_allclose(spread + reaction, 0.0, atol=1e-9)


# ---------------------------------------------------------------------------
# DEM3D external-force hook
# ---------------------------------------------------------------------------


class TestDEM3DExternalForce:
    def test_external_force_moves_particle(self):
        pos = np.array([[10.0, 10.0, 10.0]])
        vel = np.zeros((1, 3))
        sim = DEM3D(pos=pos, vel=vel, radii=np.array([2.0]), nx=20, ny=20, nz=20, gravity=0.0)
        ext = np.array([[0.5, 0.0, 0.0]])
        sim.step(5, external_forces=ext)
        assert sim.vel[0, 0] > 0.0
        assert sim.pos[0, 0] > 10.0

    def test_no_external_force_is_inert(self):
        pos = np.array([[10.0, 10.0, 10.0]])
        vel = np.zeros((1, 3))
        sim = DEM3D(pos=pos, vel=vel, radii=np.array([2.0]), nx=20, ny=20, nz=20, gravity=0.0)
        sim.step(5)
        np.testing.assert_allclose(sim.vel, 0.0, atol=1e-12)
        np.testing.assert_allclose(sim.pos, [[10.0, 10.0, 10.0]], atol=1e-12)


# ---------------------------------------------------------------------------
# Scenario 3 + stability
# ---------------------------------------------------------------------------


class TestRegressionAndStability:
    def test_2d_ibm_untouched(self):
        # The 2D IBM entry point and DEM solver remain importable and intact.
        from particulate_flow.lbm_dem import LBMDEMSolver

        assert hasattr(LBMDEMSolver, "_apply_immersed_boundary_forces")

    def test_coupled_run_is_stable(self):
        sim = _coupled_solver(pressure_drop=2e-3)
        sim.advance(300)
        rho, ux, uy, uz = sim.get_fields()
        for fld in (rho, ux, uy, uz):
            assert np.all(np.isfinite(fld))
        assert np.abs(ux).max() < 1.0
        assert np.all(np.isfinite(sim.dem.pos))
