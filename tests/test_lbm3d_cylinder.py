"""Tests for issue-19: 3D fixed finite-cylinder solid obstacles.

Acceptance scenarios:
1. A finite cylinder in 3D pressure flow blocks fluid (no-slip) and deflects it.
2. Particles collide with the fixed cylinder.

Plus mask-geometry and no-cylinder regression.
"""

from __future__ import annotations

import numpy as np

from particulate_flow.lbm3d import LBMDEMSolver3D


def _flow_solver(cylinders, *, pressure_drop=3e-3, **kwargs):
    nx, ny, nz = 40, 20, 12
    return LBMDEMSolver3D(
        nx=nx,
        ny=ny,
        nz=nz,
        Re=10.0,
        u_max=0.05,
        streamwise_boundary="pressure",
        pressure_drop=pressure_drop,
        rho_out=1.0,
        cylinders=cylinders,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Mask geometry
# ---------------------------------------------------------------------------


class TestSolidMask:
    def test_mask_marks_inside_radius_and_z_extent(self):
        # Cylinder at (20,10) r=4, z in [3,8].
        sim = _flow_solver([(20.0, 10.0, 4.0, 3.0, 8.0)])
        solid = sim.solid
        assert solid.shape == (40, 20, 12)
        # A cell on the axis within the z-extent is solid.
        assert solid[20, 10, 5]
        # Same x-y but outside the z-extent is fluid.
        assert not solid[20, 10, 0]
        assert not solid[20, 10, 11]
        # A cell well outside the radius is fluid.
        assert not solid[5, 5, 5]
        # A cell just outside the radius (dx=5 > 4) is fluid.
        assert not solid[25, 10, 5]

    def test_full_z_extent_default(self):
        sim = _flow_solver([(20.0, 10.0, 4.0)])  # no z bounds → full depth
        assert sim.solid[20, 10, 0]
        assert sim.solid[20, 10, 11]

    def test_no_cylinders_all_fluid(self):
        sim = _flow_solver([])
        assert not sim.solid.any()


# ---------------------------------------------------------------------------
# Scenario 1: no-slip + deflection
# ---------------------------------------------------------------------------


class TestFluidObstacle:
    def test_no_slip_inside_cylinder(self):
        sim = _flow_solver([(20.0, 10.0, 4.0)])
        sim.advance(400)
        _, ux, uy, uz = sim.get_fields()
        speed = np.sqrt(ux**2 + uy**2 + uz**2)
        free_stream = speed[~sim.solid].max()
        # Halfway bounce-back enforces no-slip at the staircased surface: velocity
        # inside the obstacle is strongly suppressed relative to the free stream
        # (residual is a known staircase artefact, not free flow).
        assert speed[sim.solid].max() < 0.2 * free_stream
        assert speed[sim.solid].mean() < 0.05 * free_stream
        assert np.all(np.isfinite(ux))

    def test_flow_is_deflected_around_cylinder(self):
        sim = _flow_solver([(20.0, 10.0, 4.0)])
        sim.advance(400)
        _, ux, uy, uz = sim.get_fields()
        # Transverse (y) velocity should be nonzero near the cylinder flanks
        # (flow diverts around the obstacle).
        near = uy[16:25, :, :]
        assert np.abs(near).max() > 1e-4

    def test_centreline_slows_versus_open_channel(self):
        cyl = _flow_solver([(20.0, 10.0, 4.0)])
        opn = _flow_solver([])
        cyl.advance(400)
        opn.advance(400)
        _, ux_c, _, _ = cyl.get_fields()
        _, ux_o, _, _ = opn.get_fields()
        # Just downstream of the cylinder, on the centreline, the obstacle case
        # has lower streamwise velocity than the open channel.
        assert ux_c[26, 10, 6] < ux_o[26, 10, 6]


# ---------------------------------------------------------------------------
# Scenario 2: particle collides with the cylinder
# ---------------------------------------------------------------------------


class TestParticleCollision:
    def test_particle_does_not_penetrate_cylinder(self):
        cx, cy, cr = 20.0, 10.0, 4.0
        rp = 2.0
        # Particle starts upstream on the centreline, will be advected toward the
        # cylinder by the pressure flow.
        pos = np.array([[10.0, cy, 6.0]])
        sim = _flow_solver(
            [(cx, cy, cr)],
            n_particles=1,
            particle_fluid_coupling="immersed_boundary",
            particle_positions=pos,
            particle_radius=rp,
            density_ratio=2.0,
            gravity=0.0,
        )
        sim.advance(600)
        p = sim.dem.pos[0]
        radial = np.hypot(p[0] - cx, p[1] - cy)
        # Must not penetrate inside (r_cyl + r_p), allowing a small contact overlap.
        assert radial > (cr + rp) - 0.75
        assert np.all(np.isfinite(sim.dem.pos))
