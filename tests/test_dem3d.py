"""Tests for issue-18: fresh 3D DEM contact physics (fluid-free core).

Covers the three acceptance scenarios from issue #18:
1. Two approaching spheres repel via Hertz normal contact.
2. A sphere settles under gravity onto a wall and rests.
3. Tangential friction + rolling resistance act in 3D.

Plus momentum/finiteness sanity, a cylinder-deflection case, and a guard that
importing the 3D module does not perturb the 2D DEM solver.
"""

from __future__ import annotations

import numpy as np
import pytest

from particulate_flow.dem.contact3d import DEM3D


def _two_sphere_solver(gap: float, *, vx: float = 0.0, **kwargs) -> DEM3D:
    """Two equal spheres on the x-axis separated centre-to-centre by 2r - gap.

    Positive ``gap`` means overlap (contact); negative means a clear separation.
    """
    r = 3.0
    dist = 2 * r - gap
    pos = np.array([[0.0, 10.0, 10.0], [dist, 10.0, 10.0]])
    vel = np.array([[vx, 0.0, 0.0], [-vx, 0.0, 0.0]])
    return DEM3D(
        pos=pos,
        vel=vel,
        radii=np.full(2, r),
        nx=40,
        ny=20,
        nz=20,
        gravity=0.0,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Scenario 1: two approaching spheres repel
# ---------------------------------------------------------------------------


class TestSphereSphereNormal:
    def test_overlapping_spheres_repel(self):
        sim = _two_sphere_solver(gap=0.5)  # 0.5 lattice-unit overlap
        d0 = np.linalg.norm(sim.pos[1] - sim.pos[0])
        sim.step(50)
        d1 = np.linalg.norm(sim.pos[1] - sim.pos[0])
        assert d1 > d0, "overlapping spheres should push apart"
        # Velocities should point away from each other along x.
        assert sim.vel[0, 0] < 0.0 < sim.vel[1, 0]
        assert np.all(np.isfinite(sim.pos))

    def test_separated_spheres_feel_no_contact(self):
        sim = _two_sphere_solver(gap=-2.0)  # 2 units clear gap
        v_before = sim.vel.copy()
        forces, _ = sim.compute_loads()
        np.testing.assert_allclose(forces, 0.0, atol=1e-12)
        # Stepping should not change velocity (no gravity, no contact).
        sim.step(5)
        np.testing.assert_allclose(sim.vel, v_before, atol=1e-12)

    def test_head_on_collision_conserves_momentum(self):
        sim = _two_sphere_solver(gap=0.2, vx=0.01)
        p0 = (sim.masses[:, None] * sim.vel).sum(axis=0)
        sim.step(30)
        p1 = (sim.masses[:, None] * sim.vel).sum(axis=0)
        np.testing.assert_allclose(p1, p0, atol=1e-9)

    def test_shear_contact_spins_both_spheres(self):
        # Two overlapping spheres with opposite transverse (y) velocities create
        # a shearing contact: both partners must receive a tangential torque and
        # the linear forces must obey Newton's 3rd law.
        r = 3.0
        pos = np.array([[0.0, 10.0, 10.0], [5.5, 10.0, 10.0]])  # 0.5 overlap
        vel = np.array([[0.0, 0.02, 0.0], [0.0, -0.02, 0.0]])
        sim = DEM3D(
            pos=pos,
            vel=vel,
            radii=np.full(2, r),
            nx=40,
            ny=20,
            nz=20,
            gravity=0.0,
            sliding_friction=0.5,
        )
        forces, torques = sim.compute_loads()
        assert np.linalg.norm(torques[0]) > 0.0
        assert np.linalg.norm(torques[1]) > 0.0
        np.testing.assert_allclose(forces[0] + forces[1], 0.0, atol=1e-9)


# ---------------------------------------------------------------------------
# Scenario 2: sphere settles on a wall under gravity
# ---------------------------------------------------------------------------


class TestSettleOnWall:
    def test_sphere_rests_above_wall(self):
        r = 3.0
        pos = np.array([[10.0, 10.0, 10.0]])
        vel = np.zeros((1, 3))
        sim = DEM3D(
            pos=pos,
            vel=vel,
            radii=np.array([r]),
            nx=20,
            ny=20,
            nz=20,
            gravity=1e-3,
            walls=("y_min",),
            dem_substeps=8,
        )
        for _ in range(4000):
            sim.step(1)
        # Comes to rest just above the y=0 wall, not tunnelled through it.
        assert sim.pos[0, 1] > 0.0, "must not tunnel through the wall"
        assert sim.pos[0, 1] == pytest.approx(r, abs=0.5)
        assert abs(sim.vel[0, 1]) < 1e-3, "should be ~at rest"
        assert np.all(np.isfinite(sim.pos))


# ---------------------------------------------------------------------------
# Scenario 3: tangential friction + rolling resistance
# ---------------------------------------------------------------------------


class TestFrictionAndRolling:
    def test_tangential_force_opposes_slip_at_wall(self):
        r = 3.0
        # Sphere pressed into the y=0 wall (slightly overlapping) and sliding in +x.
        pos = np.array([[10.0, r - 0.2, 10.0]])
        vel = np.array([[0.05, -0.01, 0.0]])
        sim = DEM3D(
            pos=pos,
            vel=vel,
            radii=np.array([r]),
            nx=20,
            ny=20,
            nz=20,
            gravity=0.0,
            walls=("y_min",),
            sliding_friction=0.5,
        )
        forces, torques = sim.compute_loads()
        # Tangential (x) force should oppose the +x slip → negative.
        assert forces[0, 0] < 0.0
        # Coulomb bound: |f_t| <= mu * |f_n|.
        assert abs(forces[0, 0]) <= 0.5 * abs(forces[0, 1]) + 1e-9
        # A rolling-resistance torque should be present (nonzero vector).
        assert np.linalg.norm(torques[0]) > 0.0

    def test_no_tangential_without_contact(self):
        r = 3.0
        pos = np.array([[10.0, 10.0, 10.0]])  # far from any wall
        vel = np.array([[0.05, 0.0, 0.0]])
        sim = DEM3D(
            pos=pos,
            vel=vel,
            radii=np.array([r]),
            nx=20,
            ny=20,
            nz=20,
            gravity=0.0,
            walls=("y_min",),
        )
        forces, torques = sim.compute_loads()
        np.testing.assert_allclose(forces, 0.0, atol=1e-12)
        np.testing.assert_allclose(torques, 0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# Cylinder deflection
# ---------------------------------------------------------------------------


class TestCylinderContact:
    def test_sphere_deflected_by_cylinder(self):
        r = 3.0
        # z-aligned cylinder at (x=20, y=10), radius 4, full z-extent.
        # Sphere approaches from the left moving +x toward the axis.
        pos = np.array([[14.5, 10.0, 10.0]])  # overlapping: 20-14.5=5.5 < 4+3=7
        vel = np.array([[0.02, 0.0, 0.0]])
        sim = DEM3D(
            pos=pos,
            vel=vel,
            radii=np.array([r]),
            nx=40,
            ny=20,
            nz=20,
            gravity=0.0,
            cylinders=[(20.0, 10.0, 4.0)],
        )
        forces, _ = sim.compute_loads()
        # Contact pushes the sphere away from the axis → -x component.
        assert forces[0, 0] < 0.0
        assert np.all(np.isfinite(forces))


# ---------------------------------------------------------------------------
# 2D regression guard
# ---------------------------------------------------------------------------


class TestNo2DRegression:
    def test_importing_contact3d_does_not_break_2d_solver(self):
        # Importing the 3D module must not perturb the 2D DEM API surface.
        from particulate_flow.dem.solver import DEMSolver

        assert hasattr(DEMSolver, "tangent_from_normal")
        assert hasattr(DEMSolver, "compute_loads")
