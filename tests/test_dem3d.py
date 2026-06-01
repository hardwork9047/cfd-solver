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
        # A torque should be present (here from the tangential friction force).
        assert np.linalg.norm(torques[0]) > 0.0

    def test_rolling_resistance_opposes_spin(self):
        # A spinning sphere pressed into the wall (no slip velocity) should feel
        # a rolling-resistance torque that anti-aligns with its angular velocity.
        r = 3.0
        pos = np.array([[10.0, r - 0.2, 10.0]])  # overlapping the y=0 wall
        vel = np.zeros((1, 3))
        sim = DEM3D(
            pos=pos,
            vel=vel,
            radii=np.array([r]),
            nx=20,
            ny=20,
            nz=20,
            gravity=0.0,
            walls=("y_min",),
            rolling_friction_coeff=0.1,
        )
        sim.omega[0] = np.array([0.0, 0.0, 0.5])  # spin about +z
        _, torques = sim.compute_loads()
        # No tangential slip (centre at rest, surface velocity from spin is
        # tangential and produces friction, but the rolling term must oppose ω).
        # The z-component of the total torque must oppose the +z spin.
        assert torques[0, 2] < 0.0

    def test_rolling_resistance_vanishes_without_contact(self):
        r = 3.0
        pos = np.array([[10.0, 10.0, 10.0]])  # far from the wall
        vel = np.zeros((1, 3))
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
        sim.omega[0] = np.array([0.0, 0.0, 0.5])
        _, torques = sim.compute_loads()
        np.testing.assert_allclose(torques, 0.0, atol=1e-12)

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
# Hamaker attraction / repulsion (issue #36)
# ---------------------------------------------------------------------------


class TestSphereSphereAttraction:
    """Sphere-sphere Hamaker attraction and repulsion."""

    def test_attraction_pulls_nearby_spheres_together(self):
        # Two spheres with a clear gap (no contact) but within attraction cutoff.
        # attraction force should bring them closer.
        r = 3.0
        gap = -1.5  # 1.5 lattice units clear gap (negative = not overlapping)
        dist = 2 * r - gap  # 7.5
        pos = np.array([[0.0, 10.0, 10.0], [dist, 10.0, 10.0]])
        sim = DEM3D(
            pos=pos,
            vel=np.zeros((2, 3)),
            radii=np.full(2, r),
            nx=40,
            ny=20,
            nz=20,
            gravity=0.0,
            particle_attraction=True,
            attraction_strength=0.01,
            attraction_cutoff=3.0,
        )
        forces, _ = sim.compute_loads()
        # Sphere 0 should be pulled in +x direction, sphere 1 in -x.
        assert forces[0, 0] > 0.0
        assert forces[1, 0] < 0.0
        np.testing.assert_allclose(forces[0] + forces[1], 0.0, atol=1e-12)

    def test_attraction_zero_beyond_cutoff(self):
        r = 3.0
        gap = -5.0  # 5 units clear gap, beyond cutoff=3.0
        dist = 2 * r - gap  # 11.0
        pos = np.array([[0.0, 10.0, 10.0], [dist, 10.0, 10.0]])
        sim = DEM3D(
            pos=pos,
            vel=np.zeros((2, 3)),
            radii=np.full(2, r),
            nx=40,
            ny=20,
            nz=20,
            gravity=0.0,
            particle_attraction=True,
            attraction_strength=0.01,
            attraction_cutoff=3.0,
        )
        forces, _ = sim.compute_loads()
        np.testing.assert_allclose(forces, 0.0, atol=1e-12)

    def test_repulsion_pushes_nearby_spheres_apart(self):
        r = 3.0
        gap = -1.0  # 1 unit gap, within cutoff
        dist = 2 * r - gap  # 7.0
        pos = np.array([[0.0, 10.0, 10.0], [dist, 10.0, 10.0]])
        sim = DEM3D(
            pos=pos,
            vel=np.zeros((2, 3)),
            radii=np.full(2, r),
            nx=40,
            ny=20,
            nz=20,
            gravity=0.0,
            particle_repulsion=True,
            repulsion_strength=0.01,
            repulsion_cutoff=3.0,
        )
        forces, _ = sim.compute_loads()
        # Sphere 0 pushed in -x, sphere 1 in +x.
        assert forces[0, 0] < 0.0
        assert forces[1, 0] > 0.0
        np.testing.assert_allclose(forces[0] + forces[1], 0.0, atol=1e-12)

    def test_attraction_and_repulsion_mutually_exclusive(self):
        with pytest.raises(ValueError, match="mutually exclusive"):
            DEM3D(
                pos=np.array([[0.0, 5.0, 5.0]]),
                vel=np.zeros((1, 3)),
                radii=np.array([2.0]),
                nx=20,
                ny=20,
                nz=20,
                particle_attraction=True,
                particle_repulsion=True,
            )

    def test_default_no_surface_force(self):
        # Both disabled by default: nearby spheres (gap within typical cutoff)
        # produce zero force when there is no contact.
        r = 3.0
        gap = -1.5
        dist = 2 * r - gap
        pos = np.array([[0.0, 10.0, 10.0], [dist, 10.0, 10.0]])
        sim = DEM3D(
            pos=pos,
            vel=np.zeros((2, 3)),
            radii=np.full(2, r),
            nx=40,
            ny=20,
            nz=20,
            gravity=0.0,
        )
        forces, _ = sim.compute_loads()
        np.testing.assert_allclose(forces, 0.0, atol=1e-12)

    def test_attraction_causes_aggregation(self):
        # Integration test: two non-overlapping spheres with attraction should
        # approach each other over time.
        r = 3.0
        gap = -1.5
        dist = 2 * r - gap
        pos = np.array([[0.0, 10.0, 10.0], [dist, 10.0, 10.0]])
        sim = DEM3D(
            pos=pos,
            vel=np.zeros((2, 3)),
            radii=np.full(2, r),
            nx=40,
            ny=20,
            nz=20,
            gravity=0.0,
            particle_attraction=True,
            attraction_strength=0.01,
            attraction_cutoff=3.0,
            dem_substeps=4,
        )
        d0 = np.linalg.norm(sim.pos[1] - sim.pos[0])
        sim.step(100)
        d1 = np.linalg.norm(sim.pos[1] - sim.pos[0])
        assert d1 < d0, "attracted spheres should move closer"


class TestCylinderAttraction:
    """Sphere-cylinder Hamaker attraction and repulsion."""

    def test_cylinder_attraction_pulls_sphere_inward(self):
        # Sphere just outside cylinder, within attraction cutoff.
        r_sphere = 2.0
        r_cyl = 4.0
        # Place sphere so surface gap = 1.0 (within cutoff=3.0)
        # radial dist from axis = r_sphere + r_cyl + 1.0 = 7.0
        pos = np.array([[7.0, 10.0, 10.0]])  # cylinder at (0,10)
        sim = DEM3D(
            pos=pos,
            vel=np.zeros((1, 3)),
            radii=np.array([r_sphere]),
            nx=30,
            ny=20,
            nz=20,
            gravity=0.0,
            cylinders=[(0.0, 10.0, r_cyl)],
            particle_attraction=True,
            attraction_strength=0.01,
            attraction_cutoff=3.0,
        )
        forces, _ = sim.compute_loads()
        # Sphere should be pulled toward cylinder axis (-x direction).
        assert forces[0, 0] < 0.0
        assert np.all(np.isfinite(forces))

    def test_cylinder_repulsion_pushes_sphere_outward(self):
        r_sphere = 2.0
        r_cyl = 4.0
        pos = np.array([[7.0, 10.0, 10.0]])
        sim = DEM3D(
            pos=pos,
            vel=np.zeros((1, 3)),
            radii=np.array([r_sphere]),
            nx=30,
            ny=20,
            nz=20,
            gravity=0.0,
            cylinders=[(0.0, 10.0, r_cyl)],
            particle_repulsion=True,
            repulsion_strength=0.01,
            repulsion_cutoff=3.0,
        )
        forces, _ = sim.compute_loads()
        # Sphere should be pushed away from cylinder axis (+x direction).
        assert forces[0, 0] > 0.0


class TestAttractionConfigWiring:
    """attraction/repulsion propagates through LBMDEMSolver3D and runner3d."""

    def test_lbm3d_wires_attraction_to_dem(self):
        from particulate_flow.lbm3d import LBMDEMSolver3D

        sim = LBMDEMSolver3D(
            nx=12,
            ny=8,
            nz=8,
            particle_source="left_inlet",
            particle_fluid_coupling="immersed_boundary",
            ibm_stiffness=0.3,
            pressure_drop=1e-3,
            streamwise_boundary="pressure",
            source_volume_fraction=0.1,
            particle_radius=2.0,
            particle_attraction=True,
            attraction_strength=0.005,
            attraction_cutoff=2.5,
            attraction_min_gap=0.1,
        )
        assert sim.dem.particle_attraction is True
        assert sim.dem.attraction_strength == pytest.approx(0.005)
        assert sim.dem.attraction_cutoff == pytest.approx(2.5)
        assert sim.dem.attraction_min_gap == pytest.approx(0.1)

    def test_lbm3d_raises_on_conflicting_flags(self):
        from particulate_flow.lbm3d import LBMDEMSolver3D

        with pytest.raises(ValueError, match="mutually exclusive"):
            LBMDEMSolver3D(
                nx=12,
                ny=8,
                nz=8,
                particle_source="left_inlet",
                particle_fluid_coupling="immersed_boundary",
                pressure_drop=1e-3,
                streamwise_boundary="pressure",
                source_volume_fraction=0.1,
                particle_radius=2.0,
                particle_attraction=True,
                particle_repulsion=True,
            )

    def test_runner3d_passes_attraction_from_args(self):
        import argparse
        from particulate_flow.runner3d import build_3d_solver

        args = argparse.Namespace(
            nx=12, ny=8, nz=8,
            reynolds_number=10.0, u_max=0.02,
            streamwise_boundary="pressure", pressure_drop=1e-3, rho_out=1.0,
            cylinder_spec=[],
            particle_source="left-inlet",
            particle_fluid_coupling="immersed_boundary",
            ibm_stiffness=0.3, ibm_marker_spacing=2.0,
            particle_radius=2.0, density_ratio=2.0, gravity=0.0, dem_substeps=4,
            sliding_friction=0.5, rolling_friction_coeff=0.05,
            particle_volume_fraction=0.1,
            particle_attraction=True,
            particle_repulsion=False,
            attraction_strength=0.007,
            repulsion_strength=1e-3,
            attraction_cutoff=3.0,
            repulsion_cutoff=3.0,
            attraction_min_gap=0.05,
            repulsion_min_gap=0.05,
        )
        solver = build_3d_solver(args)
        assert solver.dem.particle_attraction is True
        assert solver.dem.attraction_strength == pytest.approx(0.007)


# ---------------------------------------------------------------------------
# 2D regression guard
# ---------------------------------------------------------------------------


class TestNo2DRegression:
    def test_importing_contact3d_does_not_break_2d_solver(self):
        # Importing the 3D module must not perturb the 2D DEM API surface.
        from particulate_flow.dem.solver import DEMSolver

        assert hasattr(DEMSolver, "tangent_from_normal")
        assert hasattr(DEMSolver, "compute_loads")
