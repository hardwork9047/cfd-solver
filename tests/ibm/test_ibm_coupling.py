"""Tests for the coupled LBM-DEM solver."""

import numpy as np
import pytest

from particulate_flow import DEMSolver, FastLBMDEM, LBMDEMSolver

def test_solid_boundary_coupling_overlays_particles_on_solid_mask():
    """Solid-boundary coupling should expose moving particles to the LBM mask."""
    sim = LBMDEMSolver(
        nx=40,
        ny=24,
        Re=50.0,
        u_max=0.02,
        n_particles=1,
        particle_radius=3.0,
        gravity=0.0,
        particle_fluid_coupling="solid_boundary",
        seed=9,
    )
    sim.pos[:] = np.array([[20.0, 12.0]])
    sim.vel[:] = 0.0
    sim._update_particle_solid_mask()

    assert sim.particle_fluid_coupling == "solid_boundary"
    assert sim.solid[20, 12]
    assert not sim.fixed_solid[20, 12]
    assert sim.dynamic_solid_fraction() > 0.0


def test_immersed_boundary_coupling_applies_particle_reaction_force():
    """IBM coupling should add fluid reaction forces without solidifying particles."""
    sim = LBMDEMSolver(
        nx=40,
        ny=24,
        Re=50.0,
        u_max=0.02,
        n_particles=1,
        particle_radius=3.0,
        gravity=0.0,
        particle_fluid_coupling="immersed_boundary",
        ibm_stiffness=1.0,
        seed=10,
    )
    sim.pos[:] = np.array([[20.0, 12.0]])
    sim.vel[:] = np.array([[0.1, 0.0]])
    sim.omega_p[:] = 0.0
    _, ux, uy = sim.get_fields()
    sim.Fx[:] = sim.F_drive
    sim.Fy[:] = 0.0

    sim._apply_immersed_boundary_forces(ux, uy)

    assert not sim.solid[20, 12]
    assert sim.ibm_forces_p[0, 0] < 0.0
    assert np.any(np.abs(sim.Fx - sim.F_drive) > 0.0)


def test_immersed_boundary_marker_points_are_available_for_visualization():
    """IBM marker coordinates should expose the embedded boundary points."""
    sim = LBMDEMSolver(
        nx=40,
        ny=24,
        Re=50.0,
        u_max=0.02,
        n_particles=1,
        particle_radius=3.0,
        gravity=0.0,
        particle_fluid_coupling="immersed_boundary",
        ibm_marker_spacing=1.0,
        seed=10,
    )
    sim.pos[:] = np.array([[20.0, 12.0]])

    markers = sim.ibm_marker_points()

    assert len(markers["x"]) >= 8
    assert markers["x"].shape == markers["y"].shape == markers["owner"].shape
    assert set(markers["owner"].tolist()) == {0}
    np.testing.assert_allclose(
        np.hypot(markers["x"] - 20.0, markers["y"] - 12.0),
        3.0,
        atol=1e-12,
    )


def test_rolling_friction_resists_wall_slip_and_spins_particle():
    """Wall contact friction must oppose tangential slip and add a torque."""
    sim = LBMDEMSolver(
        nx=40,
        ny=20,
        Re=100.0,
        u_max=0.0,
        n_particles=1,
        particle_radius=2.0,
        density_ratio=2.0,
        gravity=0.0,
        rolling_friction=True,
        y_boundary="wall",
        streamwise_boundary="periodic_force",
        sliding_friction=0.5,
        tangential_damping=0.4,
        seed=2,
    )
    sim.pos[:] = np.array([[20.0, 2.2]])
    sim.vel[:] = np.array([[1.0, 0.0]])
    sim.omega_p[:] = 0.0

    forces, torques = sim._dem_loads(dt_sub=1.0)

    assert forces[0, 0] < 0.0
    assert torques[0] < 0.0


