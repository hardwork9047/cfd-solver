"""Tests for the coupled LBM-DEM solver."""

import numpy as np

from cfd_dem_lbm import LBMDEMSolver


def _two_particle_solver(
    *,
    particle_attraction: bool = False,
    particle_repulsion: bool = False,
) -> LBMDEMSolver:
    sim = LBMDEMSolver(
        nx=50,
        ny=30,
        Re=100.0,
        u_max=0.0,
        n_particles=2,
        particle_radius=2.0,
        density_ratio=2.0,
        gravity=0.0,
        particle_attraction=particle_attraction,
        particle_repulsion=particle_repulsion,
        attraction_strength=6e-3,
        repulsion_strength=6e-3,
        attraction_cutoff=3.0,
        repulsion_cutoff=3.0,
        attraction_min_gap=0.05,
        repulsion_min_gap=0.05,
        seed=1,
    )
    sim.pos[:] = np.array([[20.0, 15.0], [25.0, 15.0]])
    sim.vel[:] = 0.0
    return sim


def test_particle_attraction_is_disabled_by_default_for_separated_particles():
    """Separated particles exert no DEM force when attraction is disabled."""
    sim = _two_particle_solver(particle_attraction=False)

    forces = sim._dem_forces(dt_sub=1.0)

    np.testing.assert_allclose(forces, 0.0, atol=1e-12)


def test_particle_attraction_pulls_separated_particles_together():
    """Hamaker-like attraction pulls separated particles toward each other."""
    sim = _two_particle_solver(particle_attraction=True)

    forces = sim._dem_forces(dt_sub=1.0)

    assert forces[0, 0] > 0.0
    assert forces[1, 0] < 0.0
    np.testing.assert_allclose(forces[:, 1], 0.0, atol=1e-12)
    np.testing.assert_allclose(forces[0], -forces[1], rtol=1e-12, atol=1e-12)


def test_particle_repulsion_pushes_separated_particles_apart():
    """Hamaker-like repulsion pushes separated particles away from each other."""
    sim = _two_particle_solver(particle_repulsion=True)

    forces = sim._dem_forces(dt_sub=1.0)

    assert forces[0, 0] < 0.0
    assert forces[1, 0] > 0.0
    np.testing.assert_allclose(forces[:, 1], 0.0, atol=1e-12)
    np.testing.assert_allclose(forces[0], -forces[1], rtol=1e-12, atol=1e-12)


def test_particle_attraction_and_repulsion_are_mutually_exclusive():
    """The two non-contact surface-force modes cannot be enabled together."""
    try:
        _two_particle_solver(particle_attraction=True, particle_repulsion=True)
    except ValueError as exc:
        assert "mutually exclusive" in str(exc)
    else:
        raise AssertionError("Expected mutually exclusive surface-force modes to fail")


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


def test_left_inlet_mode_places_particles_upstream_and_deletes_right_outflow():
    """Left-inlet mode feeds particles upstream and removes right-outflow particles."""
    sim = LBMDEMSolver(
        nx=60,
        ny=30,
        Re=100.0,
        u_max=0.05,
        n_particles=6,
        particle_radius=2.0,
        density_ratio=2.0,
        gravity=0.0,
        cylinder=(25.0, 15.0, 4.0),
        particle_source="left_inlet",
        seed=3,
    )

    assert sim.n_p == 6
    assert np.all(sim.pos[:, 0] < 25.0 - 4.0)
    assert np.ptp(sim.pos[:, 1]) > 0.0

    sim.pos[0, 0] = sim.nx + sim.radii[0] + 0.1
    sim._delete_right_outflow_particles()

    assert sim.n_p == 5
    assert sim.removed_particles == 1
