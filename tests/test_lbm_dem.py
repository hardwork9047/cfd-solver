"""Tests for the coupled LBM-DEM solver."""

import numpy as np

from cfd_dem_lbm import FastLBMDEM, LBMDEMSolver


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


def test_particle_pair_candidates_use_interaction_cutoff():
    """The cell list must include separated particles inside the surface-force range."""
    sim = _two_particle_solver(particle_attraction=True)
    sim.pos[:] = np.array([[20.0, 15.0], [26.5, 15.0]])

    pairs = list(sim._particle_pair_candidates())

    assert pairs == [(0, 1)]


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


def test_multiple_cylinders_define_solid_nodes_and_particle_contacts():
    """Multiple fixed cylinders can represent pore-scale membrane geometry."""
    sim = LBMDEMSolver(
        nx=60,
        ny=30,
        Re=100.0,
        u_max=0.0,
        n_particles=1,
        particle_radius=2.0,
        density_ratio=2.0,
        gravity=0.0,
        cylinders=[(20.0, 15.0, 4.0), (35.0, 15.0, 3.0)],
        seed=4,
    )
    sim.pos[:] = np.array([[23.0, 15.0]])
    sim.vel[:] = 0.0

    forces = sim._dem_forces(dt_sub=1.0)

    assert sim.solid[20, 15]
    assert sim.solid[35, 15]
    assert forces[0, 0] > 0.0
    assert len(sim.cylinders) == 2


def test_left_inlet_mode_generates_particles_from_flux_budget_and_deletes_outflow():
    """Left-inlet mode starts empty, injects from budget, and removes right outflow."""
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
        source_volume_fraction=0.05,
        seed=3,
    )

    assert sim.n_p == 0
    assert len(sim._pending_radii) == 6

    sim.inlet_particle_area_budget = float(np.pi * sim._pending_radii[0] ** 2)
    sim._try_feed_left_inlet_particles(max_new=1)

    assert sim.n_p == 1
    assert sim.generated_particles == 1
    assert np.all(sim.pos[:, 0] < 25.0 - 4.0)

    sim.pos[0, 0] = sim.nx + sim.radii[0] + 0.1
    sim._delete_right_outflow_particles()

    assert sim.n_p == 0
    assert sim.removed_particles == 1


def test_reynolds_length_sets_viscosity_from_particle_diameter_scale():
    """The Reynolds-number length scale can be decoupled from channel height."""
    sim = LBMDEMSolver(
        nx=40,
        ny=30,
        Re=100.0,
        u_max=0.05,
        reynolds_length=4.0,
        n_particles=1,
        particle_radius=2.0,
        gravity=0.0,
        seed=5,
    )

    assert np.isclose(sim.reynolds_length, 4.0)
    assert np.isclose(sim.nu, 0.05 * 4.0 / 100.0)


def test_target_max_velocity_flow_control_relaxes_drive_force():
    """Target-max-speed control should reduce drive force when the flow is too fast."""
    sim = LBMDEMSolver(
        nx=40,
        ny=30,
        Re=100.0,
        u_max=0.05,
        reynolds_length=4.0,
        flow_control="target_max_velocity",
        flow_control_gain=0.5,
        n_particles=1,
        particle_radius=2.0,
        gravity=0.0,
        seed=6,
    )
    ux = np.full((sim.nx, sim.ny), 0.10)
    uy = np.zeros((sim.nx, sim.ny))
    old_drive = sim.F_drive

    sim._control_drive_force(ux, uy)

    assert sim.F_drive < old_drive


def test_fast_solver_matches_base_solver_fluid_fields_without_particles():
    """The production cached solver must preserve the base LBM fluid path."""
    kwargs = dict(
        nx=32,
        ny=20,
        Re=50.0,
        u_max=0.03,
        reynolds_length=20.0,
        flow_control="fixed_pressure",
        n_particles=1,
        particle_radius=2.0,
        gravity=0.0,
        seed=7,
    )
    base = LBMDEMSolver(**kwargs)
    fast = FastLBMDEM(**kwargs)
    base._delete_particles(np.ones(base.n_p, dtype=bool))
    fast._delete_particles(np.ones(fast.n_p, dtype=bool))

    base.advance(20)
    fast.advance(20)
    base_fields = base.get_fields()
    fast_fields = fast.get_fields()

    for base_field, fast_field in zip(base_fields, fast_fields):
        np.testing.assert_allclose(fast_field, base_field, rtol=1e-12, atol=1e-12)
