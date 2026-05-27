"""Tests for the coupled LBM-DEM solver."""

import numpy as np
import pytest

from particulate_flow import DEMSolver, FastLBMDEM, LBMDEMSolver

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
        streamwise_boundary="periodic_force",
        seed=3,
    )

    assert sim.n_p == 0
    assert len(sim._pending_radii) == 6

    sim.inlet_particle_area_budget = float(np.pi * sim._pending_radii[0] ** 2)
    sim._try_feed_left_inlet_particles(max_new=1)

    assert sim.n_p == 1
    assert sim.generated_particles == 1
    assert sim.cumulative_inlet_flow_area > 0.0
    assert sim.injected_particle_area > 0.0
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
        streamwise_boundary="periodic_force",
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


def test_target_max_velocity_flow_control_keeps_positive_drive_floor():
    """The controller must not erase the drive force during high-speed transients."""
    sim = LBMDEMSolver(
        nx=40,
        ny=30,
        Re=100.0,
        u_max=0.05,
        reynolds_length=4.0,
        flow_control="target_max_velocity",
        flow_control_gain=1.0,
        n_particles=1,
        particle_radius=2.0,
        gravity=0.0,
        seed=11,
    )
    ux = np.full((sim.nx, sim.ny), 10.0)
    uy = np.zeros((sim.nx, sim.ny))

    for _ in range(100):
        sim._control_drive_force(ux, uy)

    assert np.isclose(sim.F_drive, 0.05 * sim.initial_F_drive)


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


def test_numba_fluid_accelerator_matches_numpy_when_available():
    """The compiled LBM path should remain numerically close to the NumPy path."""
    kwargs = dict(
        nx=24,
        ny=16,
        Re=25.0,
        u_max=0.02,
        reynolds_length=6.0,
        flow_control="fixed_pressure",
        n_particles=1,
        particle_radius=1.5,
        gravity=0.0,
        cylinders=[(9.0, 8.0, 2.0)],
        seed=12,
    )
    numpy_sim = FastLBMDEM(**kwargs, fluid_accelerator="numpy")
    numba_sim = FastLBMDEM(**kwargs, fluid_accelerator="numba")
    if not numba_sim.uses_numba_lbm:
        return
    numpy_sim._delete_particles(np.ones(numpy_sim.n_p, dtype=bool))
    numba_sim._delete_particles(np.ones(numba_sim.n_p, dtype=bool))

    numpy_sim.advance(4)
    numba_sim.advance(4)

    for numpy_field, numba_field in zip(numpy_sim.get_fields(), numba_sim.get_fields()):
        np.testing.assert_allclose(numba_field, numpy_field, rtol=1e-12, atol=1e-12)


def test_numba_compute_accelerator_matches_numpy_when_available():
    """Compiled DEM boundary and solid-mask helpers should match NumPy helpers."""
    kwargs = dict(
        nx=44,
        ny=24,
        Re=30.0,
        u_max=0.02,
        n_particles=2,
        particle_radius=2.0,
        gravity=0.0,
        rolling_friction=True,
        cylinders=[(18.0, 12.0, 4.0)],
        particle_fluid_coupling="solid_boundary",
        seed=13,
    )
    numpy_sim = FastLBMDEM(**kwargs, compute_accelerator="numpy")
    numba_sim = FastLBMDEM(**kwargs, compute_accelerator="numba")
    if not numba_sim.uses_numba_compute:
        return
    positions = np.array([[21.0, 12.0], [30.0, 3.0]])
    velocities = np.array([[-0.2, 0.0], [0.1, -0.1]])
    for sim in (numpy_sim, numba_sim):
        sim.pos[:] = positions
        sim.vel[:] = velocities
        sim.omega_p[:] = np.array([0.3, -0.2])
        sim._update_particle_solid_mask()

    np.testing.assert_array_equal(numba_sim.particle_solid, numpy_sim.particle_solid)
    np.testing.assert_allclose(
        numba_sim._dem_forces(dt_sub=1.0),
        numpy_sim._dem_forces(dt_sub=1.0),
        rtol=1e-12,
        atol=1e-12,
    )
