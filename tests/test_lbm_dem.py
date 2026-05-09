"""Tests for the coupled LBM-DEM solver."""

import numpy as np

from cfd_dem_lbm import DEMSolver, FastLBMDEM, LBMDEMSolver


def _two_particle_solver(
    *,
    particle_attraction: bool = False,
    particle_repulsion: bool = False,
    particle_method: str = "dem-hertz",
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
        particle_method=particle_method,
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
    assert isinstance(sim.dem_solver, DEMSolver)


def test_dem_solver_directly_matches_coupled_dem_force_api():
    """Verification code can use the same DEM solver as production coupling."""
    sim = _two_particle_solver(particle_attraction=True)

    direct_forces, direct_torques = sim.dem_solver.compute_loads(dt_sub=1.0)
    api_forces = sim._dem_forces(dt_sub=1.0)

    np.testing.assert_allclose(direct_forces, api_forces, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(direct_torques, sim.torques_p, rtol=1e-12, atol=1e-12)


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


def _particle_cylinder_solver(
    *,
    particle_attraction: bool = False,
    particle_repulsion: bool = False,
    compute_accelerator: str = "numpy",
) -> LBMDEMSolver:
    """Return one particle just outside a fixed cylinder surface."""
    sim = LBMDEMSolver(
        nx=60,
        ny=32,
        Re=100.0,
        u_max=0.0,
        n_particles=1,
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
        cylinders=[(24.0, 16.0, 5.0)],
        compute_accelerator=compute_accelerator,
        seed=3,
    )
    cx, cy, cr = sim.cylinders[0]
    sim.pos[:] = np.array([[cx + cr + sim.radii[0] + 1.0, cy]])
    sim.vel[:] = 0.0
    return sim


def test_particle_cylinder_attraction_pulls_particle_toward_cylinder():
    """The Hamaker-like attraction also acts between particles and cylinders."""
    sim = _particle_cylinder_solver(particle_attraction=True)

    forces = sim._dem_forces(dt_sub=1.0)

    assert forces[0, 0] < 0.0
    np.testing.assert_allclose(forces[0, 1], 0.0, atol=1e-12)


def test_particle_cylinder_repulsion_pushes_particle_away_from_cylinder():
    """The Hamaker-like repulsion also acts between particles and cylinders."""
    sim = _particle_cylinder_solver(particle_repulsion=True)

    forces = sim._dem_forces(dt_sub=1.0)

    assert forces[0, 0] > 0.0
    np.testing.assert_allclose(forces[0, 1], 0.0, atol=1e-12)


def test_particle_cylinder_surface_force_matches_numba_backend():
    """Compiled and NumPy boundary kernels should give the same cylinder force."""
    numpy_sim = _particle_cylinder_solver(particle_attraction=True, compute_accelerator="numpy")
    numba_sim = _particle_cylinder_solver(particle_attraction=True, compute_accelerator="auto")

    numpy_forces = numpy_sim._dem_forces(dt_sub=1.0)
    numba_forces = numba_sim._dem_forces(dt_sub=1.0)

    np.testing.assert_allclose(numba_forces, numpy_forces, rtol=1e-12, atol=1e-12)


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


def test_invalid_solver_methods_are_rejected():
    """Fluid and DEM method options should fail early when misspelled."""
    try:
        LBMDEMSolver(fluid_method="not-a-fluid-method")
    except ValueError as exc:
        assert "fluid_method" in str(exc)
    else:
        raise AssertionError("Expected invalid fluid method to fail")

    try:
        LBMDEMSolver(particle_method="not-a-particle-method")
    except ValueError as exc:
        assert "particle_method" in str(exc)
    else:
        raise AssertionError("Expected invalid particle method to fail")

    try:
        LBMDEMSolver(particle_fluid_coupling="not-a-coupling")
    except ValueError as exc:
        assert "particle_fluid_coupling" in str(exc)
    else:
        raise AssertionError("Expected invalid particle-fluid coupling to fail")

    try:
        LBMDEMSolver(fluid_accelerator="not-an-accelerator")
    except ValueError as exc:
        assert "fluid_accelerator" in str(exc)
    else:
        raise AssertionError("Expected invalid fluid accelerator to fail")

    try:
        LBMDEMSolver(compute_accelerator="not-an-accelerator")
    except ValueError as exc:
        assert "compute_accelerator" in str(exc)
    else:
        raise AssertionError("Expected invalid compute accelerator to fail")

    try:
        LBMDEMSolver(particle_search="not-a-search")
    except ValueError as exc:
        assert "particle_search" in str(exc)
    else:
        raise AssertionError("Expected invalid particle search to fail")

    try:
        LBMDEMSolver(spatial_dimension=3)
    except ValueError as exc:
        assert "spatial_dimension" in str(exc)
    else:
        raise AssertionError("Expected unsupported spatial dimension to fail")

    try:
        LBMDEMSolver(particle_fluid_coupling="immersed_boundary", ibm_stiffness=0.0)
    except ValueError as exc:
        assert "ibm_stiffness" in str(exc)
    else:
        raise AssertionError("Expected invalid IBM stiffness to fail")


def test_particle_search_methods_are_selectable():
    """Cell-list search can be replaced by exhaustive pairs for debugging."""
    cell = _two_particle_solver(particle_attraction=True)
    all_pairs = _two_particle_solver(particle_attraction=True)
    all_pairs.particle_search = "all_pairs"

    assert list(cell._particle_pair_candidates()) == [(0, 1)]
    assert list(all_pairs._particle_pair_candidates()) == [(0, 1)]

    np.testing.assert_allclose(
        cell._dem_forces(dt_sub=1.0),
        all_pairs._dem_forces(dt_sub=1.0),
        rtol=1e-12,
        atol=1e-12,
    )


def test_linear_dem_contact_model_changes_normal_contact_force():
    """Linear DEM contact is selectable and differs from Hertz for small overlaps."""
    sim_hertz = _two_particle_solver()
    sim_linear = _two_particle_solver(particle_method="dem-linear")

    sim_hertz.pos[:] = np.array([[20.0, 15.0], [23.5, 15.0]])
    sim_linear.pos[:] = sim_hertz.pos.copy()

    hertz_forces = sim_hertz._dem_forces(dt_sub=1.0)
    linear_forces = sim_linear._dem_forces(dt_sub=1.0)

    assert linear_forces[1, 0] > hertz_forces[1, 0]
    assert sim_linear.dem_solver.pair_kernel is None


def test_trt_lbm_method_advances_with_finite_fields():
    """TRT LBM is selectable for the same coupled solver path."""
    sim = LBMDEMSolver(
        nx=32,
        ny=20,
        Re=50.0,
        u_max=0.03,
        n_particles=1,
        particle_radius=2.0,
        gravity=0.0,
        fluid_method="lbm-trt-guo",
        seed=8,
    )
    sim._delete_particles(np.ones(sim.n_p, dtype=bool))

    sim.advance(5)
    fields = sim.get_fields()

    assert sim.fluid_method == "lbm-trt-guo"
    for field in fields:
        assert np.all(np.isfinite(field))


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
