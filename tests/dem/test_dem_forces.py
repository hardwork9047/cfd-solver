"""Tests for the coupled LBM-DEM solver."""

import numpy as np
import pytest

from particulate_flow import DEMSolver, FastLBMDEM, LBMDEMSolver

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


