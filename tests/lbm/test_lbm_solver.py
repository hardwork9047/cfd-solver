"""Tests for the coupled LBM-DEM solver."""

import numpy as np
import pytest

from particulate_flow import DEMSolver, FastLBMDEM, LBMDEMSolver

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


def test_periodic_y_pressure_boundary_advances_with_open_top_bottom():
    """Pressure inlet/outlet with y-periodicity should run without top/bottom walls."""
    sim = FastLBMDEM(
        nx=48,
        ny=20,
        Re=40.0,
        u_max=0.02,
        n_particles=0,
        particle_radius=1.0,
        gravity=0.0,
        y_boundary="periodic",
        streamwise_boundary="pressure",
        pressure_drop=1e-5,
        fluid_accelerator="numpy",
    )

    assert not sim.fixed_solid[:, 0].any()
    assert not sim.fixed_solid[:, -1].any()
    assert np.isclose(sim.F_drive, 0.0)

    sim.advance(20)
    rho, ux, uy = sim.get_fields()

    assert np.all(np.isfinite(rho))
    assert np.all(np.isfinite(ux))
    assert np.all(np.isfinite(uy))
    assert np.mean(ux[1:-1, :]) > 0.0
    assert np.isclose(np.mean(rho[0, :]), sim.rho_in, rtol=5e-4)
    assert np.isclose(np.mean(rho[-1, :]), sim.rho_out, rtol=5e-4)


def test_periodic_y_particles_wrap_without_wall_restitution():
    """DEM particles should wrap in y when the transverse boundary is periodic."""
    sim = FastLBMDEM(
        nx=40,
        ny=20,
        Re=40.0,
        u_max=0.0,
        n_particles=1,
        particle_radius=1.0,
        gravity=0.0,
        y_boundary="periodic",
        streamwise_boundary="pressure",
        pressure_drop=0.0,
        fluid_accelerator="numpy",
    )
    sim.pos[:] = np.array([[10.0, 19.8]])
    sim.vel[:] = np.array([[0.0, 2.0]])
    sim._dem_substep(0.5)

    assert 0.0 <= sim.pos[0, 1] < sim.ny


def test_periodic_y_particles_wrap_from_bottom_to_top():
    """A particle crossing below y=0 should re-enter from the top side."""
    sim = FastLBMDEM(
        nx=40,
        ny=20,
        Re=40.0,
        u_max=0.0,
        n_particles=1,
        particle_radius=1.0,
        gravity=0.0,
        y_boundary="periodic",
        streamwise_boundary="pressure",
        pressure_drop=0.0,
        fluid_accelerator="numpy",
    )
    sim.pos[:] = np.array([[10.0, 0.2]])
    sim.vel[:] = np.array([[0.0, -2.0]])
    sim._dem_substep(0.5)

    assert 19.0 < sim.pos[0, 1] < sim.ny
    assert sim.vel[0, 1] < 0.0


def test_periodic_y_particle_fluid_interpolation_and_force_spreading_wrap():
    """Particle-fluid interpolation and back-reaction should use periodic y."""
    sim = FastLBMDEM(
        nx=12,
        ny=20,
        Re=40.0,
        u_max=0.0,
        n_particles=0,
        particle_radius=1.0,
        gravity=0.0,
        y_boundary="periodic",
        streamwise_boundary="pressure",
        pressure_drop=0.0,
        fluid_accelerator="numpy",
    )
    ux = np.zeros((sim.nx, sim.ny))
    uy = np.zeros((sim.nx, sim.ny))
    ux[:, -1] = 20.0
    ux[:, 0] = 10.0

    uf_x, uf_y = sim._interp_velocity_many(np.array([5.0]), np.array([-0.25]), ux=ux, uy=uy)

    assert np.isclose(uf_x[0], 12.5)
    assert np.isclose(uf_y[0], 0.0)

    sim.Fx[:] = 0.0
    sim.Fy[:] = 0.0
    sim._distribute_force(5.0, -0.25, 0.0, 1.0, radius=1.0)

    assert sim.Fy[5, -1] > 0.0
    assert sim.Fy[5, 0] > 0.0


def test_periodic_y_ibm_markers_wrap_across_transverse_boundary():
    """IBM marker points should wrap, not clip, when y is periodic."""
    sim = LBMDEMSolver(
        nx=40,
        ny=20,
        Re=40.0,
        u_max=0.0,
        n_particles=1,
        particle_radius=1.0,
        gravity=0.0,
        y_boundary="periodic",
        streamwise_boundary="pressure",
        pressure_drop=0.0,
        particle_fluid_coupling="immersed_boundary",
        compute_accelerator="numpy",
        seed=11,
    )
    sim.pos[:] = np.array([[10.0, 0.2]])

    markers = sim.ibm_marker_points()

    assert np.all((0.0 <= markers["y"]) & (markers["y"] < sim.ny))
    assert np.any(markers["y"] > sim.ny - 1.0)


@pytest.mark.parametrize("fluid_method", ["lbm-bgk-guo", "lbm-trt-guo"])
@pytest.mark.parametrize("particle_method", ["dem-hertz", "dem-linear"])
@pytest.mark.parametrize("coupling", ["point_force", "solid_boundary", "immersed_boundary"])
def test_periodic_pressure_pore_mode_runs_for_solver_method_combinations(
    fluid_method,
    particle_method,
    coupling,
):
    """Membrane-pore boundary mode should remain finite for all solver options."""
    sim = FastLBMDEM(
        nx=36,
        ny=24,
        Re=20.0,
        u_max=0.01,
        n_particles=1,
        particle_radius=1.0,
        gravity=0.0,
        y_boundary="periodic",
        streamwise_boundary="pressure",
        pressure_drop=1e-6,
        fluid_method=fluid_method,
        particle_method=particle_method,
        particle_fluid_coupling=coupling,
        cylinders=[(22.0, 12.0, 2.5)],
        fluid_accelerator="numpy",
        compute_accelerator="numpy",
        seed=14,
    )
    sim.pos[:] = np.array([[8.0, 0.2]])
    sim.vel[:] = np.array([[0.0, -0.4]])

    sim.advance(3)
    rho, ux, uy = sim.get_fields()

    assert 0.0 <= sim.pos[0, 1] < sim.ny
    assert np.all(np.isfinite(rho))
    assert np.all(np.isfinite(ux))
    assert np.all(np.isfinite(uy))
    assert np.all(np.isfinite(sim.pos))


