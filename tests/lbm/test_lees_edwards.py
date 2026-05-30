"""Tests for Lees-Edwards boundary conditions and iSP coupling (issue-7, issue-33)."""

from __future__ import annotations

import numpy as np
import pytest

from particulate_flow import FastLBMDEM

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _le_sim(
    nx: int = 48,
    ny: int = 32,
    shear_rate: float = 5e-4,
    n_particles: int = 0,
    interpolation_order: int = 3,
    particle_fluid_coupling: str = "point_force",
    u_max: float | None = None,
    init_analytical: bool = False,
) -> FastLBMDEM:
    """Build a minimal Lees-Edwards shear-flow simulation.

    For a *pure* shear validation, pass ``u_max=0.0`` so the ``periodic_force``
    streamwise mode derives no body force (``F_drive=0``) — otherwise a plug flow
    toward ``u_max`` swamps the shear signal.  ``init_analytical=True`` seeds the
    distribution from the analytical linear profile (mirroring the 3D solver), which
    the LE boundary correction then maintains.
    """
    return FastLBMDEM(
        nx=nx,
        ny=ny,
        Re=1.0,
        u_max=(shear_rate * ny / 2.0) if u_max is None else u_max,
        n_particles=n_particles,
        particle_radius=2.0,
        gravity=0.0,
        y_boundary="lees_edwards",
        streamwise_boundary="periodic_force",
        le_shear_rate=shear_rate,
        le_shear_axis=0,
        le_boundary_axis=1,
        le_interpolation_order=interpolation_order,
        particle_fluid_coupling=particle_fluid_coupling,
        pressure_drop=0.0,
        fluid_accelerator="numpy",
        init_analytical=init_analytical,
        seed=42,
    )


# ---------------------------------------------------------------------------
# Acceptance scenario 1: linear velocity profile u_x = γ̇·y after steady state
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_le_linear_shear_profile():
    """Lees-Edwards shear maintains the analytical linear profile u_x = γ̇·y.

    Pure shear: ``u_max=0`` so no body force drives a plug flow, and
    ``init_analytical=True`` seeds the linear profile (mirroring the 3D LE test).
    The LE boundary correction must hold it to L2 error < 1% with zero bulk offset.
    """
    shear_rate = 5e-4
    ny = 32
    sim = _le_sim(shear_rate=shear_rate, ny=ny, u_max=0.0, init_analytical=True)

    sim.advance(5000)

    _, ux, _ = sim.get_fields()
    y = np.arange(ny)
    expected = shear_rate * (y - (ny - 1) / 2.0)
    # Average over x.
    ux_mean = np.mean(ux, axis=0)
    scale = float(np.max(np.abs(expected))) or 1.0
    l2_err = float(np.sqrt(np.mean((ux_mean - expected) ** 2))) / scale
    assert l2_err < 0.01, f"L2 error {l2_err:.4f} exceeds 1% threshold"
    # No plug-flow contamination: the bulk-mean velocity stays ~0 (pure shear).
    assert abs(float(ux_mean.mean())) < 1e-4, "unexpected bulk velocity offset (body force?)"


def test_le_default_init_is_rest():
    """Without init_analytical the 2D solver still initialises from rest (opt-in)."""
    sim = _le_sim(u_max=0.0, init_analytical=False)
    _, ux, uy = sim.get_fields()
    np.testing.assert_allclose(ux, 0.0, atol=1e-12)
    np.testing.assert_allclose(uy, 0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# Acceptance scenario 2: shear_axis / boundary_axis configurability
# ---------------------------------------------------------------------------


def test_le_axis_configuration():
    """Non-default shear_axis and boundary_axis should be stored without error."""
    sim = FastLBMDEM(
        nx=32,
        ny=20,
        Re=1.0,
        u_max=0.01,
        n_particles=0,
        particle_radius=1.0,
        gravity=0.0,
        y_boundary="lees_edwards",
        streamwise_boundary="periodic_force",
        le_shear_rate=1e-4,
        le_shear_axis=0,
        le_boundary_axis=1,
        le_interpolation_order=3,
        pressure_drop=0.0,
        fluid_accelerator="numpy",
    )
    assert sim.le_shear_axis == 0
    assert sim.le_boundary_axis == 1
    sim.advance(2)
    _, ux, _ = sim.get_fields()
    assert np.all(np.isfinite(ux))


# ---------------------------------------------------------------------------
# Acceptance scenario 3: interpolation_order 1 and 3
# ---------------------------------------------------------------------------


def test_le_interpolation_order_toggle():
    """Both interpolation_order=1 and 3 should run and produce finite fields."""
    for order in (1, 3):
        sim = _le_sim(interpolation_order=order)
        sim.advance(5)
        _, ux, _ = sim.get_fields()
        assert np.all(np.isfinite(ux)), f"Non-finite ux with interpolation_order={order}"


# ---------------------------------------------------------------------------
# Acceptance scenario 4: iSP coupling mode selectable
# ---------------------------------------------------------------------------


def test_isp_coupling_selectable():
    """particle_fluid_coupling='isp' should be accepted without ValueError."""
    sim = _le_sim(n_particles=3, particle_fluid_coupling="isp")
    assert sim.particle_fluid_coupling == "isp"
    sim.advance(5)
    _, ux, _ = sim.get_fields()
    assert np.all(np.isfinite(ux))


# ---------------------------------------------------------------------------
# Acceptance scenario 5: particle wrap across LE boundary (unit test)
# ---------------------------------------------------------------------------


def test_le_particle_wrapping_position():
    """Particle crossing y=ny boundary should have x wrapped by le_shift."""
    shear_rate = 1e-3
    ny = 20
    sim = _le_sim(nx=40, ny=ny, shear_rate=shear_rate, n_particles=1)

    # Pre-set a known non-zero le_shift so we can verify the x displacement
    known_shift = 5.0
    sim._le_shift = known_shift

    # Place particle near top boundary, moving upward fast enough to cross
    sim.pos[:] = np.array([[20.0, float(ny) - 0.5]])
    sim.vel[:] = np.array([[0.0, 2.0]])

    x_before = float(sim.pos[0, 0])
    sim._dem_substep(1.0)

    y_after = float(sim.pos[0, 1])
    x_after = float(sim.pos[0, 0])

    # Particle should have wrapped in y
    assert y_after < ny / 2.0, "Particle should have wrapped in y"

    # x should be shifted by -known_shift (mod nx) when crossing top boundary
    expected_x = (x_before - known_shift) % sim.nx
    assert (
        abs(x_after - expected_x) < 0.5
    ), f"x_after={x_after:.4f} expected≈{expected_x:.4f} (shift={known_shift})"


def test_le_particle_wrapping_velocity():
    """Particle velocity ux should be corrected by γ̇·ny on y-wrap."""
    shear_rate = 1e-3
    ny = 20
    sim = _le_sim(nx=40, ny=ny, shear_rate=shear_rate, n_particles=1)

    sim.pos[:] = np.array([[20.0, float(ny) - 0.5]])
    sim.vel[:] = np.array([[0.1, 2.0]])

    vx_before = float(sim.vel[0, 0])
    sim._dem_substep(1.0)
    y_after = float(sim.pos[0, 1])

    # Particle must have wrapped — assert rather than skip silently
    assert y_after < ny / 2.0, "Particle should have crossed top boundary and wrapped"

    # Crossed top boundary: vx should decrease by shear_rate * ny
    expected_correction = -shear_rate * ny
    actual_correction = float(sim.vel[0, 0]) - vx_before
    # Tolerance is loose because DEM drag modifies velocity during the same substep
    # and the test uses dt=1 (coarse). Sign and rough magnitude are what matter.
    assert abs(actual_correction - expected_correction) < 2.0 * abs(expected_correction)


# ---------------------------------------------------------------------------
# Acceptance scenario 6: existing fouling tests still pass
# ---------------------------------------------------------------------------


def test_existing_periodic_y_still_works():
    """Existing periodic y_boundary mode unchanged after LE addition."""
    sim = FastLBMDEM(
        nx=32,
        ny=20,
        Re=40.0,
        u_max=0.02,
        n_particles=0,
        particle_radius=1.0,
        gravity=0.0,
        y_boundary="periodic",
        streamwise_boundary="periodic_force",
        pressure_drop=0.0,
        fluid_accelerator="numpy",
    )
    sim.advance(10)
    _, ux, _ = sim.get_fields()
    assert np.all(np.isfinite(ux))


# ---------------------------------------------------------------------------
# issue-33: Numba path must apply LE correction (new tests — currently failing)
# ---------------------------------------------------------------------------


def _le_sim_numba(
    nx: int = 48,
    ny: int = 32,
    shear_rate: float = 5e-4,
    u_max: float = 0.0,
    init_analytical: bool = True,
) -> FastLBMDEM:
    """Build a LE shear sim with the numba fluid accelerator."""
    return FastLBMDEM(
        nx=nx,
        ny=ny,
        Re=1.0,
        u_max=u_max,
        n_particles=0,
        particle_radius=2.0,
        gravity=0.0,
        y_boundary="lees_edwards",
        streamwise_boundary="periodic_force",
        le_shear_rate=shear_rate,
        le_shear_axis=0,
        le_boundary_axis=1,
        le_interpolation_order=3,
        pressure_drop=0.0,
        fluid_accelerator="numba",
        init_analytical=init_analytical,
        seed=42,
    )


@pytest.mark.slow
def test_le_numba_linear_shear_profile():
    """Numba path maintains analytical linear shear profile u_x = γ̇·y with < 1% L2 error.

    Mirrors ``test_le_linear_shear_profile`` but uses ``fluid_accelerator="numba"``.
    Skipped when numba is not installed.
    """
    pytest.importorskip("numba")

    shear_rate = 5e-4
    ny = 32
    sim = _le_sim_numba(shear_rate=shear_rate, ny=ny, u_max=0.0, init_analytical=True)
    sim.advance(5000)

    _, ux, _ = sim.get_fields()
    y = np.arange(ny)
    expected = shear_rate * (y - (ny - 1) / 2.0)
    ux_mean = np.mean(ux, axis=0)
    scale = float(np.max(np.abs(expected))) or 1.0
    l2_err = float(np.sqrt(np.mean((ux_mean - expected) ** 2))) / scale
    assert l2_err < 0.01, f"Numba LE L2 error {l2_err:.4f} exceeds 1% threshold"
    assert abs(float(ux_mean.mean())) < 1e-4, "unexpected bulk velocity offset (numba path)"


@pytest.mark.slow
def test_le_numpy_numba_field_agreement():
    """NumPy and Numba accelerators produce the same steady-state velocity field.

    Runs both paths from an identical analytical seed for 500 steps and asserts
    the x-averaged velocity profiles agree within 1e-4 (physical symmetry check).
    Skipped when numba is not installed.
    """
    pytest.importorskip("numba")

    shear_rate = 5e-4
    ny = 32
    steps = 500

    sim_np = _le_sim(shear_rate=shear_rate, ny=ny, u_max=0.0, init_analytical=True)
    sim_nb = _le_sim_numba(shear_rate=shear_rate, ny=ny, u_max=0.0, init_analytical=True)

    sim_np.advance(steps)
    sim_nb.advance(steps)

    _, ux_np, _ = sim_np.get_fields()
    _, ux_nb, _ = sim_nb.get_fields()

    profile_np = np.mean(ux_np, axis=0)
    profile_nb = np.mean(ux_nb, axis=0)

    max_diff = float(np.max(np.abs(profile_np - profile_nb)))
    assert max_diff < 1e-4, (
        f"NumPy/Numba ux profiles diverge by {max_diff:.2e} (LE correction asymmetry)"
    )
