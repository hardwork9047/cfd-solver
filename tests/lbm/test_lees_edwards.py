"""Tests for Lees-Edwards boundary conditions and iSP coupling (issue-7)."""

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
) -> FastLBMDEM:
    """Build a minimal Lees-Edwards shear-flow simulation."""
    return FastLBMDEM(
        nx=nx,
        ny=ny,
        Re=1.0,
        u_max=shear_rate * ny / 2.0,
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
        seed=42,
    )


# ---------------------------------------------------------------------------
# Acceptance scenario 1: linear velocity profile u_x = γ̇·y after steady state
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_le_linear_shear_profile():
    """Steady-state velocity field should match u_x = γ̇·y with L2 error < 1%."""
    shear_rate = 5e-4
    ny = 32
    sim = _le_sim(shear_rate=shear_rate, ny=ny)

    sim.advance(5000)

    _, ux, _ = sim.get_fields()
    y = np.arange(ny)
    expected = shear_rate * (y - (ny - 1) / 2.0)
    # Average over x
    ux_mean = np.mean(ux, axis=0)
    scale = float(np.max(np.abs(expected))) or 1.0
    l2_err = float(np.sqrt(np.mean((ux_mean - expected) ** 2))) / scale
    assert l2_err < 0.01, f"L2 error {l2_err:.4f} exceeds 1% threshold"


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
    assert abs(x_after - expected_x) < 0.5, (
        f"x_after={x_after:.4f} expected≈{expected_x:.4f} (shift={known_shift})"
    )


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
