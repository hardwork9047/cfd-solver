"""Analytical benchmark and grid-convergence tests (issue #27).

Verifies that the solvers not only run but converge to the correct physics:

- 2D body-force Poiseuille profile: steady-state is an exact parabola (L2 < 0.1%)
- Profile symmetry and positive centreline velocity
- Grid convergence: the parabola fit quality improves (remains perfect) as grid refines
- Stokes drag unit test: _stokes_drag implements F = 3π μ d Δu exactly
- 3D pressure-driven uniform flow: bulk ux positive, cross-sectional variation < 10%

Note on analytical formula: the LBM body-force Poiseuille steady state is an exact
parabola, but the amplitude depends on the bounce-back wall offset (mid-link BC places
effective no-slip ~0.05 nodes inside the wall node due to the Guo forcing structure).
Tests therefore verify SHAPE (parabola fit residual) and SIGN rather than the exact
amplitude, which would require accounting for boundary-condition discretisation offsets.
"""

from __future__ import annotations

import numpy as np
import pytest

from particulate_flow import FastLBMDEM
from particulate_flow.lbm3d import LBMDEMSolver3D

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _poiseuille_2d(
    ny: int = 32,
    u_max: float = 0.01,
    Re: float = 1.0,
    nx: int = 4,
) -> FastLBMDEM:
    """Build a 2D body-force Poiseuille channel (wall BCs, periodic_force)."""
    return FastLBMDEM(
        nx=nx,
        ny=ny,
        Re=Re,
        u_max=u_max,
        n_particles=0,
        particle_radius=1.0,
        gravity=0.0,
        y_boundary="wall",
        streamwise_boundary="periodic_force",
        pressure_drop=0.0,
        fluid_accelerator="numpy",
        seed=0,
    )


def _parabola_fit_l2(y: np.ndarray, u: np.ndarray) -> float:
    """Relative L2 residual after fitting a degree-2 polynomial to (y, u)."""
    p = np.polyfit(y, u, 2)
    u_fit = np.polyval(p, y)
    scale = float(np.max(np.abs(u))) or 1.0
    return float(np.sqrt(np.mean((u - u_fit) ** 2))) / scale


# ---------------------------------------------------------------------------
# Acceptance scenario 1a: 2D Poiseuille profile is an exact parabola
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_poiseuille_2d_profile_is_parabolic():
    """2D body-force Poiseuille steady state is an exact parabolic profile.

    Fits a degree-2 polynomial to the x-averaged velocity field at the fluid
    nodes and checks the relative L2 residual is < 0.1%.  Also verifies:
    - centreline velocity is positive (body force drives flow)
    - profile is symmetric (L/R halves differ by < 0.1%)
    """
    ny = 32
    sim = _poiseuille_2d(ny=ny)
    sim.advance(8000)

    _, ux, _ = sim.get_fields()
    ux_mean = np.mean(ux, axis=0)

    y = np.arange(ny, dtype=float)
    fluid = (y >= 1) & (y <= ny - 2)
    yf = y[fluid]
    uf = ux_mean[fluid]

    # Shape: must be an exact parabola
    l2_fit = _parabola_fit_l2(yf, uf)
    assert (
        l2_fit < 0.001
    ), f"Poiseuille profile is not parabolic: fit residual L2 = {l2_fit:.4f} > 0.1%"

    # Sign: centreline velocity must be positive (body force in +x)
    assert uf.max() > 0.0, "Centreline velocity must be positive (body force drives +x flow)"

    # Symmetry: left and right halves of the profile should mirror
    mid = len(uf) // 2
    left = uf[:mid]
    right = uf[len(uf) - mid :][::-1]
    scale = float(np.max(np.abs(uf))) or 1.0
    sym_err = float(np.max(np.abs(left - right))) / scale
    assert sym_err < 0.001, f"Profile asymmetry {sym_err:.4f} exceeds 0.1%"


@pytest.mark.slow
def test_poiseuille_2d_umax_matches_body_force():
    """2D Poiseuille peak velocity is consistent with the body-force formula.

    The solver sets F_drive = 8ν·u_max/ny², which predicts a peak velocity of
    u_max.  The actual LBM peak differs from this target due to the bounce-back
    wall offset (effectively ~1 node at each wall).  At ny=32, the measured
    offset is ~7%; we accept ±12% to accommodate the documented offset plus
    floating-point noise while catching large implementation errors.
    """
    ny = 32
    u_max_target = 0.01
    sim = _poiseuille_2d(ny=ny, u_max=u_max_target)
    sim.advance(8000)

    _, ux, _ = sim.get_fields()
    ux_mean = np.mean(ux, axis=0)

    y = np.arange(ny, dtype=float)
    fluid = (y >= 1) & (y <= ny - 2)
    u_max_sim = float(ux_mean[fluid].max())

    rel_err = abs(u_max_sim - u_max_target) / u_max_target
    assert rel_err < 0.12, (
        f"u_max_sim={u_max_sim:.5f} deviates from target {u_max_target:.5f} "
        f"by {rel_err*100:.1f}% (> 12%)"
    )


# ---------------------------------------------------------------------------
# Acceptance scenario 2: Stokes drag unit test
# ---------------------------------------------------------------------------


def test_stokes_drag_unit():
    """_stokes_drag computes F = 3π μ d (u_fluid - v_particle) exactly."""
    sim = _poiseuille_2d(ny=20, u_max=0.01, Re=1.0)

    radius = 2.0
    nu = sim.nu
    uf_x, uf_y = 0.05, 0.0
    vp_x, vp_y = 0.01, 0.02

    fx, fy = sim._stokes_drag(uf_x, uf_y, vp_x, vp_y, radius=radius)

    coeff = 3.0 * np.pi * nu * 2.0 * radius
    expected_fx = coeff * (uf_x - vp_x)
    expected_fy = coeff * (uf_y - vp_y)

    assert abs(fx - expected_fx) < 1e-14, f"fx={fx:.6e} expected={expected_fx:.6e}"
    assert abs(fy - expected_fy) < 1e-14, f"fy={fy:.6e} expected={expected_fy:.6e}"


def test_stokes_drag_zero_slip():
    """Zero slip velocity produces zero Stokes drag."""
    sim = _poiseuille_2d(ny=20, u_max=0.01, Re=1.0)
    fx, fy = sim._stokes_drag(0.05, 0.03, 0.05, 0.03, radius=1.5)
    assert fx == 0.0
    assert fy == 0.0


def test_stokes_drag_scales_with_radius():
    """Stokes drag scales linearly with particle radius (diameter)."""
    sim = _poiseuille_2d(ny=20, u_max=0.01, Re=1.0)
    r1, r2 = 1.0, 2.0
    fx1, _ = sim._stokes_drag(0.05, 0.0, 0.0, 0.0, radius=r1)
    fx2, _ = sim._stokes_drag(0.05, 0.0, 0.0, 0.0, radius=r2)
    assert (
        abs(fx2 / fx1 - (r2 / r1)) < 1e-12
    ), f"Expected linear scaling by r, got ratio {fx2 / fx1:.6f}"


# ---------------------------------------------------------------------------
# Acceptance scenario 3: grid convergence
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_poiseuille_2d_grid_convergence():
    """2D Poiseuille parabola fit quality and u_max convergence across grid sizes.

    Tests ny = 16, 24, 32.  The parabola fit residual should remain below 0.1%
    at all resolutions — the LBM steady state is analytically exact at all grid
    sizes.  Also verifies that the u_max amplitude error decreases monotonically
    as the grid is refined (the bounce-back wall offset shrinks relative to ny,
    producing first-order convergence in the amplitude toward the target u_max).
    """
    grids = [16, 24, 32]
    u_max_target = 0.01
    fit_residuals = []
    u_max_errors = []

    for ny in grids:
        sim = _poiseuille_2d(ny=ny, u_max=u_max_target)
        steps = max(4000, 10 * ny**2)
        sim.advance(steps)

        _, ux, _ = sim.get_fields()
        ux_mean = np.mean(ux, axis=0)

        y = np.arange(ny, dtype=float)
        fluid = (y >= 1) & (y <= ny - 2)
        uf = ux_mean[fluid]
        yf = y[fluid]

        l2_fit = _parabola_fit_l2(yf, uf)
        fit_residuals.append(l2_fit)

        u_max_err = abs(uf.max() - u_max_target) / u_max_target
        u_max_errors.append(u_max_err)

    for ny, residual in zip(grids, fit_residuals):
        assert residual < 0.001, f"ny={ny}: parabola fit residual {residual:.4f} > 0.1%"

    # u_max convergence: finer grids must have strictly smaller amplitude error.
    # Measured: ny=16 ~14%, ny=24 ~9%, ny=32 ~7% (monotone decrease).
    assert u_max_errors[1] < u_max_errors[0], (
        f"u_max error should decrease from ny=16 to ny=24: "
        f"{u_max_errors[0]*100:.1f}% → {u_max_errors[1]*100:.1f}%"
    )
    assert u_max_errors[2] < u_max_errors[1], (
        f"u_max error should decrease from ny=24 to ny=32: "
        f"{u_max_errors[1]*100:.1f}% → {u_max_errors[2]*100:.1f}%"
    )


# ---------------------------------------------------------------------------
# Acceptance scenario 1b: 3D uniform pressure-driven flow
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_3d_uniform_flow_profile():
    """3D pressure-driven flow develops a uniform streamwise profile at low Re.

    With periodic y and z boundaries and a pressure inlet/outlet, the 3D solver
    should produce a spatially uniform ux field at steady state (no walls →
    no shear, just plug flow). The bulk velocity should be positive and the
    cross-channel variation should be small (rel. std < 10%).
    """
    nx, ny, nz = 24, 16, 8
    pressure_drop = 1e-3

    sim = LBMDEMSolver3D(
        nx=nx,
        ny=ny,
        nz=nz,
        Re=5.0,
        u_max=0.02,
        le_shear_rate=0.0,
        streamwise_boundary="pressure",
        pressure_drop=pressure_drop,
        rho_out=1.0,
    )
    sim.advance(2000)

    _, ux, _, _ = sim.get_fields()
    bulk_ux = float(np.mean(ux[1:-1]))  # exclude inlet/outlet planes
    assert bulk_ux > 0.0, f"Expected positive bulk ux, got {bulk_ux:.4e}"

    # Cross-sectional uniformity: std/mean < 10% (no walls → no parabolic shear)
    interior = ux[1:-1]
    rel_std = float(np.std(interior)) / (abs(bulk_ux) + 1e-12)
    assert rel_std < 0.10, f"3D ux field non-uniform: std/mean = {rel_std:.3f} > 10%"
