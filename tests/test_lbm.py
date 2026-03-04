"""
Unit tests for the cfd_lbm (D2Q9 Lattice-Boltzmann) package.
"""

import numpy as np
import pytest

from cfd_lbm import LBM, equilibrium


class TestEquilibrium:
    """Tests for the equilibrium distribution function."""

    def test_mass_conservation(self):
        """Sum of feq over all velocities must equal rho."""
        rho = np.ones((8, 8)) * 1.0
        ux = np.zeros((8, 8))
        uy = np.zeros((8, 8))

        feq = equilibrium(rho, ux, uy)

        np.testing.assert_allclose(feq.sum(axis=0), rho, rtol=1e-12)

    def test_output_shape(self):
        """feq must have shape (Q, nx, ny) = (9, nx, ny)."""
        nx, ny = 16, 12
        rho = np.ones((nx, ny))
        ux = np.zeros((nx, ny))
        uy = np.zeros((nx, ny))

        feq = equilibrium(rho, ux, uy)

        assert feq.shape == (9, nx, ny)

    def test_non_negative(self):
        """feq must be non-negative for physically valid (low-Ma) inputs."""
        rho = np.ones((8, 8))
        ux = np.full((8, 8), 0.05)
        uy = np.full((8, 8), 0.02)

        feq = equilibrium(rho, ux, uy)

        assert np.all(feq >= 0), "Equilibrium distribution should be non-negative."

    def test_rest_state_weights(self):
        """At rest (u=0), feq[i] = W[i] * rho."""
        from cfd_lbm.lbm import W

        rho = np.full((4, 4), 1.5)
        ux = np.zeros((4, 4))
        uy = np.zeros((4, 4))

        feq = equilibrium(rho, ux, uy)

        for i in range(9):
            expected = W[i] * rho
            np.testing.assert_allclose(feq[i], expected, rtol=1e-12)


class TestLBM:
    """Tests for the LBM D2Q9 lid-driven cavity solver."""

    def test_initialization(self):
        """Test that fields are initialized with correct shapes and parameters."""
        sim = LBM(nx=32, ny=32, Re=100.0, u_lid=0.1)

        assert sim.nx == 32
        assert sim.ny == 32
        assert sim.Re == 100.0
        assert sim.u_lid == 0.1
        assert sim.f.shape == (9, 32, 32)
        assert sim.step == 0

    def test_relaxation_time(self):
        """tau must satisfy tau = nu/cs² + 0.5; tau > 0.5 for stability."""
        sim = LBM(nx=32, ny=32, Re=100.0, u_lid=0.1)

        assert sim.tau > 0.5, "tau <= 0.5 leads to numerical instability."
        assert sim.omega == pytest.approx(1.0 / sim.tau)

    def test_distribution_sum_is_rho(self):
        """After initialization, sum of f over velocities should equal ~1 everywhere."""
        sim = LBM(nx=16, ny=16, Re=100.0, u_lid=0.1)

        rho, _, _ = sim.get_fields()

        np.testing.assert_allclose(rho, 1.0, rtol=1e-10)

    def test_solid_mask(self):
        """Solid boundary nodes must include all four walls (excluding lid)."""
        sim = LBM(nx=16, ny=16, Re=100.0, u_lid=0.1)

        assert sim.solid[:, 0].all(), "Bottom wall must be solid."
        assert sim.solid[0, :].all(), "Left wall must be solid."
        assert sim.solid[-1, :].all(), "Right wall must be solid."

    def test_advance_preserves_mass(self):
        """Total mass (sum of rho) should be approximately conserved.

        Zou-He moving-lid BC injects a small amount of mass; allow 0.1%.
        """
        sim = LBM(nx=16, ny=16, Re=100.0, u_lid=0.1)

        rho_init, _, _ = sim.get_fields()
        mass_init = rho_init.sum()

        sim.advance(10)

        rho_after, _, _ = sim.get_fields()
        mass_after = rho_after.sum()

        np.testing.assert_allclose(mass_after, mass_init, rtol=1e-3)

    def test_advance_no_nan_or_inf(self):
        """Fields must not contain NaN or Inf after advancing."""
        sim = LBM(nx=16, ny=16, Re=100.0, u_lid=0.1)

        sim.advance(20)

        rho, ux, uy = sim.get_fields()
        assert not np.any(np.isnan(rho))
        assert not np.any(np.isnan(ux))
        assert not np.any(np.isnan(uy))
        assert not np.any(np.isinf(rho))
        assert not np.any(np.isinf(ux))
        assert not np.any(np.isinf(uy))

    def test_step_counter_increments(self):
        """step counter must increment by n_steps on each advance call."""
        sim = LBM(nx=16, ny=16, Re=100.0, u_lid=0.1)

        sim.advance(5)
        assert sim.step == 5

        sim.advance(3)
        assert sim.step == 8

    def test_lid_velocity_applied(self):
        """After advancing, top-wall nodes should show non-zero x-velocity."""
        sim = LBM(nx=32, ny=32, Re=100.0, u_lid=0.1)
        sim.advance(50)

        _, ux, _ = sim.get_fields()
        # Top row (y=ny-1): average ux should be > 0 due to lid
        assert np.mean(ux[:, -1]) > 0, "Lid velocity not being applied."

    def test_get_fields_shape(self):
        """get_fields() must return arrays of shape (nx, ny)."""
        sim = LBM(nx=24, ny=20, Re=200.0, u_lid=0.1)

        rho, ux, uy = sim.get_fields()

        assert rho.shape == (24, 20)
        assert ux.shape == (24, 20)
        assert uy.shape == (24, 20)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
