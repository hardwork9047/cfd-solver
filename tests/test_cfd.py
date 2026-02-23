"""
CFD Solver Test Suite
"""

import numpy as np
import pytest

from cfd import CavityFlow, CircularPoiseuille, PlanePoiseuille


class TestPlanePoiseuille:
    """Tests for Plane Poiseuille Flow."""

    def test_initialization(self):
        """Test basic initialization."""
        L = 0.01
        dp_dx = -100
        mu = 1e-3

        flow = PlanePoiseuille(L, dp_dx, mu)

        assert flow.L == L
        assert flow.dp_dx == dp_dx
        assert flow.mu == mu
        assert len(flow.y) == flow.ny

    def test_analytical_solution_shape(self):
        """Test the shape and validity of the analytical solution."""
        L = 0.01
        dp_dx = -100
        mu = 1e-3

        flow = PlanePoiseuille(L, dp_dx, mu, ny=51)
        u = flow.analytical_solution()

        assert u.shape == (51,)
        assert np.max(u) > 0  # Velocity should be positive for negative pressure gradient

    def test_max_velocity(self):
        """Test max velocity calculation."""
        L = 0.01
        dp_dx = -100
        mu = 1e-3

        flow = PlanePoiseuille(L, dp_dx, mu)
        u_max = flow.get_max_velocity()

        # Compare with theoretical formula
        expected = -dp_dx * L**2 / (8 * mu)
        assert np.isclose(u_max, expected)

    def test_flow_rate(self):
        """Test flow rate calculation."""
        L = 0.01
        dp_dx = -100
        mu = 1e-3

        flow = PlanePoiseuille(L, dp_dx, mu)
        Q = flow.get_flow_rate()

        # Compare with theoretical formula
        expected = -dp_dx * L**3 / (12 * mu)
        assert np.isclose(Q, expected)


class TestCircularPoiseuille:
    """Tests for Circular Poiseuille Flow."""

    def test_initialization(self):
        """Test initialization of pipe flow."""
        R = 0.005
        dp_dx = -100
        mu = 1e-3

        flow = CircularPoiseuille(R, dp_dx, mu)

        assert flow.R == R
        assert flow.dp_dx == dp_dx
        assert flow.mu == mu
        assert len(flow.r) == flow.nr

    def test_analytical_solution_shape(self):
        """Test shape and boundary values of pipe flow analytical solution."""
        R = 0.005
        dp_dx = -100
        mu = 1e-3

        flow = CircularPoiseuille(R, dp_dx, mu, nr=51)
        u = flow.analytical_solution()

        assert u.shape == (51,)
        assert np.max(u) > 0
        assert u[-1] == 0  # Velocity must be zero at the wall

    def test_max_velocity(self):
        """Test max velocity for circular flow."""
        R = 0.005
        dp_dx = -100
        mu = 1e-3

        flow = CircularPoiseuille(R, dp_dx, mu)
        u_max = flow.get_max_velocity()

        expected = -dp_dx * R**2 / (4 * mu)
        assert np.isclose(u_max, expected)

    def test_hagen_poiseuille_formula(self):
        """Verify against Hagen-Poiseuille formula."""
        R = 0.005
        dp_dx = -100
        mu = 1e-3

        flow = CircularPoiseuille(R, dp_dx, mu)
        Q = flow.get_flow_rate()

        expected = -np.pi * dp_dx * R**4 / (8 * mu)
        assert np.isclose(Q, expected)


class TestCavityFlow:
    """Tests for Lid-Driven Cavity Flow."""

    def test_initialization(self):
        """Test field initialization."""
        L = 1.0
        rho = 1.0
        mu = 0.01

        flow = CavityFlow(L, rho, mu)

        assert flow.L == L
        assert flow.rho == rho
        assert flow.mu == mu
        assert flow.u.shape == (flow.ny, flow.nx)
        assert flow.v.shape == (flow.ny, flow.nx)
        assert flow.p.shape == (flow.ny, flow.nx)

    def test_reynolds_number(self):
        """Test Reynolds number calculation."""
        L = 1.0
        rho = 1.0
        mu = 0.01

        flow = CavityFlow(L, rho, mu)
        Re = flow.get_Reynolds_number()

        # Re = rho * U_lid * L / mu = 1 * 1 * 1 / 0.01 = 100
        expected = 100.0
        assert np.isclose(Re, expected)

    def test_boundary_conditions(self):
        """Test boundary condition application."""
        flow = CavityFlow(1.0, 1.0, 0.01, nx=9, ny=9)
        flow.set_boundary_conditions()

        # Top lid velocity
        assert np.allclose(flow.u[-1, 1:-1], flow.U_lid)
        assert np.allclose(flow.v[-1, :], 0.0)

        # Walls (no-slip)
        assert np.allclose(flow.u[0, :], 0.0)
        assert np.allclose(flow.u[:, 0], 0.0)
        assert np.allclose(flow.u[:, -1], 0.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
