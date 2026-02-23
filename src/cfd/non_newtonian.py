"""
Non-Newtonian Flow Calculation Module

This module implements solvers for non-Newtonian fluids, specifically
focusing on the Power-Law (Ostwald-de Waele) model.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# Configuration for logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Japanese font setting
rcParams["font.family"] = "sans-serif"


class PowerLawPlanePoiseuille:
    """
    Plane Poiseuille Flow for a Power-Law fluid.
    The viscosity follows: mu_eff = K * (gamma_dot)^(n-1)
    n < 1: Shear-thinning (Pseudoplastic)
    n = 1: Newtonian
    n > 1: Shear-thickening (Dilatant)
    """

    def __init__(self, L, dp_dx, K, n, ny=101):
        """
        Parameters:
        -----------
        L : float
            Distance between plates [m]
        dp_dx : float
            Pressure gradient dp/dx [Pa/m] (should be negative)
        K : float
            Consistency index [Pa·s^n]
        n : float
            Flow behavior index
        ny : int
            Number of grid points
        """
        if L <= 0 or K <= 0 or n <= 0:
            raise ValueError("L, K, and n must be positive values.")

        self.L = L
        self.dp_dx = dp_dx
        self.K = K
        self.n = n
        self.ny = ny
        self.y = np.linspace(-L / 2, L / 2, ny)

        logger.info(f"Initialized PowerLawPlanePoiseuille: K={K}, n={n}, ny={ny}")

    def analytical_solution(self):
        """
        Calculate the analytical solution for a Power-Law fluid.
        u(y) = [n/(n+1)] * (|dp/dx|/K)^(1/n) * [(L/2)^((n+1)/n) - |y|^((n+1)/n)]
        """
        abs_dpdx = abs(self.dp_dx)
        term1 = self.n / (self.n + 1)
        term2 = (abs_dpdx / self.K) ** (1 / self.n)
        h_half = self.L / 2

        u = (
            term1
            * term2
            * (h_half ** ((self.n + 1) / self.n) - np.abs(self.y) ** ((self.n + 1) / self.n))
        )
        return u

    def numerical_solution(self, tol=1e-8, max_iter=200000, omega=0.1):
        """
        Numerical solution using an iterative FDM approach.
        We initialize with a Newtonian profile and use a small omega for non-linear stability.
        """
        # Start with Newtonian initialization (n=1 solution)
        mu_ref = self.K
        u = -(self.dp_dx / (2 * mu_ref)) * (self.L**2 / 4 - self.y**2)
        u[0] = 0
        u[-1] = 0

        dy = self.y[1] - self.y[0]

        logger.info(f"Starting non-linear solver for Power-Law fluid (n={self.n}, omega={omega})")

        for iteration in range(max_iter):
            u_old = u.copy()

            # Calculate effective viscosity at staggered points (i+1/2)
            # gamma_dot = |du/dy|
            # We use a small epsilon to avoid singularity at zero gradient
            eps = 1e-10
            du_dy = np.abs(np.diff(u) / dy)
            mu_eff = self.K * (du_dy + eps) ** (self.n - 1)

            # Update velocity at interior points
            for i in range(1, self.ny - 1):
                # Viscosity at staggered points: m_plus = mu_{i+1/2}, m_minus = mu_{i-1/2}
                m_plus = mu_eff[i]
                m_minus = mu_eff[i - 1]

                # Discrete form: (m_plus*(u[i+1]-u[i])/dy - m_minus*(u[i]-u[i-1])/dy) / dy = dp_dx
                u_gs = (m_plus * u[i + 1] + m_minus * u[i - 1] - dy**2 * self.dp_dx) / (
                    m_plus + m_minus
                )
                u[i] = u_old[i] + omega * (u_gs - u_old[i])

            # Stability check
            if np.any(np.isnan(u)):
                logger.error("Numerical instability: NaN detected.")
                raise RuntimeError("Solver diverged.")

            # Check for convergence
            residual = np.max(np.abs(u - u_old))
            if residual < tol and iteration > 100:
                logger.info(f"Converged in {iteration+1} iterations. Residual: {residual:.2e}")
                return u

        logger.warning(f"Reached max iterations. Residual: {residual:.2e}")
        return u

    def plot_comparison(self, u_numerical=None):
        """Visualize how the velocity profile changes with n."""
        u_ana = self.analytical_solution()

        plt.figure(figsize=(10, 6))
        plt.plot(u_ana, self.y, "b-", label=f"Analytical (n={self.n})")

        if u_numerical is not None:
            plt.plot(u_numerical, self.y, "ro", markersize=3, label="Numerical", alpha=0.5)

        plt.xlabel("Velocity u [m/s]")
        plt.ylabel("y [m]")
        plt.title(f"Power-Law Fluid Profile (n={self.n})")
        plt.legend()
        plt.grid(True, alpha=0.3)
        return plt.gcf()
