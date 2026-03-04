"""
Cavity Flow Calculation Module

This module simulates the Lid-Driven Cavity Flow using the fractional-step
(projection) method for the incompressible Navier-Stokes equations.
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


class CavityFlow:
    """Lid-Driven Cavity Flow Solver using the fractional-step projection method."""

    def __init__(self, L=1.0, rho=1.0, mu=0.01, nx=65, ny=65, advection_scheme="upwind"):
        """
        Parameters:
        -----------
        L : float
            Size of the square cavity [m]
        rho : float
            Fluid density [kg/m³]
        mu : float
            Dynamic viscosity [Pa·s]
        nx, ny : int
            Grid resolution
        advection_scheme : str
            "upwind" for 1st-order or "upwind2" for 2nd-order linear-upwind
        """
        self.L = L
        self.rho = rho
        self.mu = mu
        self.nu = mu / rho
        self.nx = nx
        self.ny = ny
        self.advection_scheme = advection_scheme

        self.dx = L / (nx - 1)
        self.dy = L / (ny - 1)
        self.x = np.linspace(0, L, nx)
        self.y = np.linspace(0, L, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        self.U_lid = 1.0

        self.u = np.zeros((ny, nx))
        self.v = np.zeros((ny, nx))
        self.p = np.zeros((ny, nx))

        self.dt_base = 0.01
        self.dt = self.dt_base

        self.time = 0.0
        self.iteration = 0

        logger.info(
            f"Initialized CavityFlow: Re={self.get_Reynolds_number():.1f}, "
            f"Grid={nx}x{ny}, scheme={advection_scheme}"
        )

    def get_Reynolds_number(self):
        return self.rho * self.U_lid * self.L / self.mu

    def _compute_optimal_dt(self):
        """Dynamic time-step control based on CFL and Viscous conditions."""
        u_max = np.max(np.abs(self.u)) + np.max(np.abs(self.v)) + 1e-10
        dt_viscous = 0.1 * min(self.dx**2, self.dy**2) / (self.nu + 1e-10)
        dt_convective = 0.25 * min(self.dx, self.dy) / u_max

        self.dt = min(self.dt_base, dt_viscous * 0.1, dt_convective * 0.5)

    def set_boundary_conditions(self):
        # Top Lid
        self.u[-1, :] = self.U_lid
        self.v[-1, :] = 0.0
        # Walls
        self.u[0, :] = 0.0
        self.u[:, 0] = 0.0
        self.u[:, -1] = 0.0
        self.v[0, :] = 0.0
        self.v[:, 0] = 0.0
        self.v[:, -1] = 0.0
        self.v[-1, :] = 0.0

    def _advect(self, phi, u, v):
        """Return u·∂φ/∂x + v·∂φ/∂y for interior cells [1:-1,1:-1].

        Dispatches to the scheme selected at construction time.
        Returns an array of shape (ny-2, nx-2).
        """
        ii = slice(1, -1)
        jj = slice(1, -1)
        u_ij = u[ii, jj]
        v_ij = v[ii, jj]

        # 1st-order upwind (boundary fallback and "upwind" scheme)
        dphi_dx = (
            np.maximum(u_ij, 0) * (phi[ii, jj] - phi[ii, :-2]) / self.dx
            + np.minimum(u_ij, 0) * (phi[ii, 2:] - phi[ii, jj]) / self.dx
        )
        dphi_dy = (
            np.maximum(v_ij, 0) * (phi[ii, jj] - phi[:-2, jj]) / self.dy
            + np.minimum(v_ij, 0) * (phi[2:, jj] - phi[ii, jj]) / self.dy
        )

        if self.advection_scheme == "upwind2":
            ny_phi, nx_phi = phi.shape
            if nx_phi >= 5 and ny_phi >= 5:
                u_in = u[2:-2, 2:-2]
                v_in = v[2:-2, 2:-2]
                dphi_dx[1:-1, 1:-1] = np.maximum(u_in, 0) * (
                    3 * phi[2:-2, 2:-2] - 4 * phi[2:-2, 1:-3] + phi[2:-2, :-4]
                ) / (2.0 * self.dx) + np.minimum(u_in, 0) * (
                    -phi[2:-2, 4:] + 4 * phi[2:-2, 3:-1] - 3 * phi[2:-2, 2:-2]
                ) / (
                    2.0 * self.dx
                )
                dphi_dy[1:-1, 1:-1] = np.maximum(v_in, 0) * (
                    3 * phi[2:-2, 2:-2] - 4 * phi[1:-3, 2:-2] + phi[:-4, 2:-2]
                ) / (2.0 * self.dy) + np.minimum(v_in, 0) * (
                    -phi[4:, 2:-2] + 4 * phi[3:-1, 2:-2] - 3 * phi[2:-2, 2:-2]
                ) / (
                    2.0 * self.dy
                )

        return dphi_dx + dphi_dy

    def solve_poisson_pressure(self, u=None, v=None, max_iter=50, tol=1e-4):
        """Pressure Poisson equation: ∇²p = (ρ/dt) ∇·u* with SOR.

        Parameters
        ----------
        u, v : ndarray, optional
            Intermediate (predicted) velocity fields. Defaults to self.u, self.v.
        """
        if u is None:
            u = self.u
        if v is None:
            v = self.v
        p_new = self.p.copy()

        # RHS = +rho/dt * div(u*)  — positive sign for divergence-free projection
        dudx = (u[1:-1, 2:] - u[1:-1, :-2]) / (2 * self.dx)
        dvdy = (v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * self.dy)
        rhs = self.rho * np.clip(dudx + dvdy, -10, 10) / (self.dt + 1e-10)

        omega_sor = 1.4
        coeff = 2.0 / (self.dx**2) + 2.0 / (self.dy**2)

        for _ in range(max_iter):
            p_old = p_new.copy()

            lap_x = (p_new[1:-1, 2:] + p_new[1:-1, :-2]) / (self.dx**2)
            lap_y = (p_new[2:, 1:-1] + p_new[:-2, 1:-1]) / (self.dy**2)

            res = np.clip(rhs - (lap_x + lap_y - coeff * p_new[1:-1, 1:-1]), -1e3, 1e3)
            p_new[1:-1, 1:-1] -= omega_sor * res / coeff

            # Neumann BCs (all walls)
            p_new[0, :] = p_new[1, :]
            p_new[-1, :] = p_new[-2, :]
            p_new[:, 0] = p_new[:, 1]
            p_new[:, -1] = p_new[:, -2]

            if np.any(np.isnan(p_new)):
                p_new = np.nan_to_num(p_old)
                break

            if np.max(np.abs(p_new - p_old)) < tol:
                break

        self.p = p_new

    def step(self):
        """Fractional-step (projection) method.

        Algorithm:
          1. Predict u* with advection + diffusion (no pressure).
          2. Solve ∇²p = (ρ/dt) ∇·u*
          3. Correct u = u* - (dt/ρ) ∇p → divergence-free
        """
        self._compute_optimal_dt()
        self.set_boundary_conditions()

        ii, jj = slice(1, -1), slice(1, -1)

        # Viscous terms
        lap_u = (self.u[2:, jj] - 2 * self.u[ii, jj] + self.u[:-2, jj]) / self.dy**2 + (
            self.u[ii, 2:] - 2 * self.u[ii, jj] + self.u[ii, :-2]
        ) / self.dx**2
        lap_v = (self.v[2:, jj] - 2 * self.v[ii, jj] + self.v[:-2, jj]) / self.dy**2 + (
            self.v[ii, 2:] - 2 * self.v[ii, jj] + self.v[ii, :-2]
        ) / self.dx**2

        # Advection terms
        u_adv = self._advect(self.u, self.u, self.v)
        v_adv = self._advect(self.v, self.u, self.v)

        # --- Step 1: Predict u* (advection + diffusion, no pressure) ---
        u_star = self.u.copy()
        v_star = self.v.copy()
        u_star[ii, jj] = np.clip(self.u[ii, jj] + self.dt * (-u_adv + self.nu * lap_u), -2.0, 2.0)
        v_star[ii, jj] = np.clip(self.v[ii, jj] + self.dt * (-v_adv + self.nu * lap_v), -2.0, 2.0)
        # Apply BCs to predicted velocity
        u_star[-1, :] = self.U_lid
        v_star[-1, :] = 0.0
        u_star[0, :] = 0.0
        u_star[:, 0] = 0.0
        u_star[:, -1] = 0.0
        v_star[0, :] = 0.0
        v_star[:, 0] = 0.0
        v_star[:, -1] = 0.0
        v_star[-1, :] = 0.0

        # --- Step 2: Solve ∇²p = (ρ/dt) ∇·u* ---
        self.solve_poisson_pressure(u=u_star, v=v_star)

        # --- Step 3: Correct velocity with pressure gradient ---
        dpdx = (self.p[ii, 2:] - self.p[ii, :-2]) / (2 * self.dx)
        dpdy = (self.p[2:, jj] - self.p[:-2, jj]) / (2 * self.dy)

        u_new = u_star.copy()
        v_new = v_star.copy()
        u_new[ii, jj] = np.clip(u_star[ii, jj] - self.dt / self.rho * dpdx, -2.0, 2.0)
        v_new[ii, jj] = np.clip(v_star[ii, jj] - self.dt / self.rho * dpdy, -2.0, 2.0)

        if np.any(np.isnan(u_new)) or np.any(np.isnan(v_new)):
            logger.error("Instability detected (NaN). Aborting step.")
            raise RuntimeError("Solver diverged.")

        self.u, self.v = u_new, v_new
        self.set_boundary_conditions()
        self.time += self.dt
        self.iteration += 1

    def solve_steady_state(self, max_iterations=500, verbose=True):
        for i in range(max_iterations):
            u_old = self.u.copy()
            try:
                self.step()
            except RuntimeError as e:
                logger.error(f"Aborted: {e}")
                break
            res = np.max(np.abs(self.u - u_old))
            if verbose and (i + 1) % 100 == 0:
                logger.info(f"Step {i+1}: t={self.time:.4f}, residual={res:.2e}")
            if res < 1e-6 and i > 50:
                logger.info("Steady state reached.")
                break

    def plot_velocity_field(self, figsize=(10, 8)):
        fig, ax = plt.subplots(figsize=figsize)
        mag = np.sqrt(self.u**2 + self.v**2)
        cf = ax.contourf(self.X, self.Y, mag, levels=20, cmap="viridis")
        ax.streamplot(self.x, self.y, self.u, self.v, color="white", density=1.0)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(
            f"Velocity Magnitude — Re={self.get_Reynolds_number():.0f}, t={self.time:.3f}s"
        )
        plt.colorbar(cf, label="Velocity [m/s]")
        return fig

    def plot_streamlines(self, figsize=(10, 8)):
        """Plot streamlines coloured by velocity magnitude."""
        fig, ax = plt.subplots(figsize=figsize)
        mag = np.sqrt(self.u**2 + self.v**2)
        strm = ax.streamplot(
            self.x,
            self.y,
            self.u,
            self.v,
            color=mag,
            cmap="cool",
            linewidth=1.5,
            density=1.5,
        )
        ax.set_xlim(0, self.L)
        ax.set_ylim(0, self.L)
        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(f"Streamlines — Re={self.get_Reynolds_number():.0f}, t={self.time:.3f}s")
        plt.colorbar(strm.lines, ax=ax, label="Velocity [m/s]")
        return fig

    def plot_pressure_field(self, figsize=(10, 8)):
        """Plot the pressure distribution in the cavity."""
        fig, ax = plt.subplots(figsize=figsize)
        cf = ax.contourf(self.X, self.Y, self.p, levels=20, cmap="RdBu_r")
        ax.set_xlim(0, self.L)
        ax.set_ylim(0, self.L)
        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(f"Pressure — Re={self.get_Reynolds_number():.0f}, t={self.time:.3f}s")
        plt.colorbar(cf, ax=ax, label="Pressure [Pa]")
        return fig
