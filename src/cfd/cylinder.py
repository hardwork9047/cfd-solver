"""
Cylinder Flow Calculation Module

This module simulates flow around a circular cylinder using the
incompressible Navier-Stokes equations with numerical stabilization.
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


class CylinderFlow:
    """Flow around a circular cylinder solver with numerical stabilization."""

    def __init__(self, L=2.0, H=1.0, D=0.1, U_inf=1.0, rho=1.0, mu=0.01, nx=200, ny=100):
        """
        Parameters:
        -----------
        L : float
            Length of the computational domain [m]
        H : float
            Height of the computational domain [m]
        D : float
            Diameter of the cylinder [m]
        U_inf : float
            Free-stream velocity [m/s]
        rho : float
            Fluid density [kg/m³]
        mu : float
            Dynamic viscosity [Pa·s]
        nx, ny : int
            Grid resolution
        """
        if D >= H:
            raise ValueError("Cylinder diameter D must be smaller than domain height H.")
        if D >= L / 2:
            raise ValueError("Cylinder diameter D must be smaller than half domain length L.")

        self.L = L
        self.H = H
        self.D = D
        self.R = D / 2  # Cylinder radius
        self.U_inf = U_inf
        self.rho = rho
        self.mu = mu
        self.nu = mu / rho
        self.nx = nx
        self.ny = ny

        # Grid
        self.dx = L / (nx - 1)
        self.dy = H / (ny - 1)
        self.x = np.linspace(0, L, nx)
        self.y = np.linspace(0, H, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        # Cylinder center (positioned at 1/4 of domain length)
        self.x_c = L / 4
        self.y_c = H / 2

        # Initialize flow fields
        self.u = np.ones((ny, nx)) * U_inf
        self.v = np.zeros((ny, nx))
        self.p = np.zeros((ny, nx))

        # Create mask for solid cylinder
        self.mask = self._create_cylinder_mask()

        # Time stepping
        self.dt_base = 0.001
        self.dt = self.dt_base
        self.time = 0.0
        self.iteration = 0

        # Apply initial boundary conditions
        self.set_boundary_conditions()

        Re = self.get_Reynolds_number()
        logger.info(
            f"Initialized CylinderFlow: Re={Re:.1f}, D={D}m, Grid={nx}x{ny}, "
            f"Cylinder at ({self.x_c:.2f}, {self.y_c:.2f})"
        )

    def _create_cylinder_mask(self):
        """Create a boolean mask for the cylinder interior."""
        dist = np.sqrt((self.X - self.x_c) ** 2 + (self.Y - self.y_c) ** 2)
        return dist <= self.R

    def get_Reynolds_number(self):
        """Calculate Reynolds number based on cylinder diameter."""
        return self.rho * self.U_inf * self.D / self.mu

    def _compute_optimal_dt(self):
        """Dynamic time-step control based on CFL and viscous conditions."""
        u_max = np.max(np.abs(self.u)) + np.max(np.abs(self.v)) + 1e-10
        dt_viscous = 0.1 * min(self.dx**2, self.dy**2) / (self.nu + 1e-10)
        dt_convective = 0.25 * min(self.dx, self.dy) / u_max

        self.dt = min(self.dt_base, dt_viscous * 0.1, dt_convective * 0.5)

    def set_boundary_conditions(self):
        """Apply boundary conditions to velocity fields."""
        # Inlet (left): uniform flow
        self.u[:, 0] = self.U_inf
        self.v[:, 0] = 0.0

        # Outlet (right): zero gradient (Neumann)
        self.u[:, -1] = self.u[:, -2]
        self.v[:, -1] = self.v[:, -2]

        # Top and bottom walls: slip or symmetry
        self.u[0, :] = self.u[1, :]
        self.u[-1, :] = self.u[-2, :]
        self.v[0, :] = 0.0
        self.v[-1, :] = 0.0

        # Cylinder surface: no-slip
        self.u[self.mask] = 0.0
        self.v[self.mask] = 0.0

    def solve_poisson_pressure(self, max_iter=50, tol=1e-4):
        """Pressure Poisson equation solver with SOR."""
        p_new = self.p.copy()

        # RHS = -rho/dt * (du/dx + dv/dy)
        dudx = np.zeros_like(self.u)
        dvdy = np.zeros_like(self.v)

        dudx[1:-1, 1:-1] = (self.u[1:-1, 2:] - self.u[1:-1, :-2]) / (2 * self.dx)
        dvdy[1:-1, 1:-1] = (self.v[2:, 1:-1] - self.v[:-2, 1:-1]) / (2 * self.dy)

        rhs = -self.rho * np.clip(dudx + dvdy, -10, 10) / (self.dt + 1e-10)

        omega_sor = 1.5
        coeff = 2.0 / (self.dx**2) + 2.0 / (self.dy**2)

        for it in range(max_iter):
            for i in range(1, self.ny - 1):
                for j in range(1, self.nx - 1):
                    if self.mask[i, j]:
                        continue  # Skip cylinder interior

                    p_gs = (
                        (p_new[i - 1, j] + p_new[i + 1, j]) / (self.dy**2)
                        + (p_new[i, j - 1] + p_new[i, j + 1]) / (self.dx**2)
                        - rhs[i, j]
                    ) / coeff

                    p_new[i, j] = p_new[i, j] + omega_sor * (p_gs - p_new[i, j])

            # Boundary conditions for pressure
            p_new[:, 0] = p_new[:, 1]  # Inlet
            p_new[:, -1] = 0  # Outlet (reference)
            p_new[0, :] = p_new[1, :]  # Bottom
            p_new[-1, :] = p_new[-2, :]  # Top

            # Cylinder surface: Neumann condition
            p_new[self.mask] = np.mean(
                [
                    p_new[np.roll(self.mask, 1, axis=0)],
                    p_new[np.roll(self.mask, -1, axis=0)],
                    p_new[np.roll(self.mask, 1, axis=1)],
                    p_new[np.roll(self.mask, -1, axis=1)],
                ]
            )

        self.p = p_new

    def advection_upwind(self, phi, u, v):
        """1st-order upwind scheme for advection."""
        phi_new = phi.copy()

        for i in range(1, self.ny - 1):
            for j in range(1, self.nx - 1):
                if self.mask[i, j]:
                    continue

                dphi_dx = 0.0
                dphi_dy = 0.0

                # Upwind in x-direction
                if u[i, j] > 0:
                    dphi_dx = (phi[i, j] - phi[i, j - 1]) / self.dx
                else:
                    dphi_dx = (phi[i, j + 1] - phi[i, j]) / self.dx

                # Upwind in y-direction
                if v[i, j] > 0:
                    dphi_dy = (phi[i, j] - phi[i - 1, j]) / self.dy
                else:
                    dphi_dy = (phi[i + 1, j] - phi[i, j]) / self.dy

                phi_new[i, j] = phi[i, j] - self.dt * (u[i, j] * dphi_dx + v[i, j] * dphi_dy)

        return phi_new

    def diffusion_central(self, phi):
        """2nd-order central difference for diffusion."""
        phi_new = phi.copy()

        for i in range(1, self.ny - 1):
            for j in range(1, self.nx - 1):
                if self.mask[i, j]:
                    continue

                d2phi_dx2 = (phi[i, j + 1] - 2 * phi[i, j] + phi[i, j - 1]) / (self.dx**2)
                d2phi_dy2 = (phi[i + 1, j] - 2 * phi[i, j] + phi[i - 1, j]) / (self.dy**2)

                phi_new[i, j] = phi[i, j] + self.dt * self.nu * (d2phi_dx2 + d2phi_dy2)

        return phi_new

    def step(self):
        """Perform one time step of the simulation."""
        self._compute_optimal_dt()

        # Store old velocities
        u_old = self.u.copy()
        v_old = self.v.copy()

        # Advection step (upwind)
        u_star = self.advection_upwind(u_old, u_old, v_old)
        v_star = self.advection_upwind(v_old, u_old, v_old)

        # Diffusion step (central)
        u_star = self.diffusion_central(u_star)
        v_star = self.diffusion_central(v_star)

        # Solve pressure Poisson equation
        self.u = u_star
        self.v = v_star
        self.solve_poisson_pressure()

        # Pressure correction
        dp_dx = np.zeros_like(self.p)
        dp_dy = np.zeros_like(self.p)

        dp_dx[:, 1:-1] = (self.p[:, 2:] - self.p[:, :-2]) / (2 * self.dx)
        dp_dy[1:-1, :] = (self.p[2:, :] - self.p[:-2, :]) / (2 * self.dy)

        self.u = u_star - self.dt / self.rho * dp_dx
        self.v = v_star - self.dt / self.rho * dp_dy

        # Apply boundary conditions
        self.set_boundary_conditions()

        # Velocity clipping for stability
        self.u = np.clip(self.u, -5 * self.U_inf, 5 * self.U_inf)
        self.v = np.clip(self.v, -5 * self.U_inf, 5 * self.U_inf)

        self.time += self.dt
        self.iteration += 1

    def run(self, t_end=10.0, output_interval=100):
        """Run simulation until t_end."""
        logger.info(f"Starting simulation until t={t_end}s (output every {output_interval} steps)")

        while self.time < t_end:
            self.step()

            if self.iteration % output_interval == 0:
                u_max = np.max(np.abs(self.u))
                v_max = np.max(np.abs(self.v))
                logger.info(
                    f"Iteration {self.iteration}, Time {self.time:.4f}s, "
                    f"dt={self.dt:.2e}, u_max={u_max:.3f}, v_max={v_max:.3f}"
                )

        logger.info(f"Simulation complete at t={self.time:.4f}s")

    def compute_drag_coefficient(self):
        """
        Compute drag coefficient using pressure and viscous forces.

        Note: This is a simplified calculation for educational purposes.
        """
        # Find cylinder surface points
        surface_mask = self.mask & (
            np.roll(self.mask, 1, axis=0) == False
            | np.roll(self.mask, -1, axis=0) == False
            | np.roll(self.mask, 1, axis=1) == False
            | np.roll(self.mask, -1, axis=1) == False
        )

        # Pressure drag (simplified)
        p_surface = self.p[surface_mask]
        F_p = np.sum(p_surface) * self.dx * self.dy

        # Dynamic pressure
        q_inf = 0.5 * self.rho * self.U_inf**2

        # Drag coefficient (approximate)
        C_d = F_p / (q_inf * self.D * 1.0)  # 1.0 is unit depth

        return C_d

    def plot_velocity_field(self, skip=5, figsize=(12, 5)):
        """Plot velocity field with streamlines."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # Velocity magnitude
        vel_mag = np.sqrt(self.u**2 + self.v**2)
        vel_mag_masked = np.ma.array(vel_mag, mask=self.mask)

        im1 = ax1.contourf(self.X, self.Y, vel_mag_masked, levels=20, cmap="jet")
        ax1.quiver(
            self.X[::skip, ::skip],
            self.Y[::skip, ::skip],
            self.u[::skip, ::skip],
            self.v[::skip, ::skip],
            scale=20,
            color="white",
            alpha=0.6,
        )
        circle = plt.Circle((self.x_c, self.y_c), self.R, color="black", fill=True)
        ax1.add_patch(circle)
        ax1.set_xlabel("x [m]")
        ax1.set_ylabel("y [m]")
        ax1.set_title(f"Velocity Magnitude at t={self.time:.2f}s")
        ax1.set_aspect("equal")
        plt.colorbar(im1, ax=ax1, label="Velocity [m/s]")

        # Streamlines
        p_masked = np.ma.array(self.p, mask=self.mask)
        im2 = ax2.contourf(self.X, self.Y, p_masked, levels=20, cmap="RdBu_r")
        ax2.streamplot(
            self.x,
            self.y,
            self.u.T,
            self.v.T,
            color="black",
            linewidth=0.5,
            density=1.5,
            arrowsize=0.8,
        )
        circle2 = plt.Circle((self.x_c, self.y_c), self.R, color="gray", fill=True)
        ax2.add_patch(circle2)
        ax2.set_xlabel("x [m]")
        ax2.set_ylabel("y [m]")
        ax2.set_title(f"Pressure Field and Streamlines")
        ax2.set_aspect("equal")
        plt.colorbar(im2, ax=ax2, label="Pressure [Pa]")

        plt.tight_layout()
        return fig

    def get_velocity_on_centerline(self):
        """Get velocity profile along horizontal centerline."""
        j_center = self.ny // 2
        u_centerline = self.u[j_center, :]
        x_centerline = self.x
        return x_centerline, u_centerline
