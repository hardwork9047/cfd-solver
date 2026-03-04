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

    def __init__(self, L=2.0, H=1.0, D=0.1, U_inf=1.0, rho=1.0, mu=0.01, nx=200, ny=100,
                 advection_scheme="upwind2"):
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
        self.advection_scheme = advection_scheme  # "upwind" or "upwind2"

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

    def solve_poisson_pressure(self, u=None, v=None, max_iter=100, tol=1e-4):
        """Pressure Poisson equation solver - vectorized SOR.

        Parameters
        ----------
        u, v : ndarray, optional
            Velocity fields to compute RHS from. Defaults to self.u, self.v.
            Pass the predicted (intermediate) velocity for a projection method.
        """
        if u is None:
            u = self.u
        if v is None:
            v = self.v
        p_new = self.p.copy()
        ii, jj = slice(1, -1), slice(1, -1)

        # RHS = +rho/dt * (du/dx + dv/dy)  — positive sign for divergence-free correction
        dudx = np.zeros_like(u)
        dvdy = np.zeros_like(v)
        dudx[ii, jj] = (u[ii, 2:] - u[ii, :-2]) / (2 * self.dx)
        dvdy[ii, jj] = (v[2:, jj] - v[:-2, jj]) / (2 * self.dy)
        rhs = self.rho * np.clip(dudx + dvdy, -10, 10) / (self.dt + 1e-10)

        omega = 1.5
        coeff = 2.0 / self.dx**2 + 2.0 / self.dy**2

        for _ in range(max_iter):
            p_old = p_new.copy()

            lap_x = (p_new[ii, 2:] + p_new[ii, :-2]) / self.dx**2
            lap_y = (p_new[2:, jj] + p_new[:-2, jj]) / self.dy**2

            # Residual-based SOR with clipping for stability
            res = np.clip(rhs[ii, jj] - (lap_x + lap_y - coeff * p_new[ii, jj]), -1e3, 1e3)
            p_new[ii, jj] -= omega * res / coeff

            # Boundary conditions for pressure
            p_new[:, 0] = p_new[:, 1]   # Inlet (Neumann)
            p_new[:, -1] = 0.0          # Outlet (Dirichlet reference)
            p_new[0, :] = p_new[1, :]   # Bottom (Neumann)
            p_new[-1, :] = p_new[-2, :] # Top (Neumann)

            # Cylinder: zero pressure (Dirichlet) keeps solution bounded inside obstacle
            p_new[self.mask] = 0.0

            if np.max(np.abs(p_new - p_old)) < tol:
                break

        self.p = p_new

    def advection_upwind(self, phi, u, v):
        """1st-order upwind scheme for advection - vectorized."""
        phi_new = phi.copy()
        ii, jj = slice(1, -1), slice(1, -1)

        u_ij = u[ii, jj]
        v_ij = v[ii, jj]

        dphi_dx = (
            np.maximum(u_ij, 0) * (phi[ii, jj] - phi[ii, :-2])
            + np.minimum(u_ij, 0) * (phi[ii, 2:] - phi[ii, jj])
        ) / self.dx
        dphi_dy = (
            np.maximum(v_ij, 0) * (phi[ii, jj] - phi[:-2, jj])
            + np.minimum(v_ij, 0) * (phi[2:, jj] - phi[ii, jj])
        ) / self.dy

        phi_new[ii, jj] = phi[ii, jj] - self.dt * (dphi_dx + dphi_dy)
        phi_new[self.mask] = 0.0
        return phi_new

    def advection_upwind2(self, phi, u, v):
        """2nd-order upwind (linear-upwind) advection scheme.

        3-point one-sided stencils (O(Δx²), dispersive rather than diffusive):
          u>0: dφ/dx ≈ (3φ_j - 4φ_{j-1} + φ_{j-2}) / (2Δx)
          u<0: dφ/dx ≈ (-φ_{j+2} + 4φ_{j+1} - 3φ_j) / (2Δx)
        Numerical viscosity ≈ 0, so Re_eff ≈ Re_physical.
        Boundary layer of the interior (j=1, j=nx-2) uses 1st-order fallback.
        """
        phi_new = phi.copy()
        ii = slice(1, -1)
        jj = slice(1, -1)
        u_ij = u[ii, jj]
        v_ij = v[ii, jj]

        # --- 1st-order upwind for full interior (boundary layer fallback) ---
        dphi_dx = (
            np.maximum(u_ij, 0) * (phi[ii, jj] - phi[ii, :-2]) / self.dx
            + np.minimum(u_ij, 0) * (phi[ii, 2:] - phi[ii, jj]) / self.dx
        )
        dphi_dy = (
            np.maximum(v_ij, 0) * (phi[ii, jj] - phi[:-2, jj]) / self.dy
            + np.minimum(v_ij, 0) * (phi[2:, jj] - phi[ii, jj]) / self.dy
        )

        # --- 2nd-order upwind for inner cells [2:-2, 2:-2] of the grid ---
        # (needs 2 upwind neighbours; shape checks guard small test grids)
        ny, nx = phi.shape
        if nx >= 5 and ny >= 5:
            u_in = u[2:-2, 2:-2]
            v_in = v[2:-2, 2:-2]
            # x: phi at j, j±1, j±2
            dphi_dx_in = (
                np.maximum(u_in, 0)
                * (3 * phi[2:-2, 2:-2] - 4 * phi[2:-2, 1:-3] + phi[2:-2, :-4])
                / (2.0 * self.dx)
                + np.minimum(u_in, 0)
                * (-phi[2:-2, 4:] + 4 * phi[2:-2, 3:-1] - 3 * phi[2:-2, 2:-2])
                / (2.0 * self.dx)
            )
            # y: phi at i, i±1, i±2
            dphi_dy_in = (
                np.maximum(v_in, 0)
                * (3 * phi[2:-2, 2:-2] - 4 * phi[1:-3, 2:-2] + phi[:-4, 2:-2])
                / (2.0 * self.dy)
                + np.minimum(v_in, 0)
                * (-phi[4:, 2:-2] + 4 * phi[3:-1, 2:-2] - 3 * phi[2:-2, 2:-2])
                / (2.0 * self.dy)
            )
            # Overwrite the centre of the derivative arrays
            dphi_dx[1:-1, 1:-1] = dphi_dx_in
            dphi_dy[1:-1, 1:-1] = dphi_dy_in

        phi_new[ii, jj] = phi[ii, jj] - self.dt * (dphi_dx + dphi_dy)
        phi_new[self.mask] = 0.0
        return phi_new

    def _advect(self, phi, u, v):
        """Return the advection term u·∂φ/∂x + v·∂φ/∂y for interior cells [1:-1,1:-1].

        Dispatches to the scheme selected at construction time.
        Returns an array of shape (ny-2, nx-2) — no intermediate phi_new allocation.
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
            # Override inner cells with 2nd-order one-sided stencils
            ny_phi, nx_phi = phi.shape
            if nx_phi >= 5 and ny_phi >= 5:
                u_in = u[2:-2, 2:-2]
                v_in = v[2:-2, 2:-2]
                dphi_dx[1:-1, 1:-1] = (
                    np.maximum(u_in, 0)
                    * (3 * phi[2:-2, 2:-2] - 4 * phi[2:-2, 1:-3] + phi[2:-2, :-4])
                    / (2.0 * self.dx)
                    + np.minimum(u_in, 0)
                    * (-phi[2:-2, 4:] + 4 * phi[2:-2, 3:-1] - 3 * phi[2:-2, 2:-2])
                    / (2.0 * self.dx)
                )
                dphi_dy[1:-1, 1:-1] = (
                    np.maximum(v_in, 0)
                    * (3 * phi[2:-2, 2:-2] - 4 * phi[1:-3, 2:-2] + phi[:-4, 2:-2])
                    / (2.0 * self.dy)
                    + np.minimum(v_in, 0)
                    * (-phi[4:, 2:-2] + 4 * phi[3:-1, 2:-2] - 3 * phi[2:-2, 2:-2])
                    / (2.0 * self.dy)
                )

        return dphi_dx + dphi_dy

    def diffusion_central(self, phi):
        """2nd-order central difference for diffusion - vectorized."""
        phi_new = phi.copy()
        ii, jj = slice(1, -1), slice(1, -1)

        d2phi_dx2 = (phi[ii, 2:] - 2 * phi[ii, jj] + phi[ii, :-2]) / self.dx**2
        d2phi_dy2 = (phi[2:, jj] - 2 * phi[ii, jj] + phi[:-2, jj]) / self.dy**2

        phi_new[ii, jj] = phi[ii, jj] + self.dt * self.nu * (d2phi_dx2 + d2phi_dy2)
        phi_new[self.mask] = 0.0
        return phi_new

    def step(self):
        """Perform one time step using the fractional-step (projection) method.

        Algorithm:
          1. Predict u* with advection + diffusion (no pressure term).
          2. Solve pressure Poisson: ∇²p = (ρ/dt) ∇·u*
          3. Correct velocity: u^{n+1} = u* - (dt/ρ) ∇p  → divergence-free
        """
        self._compute_optimal_dt()
        self.set_boundary_conditions()

        ii, jj = slice(1, -1), slice(1, -1)

        # Viscous (diffusion) terms
        lap_u = (
            (self.u[2:, jj] - 2 * self.u[ii, jj] + self.u[:-2, jj]) / self.dy**2
            + (self.u[ii, 2:] - 2 * self.u[ii, jj] + self.u[ii, :-2]) / self.dx**2
        )
        lap_v = (
            (self.v[2:, jj] - 2 * self.v[ii, jj] + self.v[:-2, jj]) / self.dy**2
            + (self.v[ii, 2:] - 2 * self.v[ii, jj] + self.v[ii, :-2]) / self.dx**2
        )

        # Advection — dispatched through _advect() (upwind1 or upwind2)
        u_adv = self._advect(self.u, self.u, self.v)
        v_adv = self._advect(self.v, self.u, self.v)

        # --- Step 1: Predict velocity u* (advection + diffusion, no pressure) ---
        u_star = self.u.copy()
        v_star = self.v.copy()
        u_star[ii, jj] = np.clip(
            self.u[ii, jj] + self.dt * (-u_adv + self.nu * lap_u),
            -5 * self.U_inf, 5 * self.U_inf,
        )
        v_star[ii, jj] = np.clip(
            self.v[ii, jj] + self.dt * (-v_adv + self.nu * lap_v),
            -5 * self.U_inf, 5 * self.U_inf,
        )
        u_star[self.mask] = 0.0
        v_star[self.mask] = 0.0
        # Apply BCs to predicted velocity
        u_star[:, 0] = self.U_inf
        v_star[:, 0] = 0.0
        u_star[:, -1] = u_star[:, -2]
        v_star[:, -1] = v_star[:, -2]
        u_star[0, :] = u_star[1, :]
        u_star[-1, :] = u_star[-2, :]
        v_star[0, :] = 0.0
        v_star[-1, :] = 0.0

        # --- Step 2: Solve ∇²p = (ρ/dt) ∇·u* ---
        self.solve_poisson_pressure(u=u_star, v=v_star)

        # --- Step 3: Correct velocity with pressure gradient ---
        dpdx = (self.p[ii, 2:] - self.p[ii, :-2]) / (2 * self.dx)
        dpdy = (self.p[2:, jj] - self.p[:-2, jj]) / (2 * self.dy)

        u_new = u_star.copy()
        v_new = v_star.copy()
        u_new[ii, jj] = np.clip(
            u_star[ii, jj] - self.dt / self.rho * dpdx,
            -5 * self.U_inf, 5 * self.U_inf,
        )
        v_new[ii, jj] = np.clip(
            v_star[ii, jj] - self.dt / self.rho * dpdy,
            -5 * self.U_inf, 5 * self.U_inf,
        )

        if np.any(np.isnan(u_new)) or np.any(np.isnan(v_new)):
            logger.error("Instability detected (NaN). Aborting step.")
            raise RuntimeError("Solver diverged.")

        self.u, self.v = u_new, v_new
        self.u[self.mask] = 0.0
        self.v[self.mask] = 0.0
        self.set_boundary_conditions()

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
        Compute drag coefficient using pressure forces on cylinder surface.

        Note: This is a simplified calculation for educational purposes.
        """
        # Surface cells: cylinder cells with at least one fluid neighbor
        neighbor_fluid = (
            ~np.roll(self.mask, 1, axis=0)
            | ~np.roll(self.mask, -1, axis=0)
            | ~np.roll(self.mask, 1, axis=1)
            | ~np.roll(self.mask, -1, axis=1)
        )
        surface_mask = self.mask & neighbor_fluid

        # Pressure drag: p * n_x * dA, where n_x = outward normal from cylinder into fluid
        # Use neighboring fluid cell pressure; normal direction from cylinder center to cell
        surf_indices = np.argwhere(surface_mask)
        F_p = 0.0
        for i, j in surf_indices:
            # Average pressure from fluid neighbors
            p_neighbors = []
            for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                ni, nj = i + di, j + dj
                if 0 <= ni < self.ny and 0 <= nj < self.nx and not self.mask[ni, nj]:
                    p_neighbors.append(self.p[ni, nj])
            if p_neighbors:
                p_avg = np.mean(p_neighbors)
                # x-component of outward normal (cylinder center → surface cell)
                dx_vec = self.X[i, j] - self.x_c
                r = np.sqrt((self.X[i, j] - self.x_c)**2 + (self.Y[i, j] - self.y_c)**2)
                nx = dx_vec / (r + 1e-10)
                # Drag force contribution: pressure pushes inward (-n̂), x-component gives drag
                F_p += p_avg * (-nx) * self.dx * self.dy

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

        # Pressure field and streamlines
        p_masked = np.ma.array(self.p, mask=self.mask)
        im2 = ax2.contourf(self.X, self.Y, p_masked, levels=20, cmap="RdBu_r")
        ax2.streamplot(
            self.x,
            self.y,
            self.u,
            self.v,
            color="black",
            linewidth=0.5,
            density=1.5,
            arrowsize=0.8,
        )
        circle2 = plt.Circle((self.x_c, self.y_c), self.R, color="gray", fill=True)
        ax2.add_patch(circle2)
        ax2.set_xlabel("x [m]")
        ax2.set_ylabel("y [m]")
        ax2.set_title("Pressure Field and Streamlines")
        ax2.set_aspect("equal")
        plt.colorbar(im2, ax=ax2, label="Pressure [Pa]")

        plt.tight_layout()
        return fig

    def get_velocity_on_centerline(self):
        """Get velocity profile along horizontal centerline."""
        i_center = self.ny // 2
        u_centerline = self.u[i_center, :]
        x_centerline = self.x
        return x_centerline, u_centerline
