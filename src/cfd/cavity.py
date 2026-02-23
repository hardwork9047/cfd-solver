"""
Cavity Flow Calculation Module

This module simulates the Lid-Driven Cavity Flow using a stabilized 
numerical solver for the incompressible Navier-Stokes equations.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Configuration for logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Japanese font setting
rcParams['font.family'] = 'sans-serif'


class CavityFlow:
    """Lid-Driven Cavity Flow Solver with Numerical Stabilization."""
    
    def __init__(self, L=1.0, rho=1.0, mu=0.01, nx=65, ny=65):
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
        """
        self.L = L
        self.rho = rho
        self.mu = mu
        self.nu = mu / rho
        self.nx = nx
        self.ny = ny
        
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
        
        logger.info(f"Initialized CavityFlow: Re={self.get_Reynolds_number():.1f}, Grid={nx}x{ny}")

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

    def solve_poisson_pressure(self, max_iter=50, tol=1e-4):
        """Pressure Poisson Equation solver with SOR and clipping."""
        p_new = self.p.copy()
        
        # RHS = -rho/dt * (du/dx + dv/dy)
        dudx = (self.u[1:-1, 2:] - self.u[1:-1, :-2]) / (2 * self.dx)
        dvdy = (self.v[2:, 1:-1] - self.v[:-2, 1:-1]) / (2 * self.dy)
        rhs = -self.rho * np.clip(dudx + dvdy, -10, 10) / (self.dt + 1e-10)
        
        omega_sor = 1.4
        coeff = 2.0 / (self.dx**2) + 2.0 / (self.dy**2)
        
        for _ in range(max_iter):
            p_old = p_new.copy()
            
            # Vectorized Laplace-like update (partially)
            lap_x = (p_new[1:-1, 2:] + p_new[1:-1, :-2]) / (self.dx**2)
            lap_y = (p_new[2:, 1:-1] + p_new[:-2, 1:-1]) / (self.dy**2)
            
            # Local residual clipping for stability
            res = np.clip(rhs - (lap_x + lap_y - coeff * p_new[1:-1, 1:-1]), -1e3, 1e3)
            p_new[1:-1, 1:-1] = p_new[1:-1, 1:-1] + omega_sor * res / coeff
            
            # Neumann BCs
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

    def update_velocity(self):
        """Stabilized velocity update using Upwind differencing and clipping."""
        u_new = self.u.copy()
        v_new = self.v.copy()
        ii, jj = slice(1, -1), slice(1, -1)
        
        # Viscous terms
        lap_u = (self.u[2:, jj] - 2*self.u[ii, jj] + self.u[:-2, jj]) / self.dy**2 + \
                (self.u[ii, 2:] - 2*self.u[ii, jj] + self.u[ii, :-2]) / self.dx**2
        lap_v = (self.v[2:, jj] - 2*self.v[ii, jj] + self.v[:-2, jj]) / self.dy**2 + \
                (self.v[ii, 2:] - 2*self.v[ii, jj] + self.v[ii, :-2]) / self.dx**2
                
        # Advection (Upwind)
        u_pos = np.maximum(self.u[ii, jj], 0)
        u_neg = np.minimum(self.u[ii, jj], 0)
        v_pos = np.maximum(self.v[ii, jj], 0)
        v_neg = np.minimum(self.v[ii, jj], 0)
        
        u_adv = u_pos * (self.u[ii, jj] - self.u[ii, :-2]) / self.dx + \
                u_neg * (self.u[ii, 2:] - self.u[ii, jj]) / self.dx + \
                v_pos * (self.u[ii, jj] - self.u[:-2, jj]) / self.dy + \
                v_neg * (self.u[2:, jj] - self.u[ii, jj]) / self.dy
                
        v_adv = u_pos * (self.v[ii, jj] - self.v[ii, :-2]) / self.dx + \
                u_neg * (self.v[ii, 2:] - self.v[ii, jj]) / self.dx + \
                v_pos * (self.v[ii, jj] - self.v[:-2, jj]) / self.dy + \
                v_neg * (self.v[2:, jj] - self.v[ii, jj]) / self.dy

        # Pressure Gradient
        dpdx = (self.p[ii, 2:] - self.p[ii, :-2]) / (2 * self.dx)
        dpdy = (self.p[2:, jj] - self.p[:-2, jj]) / (2 * self.dy)
        
        # Update with clipping for physical limits
        u_new[ii, jj] = np.clip(self.u[ii, jj] + self.dt * (-u_adv - dpdx/self.rho + self.nu*lap_u), -2, 2)
        v_new[ii, jj] = np.clip(self.v[ii, jj] + self.dt * (-v_adv - dpdy/self.rho + self.nu*lap_v), -2, 2)
        
        if np.any(np.isnan(u_new)):
            logger.error("Instability detected.")
            raise RuntimeError("Solver diverged.")
            
        self.u, self.v = u_new, v_new

    def step(self):
        self.set_boundary_conditions()
        self._compute_optimal_dt()
        self.update_velocity()
        self.solve_poisson_pressure()
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
            if verbose and (i+1) % 100 == 0:
                logger.info(f"Step {i+1}: t={self.time:.4f}, residual={res:.2e}")
            if res < 1e-6 and i > 50:
                logger.info("Steady state reached.")
                break

    def plot_velocity_field(self, figsize=(10, 8)):
        fig, ax = plt.subplots(figsize=figsize)
        mag = np.sqrt(self.u**2 + self.v**2)
        cf = ax.contourf(self.X, self.Y, mag, levels=20, cmap='viridis')
        ax.streamplot(self.x, self.y, self.u, self.v, color='white', density=1.0)
        plt.colorbar(cf, label='Velocity [m/s]')
        return fig
