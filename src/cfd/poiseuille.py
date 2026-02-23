"""
Poiseuille Flow Calculation Module

This module calculates Plane and Circular Poiseuille flows using both 
analytical and numerical (Finite Difference Method with SOR) solutions.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Configuration for logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Font setting for compatibility
rcParams['font.family'] = 'sans-serif'


class PlanePoiseuille:
    """Plane Poiseuille Flow between two parallel plates."""
    
    def __init__(self, L, dp_dx, mu, ny=101):
        """
        Parameters:
        -----------
        L : float
            Distance between plates [m]
        dp_dx : float
            Pressure gradient dp/dx [Pa/m] (should be negative)
        mu : float
            Dynamic viscosity [Pa·s]
        ny : int
            Number of grid points in the y-direction
        """
        if L <= 0 or mu <= 0:
            raise ValueError("L and mu must be positive values.")
        
        self.L = L
        self.dp_dx = dp_dx
        self.mu = mu
        self.ny = ny
        self.y = np.linspace(-L/2, L/2, ny)
        logger.info(f"Initialized PlanePoiseuille: L={L}, dp/dx={dp_dx}, mu={mu}, ny={ny}")
        
    def analytical_solution(self):
        """Calculate the analytical solution."""
        return -(self.dp_dx / (2 * self.mu)) * (self.L**2 / 4 - self.y**2)
    
    def numerical_solution(self, tol=1e-6, max_iter=20000, omega=1.5):
        """
        Numerical solution using Finite Difference Method with SOR.
        
        Parameters:
        -----------
        tol : float
            Convergence tolerance
        max_iter : int
            Maximum number of iterations
        omega : float
            Over-relaxation parameter (1.0 < omega < 2.0)
        """
        u = np.zeros(self.ny)
        dy = self.y[1] - self.y[0]
        
        # Boundary conditions: no-slip at walls
        u[0] = 0
        u[-1] = 0
        
        logger.info(f"Starting SOR solver for PlanePoiseuille (omega={omega}, tol={tol})")
        
        for iteration in range(max_iter):
            u_old = u.copy()
            
            for i in range(1, self.ny - 1):
                # Gauss-Seidel update value
                u_gs = 0.5 * (u[i-1] + u[i+1] - dy**2 * self.dp_dx / self.mu)
                # SOR update
                u[i] = u_old[i] + omega * (u_gs - u_old[i])
            
            # Check for instability
            if np.any(np.isnan(u)) or np.any(np.isinf(u)):
                logger.error(f"Numerical instability detected at iteration {iteration+1}!")
                raise RuntimeError("Solver diverged (NaN or Inf detected). Try reducing omega or increasing resolution.")
            
            # Check for convergence
            residual = np.max(np.abs(u - u_old))
            if residual < tol:
                logger.info(f"Converged in {iteration+1} iterations. Final residual: {residual:.2e}")
                return u
        
        logger.warning(f"Failed to converge within {max_iter} iterations.")
        return u
    
    def get_max_velocity(self):
        """Return the maximum velocity (at the center)."""
        return -self.dp_dx * self.L**2 / (8 * self.mu)
    
    def get_flow_rate(self):
        """Return the flow rate per unit width."""
        return -self.dp_dx * self.L**3 / (12 * self.mu)
    
    def plot(self):
        """Plot the velocity distribution."""
        u_analytical = self.analytical_solution()
        u_numerical = self.numerical_solution()
        
        plt.figure(figsize=(10, 6))
        plt.plot(u_analytical, self.y, 'b-', linewidth=2, label='Analytical')
        plt.plot(u_numerical, self.y, 'ro', markersize=4, 
                markerfacecolor='none', label='Numerical')
        plt.xlabel('Velocity u [m/s]', fontsize=12)
        plt.ylabel('y [m]', fontsize=12)
        plt.title('Plane Poiseuille Flow', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=11)
        plt.tight_layout()
        
        return plt.gcf()


class CircularPoiseuille:
    """Circular Poiseuille Flow (Hagen-Poiseuille Flow) in a pipe."""
    
    def __init__(self, R, dp_dx, mu, nr=101):
        """
        Parameters:
        -----------
        R : float
            Pipe radius [m]
        dp_dx : float
            Pressure gradient dp/dx [Pa/m] (should be negative)
        mu : float
            Dynamic viscosity [Pa·s]
        nr : int
            Number of grid points in the radial direction
        """
        if R <= 0 or mu <= 0:
            raise ValueError("R and mu must be positive values.")
            
        self.R = R
        self.dp_dx = dp_dx
        self.mu = mu
        self.nr = nr
        self.r = np.linspace(0, R, nr)
        logger.info(f"Initialized CircularPoiseuille: R={R}, dp/dx={dp_dx}, mu={mu}, nr={nr}")
        
    def analytical_solution(self):
        """Calculate the analytical solution."""
        return -(self.dp_dx / (4 * self.mu)) * (self.R**2 - self.r**2)
    
    def numerical_solution(self, tol=1e-6, max_iter=20000, omega=1.5):
        """
        Numerical solution using Finite Difference Method with SOR in cylindrical coordinates.
        """
        u = np.zeros(self.nr)
        dr = self.r[1] - self.r[0]
        
        # Boundary condition: no-slip at wall (r=R)
        u[-1] = 0
        
        logger.info(f"Starting SOR solver for CircularPoiseuille (omega={omega}, tol={tol})")
        
        for iteration in range(max_iter):
            u_old = u.copy()
            
            # Center point (r=0): Symmetry condition du/dr = 0
            # 2*d²u/dr² = dp/dx / mu
            u_gs_0 = u[1] - dr**2 * self.dp_dx / (2 * self.mu)
            u[0] = u_old[0] + omega * (u_gs_0 - u_old[0])
            
            # Internal points
            for i in range(1, self.nr - 1):
                r_i = self.r[i]
                u_gs_i = (u[i+1] * (r_i + dr/2) + u[i-1] * (r_i - dr/2)) / (2 * r_i)
                u_gs_i -= dr**2 * self.dp_dx / (2 * self.mu)
                u[i] = u_old[i] + omega * (u_gs_i - u_old[i])
            
            # Check for instability
            if np.any(np.isnan(u)) or np.any(np.isinf(u)):
                logger.error(f"Numerical instability detected at iteration {iteration+1}!")
                raise RuntimeError("Solver diverged (NaN or Inf detected).")
            
            # Check for convergence
            residual = np.max(np.abs(u - u_old))
            if residual < tol:
                logger.info(f"Converged in {iteration+1} iterations. Final residual: {residual:.2e}")
                return u
        
        logger.warning(f"Failed to converge within {max_iter} iterations.")
        return u
    
    def get_max_velocity(self):
        """Return the maximum velocity (at the center)."""
        return -self.dp_dx * self.R**2 / (4 * self.mu)
    
    def get_flow_rate(self):
        """Return the flow rate (Hagen-Poiseuille equation)."""
        return -np.pi * self.dp_dx * self.R**4 / (8 * self.mu)
    
    def plot(self):
        """Plot the velocity distribution."""
        u_analytical = self.analytical_solution()
        u_numerical = self.numerical_solution()
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # 1D Plot
        ax1.plot(self.r, u_analytical, 'b-', linewidth=2, label='Analytical')
        ax1.plot(self.r, u_numerical, 'ro', markersize=4, 
                markerfacecolor='none', label='Numerical')
        ax1.set_xlabel('Radial position r [m]', fontsize=12)
        ax1.set_ylabel('Velocity u [m/s]', fontsize=12)
        ax1.set_title('Velocity Profile', fontsize=13, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=11)
        
        # 2D Plot (Cross section)
        theta = np.linspace(0, 2*np.pi, 100)
        R_mesh, Theta_mesh = np.meshgrid(self.r, theta)
        U_mesh = -(self.dp_dx / (4 * self.mu)) * (self.R**2 - R_mesh**2)
        
        X = R_mesh * np.cos(Theta_mesh)
        Y = R_mesh * np.sin(Theta_mesh)
        
        contour = ax2.contourf(X, Y, U_mesh, levels=20, cmap='jet')
        ax2.set_xlabel('x [m]', fontsize=12)
        ax2.set_ylabel('y [m]', fontsize=12)
        ax2.set_title('Velocity Field (Cross Section)', fontsize=13, fontweight='bold')
        ax2.set_aspect('equal')
        plt.colorbar(contour, ax=ax2, label='Velocity [m/s]')
        
        plt.tight_layout()
        return fig
