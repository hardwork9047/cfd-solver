"""
Particle System Module for DEM

This module implements a discrete element method (DEM) simulation for
granular materials using soft-sphere collision model.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.patches import Circle

# Configuration for logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Font setting
rcParams["font.family"] = "sans-serif"


class ParticleSystem:
    """
    Discrete Element Method (DEM) Particle System.
    
    Simulates the dynamics of spherical particles using:
    - Gravity
    - Normal contact forces (Hertzian contact)
    - Damping
    - Particle-particle and particle-wall collisions
    """

    def __init__(
        self,
        n_particles=100,
        domain_size=(10.0, 10.0),
        particle_radius=0.2,
        particle_density=2500,
        gravity=9.81,
        k_n=1e5,
        damping=0.3,
        dt=1e-4,
    ):
        """
        Parameters:
        -----------
        n_particles : int
            Number of particles
        domain_size : tuple
            (width, height) of the domain [m]
        particle_radius : float
            Radius of particles [m]
        particle_density : float
            Density of particles [kg/m³]
        gravity : float
            Gravitational acceleration [m/s²]
        k_n : float
            Normal contact stiffness [N/m]
        damping : float
            Damping coefficient (0-1)
        dt : float
            Time step [s]
        """
        self.n_particles = n_particles
        self.width, self.height = domain_size
        self.radius = particle_radius
        self.density = particle_density
        self.gravity = gravity
        self.k_n = k_n
        self.damping = damping
        self.dt = dt

        # Particle mass (assuming 2D circular particles with unit thickness)
        self.mass = particle_density * np.pi * particle_radius**2
        
        # Initialize radii and masses arrays (will be uniform initially, can be modified)
        self.radii = np.full(n_particles, particle_radius)
        self.masses = np.full(n_particles, self.mass)

        # Initialize particle positions (random distribution in upper half)
        self.positions = np.random.rand(n_particles, 2) * np.array(
            [self.width, self.height * 0.8]
        ) + np.array([0, self.height * 0.2])

        # Initialize velocities (zero initially)
        self.velocities = np.zeros((n_particles, 2))

        # Simulation time
        self.time = 0.0
        self.iteration = 0
        
        # Store forces for visualization
        self.forces = np.zeros((n_particles, 2))
        
        # Spatial partitioning grid (for efficient collision detection)
        # Set cell size to 4x the particle radius to capture nearby particles
        self.cell_size = max(4.0 * particle_radius, 0.1)
        self.grid_nx = int(np.ceil(self.width / self.cell_size)) + 1
        self.grid_ny = int(np.ceil(self.height / self.cell_size)) + 1

        logger.info(
            f"Initialized DEM ParticleSystem: {n_particles} particles, "
            f"domain {domain_size[0]}x{domain_size[1]}m, "
            f"r={particle_radius}m, m={self.mass:.3e}kg, "
            f"grid={self.grid_nx}x{self.grid_ny}"
        )

    def _build_spatial_grid(self):
        """Build spatial grid for efficient collision detection."""
        # Initialize grid with empty lists
        grid = {}
        
        # Assign each particle to grid cell(s)
        for i in range(self.n_particles):
            x, y = self.positions[i]
            r = self.radii[i]
            
            # Determine which grid cells this particle overlaps
            min_cell_x = max(0, int(np.floor((x - r) / self.cell_size)))
            max_cell_x = min(self.grid_nx - 1, int(np.floor((x + r) / self.cell_size)))
            min_cell_y = max(0, int(np.floor((y - r) / self.cell_size)))
            max_cell_y = min(self.grid_ny - 1, int(np.floor((y + r) / self.cell_size)))
            
            for cell_x in range(min_cell_x, max_cell_x + 1):
                for cell_y in range(min_cell_y, max_cell_y + 1):
                    cell_key = (cell_x, cell_y)
                    if cell_key not in grid:
                        grid[cell_key] = []
                    grid[cell_key].append(i)
        
        return grid

    def compute_forces(self):
        """Compute all forces acting on particles using spatial grid optimization."""
        forces = np.zeros((self.n_particles, 2))

        # Gravity (using individual masses)
        forces[:, 1] -= self.masses * self.gravity

        # Build spatial grid for efficient collision detection
        grid = self._build_spatial_grid()
        
        # Track which particle pairs have been checked
        checked_pairs = set()

        # Particle-particle collisions (using spatial grid)
        for cell_key, particle_indices in grid.items():
            # Check collisions within this cell
            for i_idx, i in enumerate(particle_indices):
                for j in particle_indices[i_idx + 1:]:
                    pair_key = tuple(sorted([i, j]))
                    if pair_key in checked_pairs:
                        continue
                    checked_pairs.add(pair_key)
                    
                    delta_pos = self.positions[j] - self.positions[i]
                    distance = np.linalg.norm(delta_pos)
                    contact_distance = self.radii[i] + self.radii[j]
                    
                    # Quick check: only compute if particles are close enough
                    if distance >= contact_distance:
                        continue
                    
                    overlap = contact_distance - distance

                    if overlap > 0:  # Collision detected
                        # Normal direction
                        normal = delta_pos / (distance + 1e-10)

                        # Normal force (Hertzian contact model: f_n ∝ overlap^(3/2))
                        f_n = self.k_n * (overlap ** 1.5)

                        # Relative velocity
                        rel_vel = self.velocities[j] - self.velocities[i]
                        v_n = np.dot(rel_vel, normal)

                        # Damping force (using average mass)
                        avg_mass = (self.masses[i] + self.masses[j]) / 2
                        f_d = -self.damping * v_n * np.sqrt(self.k_n * avg_mass)

                        # Total force
                        f_total = (f_n + f_d) * normal

                        forces[i] -= f_total
                        forces[j] += f_total
            
            # Check collisions with neighboring cells
            cell_x, cell_y = cell_key
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    if dx == 0 and dy == 0:
                        continue  # Already checked this cell
                    
                    neighbor_key = (cell_x + dx, cell_y + dy)
                    if neighbor_key not in grid:
                        continue
                    
                    neighbor_indices = grid[neighbor_key]
                    
                    for i in particle_indices:
                        for j in neighbor_indices:
                            if i >= j:
                                continue  # Avoid duplicate checks
                            
                            pair_key = tuple(sorted([i, j]))
                            if pair_key in checked_pairs:
                                continue
                            checked_pairs.add(pair_key)
                            
                            delta_pos = self.positions[j] - self.positions[i]
                            distance = np.linalg.norm(delta_pos)
                            contact_distance = self.radii[i] + self.radii[j]
                            
                            # Quick check
                            if distance >= contact_distance:
                                continue
                            
                            overlap = contact_distance - distance

                            if overlap > 0:  # Collision detected
                                # Normal direction
                                normal = delta_pos / (distance + 1e-10)

                                # Normal force
                                f_n = self.k_n * (overlap ** 1.5)

                                # Relative velocity
                                rel_vel = self.velocities[j] - self.velocities[i]
                                v_n = np.dot(rel_vel, normal)

                                # Damping force
                                avg_mass = (self.masses[i] + self.masses[j]) / 2
                                f_d = -self.damping * v_n * np.sqrt(self.k_n * avg_mass)

                                # Total force
                                f_total = (f_n + f_d) * normal

                                forces[i] -= f_total
                                forces[j] += f_total


        # Wall collisions (using individual radii)
        for i in range(self.n_particles):
            # Bottom wall
            if self.positions[i, 1] < self.radii[i]:
                overlap = self.radii[i] - self.positions[i, 1]
                f_n = self.k_n * (overlap ** 1.5)  # Hertzian contact
                # Damping opposes normal velocity (negative velocity = moving into wall)
                v_n = self.velocities[i, 1]  # Normal velocity (negative = into wall)
                f_d = -self.damping * v_n * np.sqrt(self.k_n * self.masses[i])
                forces[i, 1] += f_n + f_d

            # Top wall
            if self.positions[i, 1] > self.height - self.radii[i]:
                overlap = self.positions[i, 1] - (self.height - self.radii[i])
                f_n = self.k_n * (overlap ** 1.5)  # Hertzian contact
                # Damping opposes normal velocity (positive velocity = moving into wall)
                v_n = self.velocities[i, 1]  # Normal velocity (positive = into wall)
                f_d = -self.damping * v_n * np.sqrt(self.k_n * self.masses[i])
                forces[i, 1] -= f_n + f_d

            # Left wall
            if self.positions[i, 0] < self.radii[i]:
                overlap = self.radii[i] - self.positions[i, 0]
                f_n = self.k_n * (overlap ** 1.5)  # Hertzian contact
                # Damping opposes normal velocity (negative velocity = moving into wall)
                v_n = self.velocities[i, 0]  # Normal velocity (negative = into wall)
                f_d = -self.damping * v_n * np.sqrt(self.k_n * self.masses[i])
                forces[i, 0] += f_n + f_d

            # Right wall
            if self.positions[i, 0] > self.width - self.radii[i]:
                overlap = self.positions[i, 0] - (self.width - self.radii[i])
                f_n = self.k_n * (overlap ** 1.5)  # Hertzian contact
                # Damping opposes normal velocity (positive velocity = moving into wall)
                v_n = self.velocities[i, 0]  # Normal velocity (positive = into wall)
                f_d = -self.damping * v_n * np.sqrt(self.k_n * self.masses[i])
                forces[i, 0] -= f_n + f_d

        return forces

    def step(self):
        """Perform one time step using Velocity Verlet integration."""
        # Compute current forces
        forces = self.compute_forces()
        accelerations = forces / self.masses[:, np.newaxis]  # Use individual masses
        
        # Store forces for visualization
        self.forces = forces.copy()

        # Update velocities (half-step)
        self.velocities += 0.5 * accelerations * self.dt

        # Update positions
        self.positions += self.velocities * self.dt

        # Compute new forces
        forces_new = self.compute_forces()
        accelerations_new = forces_new / self.masses[:, np.newaxis]  # Use individual masses

        # Update velocities (second half-step)
        self.velocities += 0.5 * accelerations_new * self.dt

        # Enforce boundary conditions (keep particles strictly inside domain)
        # Also reset velocities when particles hit walls
        for i in range(self.n_particles):
            # Left wall
            if self.positions[i, 0] < self.radii[i]:
                self.positions[i, 0] = self.radii[i]
                if self.velocities[i, 0] < 0:
                    self.velocities[i, 0] = 0.0
            
            # Right wall
            if self.positions[i, 0] > self.width - self.radii[i]:
                self.positions[i, 0] = self.width - self.radii[i]
                if self.velocities[i, 0] > 0:
                    self.velocities[i, 0] = 0.0
            
            # Bottom wall
            if self.positions[i, 1] < self.radii[i]:
                self.positions[i, 1] = self.radii[i]
                if self.velocities[i, 1] < 0:
                    self.velocities[i, 1] = 0.0
            
            # Top wall
            if self.positions[i, 1] > self.height - self.radii[i]:
                self.positions[i, 1] = self.height - self.radii[i]
                if self.velocities[i, 1] > 0:
                    self.velocities[i, 1] = 0.0

        self.time += self.dt
        self.iteration += 1

    def run(self, t_end=10.0, output_interval=1000):
        """Run simulation until t_end."""
        logger.info(f"Starting DEM simulation until t={t_end}s")

        while self.time < t_end:
            self.step()

            if self.iteration % output_interval == 0:
                ke = self.compute_kinetic_energy()
                logger.info(
                    f"Iteration {self.iteration}, Time {self.time:.4f}s, KE={ke:.2e}J, "
                    f"Avg velocity={np.mean(np.linalg.norm(self.velocities, axis=1)):.3f}m/s"
                )

        logger.info(f"Simulation complete at t={self.time:.4f}s")

    def compute_kinetic_energy(self):
        """Compute total kinetic energy of the system."""
        v_mag = np.linalg.norm(self.velocities, axis=1)
        return 0.5 * np.sum(self.masses * v_mag**2)

    def plot(self, figsize=(10, 12), save_path=None):
        """Plot current particle positions with color based on force magnitude."""
        fig, ax = plt.subplots(figsize=figsize)

        # Draw domain boundaries
        ax.plot([0, self.width], [0, 0], "k-", linewidth=2)
        ax.plot([0, self.width], [self.height, self.height], "k-", linewidth=2)
        ax.plot([0, 0], [0, self.height], "k-", linewidth=2)
        ax.plot([self.width, self.width], [0, self.height], "k-", linewidth=2)

        # Compute force magnitudes for coloring
        force_magnitudes = np.linalg.norm(self.forces, axis=1)
        
        # Normalize force magnitudes for colormap (avoid division by zero)
        max_force = np.max(force_magnitudes) if np.max(force_magnitudes) > 0 else 1.0
        normalized_forces = force_magnitudes / max_force
        
        # Use colormap (blue for low force, red for high force)
        cmap = plt.cm.coolwarm
        
        # Draw particles with color based on force (using individual radii)
        for i in range(self.n_particles):
            color = cmap(normalized_forces[i])
            circle = Circle(
                self.positions[i],
                self.radii[i],  # Use individual particle radius
                facecolor=color,
                edgecolor="black",
                linewidth=0.5,
                alpha=0.8,
            )
            ax.add_patch(circle)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_force))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, label="Force Magnitude [N]", pad=0.02)

        ax.set_xlim(-0.5, self.width + 0.5)
        ax.set_ylim(-0.5, self.height + 0.5)
        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(f"DEM Particle System at t={self.time:.2f}s")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            logger.info(f"Plot saved to: {save_path}")

        return fig

    def get_particle_distribution(self, n_bins=20):
        """Get vertical distribution of particles."""
        hist, bin_edges = np.histogram(
            self.positions[:, 1], bins=n_bins, range=(0, self.height)
        )
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        return bin_centers, hist
