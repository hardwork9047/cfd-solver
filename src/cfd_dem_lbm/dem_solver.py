"""DEM load solver used by LBM-DEM production and verification paths."""

from __future__ import annotations

from typing import Callable

import numpy as np


class DEMSolver:
    """Compute DEM forces and torques for a coupled LBM-DEM simulation.

    The solver owns the DEM contact, surface-force, wall, and cylinder-load
    logic.  It reads particle and geometry arrays from the coupled simulation
    object and asks that object for fluid drag, so production runs and
    verification can share the same DEM implementation without duplicating
    contact physics in scripts.
    """

    def __init__(self, coupled_solver, pair_kernel: Callable | None = None):
        self.sim = coupled_solver
        self.pair_kernel = pair_kernel

    @staticmethod
    def tangent_from_normal(nx_: float, ny_: float) -> np.ndarray:
        """Return a unit tangent for a 2-D contact normal."""
        return np.array([-ny_, nx_])

    def normal_contact_magnitude(self, overlap: float, v_n: float, mass: float) -> float:
        """Hertz normal force with damping, clamped to avoid artificial tension."""
        sim = self.sim
        f_n = sim.k_n * overlap**1.5
        f_damp = -sim.damping * v_n * float(np.sqrt(sim.k_n * mass))
        return max(f_n + f_damp, 0.0)

    def tangential_force_magnitude(
        self,
        v_t: float,
        normal_force: float,
        mass: float,
    ) -> float:
        """Coulomb-limited tangential force opposing slip at a contact."""
        sim = self.sim
        if not sim.rolling_friction or normal_force <= 0.0:
            return 0.0
        f_trial = -sim.tangential_damping * float(np.sqrt(sim.k_n * mass)) * v_t
        f_limit = sim.sliding_friction * normal_force
        return float(np.clip(f_trial, -f_limit, f_limit))

    def rolling_resistance_torque(
        self,
        omega: float,
        normal_force: float,
        radius: float,
        mass: float,
    ) -> float:
        """Coulomb-limited rolling resistance torque opposing angular velocity."""
        sim = self.sim
        if not sim.rolling_friction or normal_force <= 0.0:
            return 0.0
        torque_trial = (
            -sim.rolling_damping
            * float(np.sqrt(sim.k_n * mass))
            * radius**2
            * omega
        )
        torque_limit = sim.rolling_friction_coeff * normal_force * radius
        return float(np.clip(torque_trial, -torque_limit, torque_limit))

    def compute_loads(self, dt_sub: float) -> tuple[np.ndarray, np.ndarray]:
        """Compute total DEM force and torque on every active particle."""
        sim = self.sim
        forces = np.zeros((sim.n_p, 2))
        torques = np.zeros(sim.n_p)

        if sim.n_p == 0:
            return forces, torques

        # 1. Gravity (buoyancy-corrected, per-particle mass)
        buoyancy_factor = 1.0 - 1.0 / sim.density_ratio
        forces[:, 1] -= sim.masses * sim.g * buoyancy_factor

        # 2. Stokes drag from interpolated fluid velocity (per-particle radius)
        forces += sim._particle_drag_forces()

        # 3. Particle-particle Hertz contact and optional near-surface forces.
        pair_i: np.ndarray | None = None
        pair_j: np.ndarray | None = None
        if self.pair_kernel is not None and sim.n_p >= 2:
            pair_i, pair_j = sim._particle_pair_arrays()
            self.pair_kernel(
                pair_i,
                pair_j,
                sim.pos,
                sim.vel,
                sim.radii,
                sim.masses,
                sim.omega_p,
                sim.k_n,
                sim.damping,
                sim.rolling_friction,
                sim.sliding_friction,
                sim.tangential_damping,
                sim.rolling_friction_coeff,
                sim.rolling_damping,
                sim.particle_attraction,
                sim.particle_repulsion,
                sim.attraction_strength,
                sim.repulsion_strength,
                sim.attraction_cutoff,
                sim.repulsion_cutoff,
                sim.attraction_min_gap,
                sim.repulsion_min_gap,
                forces,
                torques,
            )

        if self.pair_kernel is None:
            for i, j in sim._particle_pair_candidates():
                self._apply_particle_pair_loads(i, j, forces, torques)

        self._apply_wall_loads(forces, torques)
        self._apply_cylinder_loads(forces, torques)
        return forces, torques

    def _apply_particle_pair_loads(
        self,
        i: int,
        j: int,
        forces: np.ndarray,
        torques: np.ndarray,
    ) -> None:
        sim = self.sim
        dp = sim.pos[j] - sim.pos[i]
        dist = float(np.linalg.norm(dp))
        min_dist = sim.radii[i] + sim.radii[j]
        if dist <= 1e-10:
            return

        n = dp / dist

        if dist < min_dist:
            overlap = min_dist - dist
            v_n = float(np.dot(sim.vel[j] - sim.vel[i], n))
            avg_m = (sim.masses[i] + sim.masses[j]) / 2.0
            f_mag = self.normal_contact_magnitude(overlap, v_n, avg_m)
            f_total = f_mag * n
            forces[i] -= f_total
            forces[j] += f_total

            t = self.tangent_from_normal(n[0], n[1])
            v_t = float(np.dot(sim.vel[j] - sim.vel[i], t))
            v_t -= sim.omega_p[j] * sim.radii[j] + sim.omega_p[i] * sim.radii[i]
            f_t = self.tangential_force_magnitude(v_t, f_mag, avg_m)
            f_t_vec = f_t * t
            forces[i] -= f_t_vec
            forces[j] += f_t_vec
            torques[i] -= sim.radii[i] * f_t
            torques[j] -= sim.radii[j] * f_t
            torques[i] += self.rolling_resistance_torque(
                sim.omega_p[i], f_mag, sim.radii[i], sim.masses[i]
            )
            torques[j] += self.rolling_resistance_torque(
                sim.omega_p[j], f_mag, sim.radii[j], sim.masses[j]
            )

        self._apply_surface_force(i, j, dist, min_dist, n, forces)

    def _apply_surface_force(
        self,
        i: int,
        j: int,
        dist: float,
        min_dist: float,
        n: np.ndarray,
        forces: np.ndarray,
    ) -> None:
        sim = self.sim
        surface_gap = dist - min_dist
        if sim.particle_attraction and sim.attraction_strength > 0.0:
            if surface_gap <= sim.attraction_cutoff:
                h = max(surface_gap, max(sim.attraction_min_gap, 1e-12))
                r_eff = sim.radii[i] * sim.radii[j] / min_dist
                f_attr = sim.attraction_strength * r_eff / (6.0 * h**2)
                forces[i] += f_attr * n
                forces[j] -= f_attr * n
        elif sim.particle_repulsion and sim.repulsion_strength > 0.0:
            if surface_gap <= sim.repulsion_cutoff:
                h = max(surface_gap, max(sim.repulsion_min_gap, 1e-12))
                r_eff = sim.radii[i] * sim.radii[j] / min_dist
                f_rep = sim.repulsion_strength * r_eff / (6.0 * h**2)
                forces[i] -= f_rep * n
                forces[j] += f_rep * n

    def _apply_wall_loads(self, forces: np.ndarray, torques: np.ndarray) -> None:
        sim = self.sim
        for i in range(sim.n_p):
            wall_bot = sim.radii[i] + 0.5
            wall_top = sim.ny - 1.5 - sim.radii[i]
            if sim.pos[i, 1] < wall_bot:
                overlap = wall_bot - sim.pos[i, 1]
                v_n = -sim.vel[i, 1]
                f_mag = self.normal_contact_magnitude(overlap, v_n, sim.masses[i])
                forces[i, 1] += f_mag
                n = np.array([0.0, 1.0])
                self._apply_single_body_tangential_load(i, n, f_mag, forces, torques)
            if sim.pos[i, 1] > wall_top:
                overlap = sim.pos[i, 1] - wall_top
                v_n = sim.vel[i, 1]
                f_mag = self.normal_contact_magnitude(overlap, v_n, sim.masses[i])
                forces[i, 1] -= f_mag
                n = np.array([0.0, -1.0])
                self._apply_single_body_tangential_load(i, n, f_mag, forces, torques)

    def _apply_cylinder_loads(self, forces: np.ndarray, torques: np.ndarray) -> None:
        sim = self.sim
        for cx, cy, cr in sim.cylinders:
            for i in range(sim.n_p):
                dx = sim.pos[i, 0] - cx
                dy = sim.pos[i, 1] - cy
                dist = float(np.hypot(dx, dy))
                min_dist = cr + sim.radii[i]
                if dist < min_dist and dist > 1e-10:
                    overlap = min_dist - dist
                    nx_ = dx / dist
                    ny_ = dy / dist
                    v_n = sim.vel[i, 0] * nx_ + sim.vel[i, 1] * ny_
                    f_mag = self.normal_contact_magnitude(overlap, v_n, sim.masses[i])
                    forces[i, 0] += f_mag * nx_
                    forces[i, 1] += f_mag * ny_
                    self._apply_single_body_tangential_load(
                        i, np.array([nx_, ny_]), f_mag, forces, torques
                    )

    def _apply_single_body_tangential_load(
        self,
        i: int,
        normal: np.ndarray,
        normal_force: float,
        forces: np.ndarray,
        torques: np.ndarray,
    ) -> None:
        sim = self.sim
        t = self.tangent_from_normal(float(normal[0]), float(normal[1]))
        v_t = float(np.dot(sim.vel[i], t)) - sim.omega_p[i] * sim.radii[i]
        f_t = self.tangential_force_magnitude(v_t, normal_force, sim.masses[i])
        forces[i] += f_t * t
        torques[i] -= sim.radii[i] * f_t
        torques[i] += self.rolling_resistance_torque(
            sim.omega_p[i], normal_force, sim.radii[i], sim.masses[i]
        )
