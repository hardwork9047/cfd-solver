"""DEM load solver used by LBM-DEM production and verification paths."""

from __future__ import annotations

from typing import Callable

import numpy as np

PARTICLE_METHODS = ("dem-hertz", "dem-linear")

try:  # Optional acceleration for particle-boundary contacts.
    from numba import njit
except ImportError:  # pragma: no cover - depends on the local environment
    njit = None


if njit is not None:  # pragma: no cover - exercised when numba is installed

    @njit(cache=True)
    def _wall_cylinder_loads_numba(
        pos,
        vel,
        radii,
        masses,
        omega_p,
        cylinders,
        ny,
        k_n,
        damping,
        contact_model_id,
        rolling_friction,
        sliding_friction,
        tangential_damping,
        rolling_friction_coeff,
        rolling_damping,
        particle_attraction,
        particle_repulsion,
        attraction_strength,
        repulsion_strength,
        attraction_cutoff,
        repulsion_cutoff,
        attraction_min_gap,
        repulsion_min_gap,
        forces,
        torques,
    ):
        """Apply wall and fixed-cylinder DEM loads in compiled code."""
        n_p = pos.shape[0]
        for i in range(n_p):
            wall_bot = radii[i] + 0.5
            wall_top = ny - 1.5 - radii[i]
            if pos[i, 1] < wall_bot:
                overlap = wall_bot - pos[i, 1]
                v_n = -vel[i, 1]
                if contact_model_id == 1:
                    f_n = k_n * overlap
                else:
                    f_n = k_n * overlap**1.5
                f_mag = f_n - damping * v_n * (k_n * masses[i]) ** 0.5
                if f_mag < 0.0:
                    f_mag = 0.0
                forces[i, 1] += f_mag
                if rolling_friction and f_mag > 0.0:
                    tx = -1.0
                    ty = 0.0
                    v_t = vel[i, 0] * tx + vel[i, 1] * ty - omega_p[i] * radii[i]
                    f_trial = -tangential_damping * (k_n * masses[i]) ** 0.5 * v_t
                    f_limit = sliding_friction * f_mag
                    f_t = min(max(f_trial, -f_limit), f_limit)
                    forces[i, 0] += f_t * tx
                    forces[i, 1] += f_t * ty
                    torques[i] -= radii[i] * f_t
                    torque_trial = -rolling_damping * (k_n * masses[i]) ** 0.5 * radii[i] ** 2 * omega_p[i]
                    torque_limit = rolling_friction_coeff * f_mag * radii[i]
                    torques[i] += min(max(torque_trial, -torque_limit), torque_limit)

            if pos[i, 1] > wall_top:
                overlap = pos[i, 1] - wall_top
                v_n = vel[i, 1]
                if contact_model_id == 1:
                    f_n = k_n * overlap
                else:
                    f_n = k_n * overlap**1.5
                f_mag = f_n - damping * v_n * (k_n * masses[i]) ** 0.5
                if f_mag < 0.0:
                    f_mag = 0.0
                forces[i, 1] -= f_mag
                if rolling_friction and f_mag > 0.0:
                    tx = 1.0
                    ty = 0.0
                    v_t = vel[i, 0] * tx + vel[i, 1] * ty - omega_p[i] * radii[i]
                    f_trial = -tangential_damping * (k_n * masses[i]) ** 0.5 * v_t
                    f_limit = sliding_friction * f_mag
                    f_t = min(max(f_trial, -f_limit), f_limit)
                    forces[i, 0] += f_t * tx
                    forces[i, 1] += f_t * ty
                    torques[i] -= radii[i] * f_t
                    torque_trial = -rolling_damping * (k_n * masses[i]) ** 0.5 * radii[i] ** 2 * omega_p[i]
                    torque_limit = rolling_friction_coeff * f_mag * radii[i]
                    torques[i] += min(max(torque_trial, -torque_limit), torque_limit)

            for c in range(cylinders.shape[0]):
                dx = pos[i, 0] - cylinders[c, 0]
                dy = pos[i, 1] - cylinders[c, 1]
                dist = (dx * dx + dy * dy) ** 0.5
                min_dist = cylinders[c, 2] + radii[i]
                if dist > 1e-10:
                    nx_ = dx / dist
                    ny_ = dy / dist
                    if dist < min_dist:
                        overlap = min_dist - dist
                        v_n = vel[i, 0] * nx_ + vel[i, 1] * ny_
                        if contact_model_id == 1:
                            f_n = k_n * overlap
                        else:
                            f_n = k_n * overlap**1.5
                        f_mag = f_n - damping * v_n * (k_n * masses[i]) ** 0.5
                        if f_mag < 0.0:
                            f_mag = 0.0
                        forces[i, 0] += f_mag * nx_
                        forces[i, 1] += f_mag * ny_
                        if rolling_friction and f_mag > 0.0:
                            tx = -ny_
                            ty = nx_
                            v_t = vel[i, 0] * tx + vel[i, 1] * ty - omega_p[i] * radii[i]
                            f_trial = -tangential_damping * (k_n * masses[i]) ** 0.5 * v_t
                            f_limit = sliding_friction * f_mag
                            f_t = min(max(f_trial, -f_limit), f_limit)
                            forces[i, 0] += f_t * tx
                            forces[i, 1] += f_t * ty
                            torques[i] -= radii[i] * f_t
                            torque_trial = -rolling_damping * (k_n * masses[i]) ** 0.5 * radii[i] ** 2 * omega_p[i]
                            torque_limit = rolling_friction_coeff * f_mag * radii[i]
                            torques[i] += min(max(torque_trial, -torque_limit), torque_limit)

                    surface_gap = dist - min_dist
                    r_eff = radii[i] * cylinders[c, 2] / min_dist
                    if particle_attraction and attraction_strength > 0.0:
                        if surface_gap <= attraction_cutoff:
                            h = surface_gap
                            if h < attraction_min_gap:
                                h = attraction_min_gap
                            if h < 1e-12:
                                h = 1e-12
                            f_attr = attraction_strength * r_eff / (6.0 * h**2)
                            forces[i, 0] -= f_attr * nx_
                            forces[i, 1] -= f_attr * ny_
                    elif particle_repulsion and repulsion_strength > 0.0:
                        if surface_gap <= repulsion_cutoff:
                            h = surface_gap
                            if h < repulsion_min_gap:
                                h = repulsion_min_gap
                            if h < 1e-12:
                                h = 1e-12
                            f_rep = repulsion_strength * r_eff / (6.0 * h**2)
                            forces[i, 0] += f_rep * nx_
                            forces[i, 1] += f_rep * ny_

else:
    _wall_cylinder_loads_numba = None


class DEMSolver:
    """Compute DEM forces and torques for a coupled LBM-DEM simulation.

    The solver owns the DEM contact, surface-force, wall, and cylinder-load
    logic.  It reads particle and geometry arrays from the coupled simulation
    object and asks that object for fluid drag, so production runs and
    verification can share the same DEM implementation without duplicating
    contact physics in scripts.
    """

    def __init__(
        self,
        coupled_solver,
        pair_kernel: Callable | None = None,
        contact_model: str = "dem-hertz",
    ):
        if contact_model not in PARTICLE_METHODS:
            raise ValueError(f"contact_model must be one of {PARTICLE_METHODS}")
        self.sim = coupled_solver
        self.pair_kernel = pair_kernel
        self.contact_model = contact_model

    @staticmethod
    def tangent_from_normal(nx_: float, ny_: float) -> np.ndarray:
        """Return a unit tangent for a 2-D contact normal."""
        return np.array([-ny_, nx_])

    def normal_contact_magnitude(self, overlap: float, v_n: float, mass: float) -> float:
        """Normal contact force with damping, clamped to avoid artificial tension."""
        sim = self.sim
        if self.contact_model == "dem-linear":
            f_n = sim.k_n * overlap
        else:
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

        # 2. Fluid-particle coupling load.
        if sim.particle_fluid_coupling == "immersed_boundary":
            forces += sim.ibm_forces_p
            torques += sim.ibm_torques_p
        else:
            forces += sim._particle_drag_forces()

        # 3. Particle-particle contact and optional near-surface forces.
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
                getattr(sim, "surface_roughness", 0.0),
                forces,
                torques,
            )

        if self.pair_kernel is None:
            for i, j in sim._particle_pair_candidates():
                self._apply_particle_pair_loads(i, j, forces, torques)

        if (
            sim.uses_numba_compute
            and _wall_cylinder_loads_numba is not None
            and sim.y_boundary not in ("periodic", "lees_edwards")
        ):
            contact_model_id = 1 if self.contact_model == "dem-linear" else 0
            cylinders = np.asarray(sim.cylinders, dtype=np.float64).reshape((-1, 3))
            _wall_cylinder_loads_numba(
                sim.pos,
                sim.vel,
                sim.radii,
                sim.masses,
                sim.omega_p,
                cylinders,
                sim.ny,
                sim.k_n,
                sim.damping,
                contact_model_id,
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
        else:
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
        dp[1] = sim._periodic_y_delta(float(dp[1]))
        dist = float(np.linalg.norm(dp))
        geom_min_dist = sim.radii[i] + sim.radii[j]
        h_r = getattr(sim, "surface_roughness", 0.0)
        min_dist = geom_min_dist + h_r
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

        self._apply_surface_force(i, j, dist, geom_min_dist, n, forces)

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
        if sim.y_boundary in ("periodic", "lees_edwards"):
            return
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
                dy = sim._periodic_y_delta(sim.pos[i, 1] - cy)
                dist = float(np.hypot(dx, dy))
                min_dist = cr + sim.radii[i]
                if dist > 1e-10:
                    nx_ = dx / dist
                    ny_ = dy / dist
                    if dist < min_dist:
                        overlap = min_dist - dist
                        v_n = sim.vel[i, 0] * nx_ + sim.vel[i, 1] * ny_
                        f_mag = self.normal_contact_magnitude(overlap, v_n, sim.masses[i])
                        forces[i, 0] += f_mag * nx_
                        forces[i, 1] += f_mag * ny_
                        self._apply_single_body_tangential_load(
                            i, np.array([nx_, ny_]), f_mag, forces, torques
                        )
                    self._apply_cylinder_surface_force(
                        i,
                        dist,
                        min_dist,
                        cr,
                        np.array([nx_, ny_]),
                        forces,
                    )

    def _apply_cylinder_surface_force(
        self,
        i: int,
        dist: float,
        min_dist: float,
        cylinder_radius: float,
        normal_outward: np.ndarray,
        forces: np.ndarray,
    ) -> None:
        """Apply Hamaker-like particle-cylinder attraction or repulsion."""
        sim = self.sim
        surface_gap = dist - min_dist
        r_eff = sim.radii[i] * cylinder_radius / min_dist
        if sim.particle_attraction and sim.attraction_strength > 0.0:
            if surface_gap <= sim.attraction_cutoff:
                h = max(surface_gap, max(sim.attraction_min_gap, 1e-12))
                f_attr = sim.attraction_strength * r_eff / (6.0 * h**2)
                forces[i] -= f_attr * normal_outward
        elif sim.particle_repulsion and sim.repulsion_strength > 0.0:
            if surface_gap <= sim.repulsion_cutoff:
                h = max(surface_gap, max(sim.repulsion_min_gap, 1e-12))
                f_rep = sim.repulsion_strength * r_eff / (6.0 * h**2)
                forces[i] += f_rep * normal_outward

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

    def compute_contact_stress_xy(self) -> float:
        """Return the domain-averaged collisional off-diagonal stress σ_xy^C.

        Uses the Irving-Kirkwood formula: σ_xy^C = (1/V) Σ_{pairs} F_x r_y,
        where F is the **normal** contact force on particle j and
        r = pos_j - pos_i.  Tangential (friction) contributions are excluded.

        Only particle-particle contacts contribute; wall and cylinder contacts
        are excluded as they involve fixed objects.

        Returns:
            Collisional σ_xy^C in lattice units. Zero when no contacts exist.
        """
        sim = self.sim
        if sim.n_p < 2:
            return 0.0
        volume = float(sim.nx * sim.ny)
        stress_xy = 0.0
        h_r = getattr(sim, "surface_roughness", 0.0)
        for i, j in sim._particle_pair_candidates():
            dp = sim.pos[j] - sim.pos[i]
            dp[1] = sim._periodic_y_delta(float(dp[1]))
            dist = float(np.linalg.norm(dp))
            geom_min_dist = sim.radii[i] + sim.radii[j]
            min_dist = geom_min_dist + h_r
            if dist <= 1e-10 or dist >= min_dist:
                continue
            n = dp / dist
            overlap = min_dist - dist
            v_n = float(np.dot(sim.vel[j] - sim.vel[i], n))
            avg_m = (sim.masses[i] + sim.masses[j]) / 2.0
            f_mag = self.normal_contact_magnitude(overlap, v_n, avg_m)
            f_vec = f_mag * n  # force on j from i
            # Irving-Kirkwood: F_x * r_y where r = pos_j - pos_i
            stress_xy += float(f_vec[0] * dp[1])
        return stress_xy / volume
