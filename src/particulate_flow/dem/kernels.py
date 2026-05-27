"""Optional Numba-compiled DEM kernel. Falls back to None when Numba is unavailable."""

from __future__ import annotations

import numpy as np

from ..lbm.constants import COMPUTE_ACCELERATORS  # shared accelerator enum

try:
    from numba import njit
except ImportError:
    njit = None

if njit is not None:  # pragma: no cover - exercised only when numba is installed

    @njit(cache=True)
    def _dem_pair_loads_numba(
        pair_i,
        pair_j,
        pos,
        vel,
        radii,
        masses,
        omega_p,
        k_n,
        damping,
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
        """Apply particle-pair DEM loads in compiled code."""
        for pair_idx in range(pair_i.shape[0]):
            i = pair_i[pair_idx]
            j = pair_j[pair_idx]
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dist = (dx * dx + dy * dy) ** 0.5
            min_dist = radii[i] + radii[j]
            if dist <= 1e-10:
                continue

            nx_ = dx / dist
            ny_ = dy / dist

            if dist < min_dist:
                overlap = min_dist - dist
                rvx = vel[j, 0] - vel[i, 0]
                rvy = vel[j, 1] - vel[i, 1]
                v_n = rvx * nx_ + rvy * ny_
                avg_m = 0.5 * (masses[i] + masses[j])
                f_n = k_n * overlap**1.5
                f_damp = -damping * v_n * (k_n * avg_m) ** 0.5
                f_mag = f_n + f_damp
                if f_mag < 0.0:
                    f_mag = 0.0

                fx = f_mag * nx_
                fy = f_mag * ny_
                forces[i, 0] -= fx
                forces[i, 1] -= fy
                forces[j, 0] += fx
                forces[j, 1] += fy

                if rolling_friction and f_mag > 0.0:
                    tx = -ny_
                    ty = nx_
                    v_t = rvx * tx + rvy * ty
                    v_t -= omega_p[j] * radii[j] + omega_p[i] * radii[i]
                    f_trial = -tangential_damping * (k_n * avg_m) ** 0.5 * v_t
                    f_limit = sliding_friction * f_mag
                    f_t = min(max(f_trial, -f_limit), f_limit)
                    ft_x = f_t * tx
                    ft_y = f_t * ty
                    forces[i, 0] -= ft_x
                    forces[i, 1] -= ft_y
                    forces[j, 0] += ft_x
                    forces[j, 1] += ft_y
                    torques[i] -= radii[i] * f_t
                    torques[j] -= radii[j] * f_t

                    trial_i = -rolling_damping * (k_n * masses[i]) ** 0.5 * radii[i] ** 2 * omega_p[i]
                    limit_i = rolling_friction_coeff * f_mag * radii[i]
                    torques[i] += min(max(trial_i, -limit_i), limit_i)
                    trial_j = -rolling_damping * (k_n * masses[j]) ** 0.5 * radii[j] ** 2 * omega_p[j]
                    limit_j = rolling_friction_coeff * f_mag * radii[j]
                    torques[j] += min(max(trial_j, -limit_j), limit_j)

            surface_gap = dist - min_dist
            if particle_attraction and attraction_strength > 0.0:
                if surface_gap <= attraction_cutoff:
                    h = surface_gap
                    if h < attraction_min_gap:
                        h = attraction_min_gap
                    if h < 1e-12:
                        h = 1e-12
                    r_eff = radii[i] * radii[j] / min_dist
                    f_attr = attraction_strength * r_eff / (6.0 * h**2)
                    forces[i, 0] += f_attr * nx_
                    forces[i, 1] += f_attr * ny_
                    forces[j, 0] -= f_attr * nx_
                    forces[j, 1] -= f_attr * ny_
            elif particle_repulsion and repulsion_strength > 0.0:
                if surface_gap <= repulsion_cutoff:
                    h = surface_gap
                    if h < repulsion_min_gap:
                        h = repulsion_min_gap
                    if h < 1e-12:
                        h = 1e-12
                    r_eff = radii[i] * radii[j] / min_dist
                    f_rep = repulsion_strength * r_eff / (6.0 * h**2)
                    forces[i, 0] -= f_rep * nx_
                    forces[i, 1] -= f_rep * ny_
                    forces[j, 0] += f_rep * nx_
                    forces[j, 1] += f_rep * ny_

else:
    _dem_pair_loads_numba = None
