"""Optional Numba-compiled IBM kernels. Falls back to None when Numba is unavailable."""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
except ImportError:
    njit = None

if njit is not None:  # pragma: no cover - exercised only when numba is installed

    @njit(cache=True)
    def _particle_solid_mask_numba(pos, radii, nx, ny, fixed_solid):
        """Rasterise particle discs onto the solid mask grid."""
        particle_solid = np.zeros((nx, ny), dtype=np.bool_)
        for i in range(pos.shape[0]):
            radius = radii[i]
            x_pos = pos[i, 0]
            y_pos = pos[i, 1]
            x_min = int(np.floor(x_pos - radius))
            x_max = int(np.ceil(x_pos + radius))
            y_min = max(0, int(np.floor(y_pos - radius)))
            y_max = min(ny - 1, int(np.ceil(y_pos + radius)))
            r2 = radius * radius
            for x_raw in range(x_min, x_max + 1):
                ix = x_raw % nx
                dx = x_raw - x_pos
                for iy in range(y_min, y_max + 1):
                    dy = iy - y_pos
                    if dx * dx + dy * dy <= r2 and not fixed_solid[ix, iy]:
                        particle_solid[ix, iy] = True
        return particle_solid

    @njit(cache=True)
    def _ibm_markers_numba(pos, vel, omega_p, radii, marker_spacing, ny, periodic_y):
        """Generate IBM marker points on particle surfaces."""
        total = 0
        for i in range(pos.shape[0]):
            total += max(8, int(np.ceil(2.0 * np.pi * radii[i] / marker_spacing)))
        marker_x = np.empty(total)
        marker_y = np.empty(total)
        marker_rx = np.empty(total)
        marker_ry = np.empty(total)
        marker_ds = np.empty(total)
        marker_owner = np.empty(total, dtype=np.int64)
        marker_ubx = np.empty(total)
        marker_uby = np.empty(total)
        cursor = 0
        for i in range(pos.shape[0]):
            radius = radii[i]
            n_markers = max(8, int(np.ceil(2.0 * np.pi * radius / marker_spacing)))
            ds = 2.0 * np.pi * radius / n_markers
            for m in range(n_markers):
                theta = 2.0 * np.pi * m / n_markers
                rx = radius * np.cos(theta)
                ry = radius * np.sin(theta)
                y = pos[i, 1] + ry
                if periodic_y:
                    y = y % ny
                else:
                    if y < 0.5:
                        y = 0.5
                    elif y > ny - 1.5:
                        y = ny - 1.5
                marker_x[cursor] = pos[i, 0] + rx
                marker_y[cursor] = y
                marker_rx[cursor] = rx
                marker_ry[cursor] = ry
                marker_ds[cursor] = ds
                marker_owner[cursor] = i
                marker_ubx[cursor] = vel[i, 0] - omega_p[i] * ry
                marker_uby[cursor] = vel[i, 1] + omega_p[i] * rx
                cursor += 1
        return (
            marker_x, marker_y, marker_rx, marker_ry,
            marker_ds, marker_owner, marker_ubx, marker_uby,
        )

    @njit(cache=True)
    def _ibm_particle_reaction_numba(
        marker_owner,
        marker_rx,
        marker_ry,
        marker_fx,
        marker_fy,
        n_particles,
    ):
        """Accumulate IBM marker forces/torques back onto particles."""
        forces = np.zeros((n_particles, 2))
        torques = np.zeros(n_particles)
        for k in range(marker_owner.shape[0]):
            owner = marker_owner[k]
            forces[owner, 0] -= marker_fx[k]
            forces[owner, 1] -= marker_fy[k]
            torques[owner] -= marker_rx[k] * marker_fy[k] - marker_ry[k] * marker_fx[k]
        return forces, torques

else:
    _particle_solid_mask_numba = None
    _ibm_markers_numba = None
    _ibm_particle_reaction_numba = None
