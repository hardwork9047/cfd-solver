"""Optional Numba-compiled LBM kernel. Falls back to None when Numba is unavailable."""

from __future__ import annotations

import numpy as np

from .constants import CS2

try:
    from numba import njit
except ImportError:
    njit = None

if njit is not None:  # pragma: no cover - exercised only when numba is installed

    @njit(cache=True)
    def _lbm_step_numba(
        f,
        rho,
        ux,
        uy,
        Fx,
        Fy,
        solid,
        omega,
        trt_omega_minus,
        method_id,
        c_int,
        w,
        opposite,
    ):
        """Fused D2Q9 collision, streaming, and bounce-back kernel."""
        q_count, nx, ny = f.shape
        post = np.empty_like(f)
        streamed = np.empty_like(f)
        for q in range(q_count):
            cx = c_int[q, 0]
            cy = c_int[q, 1]
            opp = opposite[q]
            ox = c_int[opp, 0]
            oy = c_int[opp, 1]
            weight = w[q]
            for x in range(nx):
                for y in range(ny):
                    u2 = ux[x, y] * ux[x, y] + uy[x, y] * uy[x, y]
                    cu = cx * ux[x, y] + cy * uy[x, y]
                    feq = weight * rho[x, y] * (
                        1.0
                        + cu / CS2
                        + 0.5 * cu * cu / (CS2 * CS2)
                        - 0.5 * u2 / CS2
                    )
                    term_x = (cx - ux[x, y]) / CS2 + cu / (CS2 * CS2) * cx
                    term_y = (cy - uy[x, y]) / CS2 + cu / (CS2 * CS2) * cy
                    forcing = (1.0 - 0.5 * omega) * weight * (
                        term_x * Fx[x, y] + term_y * Fy[x, y]
                    )
                    if method_id == 1:
                        cu_opp = ox * ux[x, y] + oy * uy[x, y]
                        feq_opp = w[opp] * rho[x, y] * (
                            1.0
                            + cu_opp / CS2
                            + 0.5 * cu_opp * cu_opp / (CS2 * CS2)
                            - 0.5 * u2 / CS2
                        )
                        f_plus = 0.5 * (f[q, x, y] + f[opp, x, y])
                        f_minus = 0.5 * (f[q, x, y] - f[opp, x, y])
                        feq_plus = 0.5 * (feq + feq_opp)
                        feq_minus = 0.5 * (feq - feq_opp)
                        post[q, x, y] = (
                            f[q, x, y]
                            - omega * (f_plus - feq_plus)
                            - trt_omega_minus * (f_minus - feq_minus)
                            + forcing
                        )
                    else:
                        post[q, x, y] = f[q, x, y] + omega * (feq - f[q, x, y]) + forcing

        for q in range(q_count):
            cx = c_int[q, 0]
            cy = c_int[q, 1]
            for x in range(nx):
                xd = (x + cx) % nx
                for y in range(ny):
                    yd = (y + cy) % ny
                    streamed[q, xd, yd] = post[q, x, y]

        out = streamed.copy()
        for q in range(q_count):
            opp = opposite[q]
            for x in range(nx):
                for y in range(ny):
                    if solid[x, y]:
                        out[q, x, y] = streamed[opp, x, y]
        return out

else:
    _lbm_step_numba = None
