"""Pure LBM operator functions: equilibrium distribution and Guo forcing term."""

from __future__ import annotations

import numpy as np

from .constants import CS2, C, Q, W


def equilibrium(rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    """Maxwell-Boltzmann equilibrium f_eq[q, x, y]."""
    u2 = ux**2 + uy**2
    feq = np.empty((Q,) + rho.shape)
    for i in range(Q):
        cu = C[i, 0] * ux + C[i, 1] * uy
        feq[i] = W[i] * rho * (1.0 + cu / CS2 + 0.5 * cu**2 / CS2**2 - 0.5 * u2 / CS2)
    return feq


def guo_forcing(
    ux: np.ndarray,
    uy: np.ndarray,
    Fx: np.ndarray,
    Fy: np.ndarray,
    omega: float,
) -> np.ndarray:
    """
    Guo et al. (2002) forcing term S[q, x, y].

    S_i = (1 - omega/2) * W_i * [ (e_i - u)/cs² + (e_i·u)/cs⁴ * e_i ] · F
    """
    S = np.empty((Q,) + ux.shape)
    for i in range(Q):
        cu = C[i, 0] * ux + C[i, 1] * uy
        term_x = (C[i, 0] - ux) / CS2 + cu / CS2**2 * C[i, 0]
        term_y = (C[i, 1] - uy) / CS2 + cu / CS2**2 * C[i, 1]
        S[i] = (1.0 - 0.5 * omega) * W[i] * (term_x * Fx + term_y * Fy)
    return S
