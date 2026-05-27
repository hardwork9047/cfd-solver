"""D2Q9 lattice constants and solver option enumerations."""

from __future__ import annotations

import numpy as np

# Velocity vectors: C[i] = (cx, cy)
C = np.array(
    [
        [0, 0],   # 0 rest
        [1, 0],   # 1 E
        [0, 1],   # 2 N
        [-1, 0],  # 3 W
        [0, -1],  # 4 S
        [1, 1],   # 5 NE
        [-1, 1],  # 6 NW
        [-1, -1], # 7 SW
        [1, -1],  # 8 SE
    ],
    dtype=float,
)

W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
OPPOSITE = [0, 3, 4, 1, 2, 7, 8, 5, 6]
Q = 9
CS2 = 1.0 / 3.0  # lattice speed of sound squared

FLUID_METHODS = ("lbm-bgk-guo", "lbm-trt-guo")
PARTICLE_FLUID_COUPLINGS = ("point_force", "solid_boundary", "immersed_boundary")
FLUID_ACCELERATORS = ("numpy", "numba", "auto")
COMPUTE_ACCELERATORS = ("numpy", "numba", "auto")
PARTICLE_SEARCH_METHODS = ("cell_list", "all_pairs")
Y_BOUNDARIES = ("wall", "periodic")
STREAMWISE_BOUNDARIES = ("periodic_force", "pressure")
