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
PARTICLE_FLUID_COUPLINGS = ("point_force", "solid_boundary", "immersed_boundary", "isp")
FLUID_ACCELERATORS = ("numpy", "numba", "auto")
COMPUTE_ACCELERATORS = ("numpy", "numba", "auto")
PARTICLE_SEARCH_METHODS = ("cell_list", "all_pairs")
Y_BOUNDARIES = ("wall", "periodic", "lees_edwards")
STREAMWISE_BOUNDARIES = ("periodic_force", "pressure")

# D2Q9 direction indices that cross the y boundaries (y=0 or y=ny-1 during streaming)
# Directions with cy > 0 cross the top boundary (y = ny-1 → y = 0)
LE_TOP_CROSS_DIRS = (2, 5, 6)   # N, NE, NW
# Directions with cy < 0 cross the bottom boundary (y = 0 → y = ny-1)
LE_BOT_CROSS_DIRS = (4, 7, 8)   # S, SW, SE
