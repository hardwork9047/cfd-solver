"""IBM coupling utilities placeholder.

The IBM coupling logic (immersed boundary force application, marker management)
resides in LBMDEMSolver. This module is a namespace anchor for the ibm
subpackage; full extraction is deferred to Phase 2.
"""

from __future__ import annotations

from .kernels import (
    _particle_solid_mask_numba,
    _ibm_markers_numba,
    _ibm_particle_reaction_numba,
)

__all__ = [
    "_particle_solid_mask_numba",
    "_ibm_markers_numba",
    "_ibm_particle_reaction_numba",
]
