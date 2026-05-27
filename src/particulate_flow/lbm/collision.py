"""LBM collision method enumerations.

The collision operators themselves (_collide_bgk_guo, _collide_trt_guo, _lbm_step)
reside in LBMDEMSolver pending extraction in Phase 2. This module exposes the
method/accelerator string enumerations used by the solver and its CLI.
"""

from __future__ import annotations

from .constants import COMPUTE_ACCELERATORS, FLUID_ACCELERATORS, FLUID_METHODS

__all__ = ["FLUID_METHODS", "FLUID_ACCELERATORS", "COMPUTE_ACCELERATORS"]
