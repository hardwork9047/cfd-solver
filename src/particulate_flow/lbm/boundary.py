"""LBM boundary condition utilities (pressure inlet/outlet)."""

from __future__ import annotations

# Boundary condition logic resides in LBMDEMSolver._apply_pressure_boundaries.
# This module exposes the boundary type enumerations for external consumers.

from .constants import Y_BOUNDARIES, STREAMWISE_BOUNDARIES

__all__ = ["Y_BOUNDARIES", "STREAMWISE_BOUNDARIES"]
