"""
cfd_dem_lbm — Coupled LBM (fluid) + DEM (particles) solver.
"""

from .lbm_dem import LBMDEMSolver, plot_fields, plot_particles

__all__ = ["LBMDEMSolver", "plot_fields", "plot_particles"]
