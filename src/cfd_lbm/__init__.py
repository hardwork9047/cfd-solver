"""
cfd_lbm — Lattice-Boltzmann Method (D2Q9) solver.

Provides BGK-collision D2Q9 LBM solvers for lid-driven cavity flow
and Poiseuille channel flow with Zou-He boundary conditions.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cfd-solver")
except PackageNotFoundError:
    __version__ = "0.3.0"

__author__ = "hardwork9047"
__license__ = "MIT"

from .lbm import LBM, equilibrium, plot_centerline, plot_results
from .poiseuille import PoiseuilleFlow, grid_refinement_study, plot_convergence, plot_profile

__all__ = [
    "LBM",
    "equilibrium",
    "plot_results",
    "plot_centerline",
    "PoiseuilleFlow",
    "plot_profile",
    "plot_convergence",
    "grid_refinement_study",
]
