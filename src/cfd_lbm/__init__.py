"""cfd_lbm — Lattice-Boltzmann Method (D2Q9) solver."""
from .lbm import LBM, equilibrium, plot_results, plot_centerline
from .poiseuille import PoiseuilleFlow, plot_profile, plot_convergence, grid_refinement_study

__all__ = [
    "LBM", "equilibrium", "plot_results", "plot_centerline",
    "PoiseuilleFlow", "plot_profile", "plot_convergence", "grid_refinement_study",
]
