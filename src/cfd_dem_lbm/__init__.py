"""
cfd_dem_lbm — Coupled LBM (fluid) + DEM (particles) solver.

Combines D2Q9 Lattice-Boltzmann fluid dynamics with Discrete Element
Method particle simulation via bilinear interpolation coupling.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cfd-solver")
except PackageNotFoundError:
    __version__ = "0.3.0"

__author__ = "hardwork9047"
__license__ = "MIT"

from .dem_solver import DEMSolver
from .fast_solver import FastLBMDEM
from .fluid_verification import VerificationResult, run_fluid_verification, write_verification_outputs
from .lbm_dem import LBMDEMSolver, plot_fields, plot_particles
from .simulation_config import SimulationConfig

__all__ = [
    "DEMSolver",
    "FastLBMDEM",
    "LBMDEMSolver",
    "SimulationConfig",
    "VerificationResult",
    "plot_fields",
    "plot_particles",
    "run_fluid_verification",
    "write_verification_outputs",
]
