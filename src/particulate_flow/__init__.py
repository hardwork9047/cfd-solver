"""Particulate-flow solvers with coupled LBM fluid and DEM particle models."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cfd-solver")
except PackageNotFoundError:
    __version__ = "0.3.0"

__author__ = "hardwork9047"
__license__ = "MIT"

from .dem_solver import DEMSolver
from .dem_packing import DEMPackingSimulation, PackingMetrics
from .fast_solver import FastLBMDEM
from .fluid_verification import VerificationResult, run_fluid_verification, write_verification_outputs
from .geometry import Cylinder, PoreGeometry
from .lbm_dem import LBMDEMSolver, plot_fields, plot_particles
from .simulation_config import SimulationConfig

__all__ = [
    "Cylinder",
    "DEMSolver",
    "DEMPackingSimulation",
    "FastLBMDEM",
    "LBMDEMSolver",
    "PoreGeometry",
    "PackingMetrics",
    "SimulationConfig",
    "VerificationResult",
    "plot_fields",
    "plot_particles",
    "run_fluid_verification",
    "write_verification_outputs",
]
