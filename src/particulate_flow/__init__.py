"""Particulate-flow solvers with coupled LBM fluid and DEM particle models."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cfd-solver")
except PackageNotFoundError:
    __version__ = "0.4.0"

__author__ = "hardwork9047"
__license__ = "MIT"

from .dem.packing import DEMPackingSimulation, PackingMetrics
from .dem.solver import DEMSolver
from .fast_solver import FastLBMDEM
from .geometry.pore import Cylinder, PoreGeometry
from .io.config import SimulationConfig
from .io.verification import VerificationResult, run_fluid_verification, write_verification_outputs
from .lbm_dem import LBMDEMSolver, plot_fields, plot_particles

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
