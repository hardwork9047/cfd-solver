"""
CFD - Computational Fluid Dynamics Library

A library for Poiseuille flows, Cavity flows, and Cylinder flows
using finite-difference Navier-Stokes solvers.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cfd-solver")
except PackageNotFoundError:
    __version__ = "0.3.0"

__author__ = "hardwork9047"
__license__ = "MIT"

from .cavity import CavityFlow
from .cylinder import CylinderFlow
from .non_newtonian import PowerLawPlanePoiseuille
from .poiseuille import CircularPoiseuille, PlanePoiseuille

__all__ = [
    "PlanePoiseuille",
    "CircularPoiseuille",
    "CavityFlow",
    "CylinderFlow",
    "PowerLawPlanePoiseuille",
]
