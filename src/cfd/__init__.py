"""
CFD - Computational Fluid Dynamics Library

A library for calculating Poiseuille flows, Cavity flows, and Cylinder flows.
"""

__version__ = "0.1.0"
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
