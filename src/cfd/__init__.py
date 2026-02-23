"""
CFD - Computational Fluid Dynamics Library

A library for calculating Poiseuille flows and Cavity flows.
"""

__version__ = "0.1.0"
__author__ = "Your Name"
__license__ = "MIT"

from .poiseuille import PlanePoiseuille, CircularPoiseuille
from .cavity import CavityFlow
from .non_newtonian import PowerLawPlanePoiseuille

__all__ = [
    "PlanePoiseuille",
    "CircularPoiseuille", 
    "CavityFlow",
    "PowerLawPlanePoiseuille",
]
