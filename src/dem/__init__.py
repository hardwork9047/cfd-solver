"""
DEM - Discrete Element Method Library

A library for simulating granular particle dynamics using the Discrete
Element Method with Hertzian contact and velocity-Verlet integration.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cfd-solver")
except PackageNotFoundError:
    __version__ = "0.3.0"

__author__ = "hardwork9047"
__license__ = "MIT"

from .particles import ParticleSystem

__all__ = [
    "ParticleSystem",
]
