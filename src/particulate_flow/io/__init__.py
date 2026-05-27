"""I/O subpackage: config, paths, verification, visualization."""

from .config import SimulationConfig
from .paths import program_results_dir

# NOTE: VerificationResult and helpers are imported lazily at the top-level
# particulate_flow package to avoid a circular import (verification -> fast_solver
# -> lbm_dem -> io.paths -> io.__init__). Import directly from io.verification
# when needed, or use the top-level particulate_flow namespace.

__all__ = [
    "SimulationConfig",
    "program_results_dir",
]
