"""Execution-path solver variants for LBM-DEM runs."""

from __future__ import annotations

from .lbm_dem import LBMDEMSolver


class FastLBMDEM(LBMDEMSolver):
    """LBM-DEM solver with per-step macroscopic-field caching.

    ``run_lbm_dem.py`` uses this class for production sweeps.  Verification
    code should use the same class when the goal is to validate the actual
    execution path rather than only the base solver equations.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._cached_rho = None
        self._cached_ux = None
        self._cached_uy = None
        self._cache_step = -1

    def _refresh_cache(self):
        self._cached_rho, self._cached_ux, self._cached_uy = super()._macroscopic()
        self._cache_step = self.step_count

    def _macroscopic(self):
        if self._cache_step != self.step_count:
            self._refresh_cache()
        return self._cached_rho, self._cached_ux, self._cached_uy

    def _invalidate_macroscopic_cache(self) -> None:
        self._cache_step = -1

    def advance(self, n_steps: int = 1):
        for _ in range(n_steps):
            self._refresh_cache()
            super().advance(1)
