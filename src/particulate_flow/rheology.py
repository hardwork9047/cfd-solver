"""Apparent viscosity evaluation for Lees-Edwards shear flow suspensions.

Computes the suspension shear viscosity η_s = -⟨σ_xy⟩ / γ̇ from two contributions:

- η^H (hydrodynamic): estimated from the LBM non-equilibrium stress tensor.
- η^C (collisional): accumulated from particle-particle contact impulses.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from particulate_flow.lbm.constants import CS2, C, W

if TYPE_CHECKING:
    from particulate_flow.lbm_dem import LBMDEMSolver


def _fluid_stress_xy(sim: "LBMDEMSolver") -> float:
    """Return the domain-averaged off-diagonal fluid stress σ_xy.

    Uses the standard LBM relation for the non-equilibrium stress:
        σ_xy = -(1 - 1/(2τ)) Σ_α c_αx c_αy (f_α - f_α^eq)

    Args:
        sim: An LBMDEMSolver instance with attributes ``f`` (shape
             ``(9, nx, ny)``) and ``tau``.  ``rho``, ``ux``, and ``uy``
             are derived internally from ``f`` and are not read from ``sim``.

    Returns:
        Volume-averaged σ_xy [lattice units].
    """
    f = sim.f          # shape (9, nx, ny)
    tau = sim.tau

    # Compute macroscopic fields from f.
    rho = f.sum(axis=0)                                     # (nx, ny)
    ux = (C[:, 0, np.newaxis, np.newaxis] * f).sum(axis=0) / rho  # (nx, ny)
    uy = (C[:, 1, np.newaxis, np.newaxis] * f).sum(axis=0) / rho  # (nx, ny)

    # Equilibrium: f_eq_α = w_α ρ [1 + (c·u)/cs² + (c·u)²/(2cs⁴) - u²/(2cs²)]
    cu = (
        C[:, 0, np.newaxis, np.newaxis] * ux
        + C[:, 1, np.newaxis, np.newaxis] * uy
    )  # (9, nx, ny)
    u2 = ux**2 + uy**2  # (nx, ny)
    f_eq = (
        W[:, np.newaxis, np.newaxis]
        * rho
        * (1.0 + cu / CS2 + 0.5 * cu**2 / CS2**2 - 0.5 * u2 / CS2)
    )  # (9, nx, ny)

    f_neq = f - f_eq  # (9, nx, ny)

    # σ_xy = -(1 - 1/(2τ)) Σ_α cx_α cy_α f_neq_α
    cx_cy = C[:, 0] * C[:, 1]  # (9,)
    prefactor = -(1.0 - 1.0 / (2.0 * tau))
    sigma_xy_field = prefactor * np.einsum("k,kij->ij", cx_cy, f_neq)
    solid = getattr(sim, "solid", None)
    if solid is not None and solid.any():
        fluid_mask = ~solid
        return float(np.mean(sigma_xy_field[fluid_mask]))
    return float(np.mean(sigma_xy_field))


class ViscosityEvaluator:
    """Accumulate and output apparent viscosity η_s for a shear suspension.

    Args:
        sim: The coupled solver instance.
        start_step: First step at which to begin accumulating data.
        viscosity_interval: Write a CSV row every this many steps (from start_step).
        average_steps: Number of trailing rows to average in :meth:`finalize`.
        out_dir: Directory for ``viscosity_timeseries.csv``.
        enabled: When ``False`` all methods are no-ops (zero overhead).
    """

    def __init__(
        self,
        sim: "LBMDEMSolver",
        *,
        start_step: int = 0,
        viscosity_interval: int = 100,
        average_steps: int = 1000,
        out_dir: Path | str,
        enabled: bool = True,
    ) -> None:
        self.sim = sim
        self.start_step = start_step
        self.viscosity_interval = viscosity_interval
        self.average_steps = average_steps
        self.out_dir = Path(out_dir)
        self.enabled = enabled

        self._rows: list[dict[str, float]] = []
        self._contact_stress_acc: float = 0.0
        self._contact_acc_count: int = 0
        self._next_flush_step: int = start_step
        self._csv_path: Path = self.out_dir / "viscosity_timeseries.csv"
        self._csv_handle = None
        self._csv_writer = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def record(self, step: int) -> None:
        """Record one time step.

        Can be called once per step or once per frame (coarse-grained).
        Flushes a CSV row when the step first reaches or passes the next
        scheduled flush point (``start_step + N * viscosity_interval``).

        Args:
            step: Current simulation step index (e.g. ``sim.step_count``).

        Note:
            Internally calls ``sim.dem_solver.compute_contact_stress_xy()``
            to accumulate collisional stress.  ``sim`` must therefore have a
            ``dem_solver`` attribute (all ``LBMDEMSolver`` instances do).
        """
        if not self.enabled:
            return
        if step < self.start_step:
            return

        sigma_xy_C = self.sim.dem_solver.compute_contact_stress_xy()
        self._contact_stress_acc += sigma_xy_C
        self._contact_acc_count += 1

        if step < self._next_flush_step:
            return
        self._next_flush_step = step + self.viscosity_interval

        shear_rate = getattr(self.sim, "le_shear_rate", 0.0)
        sigma_xy_H = _fluid_stress_xy(self.sim)

        eta_H = (-sigma_xy_H / shear_rate) if shear_rate != 0.0 else 0.0

        if self._contact_acc_count > 0:
            sigma_xy_C = self._contact_stress_acc / self._contact_acc_count
        else:
            sigma_xy_C = 0.0
        eta_C = (-sigma_xy_C / shear_rate) if shear_rate != 0.0 else 0.0

        self._contact_stress_acc = 0.0
        self._contact_acc_count = 0

        row = {"step": float(step), "eta_H": eta_H, "eta_C": eta_C}
        self._rows.append(row)
        self._write_csv_row(row)

    def finalize(self) -> dict[str, Any]:
        """Return time-averaged apparent viscosity over the last ``average_steps`` rows.

        Returns:
            Dict with keys ``eta_H_mean``, ``eta_C_mean``, ``eta_s_mean``.
        """
        if self._csv_handle is not None:
            self._csv_handle.close()
            self._csv_handle = None

        tail = self._rows[-self.average_steps :] if self._rows else []
        if tail:
            eta_H_mean = float(np.mean([r["eta_H"] for r in tail]))
            eta_C_mean = float(np.mean([r["eta_C"] for r in tail]))
        else:
            eta_H_mean = 0.0
            eta_C_mean = 0.0
        return {
            "eta_H_mean": eta_H_mean,
            "eta_C_mean": eta_C_mean,
            "eta_s_mean": eta_H_mean + eta_C_mean,
        }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _write_csv_row(self, row: dict[str, float]) -> None:
        if self._csv_writer is None:
            self.out_dir.mkdir(parents=True, exist_ok=True)
            self._csv_handle = self._csv_path.open("w", newline="", encoding="utf-8")
            fieldnames = ["step", "eta_H", "eta_C"]
            self._csv_writer = csv.DictWriter(self._csv_handle, fieldnames=fieldnames)
            self._csv_writer.writeheader()
        self._csv_writer.writerow(row)
        self._csv_handle.flush()
