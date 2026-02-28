"""
Lattice-Boltzmann Method (LBM) — D2Q9 Poiseuille (Channel) Flow
================================================================
Body-force-driven flow between two infinite parallel plates.
Driving force applied via Guo et al. (2002) forcing scheme.

Geometry (y-direction):
    y = 0        : bottom wall  (bounce-back, effective no-slip at y = 0.5)
    y = 1..ny-2  : fluid nodes
    y = ny-1     : top wall     (bounce-back, effective no-slip at y = ny-1.5)

    Effective channel height : H = ny - 1  (lattice units)

Analytical solution (parabolic Hagen-Poiseuille profile):
    u(y) = Fx / (2 ν) · ŷ · (H - ŷ)    where  ŷ = y - 0.5  (bounce-back offset)
    u_max = Fx · H² / (8 ν)

Reynolds number used here:
    Re = u_max · H / ν

References:
  - Guo, Z., Zheng, C., & Shi, B. (2002). Discrete lattice effects on the
    forcing term in the lattice Boltzmann method. Phys. Rev. E 65, 046308.
  - Krüger et al. (2017). The Lattice Boltzmann Method: Principles and Practice.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from .lbm import C, W, OPPOSITE, Q, CS2, equilibrium


# ---------------------------------------------------------------------------
# Guo forcing term
# ---------------------------------------------------------------------------

def guo_force(ux: np.ndarray, uy: np.ndarray,
              Fx: float, Fy: float,
              omega: float) -> np.ndarray:
    """
    Compute Guo et al. (2002) discrete forcing term F_i [Q, nx, ny].

        F_i = w_i (1 - ω/2) [ (e_i - u)/cs² + (e_i·u) e_i / cs⁴ ] · F

    Parameters
    ----------
    ux, uy : velocity arrays  (nx, ny)
    Fx, Fy : body force per unit mass  (scalar)
    omega  : BGK inverse relaxation time
    """
    coeff = 1.0 - 0.5 * omega
    Fi = np.empty((Q,) + ux.shape)
    for i in range(Q):
        cx, cy = C[i, 0], C[i, 1]
        eu = cx * ux + cy * uy                  # e_i · u
        term = ((cx - ux) * Fx + (cy - uy) * Fy) / CS2 \
               + eu * (cx * Fx + cy * Fy) / CS2**2
        Fi[i] = coeff * W[i] * term
    return Fi


# ---------------------------------------------------------------------------
# PoiseuilleFlow solver
# ---------------------------------------------------------------------------

class PoiseuilleFlow:
    """
    2D D2Q9 LBM solver for pressure-gradient-driven Poiseuille channel flow.

    The pressure gradient is mimicked by a uniform body force Fx applied
    to every fluid node (Guo forcing scheme).
    """

    def __init__(
        self,
        nx: int = 8,
        ny: int = 64,
        Re: float = 10.0,
        u_max: float = 0.05,
    ):
        """
        Parameters
        ----------
        nx    : channel length (periodic, any value works; 8 is fine)
        ny    : total grid height including walls (fluid = ny-2 layers)
        Re    : Reynolds number  Re = u_max · H / ν,  H = ny-1
        u_max : target maximum velocity (lattice units, keep ≤ 0.1)
        """
        self.nx    = nx
        self.ny    = ny
        self.Re    = Re
        self.u_max = u_max

        # Effective channel height (bounce-back half-node correction)
        # Bottom no-slip at y=0.5, top no-slip at y=ny-1.5  → H_eff = ny-2
        self.H = float(ny - 2)

        # Kinematic viscosity and relaxation time from Re & u_max
        self.nu    = u_max * self.H / Re
        self.tau   = self.nu / CS2 + 0.5
        self.omega = 1.0 / self.tau

        # Body force that produces u_max at steady state:  u_max = Fx H² / (8 ν)
        self.Fx = 8.0 * self.nu * u_max / self.H**2
        self.Fy = 0.0

        print("LBM D2Q9 — Poiseuille Channel Flow")
        print(f"  Grid    : {nx} × {ny}   (fluid layers: {ny-2})")
        print(f"  H (eff) : {self.H:.1f}  (no-slip at y=0.5 and y={ny-1.5})")
        print(f"  Re      : {Re:.2f}")
        print(f"  u_max   : {u_max:.5f}")
        print(f"  nu      : {self.nu:.6f}   tau : {self.tau:.4f}")
        print(f"  Fx      : {self.Fx:.3e}")

        # ---- distribution functions  f[q, x, y] ----
        rho0 = np.ones((nx, ny))
        ux0  = np.zeros((nx, ny))
        uy0  = np.zeros((nx, ny))
        self.f = equilibrium(rho0, ux0, uy0)

        # ---- wall mask ----
        self.solid = np.zeros((nx, ny), dtype=bool)
        self.solid[:, 0]  = True   # bottom wall
        self.solid[:, -1] = True   # top wall

        self.step = 0

    # ------------------------------------------------------------------
    # Core LBM steps
    # ------------------------------------------------------------------

    def _macroscopic(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Density and velocity; velocity corrected for Guo forcing."""
        rho = self.f.sum(axis=0)
        ux  = (C[:, 0, np.newaxis, np.newaxis] * self.f).sum(axis=0) / rho
        uy  = (C[:, 1, np.newaxis, np.newaxis] * self.f).sum(axis=0) / rho
        # Guo correction: u_true = u_streaming + F*Δt/(2ρ)
        ux += 0.5 * self.Fx / rho
        return rho, ux, uy

    def _collide(self, rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> None:
        """BGK collision + Guo body-force term."""
        feq = equilibrium(rho, ux, uy)
        Fi  = guo_force(ux, uy, self.Fx, self.Fy, self.omega)
        self.f += self.omega * (feq - self.f) + Fi

    def _stream(self) -> None:
        """Streaming (periodic in x via np.roll)."""
        for i in range(Q):
            self.f[i] = np.roll(self.f[i], int(C[i, 0]), axis=0)
            self.f[i] = np.roll(self.f[i], int(C[i, 1]), axis=1)

    def _boundary(self) -> None:
        """Full bounce-back on top and bottom walls."""
        f_tmp = self.f.copy()
        for i in range(Q):
            j = OPPOSITE[i]
            self.f[i, self.solid] = f_tmp[j, self.solid]

    def advance(self, n_steps: int = 1) -> None:
        """Advance the simulation by n_steps time steps."""
        for _ in range(n_steps):
            rho, ux, uy = self._macroscopic()
            self._collide(rho, ux, uy)
            self._stream()
            self._boundary()
            self.step += 1

    def get_fields(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self._macroscopic()

    # ------------------------------------------------------------------
    # Analytical solution
    # ------------------------------------------------------------------

    def analytical_u(self, y: np.ndarray) -> np.ndarray:
        """
        Parabolic Hagen-Poiseuille profile at fluid nodes.

        Parameters
        ----------
        y : node indices (integer array, e.g. np.arange(1, ny-1))

        Returns
        -------
        u_analytical : velocity at each node (bounce-back corrected)
        """
        # ξ = distance from bottom no-slip plane (y=0.5)  →  ξ ∈ [0, H_eff]
        xi = y - 0.5
        return self.Fx / (2.0 * self.nu) * xi * (self.H - xi)

    def l2_error(self) -> float:
        """L2 relative error between LBM u-profile and analytical solution."""
        _, ux, _ = self.get_fields()
        j = np.arange(1, self.ny - 1)
        u_lbm  = ux[0, j]                  # x-averaged (periodic → uniform)
        u_anal = self.analytical_u(j.astype(float))
        err = np.linalg.norm(u_lbm - u_anal) / (np.linalg.norm(u_anal) + 1e-30)
        return float(err)


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_profile(sim: PoiseuilleFlow, save_path: Path | None = None) -> None:
    """Plot LBM velocity profile vs. analytical Hagen-Poiseuille parabola."""
    _, ux, _ = sim.get_fields()

    # Fluid nodes
    j = np.arange(1, sim.ny - 1)
    u_lbm  = ux[0, j]                         # x is uniform (periodic)
    u_anal = sim.analytical_u(j.astype(float))

    # Normalised coordinates:  ŷ/H  ∈ [0, 1]
    y_norm = (j - 0.5) / sim.H

    l2 = sim.l2_error()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        f"D2Q9 LBM — Poiseuille Flow   Re={sim.Re:.1f}   step={sim.step:,}\n"
        f"ny={sim.ny}   u_max={sim.u_max:.4f}   L2 error={l2:.2e}",
        fontsize=12,
    )

    # (a) Velocity profile
    ax = axes[0]
    ax.plot(u_lbm,  y_norm, "b-",  lw=2.0, label="LBM")
    ax.plot(u_anal, y_norm, "r--", lw=1.5, label="Analytical (Hagen-Poiseuille)")
    ax.set_xlabel("u  (lattice units)")
    ax.set_ylabel("ŷ / H")
    ax.set_title("Velocity Profile u(y)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (b) Absolute error
    ax = axes[1]
    ax.plot(np.abs(u_lbm - u_anal), y_norm, "g-", lw=1.5)
    ax.set_xlabel("|u_LBM − u_analytical|")
    ax.set_ylabel("ŷ / H")
    ax.set_title(f"Absolute Error   (L2 = {l2:.2e})")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")
    plt.show()


def plot_convergence(sim: PoiseuilleFlow,
                     errors: list[float],
                     steps: list[int],
                     save_path: Path | None = None) -> None:
    """Plot L2 error vs. time to show convergence to steady state."""
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(steps, errors, "b-", lw=1.5)
    ax.set_xlabel("Time step")
    ax.set_ylabel("L2 relative error")
    ax.set_title(
        f"Convergence — Poiseuille Flow   Re={sim.Re:.1f}   ny={sim.ny}"
    )
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Grid-refinement study (optional)
# ---------------------------------------------------------------------------

def grid_refinement_study(
    ny_list: list[int] | None = None,
    Re: float = 10.0,
    u_max: float = 0.02,
    steps: int = 5000,
) -> None:
    """
    Run Poiseuille flow for several grid resolutions and plot L2 error vs. Δy.
    Demonstrates second-order spatial accuracy of LBM with bounce-back BCs.
    """
    if ny_list is None:
        ny_list = [16, 32, 64, 128, 256]

    dy_list: list[float] = []
    err_list: list[float] = []

    for ny in ny_list:
        sim = PoiseuilleFlow(nx=4, ny=ny, Re=Re, u_max=u_max)
        sim.advance(steps)
        dy = 1.0 / (ny - 1)
        err = sim.l2_error()
        dy_list.append(dy)
        err_list.append(err)
        print(f"  ny={ny:4d}  Δy={dy:.4f}  L2={err:.3e}")

    # Fit slope in log-log
    slope, intercept = np.polyfit(np.log(dy_list), np.log(err_list), 1)
    print(f"\n  Convergence order ≈ {slope:.2f}  (expected ≈ 2.0)")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(dy_list, err_list, "bo-", ms=6, lw=1.5, label="LBM L2 error")
    # Reference line O(Δy²)
    dy_ref = np.array([dy_list[0], dy_list[-1]])
    ax.loglog(dy_ref, np.exp(intercept) * dy_ref**slope, "r--",
              label=f"O(Δy^{slope:.2f})")
    ax.set_xlabel("Grid spacing Δy = 1/(ny−1)")
    ax.set_ylabel("L2 relative error")
    ax.set_title(f"Grid Refinement Study — Poiseuille Flow  Re={Re:.1f}")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="D2Q9 LBM Poiseuille Channel Flow"
    )
    parser.add_argument("--nx",    type=int,   default=8,
                        help="Channel length (default 8, periodic)")
    parser.add_argument("--ny",    type=int,   default=64,
                        help="Grid height incl. walls (default 64)")
    parser.add_argument("--Re",    type=float, default=10.0,
                        help="Reynolds number (default 10)")
    parser.add_argument("--umax",  type=float, default=0.05,
                        help="Target max velocity in lattice units (default 0.05)")
    parser.add_argument("--steps", type=int,   default=10000,
                        help="Time steps to run (default 10 000)")
    parser.add_argument("--refinement", action="store_true",
                        help="Run grid refinement study instead")
    parser.add_argument("--no-show", action="store_true",
                        help="Save figures instead of displaying")
    parser.add_argument("--out-dir", type=str, default="results")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    if args.refinement:
        grid_refinement_study(Re=args.Re, u_max=args.umax)
        return

    sim = PoiseuilleFlow(nx=args.nx, ny=args.ny, Re=args.Re, u_max=args.umax)

    # Run with convergence tracking
    errors: list[float] = []
    step_log: list[int] = []
    report = max(args.steps // 20, 1)

    print(f"\nRunning {args.steps:,} steps …")
    step = 0
    while step < args.steps:
        chunk = min(report, args.steps - step)
        sim.advance(chunk)
        step += chunk
        err = sim.l2_error()
        errors.append(err)
        step_log.append(sim.step)
        print(f"  step {sim.step:>7,}   L2 error = {err:.4e}")

    print(f"\nFinal L2 error : {errors[-1]:.4e}")

    suffix = f"Re{args.Re:.0f}_ny{args.ny}"
    sp1 = out_dir / f"poiseuille_profile_{suffix}.png" if args.no_show else None
    sp2 = out_dir / f"poiseuille_conv_{suffix}.png"    if args.no_show else None

    plot_profile(sim, save_path=sp1)
    plot_convergence(sim, errors, step_log, save_path=sp2)
    print("Done.")


if __name__ == "__main__":
    main()
