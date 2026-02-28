"""
Lattice-Boltzmann Method (LBM) — D2Q9 Lid-Driven Cavity Flow
=============================================================
D2Q9 model with BGK collision operator.
Benchmark: 2D lid-driven cavity (Re = 100..1000).

References:
  - Succi, S. (2001). The Lattice Boltzmann Equation for Fluid Dynamics and Beyond.
  - Krüger et al. (2017). The Lattice Boltzmann Method: Principles and Practice.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path


# ---------------------------------------------------------------------------
# D2Q9 lattice constants
# ---------------------------------------------------------------------------

# Velocity directions: e[i] = (cx, cy)
C = np.array([
    [0,  0],   # 0: rest
    [1,  0],   # 1: E
    [0,  1],   # 2: N
    [-1, 0],   # 3: W
    [0, -1],   # 4: S
    [1,  1],   # 5: NE
    [-1, 1],   # 6: NW
    [-1,-1],   # 7: SW
    [1, -1],   # 8: SE
], dtype=float)

# Weights
W = np.array([4/9,
              1/9, 1/9, 1/9, 1/9,
              1/36, 1/36, 1/36, 1/36])

# Opposite directions for bounce-back
OPPOSITE = [0, 3, 4, 1, 2, 7, 8, 5, 6]

Q = 9   # number of velocities
CS2 = 1.0 / 3.0  # speed of sound squared (lattice units)


# ---------------------------------------------------------------------------
# Equilibrium distribution
# ---------------------------------------------------------------------------

def equilibrium(rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    """Compute Maxwell-Boltzmann equilibrium f_eq[i, x, y]."""
    u2 = ux**2 + uy**2
    feq = np.empty((Q,) + rho.shape)
    for i in range(Q):
        cu = C[i, 0] * ux + C[i, 1] * uy
        feq[i] = W[i] * rho * (1.0 + cu / CS2 + cu**2 / (2 * CS2**2) - u2 / (2 * CS2))
    return feq


# ---------------------------------------------------------------------------
# LBM simulation class
# ---------------------------------------------------------------------------

class LBM:
    """2D D2Q9 LBM solver for lid-driven cavity flow."""

    def __init__(
        self,
        nx: int = 128,
        ny: int = 128,
        Re: float = 400.0,
        u_lid: float = 0.1,
    ):
        """
        Parameters
        ----------
        nx, ny  : grid size (lattice units)
        Re      : Reynolds number  Re = u_lid * ny / nu
        u_lid   : lid velocity (lattice units, keep < 0.3 for stability)
        """
        self.nx = nx
        self.ny = ny
        self.u_lid = u_lid
        self.Re = Re

        # Kinematic viscosity → relaxation time
        self.nu = u_lid * ny / Re
        self.tau = self.nu / CS2 + 0.5   # BGK relaxation time
        self.omega = 1.0 / self.tau       # inverse relaxation time

        print(f"LBM D2Q9 — Lid-Driven Cavity")
        print(f"  Grid    : {nx} x {ny}")
        print(f"  Re      : {Re:.1f}")
        print(f"  u_lid   : {u_lid:.4f}")
        print(f"  nu      : {self.nu:.6f}  tau : {self.tau:.4f}")

        # --- distribution functions ---
        # f[q, x, y]
        self.f = np.ones((Q, nx, ny)) / Q  # uniform initial state
        rho0 = np.ones((nx, ny))
        ux0  = np.zeros((nx, ny))
        uy0  = np.zeros((nx, ny))
        self.f = equilibrium(rho0, ux0, uy0)

        # --- solid wall mask (bounce-back nodes) ---
        self.solid = np.zeros((nx, ny), dtype=bool)
        self.solid[:, 0]  = True   # bottom wall
        self.solid[:, -1] = True   # top wall  (moving lid applied separately)
        self.solid[0, :]  = True   # left wall
        self.solid[-1, :] = True   # right wall

        self.step = 0

    # ------------------------------------------------------------------
    # Core LBM steps
    # ------------------------------------------------------------------

    def _macroscopic(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute density and velocity from f."""
        rho = self.f.sum(axis=0)
        ux  = (C[:, 0, np.newaxis, np.newaxis] * self.f).sum(axis=0) / rho
        uy  = (C[:, 1, np.newaxis, np.newaxis] * self.f).sum(axis=0) / rho
        return rho, ux, uy

    def _collide(self, rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> None:
        """BGK collision: relax toward equilibrium."""
        feq = equilibrium(rho, ux, uy)
        self.f += self.omega * (feq - self.f)

    def _stream(self) -> None:
        """Streaming: shift each population along its velocity."""
        for i in range(Q):
            self.f[i] = np.roll(self.f[i], int(C[i, 0]), axis=0)
            self.f[i] = np.roll(self.f[i], int(C[i, 1]), axis=1)

    def _boundary(self) -> None:
        """
        Bounce-back on walls, moving-lid condition on top wall (y = ny-1).
        Uses the non-equilibrium bounce-back (Zou & He, 1997) for the lid.
        """
        # --- Stationary walls: full bounce-back ---
        # bottom (y=0), left (x=0), right (x=-1)
        for wall_mask in [
            self.solid[:, 0:1],    # bottom
            self.solid[0:1, :],    # left
            self.solid[-1:, :],    # right
        ]:
            pass  # handled below via solid mask

        # Simple bounce-back on all solid nodes
        f_tmp = self.f.copy()
        for i in range(Q):
            j = OPPOSITE[i]
            self.f[i, self.solid] = f_tmp[j, self.solid]

        # --- Moving lid (top wall, y = ny-1): Zou-He velocity BC ---
        # ux = u_lid, uy = 0 prescribed
        y = self.ny - 1
        rho_lid = (self.f[0, :, y]
                   + self.f[1, :, y] + self.f[3, :, y]
                   + 2.0 * (self.f[2, :, y] + self.f[5, :, y] + self.f[6, :, y])
                   ) / (1.0 + self.u_lid)

        ru = rho_lid * self.u_lid
        self.f[4, :, y] = self.f[2, :, y] - (2.0 / 3.0) * rho_lid * 0.0  # uy=0 correction
        self.f[8, :, y] = self.f[6, :, y] + 0.5 * (self.f[1, :, y] - self.f[3, :, y]) + 0.5 * ru - (1.0 / 6.0) * rho_lid * 0.0
        self.f[7, :, y] = self.f[5, :, y] - 0.5 * (self.f[1, :, y] - self.f[3, :, y]) - 0.5 * ru - (1.0 / 6.0) * rho_lid * 0.0

        # Correct bounce-back for south-pointing populations at lid
        self.f[4, :, y] = self.f[2, :, y]
        self.f[7, :, y] = self.f[5, :, y] - ru / 6.0
        self.f[8, :, y] = self.f[6, :, y] + ru / 6.0

    def advance(self, n_steps: int = 1) -> None:
        """Advance the simulation by n_steps time steps."""
        for _ in range(n_steps):
            rho, ux, uy = self._macroscopic()
            self._collide(rho, ux, uy)
            self._stream()
            self._boundary()
            self.step += 1

    def get_fields(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (rho, ux, uy) arrays shaped (nx, ny)."""
        return self._macroscopic()


# ---------------------------------------------------------------------------
# Visualization helpers
# ---------------------------------------------------------------------------

def plot_results(sim: LBM, save_path: Path | None = None) -> None:
    """Plot velocity magnitude, streamlines, and vorticity."""
    rho, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)

    # Transpose for (y, x) imshow convention
    ux_T  = ux.T
    uy_T  = uy.T
    spd_T = speed.T

    # Vorticity: duy/dx - dux/dy
    dvdx = np.gradient(uy, axis=0)
    dudy = np.gradient(ux, axis=1)
    vort_T = (dvdx - dudy).T

    x = np.arange(sim.nx)
    y = np.arange(sim.ny)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(
        f"D2Q9 LBM — Lid-Driven Cavity  Re={sim.Re:.0f}  "
        f"step={sim.step:,}  ({sim.nx}×{sim.ny})",
        fontsize=13
    )

    # (a) Speed contour
    ax = axes[0]
    im = ax.imshow(spd_T, origin="lower", cmap="inferno",
                   extent=[0, sim.nx, 0, sim.ny], aspect="equal")
    ax.set_title("Velocity Magnitude  |u|")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, shrink=0.8)

    # (b) Streamlines
    ax = axes[1]
    lw = 2.0 * spd_T / (spd_T.max() + 1e-10)
    ax.streamplot(x, y, ux_T, uy_T, color=spd_T, cmap="cool",
                  linewidth=lw, density=1.5, arrowsize=1.0)
    ax.set_xlim(0, sim.nx); ax.set_ylim(0, sim.ny)
    ax.set_title("Streamlines")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_aspect("equal")

    # (c) Vorticity
    ax = axes[2]
    vmax = np.percentile(np.abs(vort_T), 98)
    im2 = ax.imshow(vort_T, origin="lower", cmap="RdBu_r",
                    vmin=-vmax, vmax=vmax,
                    extent=[0, sim.nx, 0, sim.ny], aspect="equal")
    ax.set_title("Vorticity  ∂v/∂x − ∂u/∂y")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.colorbar(im2, ax=ax, shrink=0.8)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")
    plt.show()


def plot_centerline(sim: LBM, save_path: Path | None = None) -> None:
    """
    Compare u-velocity along vertical centerline and v-velocity along
    horizontal centerline against Ghia et al. (1982) benchmark data.
    """
    _, ux, uy = sim.get_fields()
    nx, ny = sim.nx, sim.ny

    # Centerline slices (normalized coordinates)
    y_norm = np.linspace(0, 1, ny)
    x_norm = np.linspace(0, 1, nx)

    u_center = ux[nx // 2, :] / sim.u_lid    # u/U along x=0.5
    v_center = uy[:, ny // 2] / sim.u_lid    # v/U along y=0.5

    # Ghia et al. (1982) Re=400 benchmark
    ghia_y_u = [0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
                0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
                0.9688, 0.9766, 1.0000]
    ghia_u400 = [0.00000, -0.08186, -0.09266, -0.10338, -0.14612, -0.24299,
                 -0.32726, -0.17119, -0.11477, 0.02135, 0.16256, 0.29093,
                 0.55892, 0.61756, 0.68439, 0.75837, 1.00000]

    ghia_x_v = [0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
                0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
                0.9609, 0.9688, 1.0000]
    ghia_v400 = [0.00000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077,
                 0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914,
                 -0.10313, -0.08864, -0.07391, -0.05906, 0.00000]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Centerline Velocity — Re={sim.Re:.0f}  (Ghia 1982 benchmark: Re=400)",
                 fontsize=12)

    ax1.plot(u_center, y_norm, "b-", lw=1.5, label="LBM")
    if abs(sim.Re - 400) < 50:
        ax1.plot(ghia_u400, ghia_y_u, "ro", ms=5, label="Ghia et al. (1982)")
    ax1.axvline(0, color="k", lw=0.5, ls="--")
    ax1.set_xlabel("u / U"); ax1.set_ylabel("y / H")
    ax1.set_title("u-velocity at x = 0.5")
    ax1.legend(); ax1.grid(True, alpha=0.3)

    ax2.plot(x_norm, v_center, "b-", lw=1.5, label="LBM")
    if abs(sim.Re - 400) < 50:
        ax2.plot(ghia_x_v, ghia_v400, "ro", ms=5, label="Ghia et al. (1982)")
    ax2.axhline(0, color="k", lw=0.5, ls="--")
    ax2.set_xlabel("x / L"); ax2.set_ylabel("v / U")
    ax2.set_title("v-velocity at y = 0.5")
    ax2.legend(); ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import argparse

    parser = argparse.ArgumentParser(description="D2Q9 LBM Lid-Driven Cavity")
    parser.add_argument("--nx",    type=int,   default=128,  help="Grid width  (default 128)")
    parser.add_argument("--ny",    type=int,   default=128,  help="Grid height (default 128)")
    parser.add_argument("--Re",    type=float, default=400.0, help="Reynolds number (default 400)")
    parser.add_argument("--steps", type=int,   default=20000, help="Total time steps (default 20000)")
    parser.add_argument("--plot-every", type=int, default=0,
                        help="Plot interval (0 = only at end)")
    parser.add_argument("--no-show", action="store_true",
                        help="Save figures to files instead of displaying")
    parser.add_argument("--out-dir", type=str, default="results",
                        help="Output directory for figures")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    sim = LBM(nx=args.nx, ny=args.ny, Re=args.Re, u_lid=0.1)

    print(f"\nRunning {args.steps:,} steps …")
    step = 0
    report = max(args.steps // 10, 1)

    while step < args.steps:
        chunk = min(report, args.steps - step)
        sim.advance(chunk)
        step += chunk

        rho, ux, uy = sim.get_fields()
        speed_max = float(np.sqrt(ux**2 + uy**2).max())
        print(f"  step {sim.step:>8,}  |u|_max = {speed_max:.4f}")

        if args.plot_every and step % args.plot_every == 0:
            if args.no_show:
                plot_results(sim, save_path=out_dir / f"fields_{sim.step:07d}.png")
            else:
                plot_results(sim)

    print("\nFinal plots …")
    suffix = f"Re{int(args.Re)}_{args.nx}x{args.ny}"
    save1 = out_dir / f"fields_{suffix}.png"   if args.no_show else None
    save2 = out_dir / f"centerline_{suffix}.png" if args.no_show else None

    plot_results(sim, save_path=save1)
    plot_centerline(sim, save_path=save2)
    print("Done.")


if __name__ == "__main__":
    main()
