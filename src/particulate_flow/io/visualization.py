"""Visualization helpers: plot_fields and plot_particles for LBMDEMSolver."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    from ..lbm_dem import LBMDEMSolver


def plot_fields(sim: "LBMDEMSolver", save_path: Path | None = None) -> None:
    """Plot velocity magnitude, streamlines, and vorticity with particle positions."""
    rho, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)

    ux_T = ux.T
    uy_T = uy.T
    spd_T = speed.T
    dvdx = np.gradient(uy, axis=0)
    dudy = np.gradient(ux, axis=1)
    vort_T = (dvdx - dudy).T

    x = np.arange(sim.nx)
    y = np.arange(sim.ny)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    fig.suptitle(
        f"LBM-DEM  Re={sim.Re:.0f}  step={sim.step_count:,}  "
        f"{sim.n_p} particles  ({sim.nx}×{sim.ny})",
        fontsize=12,
    )

    ax = axes[0]
    im = ax.imshow(spd_T, origin="lower", cmap="inferno", extent=[0, sim.nx, 0, sim.ny], aspect="auto")
    _draw_cylinders(ax, sim)
    _draw_particles(ax, sim)
    ax.set_title("Velocity Magnitude |u|")
    ax.set_xlabel("x [lattice]")
    ax.set_ylabel("y [lattice]")
    fig.colorbar(im, ax=ax, shrink=0.7)

    ax = axes[1]
    lw = 1.5 * spd_T / (spd_T.max() + 1e-12)
    ax.streamplot(x, y, ux_T, uy_T, color=spd_T, cmap="cool", linewidth=lw, density=1.2, arrowsize=0.8)
    _draw_cylinders(ax, sim)
    _draw_particles(ax, sim)
    ax.set_xlim(0, sim.nx)
    ax.set_ylim(0, sim.ny)
    ax.set_title("Streamlines")
    ax.set_xlabel("x [lattice]")
    ax.set_aspect("auto")

    ax = axes[2]
    vmax = float(np.percentile(np.abs(vort_T), 98)) + 1e-12
    im2 = ax.imshow(vort_T, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                    extent=[0, sim.nx, 0, sim.ny], aspect="auto")
    _draw_cylinders(ax, sim)
    _draw_particles(ax, sim)
    ax.set_title("Vorticity  ∂v/∂x − ∂u/∂y")
    ax.set_xlabel("x [lattice]")
    fig.colorbar(im2, ax=ax, shrink=0.7)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")
    plt.show()


def plot_particles(sim: "LBMDEMSolver", save_path: Path | None = None) -> None:
    """Plot particle positions coloured by speed, overlaid on fluid speed."""
    rho, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.imshow(speed.T, origin="lower", cmap="Blues", extent=[0, sim.nx, 0, sim.ny],
              aspect="auto", alpha=0.7)
    _draw_cylinders(ax, sim)

    p_speeds = np.linalg.norm(sim.vel, axis=1)
    sc = ax.scatter(sim.pos[:, 0], sim.pos[:, 1], c=p_speeds, cmap="hot_r",
                    s=(sim.r_p * 4) ** 2, edgecolors="k", linewidths=0.5, zorder=5)
    fig.colorbar(sc, ax=ax, label="Particle speed [lattice]", shrink=0.8)
    ax.set_xlim(0, sim.nx)
    ax.set_ylim(0, sim.ny)
    ax.set_title(f"Particles (n={sim.n_p}) + Fluid speed  |  step={sim.step_count:,}")
    ax.set_xlabel("x [lattice]")
    ax.set_ylabel("y [lattice]")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")
    plt.show()


def _draw_particles(ax: plt.Axes, sim: "LBMDEMSolver") -> None:
    """Overlay particle circles on an existing axes."""
    for i in range(sim.n_p):
        circle = plt.Circle((sim.pos[i, 0], sim.pos[i, 1]), sim.r_p,
                             color="white", linewidth=0.8, fill=False, zorder=3)
        ax.add_patch(circle)


def _draw_cylinders(ax: plt.Axes, sim: "LBMDEMSolver") -> None:
    """Overlay fixed solid cylinders on an existing axes."""
    for cx, cy, cr in sim.cylinders:
        circle = plt.Circle((cx, cy), cr, color="cyan", alpha=0.55, linewidth=1.0, zorder=4)
        ax.add_patch(circle)
