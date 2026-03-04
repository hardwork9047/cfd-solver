"""
LBM-DEM Coupled Solver: Fluid Flow with Particle Dynamics
==========================================================
Couples D2Q9 Lattice-Boltzmann (fluid) with Discrete Element Method (particles)
in a 2-D horizontal channel.

Physics
-------
Fluid  : Pressure-gradient-driven channel flow (body-force Poiseuille).
         BGK collision with Guo forcing scheme.
Particles: Settle under gravity, dragged by fluid (Stokes drag),
           collide with each other and walls (Hertz soft-sphere).
Coupling : Two-way.  Fluid → particle via interpolated drag force.
           Particle → fluid via distributed body force (Newton 3rd law).

References
----------
- Guo, Z., Zheng, C., & Shi, B. (2002). Discrete lattice effects on the
  forcing term in the lattice Boltzmann method. Phys. Rev. E 65, 046308.
- Krüger et al. (2017). The Lattice Boltzmann Method: Principles and Practice.
- Cundall & Strack (1979). A discrete numerical model for granular assemblies.
  Géotechnique 29(1), 47–65.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# D2Q9 lattice constants
# ---------------------------------------------------------------------------

# Velocity vectors  e[i] = (cx, cy)
C = np.array(
    [
        [0, 0],  # 0 rest
        [1, 0],  # 1 E
        [0, 1],  # 2 N
        [-1, 0],  # 3 W
        [0, -1],  # 4 S
        [1, 1],  # 5 NE
        [-1, 1],  # 6 NW
        [-1, -1],  # 7 SW
        [1, -1],  # 8 SE
    ],
    dtype=float,
)

W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
OPPOSITE = [0, 3, 4, 1, 2, 7, 8, 5, 6]
Q = 9
CS2 = 1.0 / 3.0  # lattice speed of sound squared


# ---------------------------------------------------------------------------
# LBM helpers
# ---------------------------------------------------------------------------


def _equilibrium(rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    """Maxwell-Boltzmann equilibrium  f_eq[q, x, y]."""
    u2 = ux**2 + uy**2
    feq = np.empty((Q,) + rho.shape)
    for i in range(Q):
        cu = C[i, 0] * ux + C[i, 1] * uy
        feq[i] = W[i] * rho * (1.0 + cu / CS2 + 0.5 * cu**2 / CS2**2 - 0.5 * u2 / CS2)
    return feq


def _guo_forcing(
    ux: np.ndarray,
    uy: np.ndarray,
    Fx: np.ndarray,
    Fy: np.ndarray,
    omega: float,
) -> np.ndarray:
    """
    Guo et al. (2002) forcing term  S[q, x, y].

    S_i = (1 - omega/2) * W_i * [ (e_i - u)/cs² + (e_i·u)/cs⁴ * e_i ] · F
    """
    S = np.empty((Q,) + ux.shape)
    for i in range(Q):
        cu = C[i, 0] * ux + C[i, 1] * uy
        term_x = (C[i, 0] - ux) / CS2 + cu / CS2**2 * C[i, 0]
        term_y = (C[i, 1] - uy) / CS2 + cu / CS2**2 * C[i, 1]
        S[i] = (1.0 - 0.5 * omega) * W[i] * (term_x * Fx + term_y * Fy)
    return S


# ---------------------------------------------------------------------------
# Main coupled solver
# ---------------------------------------------------------------------------


class LBMDEMSolver:
    """
    2-D Coupled LBM (fluid) + DEM (particles) solver.

    Domain
    ------
    Rectangular channel of size ``nx × ny`` lattice nodes.
    Periodic boundary conditions in x, no-slip bounce-back walls at y=0 and y=ny-1.

    Parameters (all lattice units)
    ------------------------------
    nx, ny            : int   — grid dimensions
    Re                : float — Reynolds number  Re = u_max * ny / nu
    u_max             : float — target Poiseuille centreline velocity (<0.3 for stability)
    n_particles       : int   — number of DEM particles
    particle_radius   : float — mean particle radius [lattice nodes]
    radius_variation  : float — uniform ±fraction for per-particle radius (e.g. 0.15 = ±15%)
    density_ratio     : float — ρ_particle / ρ_fluid
    gravity           : float — gravitational acceleration [lattice units]
    k_n               : float — normal contact stiffness (DEM, Hertz)
    damping           : float — viscous damping coefficient (DEM, 0–1)
    dem_substeps      : int   — DEM sub-steps per LBM step (for stability)
    seed              : int   — RNG seed for particle initialisation
    cylinder          : tuple[float, float, float] | None
                        Fixed solid cylinder as (cx, cy, radius) in lattice units.
                        ``None`` (default) means no cylinder.
    """

    def __init__(
        self,
        nx: int = 200,
        ny: int = 80,
        Re: float = 100.0,
        u_max: float = 0.05,
        n_particles: int = 20,
        particle_radius: float = 3.0,
        radius_variation: float = 0.0,
        density_ratio: float = 2.0,
        gravity: float = 2e-5,
        k_n: float = 50.0,
        damping: float = 0.4,
        dem_substeps: int = 4,
        seed: int = 42,
        cylinder: tuple | None = None,
    ):
        self.nx = nx
        self.ny = ny
        self.Re = Re
        self.u_max = u_max
        self.n_p = n_particles
        self.r_p = particle_radius  # representative (mean) radius
        self.radius_variation = radius_variation
        self.density_ratio = density_ratio
        self.g = gravity
        self.k_n = k_n
        self.damping = damping
        self.dem_substeps = dem_substeps
        self.step_count = 0
        self.cylinder = cylinder  # (cx, cy, cr) or None

        # --- Fluid parameters ---
        self.nu = u_max * ny / Re
        self.tau = self.nu / CS2 + 0.5
        self.omega = 1.0 / self.tau
        # Body force to sustain Poiseuille flow: ∂p/∂x = -8 μ u_max / ny²
        self.F_drive = 8.0 * self.nu * u_max / ny**2

        # Per-particle radii (uniform ±radius_variation around r_p)
        rng_rv = np.random.default_rng(seed + 999)
        if radius_variation > 0.0:
            self.radii = particle_radius * (
                1.0 + rng_rv.uniform(-radius_variation, radius_variation, n_particles)
            )
        else:
            self.radii = np.full(n_particles, particle_radius)

        # Per-particle masses (2-D disc: area × density_ratio)
        self.masses = density_ratio * np.pi * self.radii**2
        self.mass_p = float(np.mean(self.masses))  # representative scalar (kept for display)

        print("LBM-DEM Coupled Solver")
        print(f"  Grid         : {nx} × {ny}")
        print(f"  Re           : {Re:.1f}   nu = {self.nu:.5f}   tau = {self.tau:.4f}")
        print(f"  u_max        : {u_max:.4f}   F_drive = {self.F_drive:.2e}")
        print(
            f"  Particles    : {n_particles}  r = {particle_radius:.1f} ± "
            f"{radius_variation*100:.0f}%  density_ratio = {density_ratio:.1f}"
        )
        print(f"  Gravity (latt): {gravity:.2e}   DEM substeps = {dem_substeps}")
        if cylinder is not None:
            cx, cy, cr = cylinder[0], cylinder[1], cylinder[2]
            print(f"  Cylinder      : center=({cx:.1f}, {cy:.1f})  r={cr:.1f}")

        # --- LBM distribution functions f[q, x, y] ---
        rho0 = np.ones((nx, ny))
        self.f = _equilibrium(rho0, np.zeros((nx, ny)), np.zeros((nx, ny)))

        # Solid nodes: top and bottom walls
        self.solid = np.zeros((nx, ny), dtype=bool)
        self.solid[:, 0] = True
        self.solid[:, -1] = True

        # Solid nodes: fixed cylinder
        if cylinder is not None:
            cx, cy, cr = cylinder
            ix = np.arange(nx)
            iy = np.arange(ny)
            XX, YY = np.meshgrid(ix, iy, indexing="ij")
            self.solid |= (XX - cx) ** 2 + (YY - cy) ** 2 <= cr**2

        # Fluid body-force arrays (reset each step; includes driving + particle feedback)
        self.Fx = np.full((nx, ny), self.F_drive)
        self.Fy = np.zeros((nx, ny))

        # --- DEM arrays ---
        self.pos = np.empty((n_particles, 2))
        self.vel = np.zeros((n_particles, 2))
        self.forces_p = np.zeros((n_particles, 2))
        self._init_particles(seed)

    # ------------------------------------------------------------------
    # Particle initialisation
    # ------------------------------------------------------------------

    def _init_particles(self, seed: int) -> None:
        """Place particles randomly in the top 60 % of the channel (no overlap)."""
        rng = np.random.default_rng(seed)
        placed = 0
        for _ in range(200_000):
            if placed == self.n_p:
                break
            r_new = self.radii[placed]
            x = rng.uniform(r_new + 1, self.nx - r_new - 1)
            y = rng.uniform(self.ny * 0.4 + r_new, self.ny - r_new - 2)
            # Reject if overlapping an already-placed particle (sum-of-radii + 10%)
            if placed > 0:
                dists = np.hypot(x - self.pos[:placed, 0], y - self.pos[:placed, 1])
                min_clearance = 1.1 * (self.radii[:placed] + r_new)
                if (dists < min_clearance).any():
                    continue
            # Reject if overlapping the cylinder
            if self.cylinder is not None:
                cx, cy, cr = self.cylinder
                if np.hypot(x - cx, y - cy) < cr + r_new * 1.1:
                    continue
            self.pos[placed] = [x, y]
            placed += 1

        if placed < self.n_p:
            print(f"  Warning: only placed {placed}/{self.n_p} particles")
            self.n_p = placed
            self.pos = self.pos[:placed]
            self.vel = self.vel[:placed]
            self.forces_p = self.forces_p[:placed]
            self.radii = self.radii[:placed]
            self.masses = self.masses[:placed]

    # ------------------------------------------------------------------
    # LBM — macroscopic fields
    # ------------------------------------------------------------------

    def _macroscopic(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Density and velocity from distribution functions.
        Velocity is corrected for the body force: u = Σ(ci fi)/ρ + F/(2ρ).
        """
        rho = self.f.sum(axis=0)
        ux = (C[:, 0, np.newaxis, np.newaxis] * self.f).sum(axis=0) / rho
        uy = (C[:, 1, np.newaxis, np.newaxis] * self.f).sum(axis=0) / rho
        ux += 0.5 * self.Fx / rho
        uy += 0.5 * self.Fy / rho
        return rho, ux, uy

    def _lbm_step(self, rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> None:
        """One LBM step: collide → stream → bounce-back."""
        # Collision with Guo forcing
        feq = _equilibrium(rho, ux, uy)
        S = _guo_forcing(ux, uy, self.Fx, self.Fy, self.omega)
        self.f += self.omega * (feq - self.f) + S

        # Streaming (periodic via np.roll)
        for i in range(Q):
            self.f[i] = np.roll(self.f[i], int(C[i, 0]), axis=0)
            self.f[i] = np.roll(self.f[i], int(C[i, 1]), axis=1)

        # Bounce-back on solid nodes
        f_tmp = self.f.copy()
        for i in range(Q):
            self.f[i, self.solid] = f_tmp[OPPOSITE[i], self.solid]

    # ------------------------------------------------------------------
    # Coupling helpers
    # ------------------------------------------------------------------

    def _interp_velocity(self, x: float, y: float) -> tuple[float, float]:
        """Bilinear interpolation of (ux, uy) at position (x, y) [lattice]."""
        _, ux, uy = self._macroscopic()
        xi = int(np.floor(x)) % self.nx
        yi = int(np.clip(np.floor(y), 0, self.ny - 2))
        xi1 = (xi + 1) % self.nx
        yi1 = min(yi + 1, self.ny - 1)
        tx = x - np.floor(x)
        ty = y - np.floor(y)
        w = np.array([(1 - tx) * (1 - ty), tx * (1 - ty), (1 - tx) * ty, tx * ty])
        pts_x = [xi, xi1, xi, xi1]
        pts_y = [yi, yi, yi1, yi1]
        uf_x = sum(w[k] * ux[pts_x[k], pts_y[k]] for k in range(4))
        uf_y = sum(w[k] * uy[pts_x[k], pts_y[k]] for k in range(4))
        return float(uf_x), float(uf_y)

    def _stokes_drag(
        self,
        uf_x: float,
        uf_y: float,
        vp_x: float,
        vp_y: float,
        radius: float | None = None,
    ) -> tuple[float, float]:
        """
        Stokes drag force on a sphere (circle in 2-D):
            F_drag = 3π μ d (u_fluid − v_particle)
        where μ = ρ ν ≈ ν (ρ≈1) and d = 2 r (diameter).
        ``radius`` defaults to the representative r_p if not given.
        """
        r = radius if radius is not None else self.r_p
        coeff = 3.0 * np.pi * self.nu * 2.0 * r
        return coeff * (uf_x - vp_x), coeff * (uf_y - vp_y)

    def _distribute_force(
        self, x: float, y: float, fx: float, fy: float, radius: float | None = None
    ) -> None:
        """
        Add a point force (fx, fy) at (x, y) to the fluid body-force field
        using bilinear (Peskin) distribution over the 4 nearest nodes.
        ``radius`` sets the normalisation area (defaults to r_p).
        """
        xi = int(np.floor(x)) % self.nx
        yi = int(np.clip(np.floor(y), 0, self.ny - 2))
        xi1 = (xi + 1) % self.nx
        yi1 = min(yi + 1, self.ny - 1)
        tx = x - np.floor(x)
        ty = y - np.floor(y)
        w = [(1 - tx) * (1 - ty), tx * (1 - ty), (1 - tx) * ty, tx * ty]
        nodes = [(xi, yi), (xi1, yi), (xi, yi1), (xi1, yi1)]
        r = radius if radius is not None else self.r_p
        area = np.pi * r**2
        for wk, (ix, iy) in zip(w, nodes):
            if not self.solid[ix, iy]:
                self.Fx[ix, iy] += wk * fx / area
                self.Fy[ix, iy] += wk * fy / area

    # ------------------------------------------------------------------
    # DEM
    # ------------------------------------------------------------------

    def _dem_forces(self, dt_sub: float) -> np.ndarray:
        """
        Compute total force on each particle:
          gravity + Stokes drag (from fluid) + Hertz contact (particle/wall).

        ``dt_sub`` is unused here but kept for future damping models.
        """
        forces = np.zeros((self.n_p, 2))

        # 1. Gravity (buoyancy-corrected, per-particle mass)
        buoyancy_factor = 1.0 - 1.0 / self.density_ratio
        forces[:, 1] -= self.masses * self.g * buoyancy_factor

        # 2. Stokes drag from interpolated fluid velocity (per-particle radius)
        for i in range(self.n_p):
            uf_x, uf_y = self._interp_velocity(self.pos[i, 0], self.pos[i, 1])
            fd_x, fd_y = self._stokes_drag(
                uf_x, uf_y, self.vel[i, 0], self.vel[i, 1], radius=self.radii[i]
            )
            forces[i, 0] += fd_x
            forces[i, 1] += fd_y

        # 3. Particle–particle Hertz contact (per-particle radii / masses)
        for i in range(self.n_p):
            for j in range(i + 1, self.n_p):
                dp = self.pos[j] - self.pos[i]
                dist = float(np.linalg.norm(dp))
                min_dist = self.radii[i] + self.radii[j]
                if dist < min_dist and dist > 1e-10:
                    overlap = min_dist - dist
                    n = dp / dist
                    f_n = self.k_n * overlap**1.5
                    v_n = float(np.dot(self.vel[j] - self.vel[i], n))
                    avg_m = (self.masses[i] + self.masses[j]) / 2.0
                    f_damp = -self.damping * v_n * float(np.sqrt(self.k_n * avg_m))
                    f_total = (f_n + f_damp) * n
                    forces[i] -= f_total
                    forces[j] += f_total

        # 4. Wall contacts (Hertz, per-particle radius / mass)
        for i in range(self.n_p):
            wall_bot = self.radii[i] + 0.5
            wall_top = self.ny - 1.5 - self.radii[i]
            if self.pos[i, 1] < wall_bot:
                overlap = wall_bot - self.pos[i, 1]
                f_n = self.k_n * overlap**1.5
                v_n = -self.vel[i, 1]
                f_damp = -self.damping * v_n * float(np.sqrt(self.k_n * self.masses[i]))
                forces[i, 1] += f_n + f_damp
            if self.pos[i, 1] > wall_top:
                overlap = self.pos[i, 1] - wall_top
                f_n = self.k_n * overlap**1.5
                v_n = self.vel[i, 1]
                f_damp = -self.damping * v_n * float(np.sqrt(self.k_n * self.masses[i]))
                forces[i, 1] -= f_n + f_damp

        # 5. Cylinder contact (Hertz, per-particle radius / mass)
        if self.cylinder is not None:
            cx, cy, cr = self.cylinder
            for i in range(self.n_p):
                dx = self.pos[i, 0] - cx
                dy = self.pos[i, 1] - cy
                dist = float(np.hypot(dx, dy))
                min_dist = cr + self.radii[i]
                if dist < min_dist and dist > 1e-10:
                    overlap = min_dist - dist
                    nx_ = dx / dist
                    ny_ = dy / dist
                    f_n = self.k_n * overlap**1.5
                    v_n = self.vel[i, 0] * nx_ + self.vel[i, 1] * ny_
                    f_damp = -self.damping * v_n * float(np.sqrt(self.k_n * self.masses[i]))
                    f_mag = f_n + f_damp
                    forces[i, 0] += f_mag * nx_
                    forces[i, 1] += f_mag * ny_

        return forces

    def _dem_substep(self, dt: float) -> None:
        """One DEM sub-step via Velocity Verlet with time step ``dt``."""
        # Half-velocity update (per-particle mass)
        forces = self._dem_forces(dt)
        acc = forces / self.masses[:, np.newaxis]
        self.vel += 0.5 * dt * acc

        # Position update + periodic x BC
        self.pos += dt * self.vel
        self.pos[:, 0] %= self.nx

        # Second force evaluation
        forces_new = self._dem_forces(dt)
        acc_new = forces_new / self.masses[:, np.newaxis]
        self.vel += 0.5 * dt * acc_new

        # Clamp positions and apply restitution at walls (per-particle radius)
        for i in range(self.n_p):
            wall_bot = self.radii[i] + 0.5
            wall_top = self.ny - 1.5 - self.radii[i]
            if self.pos[i, 1] < wall_bot:
                self.pos[i, 1] = wall_bot
                if self.vel[i, 1] < 0:
                    self.vel[i, 1] *= -0.2
            if self.pos[i, 1] > wall_top:
                self.pos[i, 1] = wall_top
                if self.vel[i, 1] > 0:
                    self.vel[i, 1] *= -0.2
            # Clamp position outside cylinder surface
            if self.cylinder is not None:
                cx, cy, cr = self.cylinder
                dx = self.pos[i, 0] - cx
                dy = self.pos[i, 1] - cy
                dist = float(np.hypot(dx, dy))
                min_dist = cr + self.radii[i]
                if dist < min_dist and dist > 1e-10:
                    nx_ = dx / dist
                    ny_ = dy / dist
                    self.pos[i, 0] = cx + min_dist * nx_
                    self.pos[i, 1] = cy + min_dist * ny_
                    # Kill inward normal velocity
                    v_n = self.vel[i, 0] * nx_ + self.vel[i, 1] * ny_
                    if v_n < 0:
                        self.vel[i, 0] -= v_n * nx_ * (1 + 0.2)
                        self.vel[i, 1] -= v_n * ny_ * (1 + 0.2)

        self.forces_p = forces_new

    # ------------------------------------------------------------------
    # Public advance
    # ------------------------------------------------------------------

    def advance(self, n_steps: int = 1) -> None:
        """
        Advance the coupled simulation by ``n_steps`` LBM steps.

        Each LBM step:
          1. Compute drag forces on particles from current fluid field.
          2. Apply Newton-3rd-law back-reaction to fluid body-force field.
          3. Advance LBM one step (collide → stream → bounce-back).
          4. Advance DEM ``dem_substeps`` sub-steps (dt_sub = 1/dem_substeps).
        """
        dt_sub = 1.0 / self.dem_substeps

        for _ in range(n_steps):
            # --- Reset body force to driving force only ---
            self.Fx[:] = self.F_drive
            self.Fy[:] = 0.0

            # --- Particle → fluid back-reaction (per-particle radius) ---
            for i in range(self.n_p):
                uf_x, uf_y = self._interp_velocity(self.pos[i, 0], self.pos[i, 1])
                fd_x, fd_y = self._stokes_drag(
                    uf_x, uf_y, self.vel[i, 0], self.vel[i, 1], radius=self.radii[i]
                )
                # Newton 3rd law: -drag acts on the fluid
                self._distribute_force(
                    self.pos[i, 0], self.pos[i, 1], -fd_x, -fd_y, radius=self.radii[i]
                )

            # --- LBM step ---
            rho, ux, uy = self._macroscopic()
            self._lbm_step(rho, ux, uy)

            # --- DEM sub-steps ---
            for _ in range(self.dem_substeps):
                self._dem_substep(dt_sub)

            self.step_count += 1

    def get_fields(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (rho, ux, uy) arrays shaped (nx, ny)."""
        return self._macroscopic()


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------


def plot_fields(sim: LBMDEMSolver, save_path: Path | None = None) -> None:
    """Plot velocity magnitude, streamlines, and vorticity with particle positions."""
    rho, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)

    # Transpose to (y, x) for imshow
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

    # (a) Speed
    ax = axes[0]
    im = ax.imshow(
        spd_T,
        origin="lower",
        cmap="inferno",
        extent=[0, sim.nx, 0, sim.ny],
        aspect="auto",
    )
    _draw_particles(ax, sim)
    ax.set_title("Velocity Magnitude |u|")
    ax.set_xlabel("x [lattice]")
    ax.set_ylabel("y [lattice]")
    fig.colorbar(im, ax=ax, shrink=0.7)

    # (b) Streamlines
    ax = axes[1]
    lw = 1.5 * spd_T / (spd_T.max() + 1e-12)
    ax.streamplot(
        x, y, ux_T, uy_T, color=spd_T, cmap="cool", linewidth=lw, density=1.2, arrowsize=0.8
    )
    _draw_particles(ax, sim)
    ax.set_xlim(0, sim.nx)
    ax.set_ylim(0, sim.ny)
    ax.set_title("Streamlines")
    ax.set_xlabel("x [lattice]")
    ax.set_aspect("auto")

    # (c) Vorticity
    ax = axes[2]
    vmax = float(np.percentile(np.abs(vort_T), 98)) + 1e-12
    im2 = ax.imshow(
        vort_T,
        origin="lower",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        extent=[0, sim.nx, 0, sim.ny],
        aspect="auto",
    )
    _draw_particles(ax, sim)
    ax.set_title("Vorticity  ∂v/∂x − ∂u/∂y")
    ax.set_xlabel("x [lattice]")
    fig.colorbar(im2, ax=ax, shrink=0.7)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")
    plt.show()


def plot_particles(sim: LBMDEMSolver, save_path: Path | None = None) -> None:
    """Plot particle positions coloured by speed, overlaid on fluid speed."""
    rho, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.imshow(
        speed.T,
        origin="lower",
        cmap="Blues",
        extent=[0, sim.nx, 0, sim.ny],
        aspect="auto",
        alpha=0.7,
    )

    p_speeds = np.linalg.norm(sim.vel, axis=1)
    sc = ax.scatter(
        sim.pos[:, 0],
        sim.pos[:, 1],
        c=p_speeds,
        cmap="hot_r",
        s=(sim.r_p * 4) ** 2,
        edgecolors="k",
        linewidths=0.5,
        zorder=5,
    )
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


def _draw_particles(ax: plt.Axes, sim: LBMDEMSolver) -> None:
    """Overlay particle circles on an existing axes."""
    for i in range(sim.n_p):
        circle = plt.Circle(
            (sim.pos[i, 0], sim.pos[i, 1]),
            sim.r_p,
            color="white",
            linewidth=0.8,
            fill=False,
            zorder=3,
        )
        ax.add_patch(circle)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="2-D LBM-DEM coupled fluid-particle simulation")
    parser.add_argument("--nx", type=int, default=200, help="Grid width (default 200)")
    parser.add_argument("--ny", type=int, default=80, help="Grid height (default 80)")
    parser.add_argument("--Re", type=float, default=100.0, help="Reynolds number (default 100)")
    parser.add_argument(
        "--u-max",
        type=float,
        default=0.05,
        help="Max (centreline) fluid velocity in lattice units (default 0.05)",
    )
    parser.add_argument(
        "--n-particles", type=int, default=20, help="Number of DEM particles (default 20)"
    )
    parser.add_argument(
        "--radius", type=float, default=3.0, help="Particle radius in lattice nodes (default 3.0)"
    )
    parser.add_argument(
        "--density-ratio", type=float, default=2.0, help="ρ_particle / ρ_fluid (default 2.0)"
    )
    parser.add_argument(
        "--gravity",
        type=float,
        default=2e-5,
        help="Gravitational acceleration in lattice units (default 2e-5)",
    )
    parser.add_argument("--steps", type=int, default=10000, help="Total LBM steps (default 10000)")
    parser.add_argument(
        "--report-every", type=int, default=1000, help="Report interval (default 1000)"
    )
    parser.add_argument(
        "--no-show", action="store_true", help="Save figures to files instead of displaying"
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="results",
        help="Output directory for figures (default: results/)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    sim = LBMDEMSolver(
        nx=args.nx,
        ny=args.ny,
        Re=args.Re,
        u_max=args.u_max,
        n_particles=args.n_particles,
        particle_radius=args.radius,
        density_ratio=args.density_ratio,
        gravity=args.gravity,
    )

    print(f"\nRunning {args.steps:,} steps …")
    report = max(args.report_every, 1)
    step = 0

    while step < args.steps:
        chunk = min(report, args.steps - step)
        sim.advance(chunk)
        step += chunk

        rho, ux, uy = sim.get_fields()
        speed_max = float(np.sqrt(ux**2 + uy**2).max())
        p_ke = 0.5 * sim.mass_p * float(np.sum(sim.vel**2))
        print(
            f"  step {sim.step_count:>7,}  |u|_max = {speed_max:.5f}" f"  particle KE = {p_ke:.3e}"
        )

    print("\nFinal plots …")
    suffix = f"Re{int(args.Re)}_{args.nx}x{args.ny}_np{sim.n_p}"
    save_fields = out_dir / f"fields_{suffix}.png" if args.no_show else None
    save_parts = out_dir / f"particles_{suffix}.png" if args.no_show else None

    plot_fields(sim, save_path=save_fields)
    plot_particles(sim, save_path=save_parts)
    print("Done.")


if __name__ == "__main__":
    main()
