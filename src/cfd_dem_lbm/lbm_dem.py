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
from cfd.result_paths import program_results_dir
from .dem_solver import DEMSolver, PARTICLE_METHODS

try:  # Optional acceleration. The solver keeps a pure-NumPy/Python fallback.
    from numba import njit
except ImportError:  # pragma: no cover - depends on the local environment
    njit = None

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
FLUID_METHODS = ("lbm-bgk-guo", "lbm-trt-guo")
PARTICLE_FLUID_COUPLINGS = ("point_force", "solid_boundary", "immersed_boundary")
FLUID_ACCELERATORS = ("numpy", "numba", "auto")
COMPUTE_ACCELERATORS = ("numpy", "numba", "auto")
PARTICLE_SEARCH_METHODS = ("cell_list", "all_pairs")


if njit is not None:  # pragma: no cover - exercised only when numba is installed

    @njit(cache=True)
    def _dem_pair_loads_numba(
        pair_i,
        pair_j,
        pos,
        vel,
        radii,
        masses,
        omega_p,
        k_n,
        damping,
        rolling_friction,
        sliding_friction,
        tangential_damping,
        rolling_friction_coeff,
        rolling_damping,
        particle_attraction,
        particle_repulsion,
        attraction_strength,
        repulsion_strength,
        attraction_cutoff,
        repulsion_cutoff,
        attraction_min_gap,
        repulsion_min_gap,
        forces,
        torques,
    ):
        """Apply particle-pair DEM loads in compiled code."""
        for pair_idx in range(pair_i.shape[0]):
            i = pair_i[pair_idx]
            j = pair_j[pair_idx]
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dist = (dx * dx + dy * dy) ** 0.5
            min_dist = radii[i] + radii[j]
            if dist <= 1e-10:
                continue

            nx_ = dx / dist
            ny_ = dy / dist

            if dist < min_dist:
                overlap = min_dist - dist
                rvx = vel[j, 0] - vel[i, 0]
                rvy = vel[j, 1] - vel[i, 1]
                v_n = rvx * nx_ + rvy * ny_
                avg_m = 0.5 * (masses[i] + masses[j])
                f_n = k_n * overlap**1.5
                f_damp = -damping * v_n * (k_n * avg_m) ** 0.5
                f_mag = f_n + f_damp
                if f_mag < 0.0:
                    f_mag = 0.0

                fx = f_mag * nx_
                fy = f_mag * ny_
                forces[i, 0] -= fx
                forces[i, 1] -= fy
                forces[j, 0] += fx
                forces[j, 1] += fy

                if rolling_friction and f_mag > 0.0:
                    tx = -ny_
                    ty = nx_
                    v_t = rvx * tx + rvy * ty
                    v_t -= omega_p[j] * radii[j] + omega_p[i] * radii[i]
                    f_trial = -tangential_damping * (k_n * avg_m) ** 0.5 * v_t
                    f_limit = sliding_friction * f_mag
                    f_t = min(max(f_trial, -f_limit), f_limit)
                    ft_x = f_t * tx
                    ft_y = f_t * ty
                    forces[i, 0] -= ft_x
                    forces[i, 1] -= ft_y
                    forces[j, 0] += ft_x
                    forces[j, 1] += ft_y
                    torques[i] -= radii[i] * f_t
                    torques[j] -= radii[j] * f_t

                    trial_i = -rolling_damping * (k_n * masses[i]) ** 0.5 * radii[i] ** 2 * omega_p[i]
                    limit_i = rolling_friction_coeff * f_mag * radii[i]
                    torques[i] += min(max(trial_i, -limit_i), limit_i)
                    trial_j = -rolling_damping * (k_n * masses[j]) ** 0.5 * radii[j] ** 2 * omega_p[j]
                    limit_j = rolling_friction_coeff * f_mag * radii[j]
                    torques[j] += min(max(trial_j, -limit_j), limit_j)

            surface_gap = dist - min_dist
            if particle_attraction and attraction_strength > 0.0:
                if surface_gap <= attraction_cutoff:
                    h = surface_gap
                    if h < attraction_min_gap:
                        h = attraction_min_gap
                    if h < 1e-12:
                        h = 1e-12
                    r_eff = radii[i] * radii[j] / min_dist
                    f_attr = attraction_strength * r_eff / (6.0 * h**2)
                    forces[i, 0] += f_attr * nx_
                    forces[i, 1] += f_attr * ny_
                    forces[j, 0] -= f_attr * nx_
                    forces[j, 1] -= f_attr * ny_
            elif particle_repulsion and repulsion_strength > 0.0:
                if surface_gap <= repulsion_cutoff:
                    h = surface_gap
                    if h < repulsion_min_gap:
                        h = repulsion_min_gap
                    if h < 1e-12:
                        h = 1e-12
                    r_eff = radii[i] * radii[j] / min_dist
                    f_rep = repulsion_strength * r_eff / (6.0 * h**2)
                    forces[i, 0] -= f_rep * nx_
                    forces[i, 1] -= f_rep * ny_
                    forces[j, 0] += f_rep * nx_
                    forces[j, 1] += f_rep * ny_

else:
    _dem_pair_loads_numba = None


if njit is not None:  # pragma: no cover - exercised only when numba is installed

    @njit(cache=True)
    def _lbm_step_numba(
        f,
        rho,
        ux,
        uy,
        Fx,
        Fy,
        solid,
        omega,
        trt_omega_minus,
        method_id,
        c_int,
        w,
        opposite,
    ):
        """Fused D2Q9 collision, streaming, and bounce-back kernel."""
        q_count, nx, ny = f.shape
        post = np.empty_like(f)
        streamed = np.empty_like(f)
        for q in range(q_count):
            cx = c_int[q, 0]
            cy = c_int[q, 1]
            opp = opposite[q]
            ox = c_int[opp, 0]
            oy = c_int[opp, 1]
            weight = w[q]
            for x in range(nx):
                for y in range(ny):
                    u2 = ux[x, y] * ux[x, y] + uy[x, y] * uy[x, y]
                    cu = cx * ux[x, y] + cy * uy[x, y]
                    feq = weight * rho[x, y] * (
                        1.0
                        + cu / CS2
                        + 0.5 * cu * cu / (CS2 * CS2)
                        - 0.5 * u2 / CS2
                    )
                    term_x = (cx - ux[x, y]) / CS2 + cu / (CS2 * CS2) * cx
                    term_y = (cy - uy[x, y]) / CS2 + cu / (CS2 * CS2) * cy
                    forcing = (1.0 - 0.5 * omega) * weight * (
                        term_x * Fx[x, y] + term_y * Fy[x, y]
                    )
                    if method_id == 1:
                        cu_opp = ox * ux[x, y] + oy * uy[x, y]
                        feq_opp = w[opp] * rho[x, y] * (
                            1.0
                            + cu_opp / CS2
                            + 0.5 * cu_opp * cu_opp / (CS2 * CS2)
                            - 0.5 * u2 / CS2
                        )
                        f_plus = 0.5 * (f[q, x, y] + f[opp, x, y])
                        f_minus = 0.5 * (f[q, x, y] - f[opp, x, y])
                        feq_plus = 0.5 * (feq + feq_opp)
                        feq_minus = 0.5 * (feq - feq_opp)
                        post[q, x, y] = (
                            f[q, x, y]
                            - omega * (f_plus - feq_plus)
                            - trt_omega_minus * (f_minus - feq_minus)
                            + forcing
                        )
                    else:
                        post[q, x, y] = f[q, x, y] + omega * (feq - f[q, x, y]) + forcing

        for q in range(q_count):
            cx = c_int[q, 0]
            cy = c_int[q, 1]
            for x in range(nx):
                xd = (x + cx) % nx
                for y in range(ny):
                    yd = (y + cy) % ny
                    streamed[q, xd, yd] = post[q, x, y]

        out = streamed.copy()
        for q in range(q_count):
            opp = opposite[q]
            for x in range(nx):
                for y in range(ny):
                    if solid[x, y]:
                        out[q, x, y] = streamed[opp, x, y]
        return out

else:
    _lbm_step_numba = None


if njit is not None:  # pragma: no cover - exercised only when numba is installed

    @njit(cache=True)
    def _particle_solid_mask_numba(pos, radii, nx, ny, fixed_solid):
        particle_solid = np.zeros((nx, ny), dtype=np.bool_)
        for i in range(pos.shape[0]):
            radius = radii[i]
            x_pos = pos[i, 0]
            y_pos = pos[i, 1]
            x_min = int(np.floor(x_pos - radius))
            x_max = int(np.ceil(x_pos + radius))
            y_min = max(0, int(np.floor(y_pos - radius)))
            y_max = min(ny - 1, int(np.ceil(y_pos + radius)))
            r2 = radius * radius
            for x_raw in range(x_min, x_max + 1):
                ix = x_raw % nx
                dx = x_raw - x_pos
                for iy in range(y_min, y_max + 1):
                    dy = iy - y_pos
                    if dx * dx + dy * dy <= r2 and not fixed_solid[ix, iy]:
                        particle_solid[ix, iy] = True
        return particle_solid

    @njit(cache=True)
    def _ibm_markers_numba(pos, vel, omega_p, radii, marker_spacing, ny):
        total = 0
        for i in range(pos.shape[0]):
            total += max(8, int(np.ceil(2.0 * np.pi * radii[i] / marker_spacing)))
        marker_x = np.empty(total)
        marker_y = np.empty(total)
        marker_rx = np.empty(total)
        marker_ry = np.empty(total)
        marker_ds = np.empty(total)
        marker_owner = np.empty(total, dtype=np.int64)
        marker_ubx = np.empty(total)
        marker_uby = np.empty(total)
        cursor = 0
        for i in range(pos.shape[0]):
            radius = radii[i]
            n_markers = max(8, int(np.ceil(2.0 * np.pi * radius / marker_spacing)))
            ds = 2.0 * np.pi * radius / n_markers
            for m in range(n_markers):
                theta = 2.0 * np.pi * m / n_markers
                rx = radius * np.cos(theta)
                ry = radius * np.sin(theta)
                y = pos[i, 1] + ry
                if y < 0.5:
                    y = 0.5
                elif y > ny - 1.5:
                    y = ny - 1.5
                marker_x[cursor] = pos[i, 0] + rx
                marker_y[cursor] = y
                marker_rx[cursor] = rx
                marker_ry[cursor] = ry
                marker_ds[cursor] = ds
                marker_owner[cursor] = i
                marker_ubx[cursor] = vel[i, 0] - omega_p[i] * ry
                marker_uby[cursor] = vel[i, 1] + omega_p[i] * rx
                cursor += 1
        return marker_x, marker_y, marker_rx, marker_ry, marker_ds, marker_owner, marker_ubx, marker_uby

    @njit(cache=True)
    def _ibm_particle_reaction_numba(
        marker_owner,
        marker_rx,
        marker_ry,
        marker_fx,
        marker_fy,
        n_particles,
    ):
        forces = np.zeros((n_particles, 2))
        torques = np.zeros(n_particles)
        for k in range(marker_owner.shape[0]):
            owner = marker_owner[k]
            forces[owner, 0] -= marker_fx[k]
            forces[owner, 1] -= marker_fy[k]
            torques[owner] -= marker_rx[k] * marker_fy[k] - marker_ry[k] * marker_fx[k]
        return forces, torques

else:
    _particle_solid_mask_numba = None
    _ibm_markers_numba = None
    _ibm_particle_reaction_numba = None


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
    reynolds_length   : float | None — characteristic length for Re.  Defaults
                        to channel height for backward compatibility.
    flow_control      : str — ``"fixed_pressure"`` keeps a constant body force;
                        ``"target_max_velocity"`` adjusts the body force so the
                        fluid maximum speed approaches ``u_max``.
    flow_control_gain : float — relaxation gain for target-max-velocity control.
    n_particles       : int   — number of DEM particles
    particle_radius   : float — mean particle radius [lattice nodes]
    radius_variation  : float — uniform ±fraction for per-particle radius (e.g. 0.15 = ±15%)
    density_ratio     : float — ρ_particle / ρ_fluid
    gravity           : float — gravitational acceleration [lattice units]
    k_n               : float — normal contact stiffness (DEM, Hertz)
    damping           : float — viscous damping coefficient (DEM, 0–1)
    dem_substeps      : int   — DEM sub-steps per LBM step (for stability)
    seed              : int   — RNG seed for particle initialisation
    rolling_friction  : bool — enable angular motion, tangential friction, and rolling resistance
    sliding_friction  : float — Coulomb limit for tangential contact force
    tangential_damping : float — viscous damping scale for tangential slip
    rolling_friction_coeff : float — rolling-resistance moment coefficient
    rolling_damping   : float — viscous damping scale for angular velocity
    particle_attraction : bool — enable Hamaker-like particle-particle attraction
    particle_repulsion : bool — enable Hamaker-like particle-particle repulsion
    particle_source    : str — ``"initial"`` or ``"left_inlet"``.
                        ``"left_inlet"`` injects particles gradually from the
                        left boundary according to the inlet flow rate and
                        deletes particles after right outflow.
    source_volume_fraction : float | None
                        Particle area fraction in the inlet feed.  In
                        ``"left_inlet"`` mode, each step adds
                        ``source_volume_fraction * sum(u_x at left boundary)``
                        to the particle-area injection budget.
    attraction_strength : float
                        Dimensionless Hamaker-like attraction strength in lattice force units.
    repulsion_strength : float
                        Dimensionless Hamaker-like repulsion strength in lattice force units.
    attraction_cutoff : float — surface gap cutoff for attraction [lattice nodes]
    repulsion_cutoff  : float — surface gap cutoff for repulsion [lattice nodes]
    attraction_min_gap : float
                        Lower bound for surface gap to regularise the singularity.
    repulsion_min_gap : float
                        Lower bound for surface gap to regularise the singularity.
    cylinder          : tuple[float, float, float] | None
                        Backward-compatible single fixed cylinder as
                        (cx, cy, radius) in lattice units.
    cylinders         : list[tuple[float, float, float]] | None
                        Fixed solid cylinders as (cx, cy, radius).  Use this
                        to represent multiple pore-scale membrane obstacles.
    fluid_accelerator : str — ``"numpy"``, ``"numba"``, or ``"auto"``.
                        ``"numba"`` uses a compiled LBM step when numba is
                        installed, otherwise it falls back to NumPy with a
                        clear runtime flag.
    compute_accelerator : str — ``"numpy"``, ``"numba"``, or ``"auto"``.
                        Controls non-LBM kernels such as wall/cylinder DEM
                        loads, particle solid masks, and IBM marker bookkeeping.
    particle_search   : str — ``"cell_list"`` or ``"all_pairs"``.
                        ``"cell_list"`` is the scalable DEM neighbourhood
                        search used for production runs.
    """

    @staticmethod
    def _normalise_cylinders(
        cylinder: tuple | None,
        cylinders: list[tuple[float, float, float]] | tuple[tuple[float, float, float], ...] | None,
    ) -> list[tuple[float, float, float]]:
        """Return validated fixed-cylinder definitions."""
        items: list[tuple[float, float, float]] = []
        if cylinder is not None:
            items.append(cylinder)
        if cylinders is not None:
            items.extend(cylinders)

        normalised: list[tuple[float, float, float]] = []
        for item in items:
            if len(item) != 3:
                raise ValueError("each cylinder must be (x, y, radius)")
            cx, cy, cr = (float(item[0]), float(item[1]), float(item[2]))
            if cr <= 0.0:
                raise ValueError("cylinder radius must be positive")
            normalised.append((cx, cy, cr))
        return normalised

    def __init__(
        self,
        nx: int = 200,
        ny: int = 80,
        Re: float = 100.0,
        u_max: float = 0.05,
        reynolds_length: float | None = None,
        flow_control: str = "fixed_pressure",
        flow_control_gain: float = 0.2,
        fluid_method: str = "lbm-bgk-guo",
        particle_method: str = "dem-hertz",
        particle_fluid_coupling: str = "point_force",
        fluid_accelerator: str = "numpy",
        compute_accelerator: str = "auto",
        particle_search: str = "cell_list",
        spatial_dimension: int = 2,
        ibm_stiffness: float = 1.0,
        ibm_marker_spacing: float = 1.0,
        n_particles: int = 20,
        particle_radius: float = 3.0,
        radius_variation: float = 0.0,
        density_ratio: float = 2.0,
        gravity: float = 2e-5,
        k_n: float = 50.0,
        damping: float = 0.4,
        dem_substeps: int = 4,
        seed: int = 42,
        rolling_friction: bool = False,
        sliding_friction: float = 0.5,
        tangential_damping: float = 0.4,
        rolling_friction_coeff: float = 0.05,
        rolling_damping: float = 0.2,
        particle_attraction: bool = False,
        particle_repulsion: bool = False,
        attraction_strength: float = 1e-3,
        repulsion_strength: float = 1e-3,
        attraction_cutoff: float = 3.0,
        repulsion_cutoff: float = 3.0,
        attraction_min_gap: float = 0.05,
        repulsion_min_gap: float = 0.05,
        porous_resistance: bool = False,
        porous_resistance_coeff: float = 0.0,
        cylinder: tuple | None = None,
        cylinders: list[tuple[float, float, float]] | tuple[tuple[float, float, float], ...] | None = None,
        particle_source: str = "initial",
        source_volume_fraction: float | None = None,
    ):
        if particle_attraction and particle_repulsion:
            raise ValueError("particle_attraction and particle_repulsion are mutually exclusive")
        if particle_source not in {"initial", "left_inlet"}:
            raise ValueError("particle_source must be 'initial' or 'left_inlet'")
        if source_volume_fraction is not None and source_volume_fraction < 0.0:
            raise ValueError("source_volume_fraction must be non-negative")
        if reynolds_length is not None and reynolds_length <= 0.0:
            raise ValueError("reynolds_length must be positive")
        flow_aliases = {
            "fixed_pressure": "fixed_pressure",
            "constant_pressure": "fixed_pressure",
            "target_max_velocity": "target_max_velocity",
            "constant_flux": "target_max_velocity",
        }
        if flow_control not in flow_aliases:
            raise ValueError(
                "flow_control must be one of 'fixed_pressure', 'constant_pressure', "
                "'target_max_velocity', or 'constant_flux'"
            )
        flow_control = flow_aliases[flow_control]
        if not (0.0 < flow_control_gain <= 1.0):
            raise ValueError("flow_control_gain must be in the range (0, 1]")
        if fluid_method not in FLUID_METHODS:
            raise ValueError(f"fluid_method must be one of {FLUID_METHODS}")
        if particle_method not in PARTICLE_METHODS:
            raise ValueError(f"particle_method must be one of {PARTICLE_METHODS}")
        if particle_fluid_coupling not in PARTICLE_FLUID_COUPLINGS:
            raise ValueError(f"particle_fluid_coupling must be one of {PARTICLE_FLUID_COUPLINGS}")
        if fluid_accelerator not in FLUID_ACCELERATORS:
            raise ValueError(f"fluid_accelerator must be one of {FLUID_ACCELERATORS}")
        if compute_accelerator not in COMPUTE_ACCELERATORS:
            raise ValueError(f"compute_accelerator must be one of {COMPUTE_ACCELERATORS}")
        if particle_search not in PARTICLE_SEARCH_METHODS:
            raise ValueError(f"particle_search must be one of {PARTICLE_SEARCH_METHODS}")
        if spatial_dimension != 2:
            raise ValueError("only spatial_dimension=2 is currently implemented")
        if ibm_stiffness <= 0.0:
            raise ValueError("ibm_stiffness must be positive")
        if ibm_marker_spacing <= 0.0:
            raise ValueError("ibm_marker_spacing must be positive")
        if porous_resistance_coeff < 0.0:
            raise ValueError("porous_resistance_coeff must be non-negative")

        self.nx = nx
        self.ny = ny
        self.Re = Re
        self.u_max = u_max
        self.reynolds_length = reynolds_length if reynolds_length is not None else float(ny)
        self.flow_control = flow_control
        self.flow_control_gain = flow_control_gain
        self.fluid_method = fluid_method
        self.particle_method = particle_method
        self.particle_fluid_coupling = particle_fluid_coupling
        self.fluid_accelerator = fluid_accelerator
        self.compute_accelerator = compute_accelerator
        self.particle_search = particle_search
        self.spatial_dimension = spatial_dimension
        self.uses_numba_lbm = (
            fluid_accelerator == "numba"
            or (fluid_accelerator == "auto" and _lbm_step_numba is not None)
        ) and _lbm_step_numba is not None
        self.uses_numba_compute = (
            compute_accelerator == "numba"
            or (compute_accelerator == "auto" and _particle_solid_mask_numba is not None)
        ) and _particle_solid_mask_numba is not None
        self.ibm_stiffness = ibm_stiffness
        self.ibm_marker_spacing = ibm_marker_spacing
        self.n_p = n_particles
        self.r_p = particle_radius  # representative (mean) radius
        self.radius_variation = radius_variation
        self.density_ratio = density_ratio
        self.g = gravity
        self.k_n = k_n
        self.damping = damping
        self.dem_substeps = dem_substeps
        self.step_count = 0
        self.cylinders = self._normalise_cylinders(cylinder, cylinders)
        self.cylinder = self.cylinders[0] if self.cylinders else None
        self.rolling_friction = rolling_friction
        self.sliding_friction = sliding_friction
        self.tangential_damping = tangential_damping
        self.rolling_friction_coeff = rolling_friction_coeff
        self.rolling_damping = rolling_damping
        self.particle_attraction = particle_attraction
        self.particle_repulsion = particle_repulsion
        self.attraction_strength = attraction_strength
        self.repulsion_strength = repulsion_strength
        self.attraction_cutoff = attraction_cutoff
        self.repulsion_cutoff = repulsion_cutoff
        self.attraction_min_gap = attraction_min_gap
        self.repulsion_min_gap = repulsion_min_gap
        self.porous_resistance = porous_resistance
        self.porous_resistance_coeff = porous_resistance_coeff
        self.particle_source = particle_source
        self.source_volume_fraction = source_volume_fraction
        self.removed_particles = 0
        self.injected_particle_area = 0.0
        self.inlet_particle_area_budget = 0.0
        self.last_inlet_flow_rate = 0.0
        self.cumulative_inlet_flow_area = 0.0

        # --- Fluid parameters ---
        self.nu = u_max * self.reynolds_length / Re if Re > 0.0 else 0.0
        self.tau = self.nu / CS2 + 0.5
        self.omega = 1.0 / self.tau
        self.trt_magic_parameter = 3.0 / 16.0
        self.trt_omega_minus = self._trt_omega_minus()
        self._c_int = C.astype(np.int64)
        self._opposite = np.asarray(OPPOSITE, dtype=np.int64)
        # Body force to sustain Poiseuille flow: ∂p/∂x = -8 μ u_max / ny²
        self.F_drive = 8.0 * self.nu * u_max / ny**2
        self.initial_F_drive = self.F_drive

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
        self.inertias = 0.5 * self.masses * self.radii**2

        print("LBM-DEM Coupled Solver")
        print(f"  Grid         : {nx} × {ny}")
        print(
            f"  Re           : {Re:.1f}   L_ref = {self.reynolds_length:.3f}   "
            f"nu = {self.nu:.5f}   tau = {self.tau:.4f}"
        )
        print(
            f"  u_max target : {u_max:.4f}   F_drive = {self.F_drive:.2e}   "
            f"flow_control={flow_control}"
        )
        print(f"  Fluid method : {fluid_method}")
        print(
            f"  LBM accel.   : {fluid_accelerator}"
            f"{' (numba active)' if self.uses_numba_lbm else ' (numpy path)'}"
        )
        print(
            f"  Compute acc. : {compute_accelerator}"
            f"{' (numba active)' if self.uses_numba_compute else ' (numpy path)'}"
        )
        print(f"  DEM method   : {particle_method}")
        print(f"  DEM search   : {particle_search}")
        print(f"  Coupling     : {particle_fluid_coupling}")
        if particle_fluid_coupling == "immersed_boundary":
            print(
                f"  IBM          : stiffness={ibm_stiffness:.3g}  "
                f"marker_spacing={ibm_marker_spacing:.3g}"
            )
        print(
            f"  Particles    : {n_particles}  r = {particle_radius:.1f} ± "
            f"{radius_variation*100:.0f}%  density_ratio = {density_ratio:.1f}"
        )
        print(f"  Gravity (latt): {gravity:.2e}   DEM substeps = {dem_substeps}")
        if rolling_friction:
            print(
                "  Rolling fric. : enabled  "
                f"mu_t={sliding_friction:.2f}  mu_r={rolling_friction_coeff:.3f}"
            )
        else:
            print("  Rolling fric. : disabled")
        if particle_attraction:
            print(
                "  Attraction    : enabled  "
                f"A*={attraction_strength:.2e}  cutoff={attraction_cutoff:.2f}  "
                f"min_gap={attraction_min_gap:.3f}"
            )
        elif particle_repulsion:
            print(
                "  Repulsion     : enabled  "
                f"A*={repulsion_strength:.2e}  cutoff={repulsion_cutoff:.2f}  "
                f"min_gap={repulsion_min_gap:.3f}"
            )
        else:
            print("  Surface force : disabled")
        if porous_resistance and porous_resistance_coeff > 0.0:
            print(f"  Porous drag   : enabled  coeff={porous_resistance_coeff:.3g}")
        if self.cylinders:
            print(f"  Cylinders     : {len(self.cylinders)} fixed solids")
            for idx, (cx, cy, cr) in enumerate(self.cylinders, start=1):
                print(f"    #{idx:<2d} center=({cx:.1f}, {cy:.1f})  r={cr:.1f}")
        print(f"  Particle src. : {particle_source}")

        # --- LBM distribution functions f[q, x, y] ---
        rho0 = np.ones((nx, ny))
        self.f = _equilibrium(rho0, np.zeros((nx, ny)), np.zeros((nx, ny)))
        self._stream_x_src = [
            (np.arange(nx) - int(C[i, 0])) % nx
            for i in range(Q)
        ]
        self._stream_y_src = [
            (np.arange(ny) - int(C[i, 1])) % ny
            for i in range(Q)
        ]

        # Solid nodes: top/bottom walls and fixed cylinders.  Moving particle
        # solids are overlaid in solid-boundary coupling mode.
        self.fixed_solid = np.zeros((nx, ny), dtype=bool)
        self.fixed_solid[:, 0] = True
        self.fixed_solid[:, -1] = True

        # Solid nodes: fixed cylinders
        if self.cylinders:
            ix = np.arange(nx)
            iy = np.arange(ny)
            XX, YY = np.meshgrid(ix, iy, indexing="ij")
            for cx, cy, cr in self.cylinders:
                self.fixed_solid |= (XX - cx) ** 2 + (YY - cy) ** 2 <= cr**2
        self.particle_solid = np.zeros((nx, ny), dtype=bool)
        self.solid = self.fixed_solid.copy()

        # Fluid body-force arrays (reset each step; includes driving + particle feedback)
        self.Fx = np.full((nx, ny), self.F_drive)
        self.Fy = np.zeros((nx, ny))
        self.fluid_area = float(np.count_nonzero(~self.solid))

        # --- DEM arrays ---
        self.pos = np.empty((n_particles, 2))
        self.vel = np.zeros((n_particles, 2))
        self.omega_p = np.zeros(n_particles)
        self.forces_p = np.zeros((n_particles, 2))
        self.torques_p = np.zeros(n_particles)
        self.ibm_forces_p = np.zeros((n_particles, 2))
        self.ibm_torques_p = np.zeros(n_particles)
        self.total_particles_requested = n_particles
        self.generated_particles = 0
        self._pending_radii = np.empty(0)
        self._pending_masses = np.empty(0)
        self._pending_inertias = np.empty(0)
        self._inlet_cursor = 0
        self._rng_inlet = np.random.default_rng(seed + 2027)
        self._init_particles(seed)
        self._update_particle_solid_mask()
        pair_kernel = _dem_pair_loads_numba if particle_method == "dem-hertz" else None
        self.dem_solver = DEMSolver(
            self,
            pair_kernel=pair_kernel,
            contact_model=particle_method,
        )
        if self.particle_source != "left_inlet":
            self.generated_particles = self.n_p
        self.particle_area = float(
            np.sum(np.pi * self.radii**2) + np.sum(np.pi * self._pending_radii**2)
        )
        active_particle_area = float(np.sum(np.pi * self.radii**2))
        area_label = "Queued p.area:" if self.particle_source == "left_inlet" else "Particle frac.:"
        print(
            f"  {area_label} "
            f"{self.particle_area / self.fluid_area:.3f} "
            f"({self.particle_area:.1f}/{self.fluid_area:.1f} lattice area)"
        )
        if self.particle_source == "left_inlet":
            if self.source_volume_fraction is not None:
                print(f"  Source phi    : {self.source_volume_fraction:.3f}")
            print(
                "  Active p.frac.: "
                f"{active_particle_area / self.fluid_area:.3f} "
                f"({self.n_p} active, {len(self._pending_radii)} queued)"
            )

    # ------------------------------------------------------------------
    # Particle initialisation
    # ------------------------------------------------------------------

    def _overlaps_cylinder(self, x: float, y: float, radius: float, clearance: float) -> bool:
        """Return whether a particle disc overlaps any fixed cylinder."""
        for cx, cy, cr in self.cylinders:
            if np.hypot(x - cx, y - cy) < cr + radius * clearance:
                return True
        return False

    def _init_particles(self, seed: int) -> None:
        """Place particles without overlap, with a grid fallback for dense cases."""
        if self.particle_source == "left_inlet":
            self._init_particles_from_left_inlet()
            return

        rng = np.random.default_rng(seed)
        placed = 0

        def can_place(idx: int, x: float, y: float, clearance: float) -> bool:
            r_new = self.radii[idx]
            if y < r_new + 0.5 or y > self.ny - 1.5 - r_new:
                return False
            if placed > 0:
                dists = np.hypot(x - self.pos[:placed, 0], y - self.pos[:placed, 1])
                min_clearance = clearance * (self.radii[:placed] + r_new)
                if (dists < min_clearance).any():
                    return False
            if self._overlaps_cylinder(x, y, r_new, clearance):
                return False
            return True

        max_random_attempts = 200_000 if self.n_p <= 60 else 0
        for _ in range(max_random_attempts):
            if placed == self.n_p:
                break
            r_new = self.radii[placed]
            x = rng.uniform(r_new + 1, self.nx - r_new - 1)
            y = rng.uniform(self.ny * 0.4 + r_new, self.ny - r_new - 2)
            if not can_place(placed, x, y, clearance=1.1):
                continue
            self.pos[placed] = [x, y]
            placed += 1

        if placed < self.n_p:
            max_r = float(np.max(self.radii)) if self.n_p else self.r_p
            dx = 2.05 * max_r
            dy = np.sqrt(3.0) * max_r
            row = 0
            y = self.ny - 1.5 - max_r
            while placed < self.n_p and y >= max_r + 0.5:
                offset = (row % 2) * 0.5 * dx
                x = max_r + 0.5 + offset
                while placed < self.n_p and x <= self.nx - max_r - 0.5:
                    if can_place(placed, x, y, clearance=1.02):
                        self.pos[placed] = [x, y]
                        placed += 1
                    x += dx
                row += 1
                y -= dy

        if placed < self.n_p:
            print(f"  Warning: only placed {placed}/{self.n_p} particles")
            self.n_p = placed
            self.pos = self.pos[:placed]
            self.vel = self.vel[:placed]
            self.omega_p = self.omega_p[:placed]
            self.forces_p = self.forces_p[:placed]
            self.torques_p = self.torques_p[:placed]
            self.ibm_forces_p = self.ibm_forces_p[:placed]
            self.ibm_torques_p = self.ibm_torques_p[:placed]
            self.radii = self.radii[:placed]
            self.masses = self.masses[:placed]
            self.inertias = self.inertias[:placed]

    def _init_particles_from_left_inlet(self) -> None:
        """
        Queue all particles for gradual left-inlet injection.

        This mode is intended for filtration-style calculations: particles are
        supplied from the upstream boundary according to the inlet flow rate,
        rather than being pre-packed in the channel at time zero.
        """
        queued = self.n_p
        self._pending_radii = self.radii.copy()
        self._pending_masses = self.masses.copy()
        self._pending_inertias = self.inertias.copy()
        self.n_p = 0
        self.pos = self.pos[:0]
        self.vel = self.vel[:0]
        self.omega_p = self.omega_p[:0]
        self.forces_p = self.forces_p[:0]
        self.torques_p = self.torques_p[:0]
        self.ibm_forces_p = self.ibm_forces_p[:0]
        self.ibm_torques_p = self.ibm_torques_p[:0]
        self.radii = self.radii[:0]
        self.masses = self.masses[:0]
        self.inertias = self.inertias[:0]
        self.generated_particles = 0
        print(f"  Inlet queue   : 0 active at t=0; {queued} particles queued")

    def _inlet_velocity_at_y(self, y: float) -> np.ndarray:
        """Poiseuille-like inlet particle velocity at vertical coordinate ``y``."""
        channel_height = max(self.ny - 2.0, 1.0)
        eta = float(np.clip((y - 0.5) / channel_height, 0.0, 1.0))
        return np.array([4.0 * self.u_max * eta * (1.0 - eta), 0.0])

    def _left_boundary_flow_rate(self, radius: float | None = None) -> float:
        """Return the positive volumetric flow rate through the left boundary."""
        _, ux, _ = self._macroscopic()
        if radius is None:
            y_min = 1
            y_max = self.ny - 2
        else:
            y_min = int(np.ceil(radius + 0.5))
            y_max = int(np.floor(self.ny - 1.5 - radius))
        if y_max < y_min:
            return 0.0
        inlet_speed = np.maximum(ux[0, y_min : y_max + 1], 0.0)
        return float(np.sum(inlet_speed))

    def _sample_inlet_y(self, radius: float) -> float | None:
        """
        Sample an inlet y-coordinate from the local flow flux.

        A uniform particle concentration in the incoming fluid gives a
        particle-feed probability proportional to the local positive ux.
        If the flow has not developed yet, fall back to a deterministic uniform
        sweep so the source remains well behaved at startup.
        """
        y_min = int(np.ceil(radius + 0.5))
        y_max = int(np.floor(self.ny - 1.5 - radius))
        if y_max < y_min:
            return None

        _, ux, _ = self._macroscopic()
        y_nodes = np.arange(y_min, y_max + 1)
        weights = np.maximum(ux[0, y_nodes], 0.0)
        weight_sum = float(np.sum(weights))
        if weight_sum > 1e-12:
            idx = int(self._rng_inlet.choice(len(y_nodes), p=weights / weight_sum))
            y = float(y_nodes[idx]) + self._rng_inlet.uniform(-0.45, 0.45)
        else:
            idx = self._inlet_cursor % len(y_nodes)
            self._inlet_cursor += 1
            y = float(y_nodes[idx])
        return float(np.clip(y, radius + 0.5, self.ny - 1.5 - radius))

    def _can_place_inlet_particle(
        self,
        x: float,
        y: float,
        radius: float,
        n_existing: int | None = None,
    ) -> bool:
        """Return whether an inlet particle can be placed without overlap."""
        if y < radius + 0.5 or y > self.ny - 1.5 - radius:
            return False
        if self._overlaps_cylinder(x, y, radius, 1.05):
            return False
        n_check = self.n_p if n_existing is None else n_existing
        if n_check > 0:
            dists = np.hypot(x - self.pos[:n_check, 0], y - self.pos[:n_check, 1])
            if (dists < 1.02 * (self.radii[:n_check] + radius)).any():
                return False
        return True

    def _try_feed_left_inlet_particles(self, max_new: int = 8) -> None:
        """Inject queued particles according to the inlet-flow concentration budget."""
        if self.particle_source != "left_inlet":
            return
        phi = 0.0 if self.source_volume_fraction is None else self.source_volume_fraction
        if phi <= 0.0:
            return

        self.last_inlet_flow_rate = self._left_boundary_flow_rate()
        self.cumulative_inlet_flow_area += self.last_inlet_flow_rate
        if len(self._pending_radii) == 0:
            return
        self.inlet_particle_area_budget += phi * self.last_inlet_flow_rate

        added = 0
        attempts = 0
        max_attempts = max(12, 4 * max_new)
        while len(self._pending_radii) and added < max_new and attempts < max_attempts:
            radius = float(self._pending_radii[0])
            particle_area = float(np.pi * radius**2)
            if self.inlet_particle_area_budget < particle_area:
                break

            y = self._sample_inlet_y(radius)
            if y is None:
                break
            x = radius + 0.75
            attempts += 1
            if not self._can_place_inlet_particle(x, y, radius):
                continue
            self.pos = np.vstack([self.pos, np.array([[x, y]])])
            self.vel = np.vstack([self.vel, self._inlet_velocity_at_y(y)[np.newaxis, :]])
            self.omega_p = np.append(self.omega_p, 0.0)
            self.forces_p = np.vstack([self.forces_p, np.zeros((1, 2))])
            self.torques_p = np.append(self.torques_p, 0.0)
            self.ibm_forces_p = np.vstack([self.ibm_forces_p, np.zeros((1, 2))])
            self.ibm_torques_p = np.append(self.ibm_torques_p, 0.0)
            self.radii = np.append(self.radii, self._pending_radii[0])
            self.masses = np.append(self.masses, self._pending_masses[0])
            self.inertias = np.append(self.inertias, self._pending_inertias[0])
            self._pending_radii = self._pending_radii[1:]
            self._pending_masses = self._pending_masses[1:]
            self._pending_inertias = self._pending_inertias[1:]
            self.n_p += 1
            self.generated_particles += 1
            self.injected_particle_area += particle_area
            self.inlet_particle_area_budget -= particle_area
            added += 1

    def _delete_particles(self, mask: np.ndarray) -> None:
        """Delete particles selected by ``mask`` from all DEM arrays."""
        if mask.size == 0 or not bool(np.any(mask)):
            return
        keep = ~mask
        removed = int(np.count_nonzero(mask))
        self.pos = self.pos[keep]
        self.vel = self.vel[keep]
        self.omega_p = self.omega_p[keep]
        self.forces_p = self.forces_p[keep]
        self.torques_p = self.torques_p[keep]
        self.ibm_forces_p = self.ibm_forces_p[keep]
        self.ibm_torques_p = self.ibm_torques_p[keep]
        self.radii = self.radii[keep]
        self.masses = self.masses[keep]
        self.inertias = self.inertias[keep]
        self.n_p = int(self.pos.shape[0])
        self.removed_particles += removed
        self._update_particle_solid_mask()

    def _delete_right_outflow_particles(self) -> None:
        """Remove particles whose full disc has left the right outlet."""
        if self.particle_source != "left_inlet" or self.n_p == 0:
            return
        self._delete_particles(self.pos[:, 0] - self.radii > self.nx)

    def _update_particle_solid_mask(self) -> None:
        """Overlay active DEM particles onto the LBM solid mask when requested."""
        if self.particle_fluid_coupling != "solid_boundary" or self.n_p == 0:
            self.particle_solid = np.zeros_like(self.fixed_solid)
            self.solid = self.fixed_solid.copy()
            self._invalidate_macroscopic_cache()
            return

        particle_solid = self._particle_occupancy_mask()
        self.particle_solid = particle_solid
        self.solid = self.fixed_solid | self.particle_solid
        self._invalidate_macroscopic_cache()

    def _particle_occupancy_mask(self) -> np.ndarray:
        """Return lattice cells covered by DEM particles, excluding fixed solids."""
        if self.n_p == 0:
            return np.zeros_like(self.fixed_solid)
        if self.uses_numba_compute and _particle_solid_mask_numba is not None:
            particle_solid = _particle_solid_mask_numba(
                self.pos,
                self.radii,
                self.nx,
                self.ny,
                self.fixed_solid,
            )
        else:
            particle_solid = np.zeros_like(self.fixed_solid)
            for (x_pos, y_pos), radius in zip(self.pos, self.radii):
                x_min = int(np.floor(x_pos - radius))
                x_max = int(np.ceil(x_pos + radius))
                y_min = max(0, int(np.floor(y_pos - radius)))
                y_max = min(self.ny - 1, int(np.ceil(y_pos + radius)))
                if y_max < y_min:
                    continue
                x_indices = np.arange(x_min, x_max + 1)
                y_indices = np.arange(y_min, y_max + 1)
                xx = x_indices[:, np.newaxis] % self.nx
                yy = y_indices[np.newaxis, :]
                local = (x_indices[:, np.newaxis] - x_pos) ** 2 + (y_indices[np.newaxis, :] - y_pos) ** 2 <= radius**2
                particle_solid[xx, yy] |= local
            particle_solid &= ~self.fixed_solid
        return particle_solid

    def _apply_porous_resistance(self, ux: np.ndarray, uy: np.ndarray) -> None:
        """Apply Brinkman-like drag in cells occupied by retained particles."""
        if not self.porous_resistance or self.porous_resistance_coeff <= 0.0 or self.n_p == 0:
            return
        mask = self._particle_occupancy_mask()
        if not np.any(mask):
            return
        coeff = self.porous_resistance_coeff
        self.Fx[mask] -= coeff * ux[mask]
        self.Fy[mask] -= coeff * uy[mask]
        self._invalidate_macroscopic_cache()

    def dynamic_solid_fraction(self) -> float:
        """Return the current fraction of lattice nodes blocked by moving particles."""
        return float(np.count_nonzero(self.particle_solid) / self.particle_solid.size)

    def porous_resistance_fraction(self) -> float:
        """Return the fraction of lattice cells carrying particle-occupancy porous drag."""
        if not self.porous_resistance or self.n_p == 0:
            return 0.0
        mask = self._particle_occupancy_mask()
        return float(np.count_nonzero(mask) / mask.size)

    def dimensionless_groups(self, max_speed: float | None = None) -> dict[str, float]:
        """Return commonly tracked non-dimensional groups for fouling analysis."""
        if max_speed is None:
            _, ux, uy = self._macroscopic()
            max_speed = self._fluid_max_speed(ux, uy)
        particle_diameter = 2.0 * self.r_p
        nu = max(self.nu, 1e-12)
        l_ref = max(self.reynolds_length, 1e-12)
        u_ref = max(float(max_speed), 1e-12)
        return {
            "observed_reynolds_number": float(u_ref * l_ref / nu),
            "particle_reynolds_number": float(u_ref * particle_diameter / nu),
            "stokes_number_estimate": float(self.density_ratio * particle_diameter**2 * u_ref / (18.0 * nu * l_ref)),
            "particle_to_length_ratio": float(particle_diameter / l_ref),
            "brinkman_resistance_number": float(self.porous_resistance_coeff * l_ref**2 / nu) if self.porous_resistance else 0.0,
        }

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

    def _invalidate_macroscopic_cache(self) -> None:
        """Hook for subclasses that cache macroscopic fields."""
        return None

    def _fluid_max_speed(self, ux: np.ndarray, uy: np.ndarray) -> float:
        """Return the maximum fluid speed over non-solid nodes."""
        speed = np.sqrt(ux**2 + uy**2)
        if self.fluid_area <= 0.0:
            return 0.0
        return float(np.max(speed[~self.solid]))

    def _control_drive_force(self, ux: np.ndarray, uy: np.ndarray) -> None:
        """Adjust the drive force to keep the maximum fluid speed near ``u_max``."""
        if self.flow_control != "target_max_velocity" or self.u_max <= 0.0:
            return
        current_max = self._fluid_max_speed(ux, uy)
        if current_max <= 1e-12:
            return
        ratio = np.clip(self.u_max / current_max, 0.5, 2.0)
        gain = self.flow_control_gain
        self.F_drive *= (1.0 - gain) + gain * ratio
        min_drive = 0.05 * self.initial_F_drive
        if min_drive > 0.0:
            self.F_drive = max(self.F_drive, min_drive)

    def _trt_omega_minus(self) -> float:
        """Return the odd-mode TRT relaxation rate from the magic parameter."""
        tau_plus_offset = max(1.0 / self.omega - 0.5, 1e-12)
        tau_minus = 0.5 + self.trt_magic_parameter / tau_plus_offset
        return 1.0 / tau_minus

    def _collide_bgk_guo(self, rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> None:
        """Single-relaxation-time BGK collision with Guo forcing."""
        feq = _equilibrium(rho, ux, uy)
        S = _guo_forcing(ux, uy, self.Fx, self.Fy, self.omega)
        self.f += self.omega * (feq - self.f) + S

    def _collide_trt_guo(self, rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> None:
        """Two-relaxation-time collision using even/odd D2Q9 populations."""
        feq = _equilibrium(rho, ux, uy)
        S = _guo_forcing(ux, uy, self.Fx, self.Fy, self.omega)
        f_old = self.f.copy()
        f_new = np.empty_like(f_old)
        omega_plus = self.omega
        omega_minus = self.trt_omega_minus

        for i in range(Q):
            j = OPPOSITE[i]
            f_plus = 0.5 * (f_old[i] + f_old[j])
            f_minus = 0.5 * (f_old[i] - f_old[j])
            feq_plus = 0.5 * (feq[i] + feq[j])
            feq_minus = 0.5 * (feq[i] - feq[j])
            f_new[i] = (
                f_old[i]
                - omega_plus * (f_plus - feq_plus)
                - omega_minus * (f_minus - feq_minus)
                + S[i]
            )
        self.f = f_new

    def _lbm_step(self, rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> None:
        """One LBM step: collide → stream → bounce-back."""
        if self.uses_numba_lbm:
            method_id = 1 if self.fluid_method == "lbm-trt-guo" else 0
            self.f = _lbm_step_numba(
                self.f,
                rho,
                ux,
                uy,
                self.Fx,
                self.Fy,
                self.solid,
                self.omega,
                self.trt_omega_minus,
                method_id,
                self._c_int,
                W,
                self._opposite,
            )
            return

        if self.fluid_method == "lbm-trt-guo":
            self._collide_trt_guo(rho, ux, uy)
        else:
            self._collide_bgk_guo(rho, ux, uy)

        # Streaming (periodic source indices are precomputed at initialization)
        streamed = np.empty_like(self.f)
        for i in range(Q):
            streamed[i] = self.f[i][np.ix_(self._stream_x_src[i], self._stream_y_src[i])]
        self.f = streamed

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
        uf_x, uf_y = self._interp_velocity_many(np.array([x]), np.array([y]), ux=ux, uy=uy)
        return float(uf_x[0]), float(uf_y[0])

    def _interp_velocity_many(
        self,
        x: np.ndarray,
        y: np.ndarray,
        ux: np.ndarray | None = None,
        uy: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Vectorised bilinear interpolation of (ux, uy) for many particles."""
        if ux is None or uy is None:
            _, ux, uy = self._macroscopic()
        if len(x) == 0:
            return np.empty(0), np.empty(0)

        x_floor = np.floor(x)
        y_floor = np.floor(y)
        xi = x_floor.astype(int) % self.nx
        yi = np.clip(y_floor.astype(int), 0, self.ny - 2)
        xi1 = (xi + 1) % self.nx
        yi1 = np.minimum(yi + 1, self.ny - 1)
        tx = x - x_floor
        ty = y - y_floor

        w00 = (1.0 - tx) * (1.0 - ty)
        w10 = tx * (1.0 - ty)
        w01 = (1.0 - tx) * ty
        w11 = tx * ty
        uf_x = w00 * ux[xi, yi] + w10 * ux[xi1, yi] + w01 * ux[xi, yi1] + w11 * ux[xi1, yi1]
        uf_y = w00 * uy[xi, yi] + w10 * uy[xi1, yi] + w01 * uy[xi, yi1] + w11 * uy[xi1, yi1]
        return uf_x, uf_y

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

    def _stokes_drag_many(
        self,
        uf_x: np.ndarray,
        uf_y: np.ndarray,
        vel: np.ndarray,
        radii: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Vectorised Stokes drag for all active particles."""
        coeff = 3.0 * np.pi * self.nu * 2.0 * radii
        return coeff * (uf_x - vel[:, 0]), coeff * (uf_y - vel[:, 1])

    def _particle_drag_forces(
        self,
        ux: np.ndarray | None = None,
        uy: np.ndarray | None = None,
    ) -> np.ndarray:
        """Return Stokes drag forces for all active particles."""
        if self.n_p == 0:
            return np.zeros((0, 2))
        uf_x, uf_y = self._interp_velocity_many(self.pos[:, 0], self.pos[:, 1], ux=ux, uy=uy)
        fd_x, fd_y = self._stokes_drag_many(uf_x, uf_y, self.vel, self.radii)
        return np.column_stack((fd_x, fd_y))

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

    def _distribute_forces_many(
        self,
        x: np.ndarray,
        y: np.ndarray,
        fx: np.ndarray,
        fy: np.ndarray,
        radii: np.ndarray,
    ) -> None:
        """Vectorised four-node force distribution for particle back-reaction."""
        if len(x) == 0:
            return

        x_floor = np.floor(x)
        y_floor = np.floor(y)
        xi = x_floor.astype(int) % self.nx
        yi = np.clip(y_floor.astype(int), 0, self.ny - 2)
        xi1 = (xi + 1) % self.nx
        yi1 = np.minimum(yi + 1, self.ny - 1)
        tx = x - x_floor
        ty = y - y_floor
        weights = (
            (1.0 - tx) * (1.0 - ty),
            tx * (1.0 - ty),
            (1.0 - tx) * ty,
            tx * ty,
        )
        node_x = (xi, xi1, xi, xi1)
        node_y = (yi, yi, yi1, yi1)
        area = np.pi * radii**2
        scaled_fx = fx / area
        scaled_fy = fy / area
        for wk, ix, iy in zip(weights, node_x, node_y):
            active = ~self.solid[ix, iy]
            np.add.at(self.Fx, (ix[active], iy[active]), wk[active] * scaled_fx[active])
            np.add.at(self.Fy, (ix[active], iy[active]), wk[active] * scaled_fy[active])

    def _distribute_lagrangian_forces_many(
        self,
        x: np.ndarray,
        y: np.ndarray,
        fx: np.ndarray,
        fy: np.ndarray,
    ) -> None:
        """Spread IBM marker forces to the Eulerian lattice with bilinear weights."""
        if len(x) == 0:
            return

        x_floor = np.floor(x)
        y_floor = np.floor(y)
        xi = x_floor.astype(int) % self.nx
        yi = np.clip(y_floor.astype(int), 0, self.ny - 2)
        xi1 = (xi + 1) % self.nx
        yi1 = np.minimum(yi + 1, self.ny - 1)
        tx = x - x_floor
        ty = y - y_floor
        weights = (
            (1.0 - tx) * (1.0 - ty),
            tx * (1.0 - ty),
            (1.0 - tx) * ty,
            tx * ty,
        )
        node_x = (xi, xi1, xi, xi1)
        node_y = (yi, yi, yi1, yi1)
        for wk, ix, iy in zip(weights, node_x, node_y):
            active = ~self.solid[ix, iy]
            np.add.at(self.Fx, (ix[active], iy[active]), wk[active] * fx[active])
            np.add.at(self.Fy, (ix[active], iy[active]), wk[active] * fy[active])

    def _apply_immersed_boundary_forces(self, ux: np.ndarray, uy: np.ndarray) -> None:
        """Direct-forcing IBM for circular DEM particles."""
        self.ibm_forces_p = np.zeros((self.n_p, 2))
        self.ibm_torques_p = np.zeros(self.n_p)
        if self.n_p == 0 or self.particle_fluid_coupling != "immersed_boundary":
            return

        if self.uses_numba_compute and _ibm_markers_numba is not None:
            (
                marker_x,
                marker_y,
                marker_rx,
                marker_ry,
                marker_ds,
                marker_owner,
                ub_x,
                ub_y,
            ) = _ibm_markers_numba(
                self.pos,
                self.vel,
                self.omega_p,
                self.radii,
                self.ibm_marker_spacing,
                self.ny,
            )
        else:
            marker_x_parts: list[np.ndarray] = []
            marker_y_parts: list[np.ndarray] = []
            marker_rx_parts: list[np.ndarray] = []
            marker_ry_parts: list[np.ndarray] = []
            marker_ds_parts: list[np.ndarray] = []
            marker_owner_parts: list[np.ndarray] = []

            for i in range(self.n_p):
                radius = float(self.radii[i])
                n_markers = max(8, int(np.ceil(2.0 * np.pi * radius / self.ibm_marker_spacing)))
                ds = 2.0 * np.pi * radius / n_markers
                theta = 2.0 * np.pi * np.arange(n_markers) / n_markers
                rx = radius * np.cos(theta)
                ry = radius * np.sin(theta)
                x = self.pos[i, 0] + rx
                y = np.clip(self.pos[i, 1] + ry, 0.5, self.ny - 1.5)
                marker_x_parts.append(x)
                marker_y_parts.append(y)
                marker_rx_parts.append(rx)
                marker_ry_parts.append(ry)
                marker_ds_parts.append(np.full(n_markers, ds))
                marker_owner_parts.append(np.full(n_markers, i, dtype=np.int64))

            marker_x = np.concatenate(marker_x_parts)
            marker_y = np.concatenate(marker_y_parts)
            marker_rx = np.concatenate(marker_rx_parts)
            marker_ry = np.concatenate(marker_ry_parts)
            marker_ds = np.concatenate(marker_ds_parts)
            marker_owner = np.concatenate(marker_owner_parts)
            owner_vel = self.vel[marker_owner]
            owner_omega = self.omega_p[marker_owner]
            ub_x = owner_vel[:, 0] - owner_omega * marker_ry
            ub_y = owner_vel[:, 1] + owner_omega * marker_rx

        uf_x, uf_y = self._interp_velocity_many(marker_x, marker_y, ux=ux, uy=uy)
        marker_fx = self.ibm_stiffness * (ub_x - uf_x) * marker_ds
        marker_fy = self.ibm_stiffness * (ub_y - uf_y) * marker_ds

        if self.uses_numba_compute and _ibm_particle_reaction_numba is not None:
            self.ibm_forces_p, self.ibm_torques_p = _ibm_particle_reaction_numba(
                marker_owner,
                marker_rx,
                marker_ry,
                marker_fx,
                marker_fy,
                self.n_p,
            )
        else:
            np.add.at(self.ibm_forces_p[:, 0], marker_owner, -marker_fx)
            np.add.at(self.ibm_forces_p[:, 1], marker_owner, -marker_fy)
            marker_torque = marker_rx * marker_fy - marker_ry * marker_fx
            np.add.at(self.ibm_torques_p, marker_owner, -marker_torque)

        self._distribute_lagrangian_forces_many(
            marker_x,
            marker_y,
            marker_fx,
            marker_fy,
        )

    def ibm_marker_points(self) -> dict[str, np.ndarray]:
        """Return current IBM marker coordinates for visualization/debug output."""
        if self.n_p == 0 or self.particle_fluid_coupling != "immersed_boundary":
            empty_float = np.empty(0, dtype=float)
            return {
                "x": empty_float,
                "y": empty_float,
                "owner": np.empty(0, dtype=np.int64),
                "rx": empty_float,
                "ry": empty_float,
                "ds": empty_float,
            }

        if self.uses_numba_compute and _ibm_markers_numba is not None:
            marker_x, marker_y, marker_rx, marker_ry, marker_ds, marker_owner, _, _ = (
                _ibm_markers_numba(
                    self.pos,
                    self.vel,
                    self.omega_p,
                    self.radii,
                    self.ibm_marker_spacing,
                    self.ny,
                )
            )
            return {
                "x": marker_x,
                "y": marker_y,
                "owner": marker_owner,
                "rx": marker_rx,
                "ry": marker_ry,
                "ds": marker_ds,
            }

        marker_x_parts: list[np.ndarray] = []
        marker_y_parts: list[np.ndarray] = []
        marker_rx_parts: list[np.ndarray] = []
        marker_ry_parts: list[np.ndarray] = []
        marker_ds_parts: list[np.ndarray] = []
        marker_owner_parts: list[np.ndarray] = []
        for i in range(self.n_p):
            radius = float(self.radii[i])
            n_markers = max(8, int(np.ceil(2.0 * np.pi * radius / self.ibm_marker_spacing)))
            ds = 2.0 * np.pi * radius / n_markers
            theta = 2.0 * np.pi * np.arange(n_markers) / n_markers
            rx = radius * np.cos(theta)
            ry = radius * np.sin(theta)
            marker_x_parts.append(self.pos[i, 0] + rx)
            marker_y_parts.append(np.clip(self.pos[i, 1] + ry, 0.5, self.ny - 1.5))
            marker_rx_parts.append(rx)
            marker_ry_parts.append(ry)
            marker_ds_parts.append(np.full(n_markers, ds))
            marker_owner_parts.append(np.full(n_markers, i, dtype=np.int64))

        return {
            "x": np.concatenate(marker_x_parts),
            "y": np.concatenate(marker_y_parts),
            "owner": np.concatenate(marker_owner_parts),
            "rx": np.concatenate(marker_rx_parts),
            "ry": np.concatenate(marker_ry_parts),
            "ds": np.concatenate(marker_ds_parts),
        }

    # ------------------------------------------------------------------
    # DEM
    # ------------------------------------------------------------------

    @staticmethod
    def _tangent_from_normal(nx_: float, ny_: float) -> np.ndarray:
        """Return a unit tangent for a 2-D contact normal."""
        return np.array([-ny_, nx_])

    def _normal_contact_magnitude(self, overlap: float, v_n: float, mass: float) -> float:
        """Hertz normal force with damping, clamped to avoid artificial tension."""
        f_n = self.k_n * overlap**1.5
        f_damp = -self.damping * v_n * float(np.sqrt(self.k_n * mass))
        return max(f_n + f_damp, 0.0)

    def _tangential_force_magnitude(
        self,
        v_t: float,
        normal_force: float,
        mass: float,
    ) -> float:
        """Coulomb-limited tangential force opposing slip at a contact."""
        if not self.rolling_friction or normal_force <= 0.0:
            return 0.0
        f_trial = -self.tangential_damping * float(np.sqrt(self.k_n * mass)) * v_t
        f_limit = self.sliding_friction * normal_force
        return float(np.clip(f_trial, -f_limit, f_limit))

    def _rolling_resistance_torque(
        self,
        omega: float,
        normal_force: float,
        radius: float,
        mass: float,
    ) -> float:
        """Coulomb-limited rolling resistance torque opposing angular velocity."""
        if not self.rolling_friction or normal_force <= 0.0:
            return 0.0
        torque_trial = (
            -self.rolling_damping
            * float(np.sqrt(self.k_n * mass))
            * radius**2
            * omega
        )
        torque_limit = self.rolling_friction_coeff * normal_force * radius
        return float(np.clip(torque_trial, -torque_limit, torque_limit))

    def _particle_interaction_search_radius(self) -> float:
        """Maximum centre-to-centre distance that can produce pair forces."""
        interaction_cutoff = 0.0
        if self.particle_attraction:
            interaction_cutoff = max(interaction_cutoff, self.attraction_cutoff)
        if self.particle_repulsion:
            interaction_cutoff = max(interaction_cutoff, self.repulsion_cutoff)
        return 2.0 * float(np.max(self.radii)) + interaction_cutoff

    def _particle_cell_size(self) -> float:
        """Cell-list size for contact and near-surface particle interactions."""
        if self.n_p == 0:
            return 1.0
        return max(self._particle_interaction_search_radius(), 1.0)

    def _build_particle_cell_list(self, cell_size: float) -> dict[tuple[int, int], list[int]]:
        """Assign particle centres to spatial cells."""
        cells: dict[tuple[int, int], list[int]] = {}
        for i in range(self.n_p):
            key = (
                int(np.floor(self.pos[i, 0] / cell_size)),
                int(np.floor(self.pos[i, 1] / cell_size)),
            )
            cells.setdefault(key, []).append(i)
        return cells

    def _particle_pair_candidates(self):
        """Yield nearby particle pairs using a DEM-style cell list.

        The cell size is at least the largest centre-to-centre interaction
        range.  Therefore each particle only needs its own cell and the
        adjacent cells.  Only the forward half-neighbourhood is traversed so
        pairs are generated once without a ``seen`` set.
        """
        if self.n_p < 2:
            return

        if self.particle_search == "all_pairs":
            for i in range(self.n_p):
                for j in range(i + 1, self.n_p):
                    yield i, j
            return

        cell_size = self._particle_cell_size()
        cells = self._build_particle_cell_list(cell_size)
        neighbour_offsets = ((0, 0), (1, -1), (1, 0), (1, 1), (0, 1))

        for (cx, cy), ids in cells.items():
            for dx, dy in neighbour_offsets:
                other = cells.get((cx + dx, cy + dy))
                if other is None:
                    continue
                if dx == 0 and dy == 0:
                    for i_idx, i in enumerate(ids):
                        for j in ids[i_idx + 1 :]:
                            yield i, j
                else:
                    for i in ids:
                        for j in other:
                            yield i, j

    def _particle_pair_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        """Return nearby particle-pair indices as integer arrays."""
        pairs = list(self._particle_pair_candidates())
        if not pairs:
            return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
        pair_array = np.asarray(pairs, dtype=np.int64)
        return pair_array[:, 0], pair_array[:, 1]

    def _dem_loads(self, dt_sub: float) -> tuple[np.ndarray, np.ndarray]:
        """Compute total DEM force and torque through the shared DEM solver."""
        return self.dem_solver.compute_loads(dt_sub)

    def _dem_forces(self, dt_sub: float) -> np.ndarray:
        """Return total DEM forces, preserving the historical test/helper API."""
        forces, torques = self._dem_loads(dt_sub)
        self.torques_p = torques
        return forces

    def _dem_substep(self, dt: float) -> None:
        """One DEM sub-step via Velocity Verlet with time step ``dt``."""
        # Half-velocity update (per-particle mass)
        forces, torques = self._dem_loads(dt)
        acc = forces / self.masses[:, np.newaxis]
        alpha = torques / self.inertias
        self.vel += 0.5 * dt * acc
        self.omega_p += 0.5 * dt * alpha
        # Position update.  The historical mode uses x-periodicity; inlet mode
        # lets particles leave the right outlet and deletes them from the DEM set.
        self.pos += dt * self.vel
        if self.particle_source == "left_inlet":
            self._delete_right_outflow_particles()
            if self.n_p == 0:
                self.forces_p = np.zeros((0, 2))
                self.torques_p = np.zeros(0)
                self.ibm_forces_p = np.zeros((0, 2))
                self.ibm_torques_p = np.zeros(0)
                return
        else:
            self.pos[:, 0] %= self.nx

        # Second force evaluation
        forces_new, torques_new = self._dem_loads(dt)
        acc_new = forces_new / self.masses[:, np.newaxis]
        alpha_new = torques_new / self.inertias
        self.vel += 0.5 * dt * acc_new
        self.omega_p += 0.5 * dt * alpha_new
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
            # Clamp position outside fixed cylinder surfaces
            for cx, cy, cr in self.cylinders:
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
        self.torques_p = torques_new

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
            self._try_feed_left_inlet_particles()
            self._update_particle_solid_mask()

            _, ux_ctrl, uy_ctrl = self._macroscopic()
            self._control_drive_force(ux_ctrl, uy_ctrl)

            # --- Reset body force to driving force only ---
            self.Fx[:] = self.F_drive
            self.Fy[:] = 0.0
            self._invalidate_macroscopic_cache()

            if self.porous_resistance and self.porous_resistance_coeff > 0.0:
                _, ux_res, uy_res = self._macroscopic()
                self._apply_porous_resistance(ux_res, uy_res)

            # --- Particle → fluid back-reaction (per-particle radius) ---
            if self.n_p and self.particle_fluid_coupling == "point_force":
                _, ux_drag, uy_drag = self._macroscopic()
                drag = self._particle_drag_forces(ux=ux_drag, uy=uy_drag)
                # Newton 3rd law: -drag acts on the fluid
                self._distribute_forces_many(
                    self.pos[:, 0],
                    self.pos[:, 1],
                    -drag[:, 0],
                    -drag[:, 1],
                    self.radii,
                )
                self._invalidate_macroscopic_cache()
            elif self.n_p and self.particle_fluid_coupling == "immersed_boundary":
                _, ux_ibm, uy_ibm = self._macroscopic()
                self._apply_immersed_boundary_forces(ux_ibm, uy_ibm)
                self._invalidate_macroscopic_cache()

            # --- LBM step ---
            rho, ux, uy = self._macroscopic()
            self._lbm_step(rho, ux, uy)

            # --- DEM sub-steps ---
            for _ in range(self.dem_substeps):
                self._dem_substep(dt_sub)
            self._update_particle_solid_mask()

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
    _draw_cylinders(ax, sim)
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
    _draw_cylinders(ax, sim)
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
    _draw_cylinders(ax, sim)

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


def _draw_cylinders(ax: plt.Axes, sim: LBMDEMSolver) -> None:
    """Overlay fixed solid cylinders on an existing axes."""
    for cx, cy, cr in sim.cylinders:
        circle = plt.Circle(
            (cx, cy),
            cr,
            color="cyan",
            alpha=0.55,
            linewidth=1.0,
            zorder=4,
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
        "--fluid-method",
        choices=FLUID_METHODS,
        default="lbm-bgk-guo",
        help="LBM collision/forcing method (default lbm-bgk-guo)",
    )
    parser.add_argument(
        "--fluid-accelerator",
        choices=FLUID_ACCELERATORS,
        default="numpy",
        help="LBM execution backend (default numpy)",
    )
    parser.add_argument(
        "--compute-accelerator",
        choices=COMPUTE_ACCELERATORS,
        default="auto",
        help="Non-LBM compute backend (default auto)",
    )
    parser.add_argument(
        "--particle-method",
        choices=PARTICLE_METHODS,
        default="dem-hertz",
        help="DEM normal contact model (default dem-hertz)",
    )
    parser.add_argument(
        "--particle-search",
        choices=PARTICLE_SEARCH_METHODS,
        default="cell_list",
        help="Particle-pair neighbour search (default cell_list)",
    )
    parser.add_argument(
        "--particle-fluid-coupling",
        choices=PARTICLE_FLUID_COUPLINGS,
        default="point_force",
        help="Fluid-particle coupling mode (default point_force)",
    )
    parser.add_argument(
        "--ibm-stiffness",
        type=float,
        default=1.0,
        help="Direct-forcing immersed-boundary stiffness (default 1.0)",
    )
    parser.add_argument(
        "--ibm-marker-spacing",
        type=float,
        default=1.0,
        help="Approximate immersed-boundary marker spacing in lattice units (default 1.0)",
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
    parser.add_argument(
        "--rolling-friction",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Enable angular motion, tangential friction, and rolling resistance",
    )
    parser.add_argument(
        "--sliding-friction",
        type=float,
        default=0.5,
        help="Coulomb limit for tangential contact force (default 0.5)",
    )
    parser.add_argument(
        "--tangential-damping",
        type=float,
        default=0.4,
        help="Tangential slip damping scale (default 0.4)",
    )
    parser.add_argument(
        "--rolling-friction-coeff",
        type=float,
        default=0.05,
        help="Rolling-resistance moment coefficient (default 0.05)",
    )
    parser.add_argument(
        "--rolling-damping",
        type=float,
        default=0.2,
        help="Angular velocity damping scale for rolling resistance (default 0.2)",
    )
    parser.add_argument(
        "--particle-attraction",
        action="store_true",
        help="Enable Hamaker-like particle-particle attraction (default off)",
    )
    parser.add_argument(
        "--particle-repulsion",
        action="store_true",
        help="Enable Hamaker-like particle-particle repulsion (default off)",
    )
    parser.add_argument(
        "--attraction-strength",
        type=float,
        default=1e-3,
        help="Dimensionless Hamaker-like attraction strength (default 1e-3)",
    )
    parser.add_argument(
        "--repulsion-strength",
        type=float,
        default=1e-3,
        help="Dimensionless Hamaker-like repulsion strength (default 1e-3)",
    )
    parser.add_argument(
        "--attraction-cutoff",
        type=float,
        default=3.0,
        help="Surface gap cutoff for attraction in lattice nodes (default 3.0)",
    )
    parser.add_argument(
        "--repulsion-cutoff",
        type=float,
        default=3.0,
        help="Surface gap cutoff for repulsion in lattice nodes (default 3.0)",
    )
    parser.add_argument(
        "--attraction-min-gap",
        type=float,
        default=0.05,
        help="Minimum surface gap used to regularise attraction (default 0.05)",
    )
    parser.add_argument(
        "--repulsion-min-gap",
        type=float,
        default=0.05,
        help="Minimum surface gap used to regularise repulsion (default 0.05)",
    )
    parser.add_argument(
        "--cylinder",
        type=float,
        nargs=3,
        action="append",
        metavar=("X", "Y", "R"),
        default=None,
        help="Add a fixed cylinder. Can be repeated: --cylinder X Y R",
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
        default=None,
        help="Output directory for figures (default: src/results/lbm_dem/)",
    )
    args = parser.parse_args()
    if args.particle_attraction and args.particle_repulsion:
        parser.error("--particle-attraction and --particle-repulsion are mutually exclusive")

    out_dir = Path(args.out_dir) if args.out_dir else program_results_dir(__file__)
    out_dir.mkdir(parents=True, exist_ok=True)

    sim = LBMDEMSolver(
        nx=args.nx,
        ny=args.ny,
        Re=args.Re,
        u_max=args.u_max,
        fluid_method=args.fluid_method,
        fluid_accelerator=args.fluid_accelerator,
        compute_accelerator=args.compute_accelerator,
        particle_method=args.particle_method,
        particle_search=args.particle_search,
        particle_fluid_coupling=args.particle_fluid_coupling,
        ibm_stiffness=args.ibm_stiffness,
        ibm_marker_spacing=args.ibm_marker_spacing,
        n_particles=args.n_particles,
        particle_radius=args.radius,
        density_ratio=args.density_ratio,
        gravity=args.gravity,
        rolling_friction=args.rolling_friction,
        sliding_friction=args.sliding_friction,
        tangential_damping=args.tangential_damping,
        rolling_friction_coeff=args.rolling_friction_coeff,
        rolling_damping=args.rolling_damping,
        particle_attraction=args.particle_attraction,
        particle_repulsion=args.particle_repulsion,
        attraction_strength=args.attraction_strength,
        repulsion_strength=args.repulsion_strength,
        attraction_cutoff=args.attraction_cutoff,
        repulsion_cutoff=args.repulsion_cutoff,
        attraction_min_gap=args.attraction_min_gap,
        repulsion_min_gap=args.repulsion_min_gap,
        cylinders=args.cylinder,
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
    if args.particle_attraction:
        surface_force_mode = "attr"
    elif args.particle_repulsion:
        surface_force_mode = "rep"
    else:
        surface_force_mode = "noforce"
    rolling_mode = "rollfric" if args.rolling_friction else "freeroll"
    suffix = f"Re{int(args.Re)}_{args.nx}x{args.ny}_np{sim.n_p}_{surface_force_mode}_{rolling_mode}"
    save_fields = out_dir / f"fields_{suffix}.png"
    save_parts = out_dir / f"particles_{suffix}.png"

    plot_fields(sim, save_path=save_fields)
    plot_particles(sim, save_path=save_parts)
    print("Done.")


if __name__ == "__main__":
    main()
