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

import numpy as np

from .dem.kernels import _dem_pair_loads_numba
from .dem.solver import PARTICLE_METHODS, DEMSolver
from .geometry.pore import PoreGeometry
from .ibm.kernels import (
    _ibm_markers_numba,
    _ibm_particle_reaction_numba,
    _particle_solid_mask_numba,
)
from .io.paths import program_results_dir
from .io.visualization import plot_fields, plot_particles
from .lbm.constants import (
    COMPUTE_ACCELERATORS,
    CS2,
    FLUID_ACCELERATORS,
    FLUID_METHODS,
    LE_BOT_CROSS_DIRS,
    LE_TOP_CROSS_DIRS,
    OPPOSITE,
    PARTICLE_FLUID_COUPLINGS,
    PARTICLE_SEARCH_METHODS,
    STREAMWISE_BOUNDARIES,
    Y_BOUNDARIES,
    C,
    Q,
    W,
)
from .lbm.kernels import _lbm_step_numba
from .lbm.operators import equilibrium as _equilibrium
from .lbm.operators import guo_forcing as _guo_forcing

# D2Q9 constants and solver enumerations are defined in lbm/constants.py
# and re-imported above.


# Main coupled solver
# ---------------------------------------------------------------------------


class LBMDEMSolver:
    """
    2-D Coupled LBM (fluid) + DEM (particles) solver.

    Domain
    ------
    Rectangular channel of size ``nx × ny`` lattice nodes.
    The membrane-pore setup uses periodic transverse boundaries in y and
    pressure inlet/outlet boundaries in x.  A no-slip top/bottom wall mode is
    still available for channel-style checks.

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
        y_boundary: str = "periodic",
        streamwise_boundary: str = "pressure",
        pressure_drop: float = 1e-4,
        rho_out: float = 1.0,
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
        geometry: PoreGeometry | None = None,
        cylinder: tuple | None = None,
        cylinders: list[tuple[float, float, float]] | tuple[tuple[float, float, float], ...] | None = None,
        particle_source: str = "initial",
        source_volume_fraction: float | None = None,
        le_shear_rate: float = 0.0,
        le_shear_axis: int = 0,
        le_boundary_axis: int = 1,
        le_interpolation_order: int = 3,
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
        if y_boundary not in Y_BOUNDARIES:
            raise ValueError(f"y_boundary must be one of {Y_BOUNDARIES}")
        if streamwise_boundary not in STREAMWISE_BOUNDARIES:
            raise ValueError(f"streamwise_boundary must be one of {STREAMWISE_BOUNDARIES}")
        if rho_out <= 0.0:
            raise ValueError("rho_out must be positive")
        if pressure_drop < 0.0:
            raise ValueError("pressure_drop must be non-negative")
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
        self.y_boundary = y_boundary
        self.streamwise_boundary = streamwise_boundary
        self.rho_out = rho_out
        self.rho_in = rho_out + pressure_drop / CS2
        self.pressure_drop = pressure_drop
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
        if geometry is not None and (cylinder is not None or cylinders is not None):
            raise ValueError("geometry cannot be combined with legacy cylinder/cylinders inputs")
        self.geometry = geometry or PoreGeometry.from_legacy_inputs(
            cylinder=cylinder,
            cylinders=cylinders,
        )
        self.cylinders = self.geometry.as_tuples()
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
        # Lees-Edwards shear boundary state
        self.le_shear_rate = le_shear_rate
        self.le_shear_axis = le_shear_axis
        self.le_boundary_axis = le_boundary_axis
        self.le_interpolation_order = le_interpolation_order
        self._le_shift: float = 0.0  # accumulated fractional x-shift [lattice units]
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
        self.F_drive = 0.0 if streamwise_boundary == "pressure" else 8.0 * self.nu * u_max / ny**2
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
        self.mass_p = float(np.mean(self.masses)) if n_particles else 0.0  # representative scalar
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
        print(f"  Boundaries   : y={y_boundary}, x={streamwise_boundary}")
        if streamwise_boundary == "pressure":
            print(f"  Pressure BC  : rho_in={self.rho_in:.6g}, rho_out={self.rho_out:.6g}")
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

        # Solid nodes: optional top/bottom walls and fixed cylinders.  Moving particle
        # solids are overlaid in solid-boundary coupling mode.
        self.fixed_solid = np.zeros((nx, ny), dtype=bool)
        if self.y_boundary == "wall":
            self.fixed_solid[:, 0] = True
            self.fixed_solid[:, -1] = True

        # Solid nodes: fixed pore geometry.
        self.fixed_solid |= self.geometry.cylinder_solid_mask(
            nx,
            ny,
            y_boundary=self.y_boundary,
        )
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
        pair_kernel = (
            _dem_pair_loads_numba
            if particle_method == "dem-hertz" and self.y_boundary != "periodic"
            else None
        )
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
            if np.hypot(x - cx, self._periodic_y_delta(y - cy)) < cr + radius * clearance:
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
        if self.y_boundary == "periodic":
            return np.array([max(self.u_max, 0.0), 0.0])
        channel_height = max(self.ny - 2.0, 1.0)
        eta = float(np.clip((y - 0.5) / channel_height, 0.0, 1.0))
        return np.array([4.0 * self.u_max * eta * (1.0 - eta), 0.0])

    def _left_boundary_flow_rate(self, radius: float | None = None) -> float:
        """Return the positive volumetric flow rate through the left boundary."""
        _, ux, _ = self._macroscopic()
        if radius is None:
            y_min = 0 if self.y_boundary == "periodic" else 1
            y_max = self.ny - 1 if self.y_boundary == "periodic" else self.ny - 2
        else:
            if self.y_boundary == "periodic":
                y_min = 0
                y_max = self.ny - 1
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
        if self.y_boundary == "periodic":
            y_min = 0
            y_max = self.ny - 1
        else:
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
        if self.y_boundary == "periodic":
            return float(y % self.ny)
        return float(np.clip(y, radius + 0.5, self.ny - 1.5 - radius))

    def _can_place_inlet_particle(
        self,
        x: float,
        y: float,
        radius: float,
        n_existing: int | None = None,
    ) -> bool:
        """Return whether an inlet particle can be placed without overlap."""
        if self.y_boundary == "wall" and (y < radius + 0.5 or y > self.ny - 1.5 - radius):
            return False
        if self._overlaps_cylinder(x, y, radius, 1.05):
            return False
        n_check = self.n_p if n_existing is None else n_existing
        if n_check > 0:
            dy = y - self.pos[:n_check, 1]
            if self.y_boundary == "periodic":
                dy = dy - self.ny * np.rint(dy / self.ny)
            dists = np.hypot(x - self.pos[:n_check, 0], dy)
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
        if self.uses_numba_compute and _particle_solid_mask_numba is not None and self.y_boundary not in ("periodic", "lees_edwards"):
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
                if self.y_boundary == "periodic":
                    y_min = int(np.floor(y_pos - radius))
                    y_max = int(np.ceil(y_pos + radius))
                else:
                    y_min = max(0, int(np.floor(y_pos - radius)))
                    y_max = min(self.ny - 1, int(np.ceil(y_pos + radius)))
                if y_max < y_min:
                    continue
                x_indices = np.arange(x_min, x_max + 1)
                y_indices = np.arange(y_min, y_max + 1)
                xx = x_indices[:, np.newaxis] % self.nx
                yy = y_indices[np.newaxis, :] % self.ny
                dy = y_indices[np.newaxis, :] - y_pos
                if self.y_boundary == "periodic":
                    dy = dy - self.ny * np.rint(dy / self.ny)
                local = (x_indices[:, np.newaxis] - x_pos) ** 2 + dy**2 <= radius**2
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
        if (
            self.streamwise_boundary == "pressure"
            or self.flow_control != "target_max_velocity"
            or self.u_max <= 0.0
        ):
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

    def _apply_pressure_boundaries(self) -> None:
        """Apply Zou-He pressure inlet/outlet on the left/right boundaries."""
        if self.streamwise_boundary != "pressure":
            return
        rho_l = self.rho_in
        rho_r = self.rho_out

        y_indices = np.arange(self.ny)
        if self.y_boundary == "wall":
            y_indices = y_indices[1:-1]

        x = 0
        f0 = self.f[0, x, y_indices]
        f2 = self.f[2, x, y_indices]
        f4 = self.f[4, x, y_indices]
        f3 = self.f[3, x, y_indices]
        f6 = self.f[6, x, y_indices]
        f7 = self.f[7, x, y_indices]
        ux_l = 1.0 - (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / rho_l
        self.f[1, x, y_indices] = f3 + (2.0 / 3.0) * rho_l * ux_l
        self.f[5, x, y_indices] = f7 + 0.5 * (f4 - f2) + (1.0 / 6.0) * rho_l * ux_l
        self.f[8, x, y_indices] = f6 + 0.5 * (f2 - f4) + (1.0 / 6.0) * rho_l * ux_l

        x = self.nx - 1
        f0 = self.f[0, x, y_indices]
        f2 = self.f[2, x, y_indices]
        f4 = self.f[4, x, y_indices]
        f1 = self.f[1, x, y_indices]
        f5 = self.f[5, x, y_indices]
        f8 = self.f[8, x, y_indices]
        ux_r = -1.0 + (f0 + f2 + f4 + 2.0 * (f1 + f5 + f8)) / rho_r
        self.f[3, x, y_indices] = f1 - (2.0 / 3.0) * rho_r * ux_r
        self.f[7, x, y_indices] = f5 + 0.5 * (f2 - f4) - (1.0 / 6.0) * rho_r * ux_r
        self.f[6, x, y_indices] = f8 + 0.5 * (f4 - f2) - (1.0 / 6.0) * rho_r * ux_r
        self._invalidate_macroscopic_cache()

    def _fractional_roll_x(self, row: np.ndarray, shift: float) -> np.ndarray:
        """Shift a 1-D x-row by a fractional amount using linear or cubic interpolation."""
        nx = row.shape[0]
        int_shift = int(np.floor(shift)) % nx
        frac = shift - np.floor(shift)
        rolled = np.roll(row, int_shift)
        if abs(frac) < 1e-12:
            return rolled
        if self.le_interpolation_order == 1:
            return (1.0 - frac) * rolled + frac * np.roll(rolled, -1)
        # Third-order cubic interpolation with periodic boundary via scipy
        try:
            from scipy.ndimage import map_coordinates
            x_dst = (np.arange(nx, dtype=float) - frac) % nx
            return map_coordinates(rolled, [x_dst], order=3, mode="wrap")
        except ImportError:
            # Fall back to linear if scipy unavailable
            return (1.0 - frac) * rolled + frac * np.roll(rolled, -1)

    def _apply_le_particle_wrap(self) -> None:
        """Wrap particles that cross LE y boundaries, shifting x and vx accordingly."""
        ny = self.ny
        nx = self.nx
        shift = self._le_shift
        dv = self.le_shear_rate * ny  # velocity jump across the full domain height

        crossed_top = self.pos[:, 1] >= ny
        crossed_bot = self.pos[:, 1] < 0.0

        if np.any(crossed_top):
            self.pos[crossed_top, 1] -= ny
            self.pos[crossed_top, 0] = (self.pos[crossed_top, 0] - shift) % nx
            self.vel[crossed_top, 0] -= dv

        if np.any(crossed_bot):
            self.pos[crossed_bot, 1] += ny
            self.pos[crossed_bot, 0] = (self.pos[crossed_bot, 0] + shift) % nx
            self.vel[crossed_bot, 0] += dv

    def _apply_le_streaming_correction(self) -> None:
        """Apply Lees-Edwards x-shift to populations that just crossed the y boundary."""
        shift = self._le_shift
        # Directions crossing top boundary (cy > 0) land in y=0 after streaming
        for i in LE_TOP_CROSS_DIRS:
            self.f[i, :, 0] = self._fractional_roll_x(self.f[i, :, 0], shift)
        # Directions crossing bottom boundary (cy < 0) land in y=ny-1 after streaming
        for i in LE_BOT_CROSS_DIRS:
            self.f[i, :, self.ny - 1] = self._fractional_roll_x(self.f[i, :, self.ny - 1], -shift)

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
            self._apply_pressure_boundaries()
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

        if self.y_boundary == "lees_edwards":
            self._apply_le_streaming_correction()

        # Bounce-back on solid nodes
        f_tmp = self.f.copy()
        for i in range(Q):
            self.f[i, self.solid] = f_tmp[OPPOSITE[i], self.solid]
        self._apply_pressure_boundaries()

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
        y_eval = np.mod(y, self.ny) if self.y_boundary in ("periodic", "lees_edwards") else y
        y_floor = np.floor(y_eval)
        xi = x_floor.astype(int) % self.nx
        if self.y_boundary in ("periodic", "lees_edwards"):
            yi = y_floor.astype(int) % self.ny
            yi1 = (yi + 1) % self.ny
        else:
            yi = np.clip(y_floor.astype(int), 0, self.ny - 2)
            yi1 = np.minimum(yi + 1, self.ny - 1)
        xi1 = (xi + 1) % self.nx
        tx = x - x_floor
        ty = y_eval - y_floor

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

    def _interp_velocity_cubic(
        self,
        x: np.ndarray,
        y: np.ndarray,
        ux: np.ndarray,
        uy: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Cubic (third-order) periodic velocity interpolation for iSP coupling."""
        try:
            from scipy.ndimage import map_coordinates
        except ImportError:
            return self._interp_velocity_many(x, y, ux=ux, uy=uy)
        y_eval = np.mod(y, self.ny)
        coords = np.vstack([x % self.nx, y_eval])
        mode = "wrap" if self.y_boundary in ("periodic", "lees_edwards") else "nearest"
        uf_x = map_coordinates(ux, coords, order=3, mode=mode)
        uf_y = map_coordinates(uy, coords, order=3, mode=mode)
        return uf_x, uf_y

    def _particle_drag_forces(
        self,
        ux: np.ndarray | None = None,
        uy: np.ndarray | None = None,
    ) -> np.ndarray:
        """Return Stokes drag forces for all active particles."""
        if self.n_p == 0:
            return np.zeros((0, 2))
        if self.particle_fluid_coupling == "isp":
            if ux is None or uy is None:
                _, ux, uy = self._macroscopic()
            uf_x, uf_y = self._interp_velocity_cubic(self.pos[:, 0], self.pos[:, 1], ux, uy)
        else:
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
        y_eval = y % self.ny if self.y_boundary in ("periodic", "lees_edwards") else y
        y_floor = np.floor(y_eval)
        if self.y_boundary in ("periodic", "lees_edwards"):
            yi = int(y_floor) % self.ny
            yi1 = (yi + 1) % self.ny
        else:
            yi = int(np.clip(y_floor, 0, self.ny - 2))
            yi1 = min(yi + 1, self.ny - 1)
        xi1 = (xi + 1) % self.nx
        tx = x - np.floor(x)
        ty = y_eval - y_floor
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
        y_eval = np.mod(y, self.ny) if self.y_boundary in ("periodic", "lees_edwards") else y
        y_floor = np.floor(y_eval)
        xi = x_floor.astype(int) % self.nx
        if self.y_boundary in ("periodic", "lees_edwards"):
            yi = y_floor.astype(int) % self.ny
            yi1 = (yi + 1) % self.ny
        else:
            yi = np.clip(y_floor.astype(int), 0, self.ny - 2)
            yi1 = np.minimum(yi + 1, self.ny - 1)
        xi1 = (xi + 1) % self.nx
        tx = x - x_floor
        ty = y_eval - y_floor
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
        y_eval = np.mod(y, self.ny) if self.y_boundary in ("periodic", "lees_edwards") else y
        y_floor = np.floor(y_eval)
        xi = x_floor.astype(int) % self.nx
        if self.y_boundary in ("periodic", "lees_edwards"):
            yi = y_floor.astype(int) % self.ny
            yi1 = (yi + 1) % self.ny
        else:
            yi = np.clip(y_floor.astype(int), 0, self.ny - 2)
            yi1 = np.minimum(yi + 1, self.ny - 1)
        xi1 = (xi + 1) % self.nx
        tx = x - x_floor
        ty = y_eval - y_floor
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
                self.y_boundary in ("periodic", "lees_edwards"),
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
                if self.y_boundary in ("periodic", "lees_edwards"):
                    y = np.mod(self.pos[i, 1] + ry, self.ny)
                else:
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
                    self.y_boundary in ("periodic", "lees_edwards"),
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
            if self.y_boundary in ("periodic", "lees_edwards"):
                marker_y_parts.append(np.mod(self.pos[i, 1] + ry, self.ny))
            else:
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

    def _periodic_y_delta(self, dy: float) -> float:
        """Return the minimum-image y displacement when y is periodic or lees_edwards."""
        if self.y_boundary not in ("periodic", "lees_edwards"):
            return dy
        return dy - self.ny * float(np.rint(dy / self.ny))

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
        ncy = max(1, int(np.ceil(self.ny / cell_size)))
        for i in range(self.n_p):
            cy = int(np.floor(self.pos[i, 1] / cell_size))
            if self.y_boundary == "periodic":
                cy %= ncy
            key = (
                int(np.floor(self.pos[i, 0] / cell_size)),
                cy,
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

        if self.particle_search == "all_pairs" or self.y_boundary in ("periodic", "lees_edwards"):
            for i in range(self.n_p):
                for j in range(i + 1, self.n_p):
                    yield i, j
            return

        cell_size = self._particle_cell_size()
        cells = self._build_particle_cell_list(cell_size)
        ncy = max(1, int(np.ceil(self.ny / cell_size)))
        neighbour_offsets = ((0, 0), (1, -1), (1, 0), (1, 1), (0, 1))

        for (cx, cy), ids in cells.items():
            for dx, dy in neighbour_offsets:
                other_cy = cy + dy
                if self.y_boundary in ("periodic", "lees_edwards"):
                    other_cy %= ncy
                other = cells.get((cx + dx, other_cy))
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
        if self.y_boundary == "periodic" and self.n_p:
            self.pos[:, 1] %= self.ny
        elif self.y_boundary == "lees_edwards" and self.n_p:
            self._apply_le_particle_wrap()
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
            if self.y_boundary == "wall":
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
                dy = self._periodic_y_delta(self.pos[i, 1] - cy)
                dist = float(np.hypot(dx, dy))
                min_dist = cr + self.radii[i]
                if dist < min_dist and dist > 1e-10:
                    nx_ = dx / dist
                    ny_ = dy / dist
                    self.pos[i, 0] = cx + min_dist * nx_
                    self.pos[i, 1] = cy + min_dist * ny_
                    if self.y_boundary in ("periodic", "lees_edwards"):
                        self.pos[i, 1] %= self.ny
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

            self.Fx[:] = self.F_drive
            self.Fy[:] = 0.0
            self._invalidate_macroscopic_cache()

            if self.porous_resistance and self.porous_resistance_coeff > 0.0:
                _, ux_res, uy_res = self._macroscopic()
                self._apply_porous_resistance(ux_res, uy_res)

            if self.n_p and self.particle_fluid_coupling in ("point_force", "isp"):
                _, ux_drag, uy_drag = self._macroscopic()
                drag = self._particle_drag_forces(ux=ux_drag, uy=uy_drag)
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

            if self.y_boundary == "lees_edwards":
                self._le_shift = (self._le_shift + self.le_shear_rate * self.ny) % self.nx

            rho, ux, uy = self._macroscopic()
            self._lbm_step(rho, ux, uy)

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
    parser.add_argument("--y-boundary", choices=Y_BOUNDARIES, default="periodic")
    parser.add_argument(
        "--streamwise-boundary",
        choices=("periodic-force", "pressure"),
        default="pressure",
    )
    parser.add_argument("--pressure-drop", type=float, default=1e-4)
    parser.add_argument("--rho-out", type=float, default=1.0)
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
        y_boundary=args.y_boundary,
        streamwise_boundary=args.streamwise_boundary.replace("-", "_"),
        pressure_drop=args.pressure_drop,
        rho_out=args.rho_out,
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
