"""3D LBM solver on the D3Q15 lattice.

Implements a BGK-collision fluid solver with three boundary modes: Lees-Edwards
shear (periodic), Zou-He pressure inlet/outlet (issue #14), and immersed-boundary
particle coupling driven by an internal ``DEM3D`` contact solver (issues #17/#18).
"""

from __future__ import annotations

import math

import numpy as np

# ---------------------------------------------------------------------------
# D3Q15 lattice constants
# ---------------------------------------------------------------------------

Q3 = 15
CS2_3 = 1.0 / 3.0  # speed of sound squared

# Velocity vectors shape (15, 3): [cx, cy, cz]
C3 = np.array(
    [
        [0, 0, 0],  # 0  rest
        [1, 0, 0],  # 1  +x
        [-1, 0, 0],  # 2  -x
        [0, 1, 0],  # 3  +y
        [0, -1, 0],  # 4  -y
        [0, 0, 1],  # 5  +z
        [0, 0, -1],  # 6  -z
        [1, 1, 1],  # 7
        [-1, 1, 1],  # 8
        [1, -1, 1],  # 9
        [-1, -1, 1],  # 10
        [1, 1, -1],  # 11
        [-1, 1, -1],  # 12
        [1, -1, -1],  # 13
        [-1, -1, -1],  # 14
    ],
    dtype=float,
)

W3 = np.array(
    [
        2.0 / 9.0,  # 0  rest
        1.0 / 9.0,  # 1
        1.0 / 9.0,  # 2
        1.0 / 9.0,  # 3
        1.0 / 9.0,  # 4
        1.0 / 9.0,  # 5
        1.0 / 9.0,  # 6
        1.0 / 72.0,  # 7
        1.0 / 72.0,  # 8
        1.0 / 72.0,  # 9
        1.0 / 72.0,  # 10
        1.0 / 72.0,  # 11
        1.0 / 72.0,  # 12
        1.0 / 72.0,  # 13
        1.0 / 72.0,  # 14
    ]
)

# Opposite direction indices: C3[OPPOSITE3[i]] = -C3[i]
OPPOSITE3 = [0, 2, 1, 4, 3, 6, 5, 14, 13, 12, 11, 10, 9, 8, 7]

# Directions with cy > 0 (cross top boundary y=ny-1 → y=0 after streaming)
_LE_TOP_DIRS = [3, 7, 8, 11, 12]
# Directions with cy < 0 (cross bottom boundary y=0 → y=ny-1 after streaming)
_LE_BOT_DIRS = [4, 9, 10, 13, 14]

# Direction groups along x (streamwise) for pressure inlet/outlet BC.
# At the inlet plane x=0 the unknown (incoming) populations are those with cx > 0;
# at the outlet plane x=nx-1 the unknowns are those with cx < 0.
_X_POS_DIRS = [i for i in range(Q3) if C3[i, 0] > 0]  # cx > 0
_X_NEG_DIRS = [i for i in range(Q3) if C3[i, 0] < 0]  # cx < 0
_X_ZERO_DIRS = [i for i in range(Q3) if C3[i, 0] == 0]  # cx == 0


# ---------------------------------------------------------------------------
# Vectorised equilibrium
# ---------------------------------------------------------------------------


def _equilibrium_3d(
    rho: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
) -> np.ndarray:
    """Compute D3Q15 equilibrium distribution.

    Args:
        rho: Density field, shape (nx, ny, nz).
        ux, uy, uz: Velocity components, shape (nx, ny, nz).

    Returns:
        Equilibrium array of shape (15, nx, ny, nz).
    """
    u2 = ux**2 + uy**2 + uz**2  # (nx, ny, nz)
    feq = np.empty((Q3,) + rho.shape)
    for i in range(Q3):
        cx, cy, cz = C3[i]
        cu = cx * ux + cy * uy + cz * uz
        feq[i] = W3[i] * rho * (1.0 + cu / CS2_3 + 0.5 * cu**2 / CS2_3**2 - 0.5 * u2 / CS2_3)
    return feq


# ---------------------------------------------------------------------------
# 3D solver
# ---------------------------------------------------------------------------


class LBMDEMSolver3D:
    """3D Lattice-Boltzmann solver on the D3Q15 lattice.

    Implements BGK collision with two selectable streamwise (x) boundary modes:

    * ``streamwise_boundary="periodic"`` (default) — fully periodic streaming
      with the Lees-Edwards shear boundary correction in y.  This is the
      original wall-less shear-flow mode.
    * ``streamwise_boundary="pressure"`` — Zou-He density (pressure) inlet at
      x=0 and outlet at x=nx-1, with y and z periodic.  Drives a
      pressure-gradient flow; ``le_shear_rate`` must be 0 in this mode.

    Particles are coupled to the fluid by the immersed-boundary method (issue
    #17): set ``n_particles > 0`` with ``particle_fluid_coupling=
    "immersed_boundary"`` and supply ``particle_positions``.  An internal
    :class:`particulate_flow.dem.contact3d.DEM3D` (issue #18) provides the
    contact dynamics, driven each step by the IBM reaction force.  Fixed
    z-aligned cylinder obstacles (#19) are supported via ``cylinders`` (no-slip
    bounce-back for the fluid, collision bodies for the particles).  A
    pressure-driven ``particle_source="left_inlet"`` (#20) injects particles at
    the inlet over time per a volume-flux budget and removes them past the outlet.

    Args:
        nx: Grid size in x (streamwise direction).
        ny: Grid size in y (boundary-normal / gradient direction).
        nz: Grid size in z (vorticity direction).
        Re: Reynolds number (used to set kinematic viscosity via nu = u_ref*L/Re).
        u_max: Reference velocity scale for Re.  Ignored when ``le_shear_rate``
               drives the flow (set to 0.0 for pure shear).
        le_shear_rate: Lees-Edwards shear rate γ̇ in lattice units.  The
                       velocity jump across the y domain is γ̇ * ny.
        le_shear_axis: Index of the velocity component driven by shear (default 0 = x).
        le_boundary_axis: Index of the gradient direction (default 1 = y).
        le_interpolation_order: Polynomial order for the fractional-roll interpolation
                                used in the LE boundary correction (1 = linear, 3 = cubic).
        tau: Relaxation time.  When provided, overrides Re / u_max / ny derivation.
        init_analytical: When ``True`` initialise ``f`` from the analytical linear
                         shear velocity profile ux = γ̇ · (y − ny/2).  This reduces
                         spin-up time for validation runs.
        streamwise_boundary: ``"periodic"`` or ``"pressure"`` (see above).
        pressure_drop: Density-equivalent pressure drop driving the flow in
                       ``"pressure"`` mode; sets ``rho_in = rho_out + pressure_drop / cs²``.
        rho_out: Outlet density (must be positive).
        n_particles: Number of DEM particles (0 disables the particle path).
        cylinders: Fixed z-aligned finite cylinder obstacles (#19) as
                   ``(cx, cy, r)`` or ``(cx, cy, r, z_lo, z_hi)``.  They form a
                   no-slip solid mask for the fluid and collision bodies for the
                   particles.
        particle_source: ``"none"`` for a fixed population, or ``"left_inlet"``
                         to inject particles at the inlet over time (requires
                         ``streamwise_boundary="pressure"`` and
                         ``particle_fluid_coupling="immersed_boundary"``).
        particle_fluid_coupling: ``"none"`` or ``"immersed_boundary"``.  Particles
                                 require ``"immersed_boundary"``.
        particle_positions: ``(n_particles, 3)`` initial particle centres; required
                            when ``n_particles > 0`` (ignored for ``left_inlet``,
                            which starts empty).
        particle_radius: Radius applied to every particle (lattice units).
        density_ratio: ρ_particle / ρ_fluid for particle mass and buoyancy.
        gravity: Gravitational acceleration magnitude for the particles.
        ibm_stiffness: Direct-forcing IBM stiffness (force per unit velocity error).
        ibm_marker_spacing: Target arc spacing between surface markers.
        dem_substeps: DEM sub-steps per LBM step.
        sliding_friction: Coulomb limit μ for tangential (sliding) contact force.
        rolling_friction_coeff: Rolling-resistance moment coefficient.
        particle_attraction: Enable Hamaker-like short-range particle-particle and
            particle-cylinder attraction.  Mutually exclusive with
            ``particle_repulsion``.
        particle_repulsion: Enable Hamaker-like short-range repulsion.  Mutually
            exclusive with ``particle_attraction``.
        attraction_strength: Dimensionless Hamaker prefactor A* for attraction.
            Force magnitude: ``f = A* * r_eff / (6 * h²)``.
        repulsion_strength: Dimensionless Hamaker prefactor A* for repulsion.
        attraction_cutoff: Surface-gap distance (lattice units) beyond which
            attraction is zero.
        repulsion_cutoff: Surface-gap distance (lattice units) beyond which
            repulsion is zero.
        attraction_min_gap: Minimum surface gap used in the attraction formula
            (prevents divergence at contact).
        repulsion_min_gap: Minimum surface gap used in the repulsion formula.
        source_volume_fraction: Target injected solid-volume fraction of the inlet
                                flux for ``left_inlet`` (budget ``+= phi · inlet_flux``).
        seed: RNG seed for the inlet (y, z) sampling.
    """

    def __init__(
        self,
        nx: int,
        ny: int,
        nz: int,
        Re: float = 10.0,
        u_max: float = 0.05,
        le_shear_rate: float = 0.0,
        le_shear_axis: int = 0,
        le_boundary_axis: int = 1,
        le_interpolation_order: int = 3,
        tau: float | None = None,
        init_analytical: bool = False,
        streamwise_boundary: str = "periodic",
        pressure_drop: float = 0.0,
        rho_out: float = 1.0,
        n_particles: int = 0,
        cylinders: list | tuple | None = None,
        particle_source: str = "none",
        particle_fluid_coupling: str = "none",
        particle_positions: np.ndarray | None = None,
        particle_radius: float = 2.0,
        density_ratio: float = 2.0,
        gravity: float = 0.0,
        ibm_stiffness: float = 1.0,
        ibm_marker_spacing: float = 1.0,
        dem_substeps: int = 4,
        sliding_friction: float = 0.5,
        rolling_friction_coeff: float = 0.05,
        particle_attraction: bool = False,
        particle_repulsion: bool = False,
        attraction_strength: float = 1e-3,
        repulsion_strength: float = 1e-3,
        attraction_cutoff: float = 3.0,
        repulsion_cutoff: float = 3.0,
        attraction_min_gap: float = 0.05,
        repulsion_min_gap: float = 0.05,
        source_volume_fraction: float | None = None,
        seed: int = 42,
    ) -> None:
        if particle_attraction and particle_repulsion:
            raise ValueError("particle_attraction and particle_repulsion are mutually exclusive")
        if particle_source not in ("none", None, "left_inlet"):
            raise ValueError(
                f"particle_source must be 'none' or 'left_inlet', got {particle_source!r}"
            )
        if particle_fluid_coupling not in ("none", None, "immersed_boundary"):
            raise NotImplementedError(
                f"3D particle-fluid coupling {particle_fluid_coupling!r} is not "
                "supported; only 'immersed_boundary' (issue #17) or 'none'."
            )
        if n_particles > 0 and particle_fluid_coupling != "immersed_boundary":
            raise NotImplementedError(
                "3D particles require particle_fluid_coupling='immersed_boundary' " "(issue #17)."
            )
        if particle_source == "left_inlet":
            if streamwise_boundary != "pressure":
                raise ValueError(
                    "particle_source='left_inlet' requires "
                    "streamwise_boundary='pressure' (inlet injection needs a "
                    "pressure-driven inflow)."
                )
            if particle_fluid_coupling != "immersed_boundary":
                raise ValueError(
                    "particle_source='left_inlet' requires "
                    "particle_fluid_coupling='immersed_boundary'."
                )
            if source_volume_fraction is None or source_volume_fraction <= 0.0:
                raise ValueError(
                    "particle_source='left_inlet' requires a positive "
                    "source_volume_fraction (injection budget = phi * inlet_flux)."
                )

        if streamwise_boundary not in ("periodic", "pressure"):
            raise ValueError(
                "streamwise_boundary must be 'periodic' or 'pressure', "
                f"got {streamwise_boundary!r}"
            )
        if streamwise_boundary == "pressure" and le_shear_rate != 0.0:
            raise ValueError(
                "pressure streamwise_boundary and Lees-Edwards shear are mutually "
                "exclusive (pressure-driven flow requires le_shear_rate=0)"
            )
        if rho_out <= 0.0:
            raise ValueError("rho_out must be positive")

        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.le_shear_rate = le_shear_rate
        self.le_shear_axis = le_shear_axis
        self.le_boundary_axis = le_boundary_axis
        self.le_interpolation_order = le_interpolation_order
        self.step_count = 0
        self._le_shift: float = 0.0

        # Pressure-driven boundary configuration (Zou-He on x inlet/outlet).
        self.streamwise_boundary = streamwise_boundary
        self.rho_out = float(rho_out)
        self.rho_in = float(rho_out) + float(pressure_drop) / CS2_3

        # Derive tau from Re if not given explicitly.
        if tau is not None:
            self.tau = tau
        else:
            l_ref = float(ny)
            u_ref = float(u_max) if u_max > 0.0 else le_shear_rate * ny
            nu = u_ref * l_ref / max(Re, 1e-12)
            self.tau = nu / CS2_3 + 0.5
        self.omega = 1.0 / self.tau

        # Initialise distribution function.
        rho0 = np.ones((nx, ny, nz))
        if init_analytical:
            y = np.arange(ny, dtype=float)
            ux0 = le_shear_rate * (y - (ny - 1) / 2.0)
            ux_field = np.broadcast_to(ux0[np.newaxis, :, np.newaxis], (nx, ny, nz)).copy()
        else:
            ux_field = np.zeros((nx, ny, nz))
        self.f = _equilibrium_3d(rho0, ux_field, np.zeros_like(ux_field), np.zeros_like(ux_field))

        # Eulerian body-force field (Guo forcing), populated by IBM each step.
        self.Fx = np.zeros((nx, ny, nz))
        self.Fy = np.zeros((nx, ny, nz))
        self.Fz = np.zeros((nx, ny, nz))

        # --- Fixed finite z-aligned cylinder obstacles (issue #19) ---------
        self.cylinders = [tuple(c) for c in (cylinders or [])]
        self.solid = self._build_cylinder_solid(self.cylinders)

        # --- IBM particle coupling (issue #17) -----------------------------
        self.particle_fluid_coupling = particle_fluid_coupling
        self.ibm_stiffness = float(ibm_stiffness)
        self.ibm_marker_spacing = float(ibm_marker_spacing)

        # --- Left-inlet particle source (issue #20) ------------------------
        self.particle_source = particle_source if particle_source else "none"
        self.particle_radius = float(particle_radius)
        self.source_volume_fraction = source_volume_fraction
        self._inlet_volume_budget = 0.0
        self._inlet_cursor = 0
        self._rng_inlet = np.random.default_rng(seed)
        self.injected_particle_volume = 0.0
        self.cumulative_inlet_flow_volume = 0.0
        self.removed_particles = 0
        self.generated_particles = 0

        self.dem = None
        if n_particles > 0 or self.particle_source == "left_inlet":
            from .dem.contact3d import DEM3D

            if self.particle_source == "left_inlet":
                pos = np.empty((0, 3))  # start empty; particles injected over time
            else:
                if particle_positions is None:
                    raise ValueError("particle_positions (n,3) is required when n_particles > 0")
                pos = np.asarray(particle_positions, dtype=float).reshape((-1, 3))
                if pos.shape[0] != n_particles:
                    raise ValueError(
                        f"particle_positions has {pos.shape[0]} rows but "
                        f"n_particles={n_particles}"
                    )
            radii = np.full(pos.shape[0], float(particle_radius))
            self.dem = DEM3D(
                pos=pos,
                vel=np.zeros((pos.shape[0], 3)),
                radii=radii,
                nx=nx,
                ny=ny,
                nz=nz,
                density_ratio=density_ratio,
                gravity=gravity,
                dem_substeps=dem_substeps,
                cylinders=self.cylinders,
                sliding_friction=sliding_friction,
                rolling_friction_coeff=rolling_friction_coeff,
                particle_attraction=particle_attraction,
                particle_repulsion=particle_repulsion,
                attraction_strength=attraction_strength,
                repulsion_strength=repulsion_strength,
                attraction_cutoff=attraction_cutoff,
                repulsion_cutoff=repulsion_cutoff,
                attraction_min_gap=attraction_min_gap,
                repulsion_min_gap=repulsion_min_gap,
            )

    def _build_cylinder_solid(self, cylinders: list) -> np.ndarray:
        """Return the solid-node mask for z-aligned finite cylinders.

        A node is solid when its x-y radial distance to a cylinder axis is within
        the cylinder radius and its z lies within the cylinder's z-extent.

        Args:
            cylinders: List of ``(cx, cy, r)`` or ``(cx, cy, r, z_lo, z_hi)``.
                ``z_lo``/``z_hi`` are inclusive bounds in lattice-index units;
                without them the cylinder spans the full depth (``z_hi = nz``).

        Returns:
            Boolean array of shape (nx, ny, nz); ``True`` marks solid nodes.

        Raises:
            ValueError: in ``"pressure"`` mode, if a cylinder's x-extent reaches
                the inlet (x=0) or outlet (x=nx-1) plane, where the Zou-He
                density BC would overwrite the bounce-back result.
        """
        solid = np.zeros((self.nx, self.ny, self.nz), dtype=bool)
        if not cylinders:
            return solid
        ix = np.arange(self.nx)
        iy = np.arange(self.ny)
        iz = np.arange(self.nz)
        xx, yy, zz = np.meshgrid(ix, iy, iz, indexing="ij")
        for cyl in cylinders:
            cx, cy, cr = cyl[0], cyl[1], cyl[2]
            if self.streamwise_boundary == "pressure" and (
                cx - cr <= 0.0 or cx + cr >= self.nx - 1
            ):
                raise ValueError(
                    f"cylinder at x={cx} r={cr} reaches the pressure inlet/outlet "
                    "plane; the Zou-He BC would overwrite its bounce-back. Move it "
                    "into the interior (cx - r > 0 and cx + r < nx-1)."
                )
            z_lo = cyl[3] if len(cyl) > 3 else 0.0
            z_hi = cyl[4] if len(cyl) > 4 else float(self.nz)
            radial2 = (xx - cx) ** 2 + (yy - cy) ** 2
            in_z = (zz >= z_lo) & (zz <= z_hi)
            solid |= (radial2 <= cr**2) & in_z
        return solid

    # ------------------------------------------------------------------
    # Macroscopic fields
    # ------------------------------------------------------------------

    def _macroscopic(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return (rho, ux, uy, uz) from the distribution function.

        When a body force is present (Guo forcing), the velocity includes the
        half-force correction ``u = (Σ c_i f_i + F/2) / ρ`` so the reported
        velocity is the physical fluid velocity.

        Returns:
            Tuple of four arrays each with shape (nx, ny, nz).
        """
        rho = self.f.sum(axis=0)
        ux = (C3[:, 0, np.newaxis, np.newaxis, np.newaxis] * self.f).sum(axis=0)
        uy = (C3[:, 1, np.newaxis, np.newaxis, np.newaxis] * self.f).sum(axis=0)
        uz = (C3[:, 2, np.newaxis, np.newaxis, np.newaxis] * self.f).sum(axis=0)
        ux = (ux + 0.5 * self.Fx) / rho
        uy = (uy + 0.5 * self.Fy) / rho
        uz = (uz + 0.5 * self.Fz) / rho
        return rho, ux, uy, uz

    # ------------------------------------------------------------------
    # Collision
    # ------------------------------------------------------------------

    def _collide_bgk(
        self,
        rho: np.ndarray,
        ux: np.ndarray,
        uy: np.ndarray,
        uz: np.ndarray,
    ) -> None:
        """BGK collision step with Guo forcing (in-place update of ``self.f``).

        Adds the Guo (2002) discrete forcing source term
        ``S_i = (1 - 1/(2τ)) w_i [ (c_i-u)/cs² + (c_i·u)/cs⁴ c_i ] · F`` so the
        body force ``(Fx, Fy, Fz)`` drives the fluid consistently to second order.
        With zero force this reduces to plain BGK.

        Args:
            rho, ux, uy, uz: Current macroscopic fields, each shape (nx, ny, nz).
        """
        feq = _equilibrium_3d(rho, ux, uy, uz)
        self.f += self.omega * (feq - self.f)

        if not (self.Fx.any() or self.Fy.any() or self.Fz.any()):
            return
        prefactor = 1.0 - 0.5 * self.omega
        for i in range(Q3):
            cx, cy, cz = C3[i]
            cu = cx * ux + cy * uy + cz * uz  # c_i · u
            # (c_i - u)/cs² + (c_i·u) c_i / cs⁴, dotted with F.
            term_x = (cx - ux) / CS2_3 + cu * cx / CS2_3**2
            term_y = (cy - uy) / CS2_3 + cu * cy / CS2_3**2
            term_z = (cz - uz) / CS2_3 + cu * cz / CS2_3**2
            s_i = prefactor * W3[i] * (term_x * self.Fx + term_y * self.Fy + term_z * self.Fz)
            self.f[i] += s_i

    # ------------------------------------------------------------------
    # Streaming
    # ------------------------------------------------------------------

    def _stream(self) -> None:
        """Periodic streaming step using np.roll on all three axes."""
        streamed = np.empty_like(self.f)
        for i in range(Q3):
            cx, cy, cz = C3[i].astype(int)
            tmp = np.roll(self.f[i], cx, axis=0)
            tmp = np.roll(tmp, cy, axis=1)
            tmp = np.roll(tmp, cz, axis=2)
            streamed[i] = tmp
        self.f = streamed

    # ------------------------------------------------------------------
    # Lees-Edwards boundary correction
    # ------------------------------------------------------------------

    def _fractional_roll_x(self, row: np.ndarray, shift: float) -> np.ndarray:
        """Apply a sub-integer x-shift to a 1D or 2D slice using interpolation.

        Args:
            row: Array of shape (nx,) or (nx, nz).
            shift: Fractional shift in lattice units (positive = shift left,
               matching ``np.roll(row, -shift)`` semantics).

        Returns:
            Shifted array of the same shape as ``row``.
        """
        nx = self.nx
        int_shift = math.floor(shift) % nx
        frac = shift - math.floor(shift)
        rolled = np.roll(row, -int_shift, axis=0)
        if abs(frac) < 1e-12:
            return rolled
        if self.le_interpolation_order >= 3:
            try:
                from scipy.ndimage import map_coordinates

                x_dst = (np.arange(nx, dtype=float) - frac) % nx
                if row.ndim == 1:
                    return map_coordinates(rolled, [x_dst], order=3, mode="wrap")
                # 2D case: apply per-z column
                result = np.empty_like(rolled)
                for k in range(rolled.shape[1]):
                    result[:, k] = map_coordinates(rolled[:, k], [x_dst], order=3, mode="wrap")
                return result
            except ImportError:
                pass
        # Linear fallback
        return (1.0 - frac) * rolled + frac * np.roll(rolled, -1, axis=0)

    def _apply_le_bc(self) -> None:
        """Apply Lees-Edwards boundary correction after streaming.

        Two corrections are applied to populations that crossed the y-periodic
        boundary:

        1. **Fractional x-shift** — populations are rolled by ±``_le_shift``
           lattice units to account for the accumulated relative displacement
           of the image boxes.
        2. **Velocity boost** — the equilibrium contribution is updated to
           include the shear velocity jump ``dv = γ̇ · ny``.  Top-crossers
           entered from the upper image box (frame velocity +dv), so they
           receive a boost of **-dv** when entering the main box; bottom-crossers
           (from the lower image box at -dv) receive **+dv**.  The correction
           uses the first-order expansion:
           ``Δf_i = w_i · ρ · c_ix · (±dv) / cs²``.
        """
        shift = self._le_shift
        dv = self.le_shear_rate * self.ny  # velocity jump across the domain

        # --- Top crossers: landed at y=0 after streaming from y=ny-1 ---
        # These populations entered from the top image box (offset +dv in x).
        # Transforming to the main-box frame applies a boost of -dv.
        y0 = 0
        rho_top = self.f[:, :, y0, :].sum(axis=0)  # (nx, nz)
        for i in _LE_TOP_DIRS:
            # Spatial x-shift (top image is displaced +shift lattice units)
            self.f[i, :, y0, :] = self._fractional_roll_x(self.f[i, :, y0, :], shift)
            # Velocity boost: entering from +dv frame → subtract dv
            cx = C3[i, 0]
            self.f[i, :, y0, :] += W3[i] * rho_top * cx * (-dv) / CS2_3

        # --- Bottom crossers: landed at y=ny-1 after streaming from y=0 ---
        # These populations entered from the bottom image box (offset -dv in x).
        # Transforming to the main-box frame applies a boost of +dv.
        yn = self.ny - 1
        rho_bot = self.f[:, :, yn, :].sum(axis=0)  # (nx, nz)
        for i in _LE_BOT_DIRS:
            # Spatial x-shift (bottom image is displaced -shift lattice units)
            self.f[i, :, yn, :] = self._fractional_roll_x(self.f[i, :, yn, :], -shift)
            # Velocity boost: entering from -dv frame → add dv
            cx = C3[i, 0]
            self.f[i, :, yn, :] += W3[i] * rho_bot * cx * dv / CS2_3

    # ------------------------------------------------------------------
    # Pressure inlet/outlet boundary (Zou-He)
    # ------------------------------------------------------------------

    def _apply_pressure_bc(self) -> None:
        """Apply Zou-He density (pressure) BC on the x inlet/outlet planes.

        Imposes ``rho = rho_in`` at the inlet plane x=0 and ``rho = rho_out`` at
        the outlet plane x=nx-1.  The streamwise velocity is solved exactly from
        the density constraint, then the unknown incoming populations are
        reconstructed by the Zou-He bounce-back-of-non-equilibrium rule on the
        D3Q15 lattice.

        Only ``rho`` and the streamwise velocity ``ux`` are enforced exactly:
        the equilibrium used in the reconstruction assumes zero transverse
        velocity, which biases ``uy``/``uz`` toward zero but does not pin them.
        For the intended duct (periodic y, z) this is the desired behaviour.

        This mirrors the 2-D ``LBMDEMSolver`` Zou-He inlet/outlet contract
        (``rho_in = rho_out + pressure_drop / cs^2``).  Called after streaming;
        y and z remain periodic.  ``plane`` slices are NumPy views, so the
        in-place reconstruction writes straight through to ``self.f``.
        """
        f = self.f

        # --- Inlet plane x = 0 : impose rho = rho_in ----------------------
        plane = f[:, 0, :, :]  # (Q3, ny, nz) view into self.f
        # Known populations: cx <= 0 (rest + outgoing/tangential already streamed).
        sum_zero = plane[_X_ZERO_DIRS].sum(axis=0)
        sum_neg = plane[_X_NEG_DIRS].sum(axis=0)
        # rho = sum_zero + sum_neg + sum_pos, and momentum_x = sum_pos - sum_neg
        # = rho * ux  ->  ux = 1 - (sum_zero + 2*sum_neg) / rho.
        ux_in = 1.0 - (sum_zero + 2.0 * sum_neg) / self.rho_in
        self._zou_he_reconstruct(
            plane=plane,
            unknown_dirs=_X_POS_DIRS,
            rho=self.rho_in,
            u_normal=ux_in,
            axis=0,
        )

        # --- Outlet plane x = nx-1 : impose rho = rho_out -----------------
        plane = f[:, -1, :, :]  # view into self.f
        sum_zero = plane[_X_ZERO_DIRS].sum(axis=0)
        sum_pos = plane[_X_POS_DIRS].sum(axis=0)
        # ux = -(1 - (sum_zero + 2*sum_pos) / rho)  (flow leaves in +x)
        ux_out = -1.0 + (sum_zero + 2.0 * sum_pos) / self.rho_out
        self._zou_he_reconstruct(
            plane=plane,
            unknown_dirs=_X_NEG_DIRS,
            rho=self.rho_out,
            u_normal=ux_out,
            axis=0,
        )

    @staticmethod
    def _zou_he_reconstruct(
        plane: np.ndarray,
        unknown_dirs: list[int],
        rho: float,
        u_normal: np.ndarray,
        axis: int,
    ) -> None:
        """Reconstruct unknown populations on a boundary plane (Zou-He, D3Q15).

        For each unknown direction ``i`` the opposite (known) population is
        ``OPPOSITE3[i]``; the rule sets ``f_i = feq_i + (f_opp - feq_opp)``,
        i.e. bounce-back of the non-equilibrium part, where ``feq`` is evaluated
        with only the wall-normal velocity non-zero.

        Args:
            plane: Population slice at the boundary, shape (Q3, n1, n2). Modified
                in place.
            unknown_dirs: Indices of the incoming (unknown) populations.
            rho: Imposed density on the plane (scalar).
            u_normal: Imposed wall-normal velocity, shape (n1, n2).
            axis: Lattice velocity component index of the wall normal (0=x).
        """
        # Equilibrium with only the wall-normal velocity non-zero (uy=uz=0).
        # feq_i = w_i * rho * (1 + (c.u)/cs2 + (c.u)^2/(2 cs2^2) - u^2/(2 cs2))
        u2 = u_normal**2
        for i in unknown_dirs:
            opp = OPPOSITE3[i]
            cu = C3[i, axis] * u_normal
            feq_i = W3[i] * rho * (1.0 + cu / CS2_3 + 0.5 * cu**2 / CS2_3**2 - 0.5 * u2 / CS2_3)
            cu_o = C3[opp, axis] * u_normal
            feq_opp = (
                W3[opp] * rho * (1.0 + cu_o / CS2_3 + 0.5 * cu_o**2 / CS2_3**2 - 0.5 * u2 / CS2_3)
            )
            # Bounce-back of non-equilibrium: f_i = feq_i + (f_opp - feq_opp).
            plane[i] = feq_i + (plane[opp] - feq_opp)

    # ------------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Immersed-boundary particle coupling (issue #17)
    # ------------------------------------------------------------------

    @staticmethod
    def _sphere_markers(radius: float, spacing: float) -> np.ndarray:
        """Return surface marker offsets on a sphere via a Fibonacci lattice.

        Args:
            radius: Sphere radius.
            spacing: Target arc spacing between markers (lattice units).

        Returns:
            ``(m, 3)`` array of marker offsets from the sphere centre, roughly
            evenly distributed over the surface, with ``m`` chosen so each marker
            owns about ``spacing²`` of surface area.
        """
        area = 4.0 * np.pi * radius**2
        m = max(12, int(np.ceil(area / max(spacing, 1e-6) ** 2)))
        k = np.arange(m, dtype=float)
        # Spherical Fibonacci lattice.
        phi = np.arccos(1.0 - 2.0 * (k + 0.5) / m)  # polar angle
        golden = np.pi * (1.0 + 5.0**0.5)
        theta = golden * k  # azimuth
        sin_phi = np.sin(phi)
        offsets = radius * np.column_stack(
            (sin_phi * np.cos(theta), sin_phi * np.sin(theta), np.cos(phi))
        )
        return offsets

    def _interp_velocity_3d(
        self, pts: np.ndarray, ux: np.ndarray, uy: np.ndarray, uz: np.ndarray
    ) -> np.ndarray:
        """Trilinearly interpolate the fluid velocity at marker points.

        Periodic in all three axes (matches the periodic streaming); the inlet/
        outlet planes are sampled with periodic wrap, which is acceptable because
        markers sit on particle surfaces in the interior.

        Args:
            pts: ``(m, 3)`` marker coordinates.
            ux, uy, uz: Fluid velocity fields, shape (nx, ny, nz).

        Returns:
            ``(m, 3)`` interpolated velocities.
        """
        x0 = np.floor(pts[:, 0]).astype(int)
        y0 = np.floor(pts[:, 1]).astype(int)
        z0 = np.floor(pts[:, 2]).astype(int)
        tx = pts[:, 0] - x0
        ty = pts[:, 1] - y0
        tz = pts[:, 2] - z0
        x0 %= self.nx
        y0 %= self.ny
        z0 %= self.nz
        x1 = (x0 + 1) % self.nx
        y1 = (y0 + 1) % self.ny
        z1 = (z0 + 1) % self.nz
        out = np.zeros((pts.shape[0], 3))
        for fld, comp in ((ux, 0), (uy, 1), (uz, 2)):
            c000 = fld[x0, y0, z0]
            c100 = fld[x1, y0, z0]
            c010 = fld[x0, y1, z0]
            c110 = fld[x1, y1, z0]
            c001 = fld[x0, y0, z1]
            c101 = fld[x1, y0, z1]
            c011 = fld[x0, y1, z1]
            c111 = fld[x1, y1, z1]
            c00 = c000 * (1 - tx) + c100 * tx
            c10 = c010 * (1 - tx) + c110 * tx
            c01 = c001 * (1 - tx) + c101 * tx
            c11 = c011 * (1 - tx) + c111 * tx
            c0 = c00 * (1 - ty) + c10 * ty
            c1 = c01 * (1 - ty) + c11 * ty
            out[:, comp] = c0 * (1 - tz) + c1 * tz
        return out

    def _spread_forces_3d(self, pts: np.ndarray, fmarker: np.ndarray) -> None:
        """Spread marker forces onto the 8 surrounding lattice nodes (trilinear).

        Adds the spread force into ``self.Fx/Fy/Fz`` in place, periodic in all
        axes.  Solid (obstacle) nodes are skipped so reaction force is never
        deposited inside a fixed cylinder, mirroring the 2-D ``active = ~solid``
        spreading guard.

        Args:
            pts: ``(m, 3)`` marker coordinates.
            fmarker: ``(m, 3)`` force per marker to deposit on the fluid.
        """
        x0 = np.floor(pts[:, 0]).astype(int)
        y0 = np.floor(pts[:, 1]).astype(int)
        z0 = np.floor(pts[:, 2]).astype(int)
        tx = pts[:, 0] - x0
        ty = pts[:, 1] - y0
        tz = pts[:, 2] - z0
        x0 %= self.nx
        y0 %= self.ny
        z0 %= self.nz
        x1 = (x0 + 1) % self.nx
        y1 = (y0 + 1) % self.ny
        z1 = (z0 + 1) % self.nz
        has_solid = self.solid.any()
        for comp, F in ((0, self.Fx), (1, self.Fy), (2, self.Fz)):
            fc = fmarker[:, comp]
            for ix, wx in ((x0, 1 - tx), (x1, tx)):
                for iy, wy in ((y0, 1 - ty), (y1, ty)):
                    for iz, wz in ((z0, 1 - tz), (z1, tz)):
                        contrib = wx * wy * wz * fc
                        if has_solid:
                            contrib = np.where(self.solid[ix, iy, iz], 0.0, contrib)
                        np.add.at(F, (ix, iy, iz), contrib)

    def _apply_ibm(self, ux: np.ndarray, uy: np.ndarray, uz: np.ndarray) -> np.ndarray:
        """Direct-forcing IBM exchange for one step; returns per-particle reaction.

        Builds spherical markers for every particle, interpolates the fluid
        velocity there, computes the direct-forcing marker force
        ``f = stiffness · (u_body - u_fluid) · ds``, spreads ``+f`` back onto the
        fluid body-force field, and accumulates the particle reaction ``-Σ f`` and
        torque ``-Σ r × f``.

        Args:
            ux, uy, uz: Current fluid velocity fields, shape (nx, ny, nz).

        Returns:
            ``(n, 3)`` reaction force per particle (also stored as
            ``self._ibm_reaction``; the torque is stored as
            ``self._ibm_reaction_torque``).
        """
        dem = self.dem
        n = dem.n_p
        reaction = np.zeros((n, 3))
        reaction_torque = np.zeros((n, 3))
        all_pts: list[np.ndarray] = []
        all_f: list[np.ndarray] = []
        for i in range(n):
            offsets = self._sphere_markers(float(dem.radii[i]), self.ibm_marker_spacing)
            m = offsets.shape[0]
            ds = 4.0 * np.pi * dem.radii[i] ** 2 / m  # area per marker
            pts = dem.pos[i] + offsets
            # Rigid-body surface velocity: v + ω × r.
            u_body = dem.vel[i] + np.cross(dem.omega[i], offsets)
            u_fluid = self._interp_velocity_3d(pts, ux, uy, uz)
            fmarker = self.ibm_stiffness * (u_body - u_fluid) * ds
            reaction[i] = -fmarker.sum(axis=0)
            reaction_torque[i] = -np.cross(offsets, fmarker).sum(axis=0)
            all_pts.append(pts)
            all_f.append(fmarker)
        if all_pts:
            self._spread_forces_3d(np.vstack(all_pts), np.vstack(all_f))
        self._ibm_reaction = reaction
        self._ibm_reaction_torque = reaction_torque
        return reaction

    def _ibm_force_audit(self) -> tuple[np.ndarray, np.ndarray]:
        """Return (total spread force on fluid, total reaction on particles).

        Runs a single fresh IBM exchange from the current state without advancing
        time, for momentum-conservation checks.  With no solid obstacles the two
        totals sum to zero (Newton's 3rd law).  When cylinders are present the
        spread force that would land on solid nodes is suppressed (see
        :meth:`_spread_forces_3d`), so for a particle whose marker cloud overlaps
        an obstacle the totals no longer exactly cancel — the dropped momentum is
        absorbed by the (immovable) obstacle.

        Returns:
            ``(spread_total, reaction_total)`` each a length-3 vector.
        """
        self.Fx[:] = 0.0
        self.Fy[:] = 0.0
        self.Fz[:] = 0.0
        rho, ux, uy, uz = self._macroscopic()
        reaction = self._apply_ibm(ux, uy, uz)
        spread_total = np.array([self.Fx.sum(), self.Fy.sum(), self.Fz.sum()])
        return spread_total, reaction.sum(axis=0)

    def _apply_bounce_back(self) -> None:
        """Halfway bounce-back on fixed-cylinder solid nodes (no-slip).

        After streaming, every population on a solid node is replaced by the
        population that streamed from the opposite direction, enforcing a no-slip
        wall on the staircased obstacle surface.  Mirrors the 2-D solver's
        ``f[i, solid] = f[opposite[i], solid]``.  A no-op when there are no
        obstacles.
        """
        if not self.solid.any():
            return
        f_tmp = self.f.copy()
        for i in range(Q3):
            self.f[i, self.solid] = f_tmp[OPPOSITE3[i], self.solid]

    # ------------------------------------------------------------------
    # Left-inlet particle source (issue #20)
    # ------------------------------------------------------------------

    def _inlet_flow_rate(self, ux_inlet: np.ndarray) -> float:
        """Return the positive volumetric inflow through the inlet (x=0) plane.

        Sums ``max(ux, 0)`` over the whole y-z inlet cross-section, the 3D
        analogue of the 2-D ``_left_boundary_flow_rate``.

        Args:
            ux_inlet: Inlet-plane streamwise velocity ``ux[0]``, shape (ny, nz).

        Returns:
            Non-negative inlet flux in lattice volume per step.
        """
        return float(np.maximum(ux_inlet, 0.0).sum())

    def _sample_inlet_point(self, radius: float, ux_inlet: np.ndarray) -> tuple[float, float]:
        """Sample a (y, z) injection point on the inlet plane weighted by flux.

        Cells with higher inlet ``ux`` are more likely to receive a particle
        (uniform incoming concentration).  Before the flow develops, fall back to
        a deterministic sweep so startup is well behaved.

        Args:
            radius: Particle radius (keeps the centre clear of the y/z edges).
            ux_inlet: Inlet-plane streamwise velocity, shape (ny, nz).

        Returns:
            ``(y, z)`` injection coordinates.
        """
        weights = np.maximum(ux_inlet, 0.0).ravel()
        total = float(weights.sum())
        ny, nz = ux_inlet.shape
        if total > 1e-12:
            flat = int(self._rng_inlet.choice(weights.size, p=weights / total))
            jitter_y = self._rng_inlet.uniform(-0.45, 0.45)
            jitter_z = self._rng_inlet.uniform(-0.45, 0.45)
        else:
            flat = self._inlet_cursor % weights.size
            self._inlet_cursor += 1
            jitter_y = jitter_z = 0.0
        yc, zc = divmod(flat, nz)
        y = float(np.clip(yc + jitter_y, radius, ny - 1 - radius))
        z = float(np.clip(zc + jitter_z, radius, nz - 1 - radius))
        return y, z

    def _feed_inlet_particles(self, max_new: int = 8) -> None:
        """Inject queued particles at the inlet per the volume-flux budget.

        Each call accumulates ``phi · inlet_flow_rate`` into the volume budget and
        places spheres (cost ``(4/3)π r³`` each) just inside the inlet at
        flux-weighted, non-overlapping (y, z) points until the budget is spent or
        ``max_new`` is reached.  A no-op outside ``left_inlet`` mode.

        Args:
            max_new: Maximum particles to inject in a single step.
        """
        if self.particle_source != "left_inlet" or self.source_volume_fraction is None:
            return
        phi = self.source_volume_fraction
        if phi <= 0.0:
            return
        _, ux, _, _ = self._macroscopic()
        ux_inlet = ux[0]
        flow = self._inlet_flow_rate(ux_inlet)
        self.cumulative_inlet_flow_volume += flow
        self._inlet_volume_budget += phi * flow

        r = self.particle_radius
        particle_volume = (4.0 / 3.0) * np.pi * r**3
        added = 0
        attempts = 0
        max_attempts = max(12, 4 * max_new)
        while self._inlet_volume_budget >= particle_volume and added < max_new:
            if attempts >= max_attempts:
                break
            attempts += 1
            y, z = self._sample_inlet_point(r, ux_inlet)
            x = r + 0.75
            if not self._can_place_inlet(x, y, z, r):
                continue
            u_seed = float(max(ux_inlet[int(round(y)) % self.ny, int(round(z)) % self.nz], 0.0))
            self.dem.add_particles(
                np.array([[x, y, z]]),
                np.array([[u_seed, 0.0, 0.0]]),
                np.array([r]),
            )
            self._inlet_volume_budget -= particle_volume
            self.injected_particle_volume += particle_volume
            self.generated_particles += 1
            added += 1

    def _can_place_inlet(self, x: float, y: float, z: float, radius: float) -> bool:
        """Return whether an inlet particle fits without overlapping existing bodies.

        Rejects placement that would overlap an active particle or a fixed
        cylinder obstacle, mirroring the 2-D ``_can_place_inlet_particle`` which
        checks both particle and cylinder overlap.

        Args:
            x, y, z: Candidate particle centre.
            radius: Candidate particle radius.

        Returns:
            ``True`` when no active particle and no cylinder overlaps the candidate.
        """
        # Reject overlap with a fixed z-aligned cylinder (within its z-extent).
        for cyl in self.cylinders:
            cx, cy, cr = cyl[0], cyl[1], cyl[2]
            z_lo = cyl[3] if len(cyl) > 3 else 0.0
            z_hi = cyl[4] if len(cyl) > 4 else float(self.nz)
            if z_lo <= z <= z_hi and np.hypot(x - cx, y - cy) < 1.02 * (cr + radius):
                return False
        if self.dem is None or self.dem.n_p == 0:
            return True
        d = self.dem.pos - np.array([x, y, z])
        dist = np.sqrt((d**2).sum(axis=1))
        return bool(np.all(dist >= 1.02 * (self.dem.radii + radius)))

    def _remove_outflow_particles(self) -> None:
        """Delete particles whose sphere has fully passed the outlet.

        Removes any particle whose trailing edge is past the outlet, i.e.
        ``pos_x - r > nx`` (one lattice unit beyond the last node, matching the
        2-D convention).  A no-op outside ``left_inlet`` mode or with no particles.
        """
        if self.particle_source != "left_inlet" or self.dem is None or self.dem.n_p == 0:
            return
        mask = self.dem.pos[:, 0] - self.dem.radii > self.nx
        self.removed_particles += self.dem.remove_particles(mask)

    def advance(self, n_steps: int = 1) -> None:
        """Advance the solver by ``n_steps`` LBM steps.

        Each step: (optional inlet feed) → (optional IBM exchange) → collide (with
        Guo forcing) → stream → bounce-back → boundary correction → (optional DEM
        sub-steps) → (optional outflow removal).  In the default
        ``streamwise_boundary="periodic"`` mode the correction is the
        Lees-Edwards shear BC.  In ``"pressure"`` mode a Zou-He density BC is
        applied on the x inlet/outlet planes instead (y, z stay periodic).  When
        particles are present (``immersed_boundary`` coupling) the IBM reaction
        force drives an internal :class:`DEM3D` each step.  With
        ``particle_source="left_inlet"`` particles are injected at the inlet
        before the step and removed after they pass the outlet.

        Args:
            n_steps: Number of LBM time steps to execute.
        """
        for _ in range(n_steps):
            # Inject inlet particles first so they take part in this step's coupling.
            self._feed_inlet_particles()
            coupled = self.dem is not None and self.dem.n_p > 0

            if coupled:
                # Compute the IBM exchange against the un-forced velocity, then
                # re-read the macroscopic fields so the collision sees the new
                # body force's half-correction (consistent Guo scheme, matching
                # the 2-D path which recomputes velocity after IBM).
                self.Fx[:] = 0.0
                self.Fy[:] = 0.0
                self.Fz[:] = 0.0
                _, ux, uy, uz = self._macroscopic()
                self._apply_ibm(ux, uy, uz)

            rho, ux, uy, uz = self._macroscopic()
            self._collide_bgk(rho, ux, uy, uz)
            self._stream()
            self._apply_bounce_back()
            if self.streamwise_boundary == "pressure":
                self._apply_pressure_bc()
            else:
                self._le_shift = (self._le_shift + self.le_shear_rate * self.ny) % self.nx
                self._apply_le_bc()

            if coupled:
                self.dem.step(
                    1,
                    external_forces=self._ibm_reaction,
                    external_torques=self._ibm_reaction_torque,
                )

            self._remove_outflow_particles()
            self.step_count += 1

    def get_fields(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return macroscopic fields (rho, ux, uy, uz), each of shape (nx, ny, nz).

        Returns:
            4-tuple of arrays with shape matching the grid dimensions.
        """
        return self._macroscopic()
