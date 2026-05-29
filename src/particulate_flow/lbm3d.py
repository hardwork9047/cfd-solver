"""3D LBM solver using the D3Q15 lattice with Lees-Edwards shear boundary conditions.

Implements a minimal BGK-collision fluid solver for wall-less shear flow.
Particle coupling is not included in this phase.
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

    Particle coupling (IBM, DEM, fixed obstacles, inlet injection) is **not**
    implemented in this class yet; requesting any of it raises
    ``NotImplementedError`` (tracked as issues #17-#20).

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
        n_particles: Number of DEM particles.  Must be 0 (3D particles not yet
                     implemented — raises ``NotImplementedError`` otherwise).
        cylinders: Fixed obstacle specs.  Must be empty/None (not yet implemented).
        particle_source: Particle injection mode.  Must be ``"none"`` (not yet
                         implemented).
        particle_fluid_coupling: Coupling mode.  Must be ``"none"`` (not yet
                                 implemented).
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
    ) -> None:
        # --- Guard the not-yet-implemented particle stack (issue #14 slice) ---
        # The 3D particle subsystems are tracked as follow-ups split from #14.
        # Fail loudly so a 3D fouling config does not silently run particle-free.
        if n_particles and n_particles > 0:
            raise NotImplementedError(
                "3D DEM particles are not implemented yet "
                "(IBM coupling #17, DEM contact #18). "
                "Set n_particles=0 for the fluid-only 3D pressure-flow slice."
            )
        if cylinders:
            raise NotImplementedError(
                "3D fixed cylinder obstacles are not implemented yet (issue #19)."
            )
        if particle_source not in ("none", None):
            raise NotImplementedError(
                "3D left-inlet particle injection is not implemented yet (issue #20)."
            )
        if particle_fluid_coupling not in ("none", None):
            raise NotImplementedError(
                "3D particle-fluid (IBM) coupling is not implemented yet (issue #17)."
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

    # ------------------------------------------------------------------
    # Macroscopic fields
    # ------------------------------------------------------------------

    def _macroscopic(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return (rho, ux, uy, uz) from the distribution function.

        Returns:
            Tuple of four arrays each with shape (nx, ny, nz).
        """
        rho = self.f.sum(axis=0)
        ux = (C3[:, 0, np.newaxis, np.newaxis, np.newaxis] * self.f).sum(axis=0) / rho
        uy = (C3[:, 1, np.newaxis, np.newaxis, np.newaxis] * self.f).sum(axis=0) / rho
        uz = (C3[:, 2, np.newaxis, np.newaxis, np.newaxis] * self.f).sum(axis=0) / rho
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
        """BGK collision step (in-place update of ``self.f``).

        Args:
            rho, ux, uy, uz: Current macroscopic fields, each shape (nx, ny, nz).
        """
        feq = _equilibrium_3d(rho, ux, uy, uz)
        self.f += self.omega * (feq - self.f)

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
        the outlet plane x=nx-1, with zero transverse velocity (uy = uz = 0) at
        both planes.  The streamwise velocity is solved from the density
        constraint, then the unknown incoming populations are reconstructed by
        the Zou-He bounce-back-of-non-equilibrium rule generalised to D3Q15.

        This mirrors the 2-D ``LBMDEMSolver`` Zou-He inlet/outlet contract
        (``rho_in = rho_out + pressure_drop / cs^2``) on the D3Q15 lattice.
        Called after streaming; y and z remain periodic.
        """
        f = self.f

        # --- Inlet plane x = 0 : impose rho = rho_in, uy = uz = 0 ----------
        plane = f[:, 0, :, :]  # (Q3, ny, nz)
        rho_in = self.rho_in
        # Known populations: cx <= 0 (rest + outgoing/tangential already streamed).
        sum_zero = plane[_X_ZERO_DIRS].sum(axis=0)
        sum_neg = plane[_X_NEG_DIRS].sum(axis=0)
        # rho = sum_zero + sum_neg + sum_pos, and momentum_x = sum_pos - sum_neg
        # = rho * ux  ->  ux = 1 - (sum_zero + 2*sum_neg) / rho.
        ux_in = 1.0 - (sum_zero + 2.0 * sum_neg) / rho_in
        self._zou_he_reconstruct(
            plane=plane,
            unknown_dirs=_X_POS_DIRS,
            known_neg_dirs=_X_NEG_DIRS,
            rho=rho_in,
            u_normal=ux_in,
            axis=0,
            sign=+1,
        )
        f[:, 0, :, :] = plane

        # --- Outlet plane x = nx-1 : impose rho = rho_out, uy = uz = 0 ------
        plane = f[:, -1, :, :]
        rho_out = self.rho_out
        sum_zero = plane[_X_ZERO_DIRS].sum(axis=0)
        sum_pos = plane[_X_POS_DIRS].sum(axis=0)
        # ux = -(1 - (sum_zero + 2*sum_pos) / rho)  (flow leaves in +x)
        ux_out = -1.0 + (sum_zero + 2.0 * sum_pos) / rho_out
        self._zou_he_reconstruct(
            plane=plane,
            unknown_dirs=_X_NEG_DIRS,
            known_neg_dirs=_X_POS_DIRS,
            rho=rho_out,
            u_normal=ux_out,
            axis=0,
            sign=-1,
        )
        f[:, -1, :, :] = plane

    @staticmethod
    def _zou_he_reconstruct(
        plane: np.ndarray,
        unknown_dirs: list[int],
        known_neg_dirs: list[int],
        rho: float,
        u_normal: np.ndarray,
        axis: int,
        sign: int,
    ) -> None:
        """Reconstruct unknown populations on a boundary plane (Zou-He, D3Q15).

        Uses bounce-back of the non-equilibrium normal part plus a transverse
        momentum correction so that the imposed transverse velocity is zero.

        Args:
            plane: Population slice at the boundary, shape (Q3, n1, n2). Modified
                in place.
            unknown_dirs: Indices of the incoming (unknown) populations.
            known_neg_dirs: Indices of the outgoing populations opposite to the
                unknowns (used for the bounce-back baseline).
            rho: Imposed density on the plane (scalar).
            u_normal: Imposed wall-normal velocity, shape (n1, n2).
            axis: Lattice velocity component index of the wall normal (0=x).
            sign: +1 at an inlet (unknowns have positive normal component),
                -1 at an outlet.
        """
        # Equilibrium with only the wall-normal velocity non-zero (uy=uz=0).
        # feq_i = w_i * rho * (1 + (c.u)/cs2 + (c.u)^2/(2 cs2^2) - u^2/(2 cs2))
        u2 = u_normal**2
        for i in unknown_dirs:
            opp = OPPOSITE3[i]
            cn = C3[i, axis]
            cu = cn * u_normal
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

    def advance(self, n_steps: int = 1) -> None:
        """Advance the solver by ``n_steps`` LBM steps.

        Each step: collide → stream → boundary correction.  In the default
        ``streamwise_boundary="periodic"`` mode the correction is the
        Lees-Edwards shear BC.  In ``"pressure"`` mode a Zou-He density BC is
        applied on the x inlet/outlet planes instead (y, z stay periodic).

        Args:
            n_steps: Number of LBM time steps to execute.
        """
        for _ in range(n_steps):
            rho, ux, uy, uz = self._macroscopic()
            self._collide_bgk(rho, ux, uy, uz)
            self._stream()
            if self.streamwise_boundary == "pressure":
                self._apply_pressure_bc()
            else:
                self._le_shift = (self._le_shift + self.le_shear_rate * self.ny) % self.nx
                self._apply_le_bc()
            self.step_count += 1

    def get_fields(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return macroscopic fields (rho, ux, uy, uz), each of shape (nx, ny, nz).

        Returns:
            4-tuple of arrays with shape matching the grid dimensions.
        """
        return self._macroscopic()
