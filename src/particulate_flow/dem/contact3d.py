"""Fresh 3D DEM contact physics (fluid-free core) for issue #18.

The 2D :class:`particulate_flow.dem.solver.DEMSolver` is pervasively 2D — scalar
angular velocity, a fixed ``tangent_from_normal``, ``forces[:, 1]`` gravity, and
two-column numba kernels.  This module implements the same Hertz contact model
lifted to full 3D vectors (angular velocity and torque are 3-vectors, the
tangential direction is the in-plane slip direction) without touching the 2D
path.

Scope: sphere-sphere, sphere-wall, and sphere-(z-aligned finite) cylinder normal
+ tangential Hertz contact, rolling resistance, gravity, and semi-implicit Euler
time integration with sub-stepping.  Fluid coupling (drag + IBM back-reaction) is
deliberately out of scope here — that is issue #17.
"""

from __future__ import annotations

import numpy as np

PARTICLE_METHODS = ("dem-hertz", "dem-linear")

# Valid wall identifiers and their (axis, plane-position, inward-normal-sign).
# inward sign points from the wall into the fluid domain.
_WALL_SPECS = {
    "y_min": (1, "min", +1.0),
    "y_max": (1, "max", -1.0),
    "x_min": (0, "min", +1.0),
    "x_max": (0, "max", -1.0),
    "z_min": (2, "min", +1.0),
    "z_max": (2, "max", -1.0),
}


class DEM3D:
    """3D discrete-element contact solver for spheres.

    Owns the 3D particle state and computes contact forces/torques for
    sphere-sphere, sphere-wall, and sphere-cylinder interactions using a Hertz
    (or linear) normal model with velocity damping, Coulomb-limited tangential
    friction, and Coulomb-limited rolling resistance.  Time integration is
    semi-implicit Euler with ``dem_substeps`` sub-steps per :meth:`step` call.

    Args:
        pos: Particle centres, shape (n, 3).
        vel: Particle velocities, shape (n, 3).
        radii: Particle radii, shape (n,).
        nx, ny, nz: Domain extent in lattice units (used for wall planes).
        density_ratio: ρ_particle / ρ_fluid; sets mass and buoyancy.
        gravity: Gravitational acceleration magnitude (lattice units); applied
            along ``gravity_dir``.
        gravity_dir: Unit direction of gravity (default -y).
        contact_model: ``"dem-hertz"`` (f_n ∝ overlap^1.5) or ``"dem-linear"``.
        k_n: Normal contact stiffness.
        damping: Normal viscous-damping coefficient (0-1 scale).
        sliding_friction: Coulomb limit μ for tangential force.
        tangential_damping: Viscous damping scale for tangential slip.
        rolling_friction: Enable tangential friction + rolling resistance.
        rolling_friction_coeff: Rolling-resistance moment coefficient.
        rolling_damping: Viscous damping scale for angular velocity.
        dem_substeps: Sub-steps per :meth:`step` call.
        walls: Iterable of wall ids from ``y_min, y_max, x_min, x_max, z_min,
            z_max`` that act as collision planes.
        cylinders: Z-aligned finite cylinders as ``(cx, cy, r)`` or
            ``(cx, cy, r, z_min, z_max)``; collide in the x-y plane within the
            z-extent (full z by default).
    """

    def __init__(
        self,
        pos: np.ndarray,
        vel: np.ndarray,
        radii: np.ndarray,
        nx: int,
        ny: int,
        nz: int,
        density_ratio: float = 2.0,
        gravity: float = 0.0,
        gravity_dir: tuple[float, float, float] = (0.0, -1.0, 0.0),
        contact_model: str = "dem-hertz",
        k_n: float = 50.0,
        damping: float = 0.4,
        sliding_friction: float = 0.5,
        tangential_damping: float = 0.4,
        rolling_friction: bool = True,
        rolling_friction_coeff: float = 0.05,
        rolling_damping: float = 0.35,
        dem_substeps: int = 4,
        walls: tuple[str, ...] = (),
        cylinders: list | tuple | None = None,
    ) -> None:
        if contact_model not in PARTICLE_METHODS:
            raise ValueError(f"contact_model must be one of {PARTICLE_METHODS}")
        for w in walls:
            if w not in _WALL_SPECS:
                raise ValueError(f"unknown wall {w!r}; valid: {sorted(_WALL_SPECS)}")

        self.pos = np.asarray(pos, dtype=float).reshape((-1, 3))
        self.vel = np.asarray(vel, dtype=float).reshape((-1, 3))
        self.radii = np.asarray(radii, dtype=float).reshape(-1)
        self.n_p = self.pos.shape[0]

        self.nx, self.ny, self.nz = int(nx), int(ny), int(nz)
        self.density_ratio = float(density_ratio)
        self.g = float(gravity)
        gdir = np.asarray(gravity_dir, dtype=float)
        self.gravity_dir = gdir / (np.linalg.norm(gdir) or 1.0)

        self.contact_model = contact_model
        self.k_n = float(k_n)
        self.damping = float(damping)
        self.sliding_friction = float(sliding_friction)
        self.tangential_damping = float(tangential_damping)
        self.rolling_friction = bool(rolling_friction)
        self.rolling_friction_coeff = float(rolling_friction_coeff)
        self.rolling_damping = float(rolling_damping)
        self.dem_substeps = int(dem_substeps)
        self.walls = tuple(walls)
        self.cylinders = [tuple(c) for c in (cylinders or [])]

        # Angular velocity as a 3-vector per particle.
        self.omega = np.zeros((self.n_p, 3))
        # Solid-sphere mass and moment of inertia.
        self.masses = self.density_ratio * (4.0 / 3.0) * np.pi * self.radii**3
        self.inertias = 0.4 * self.masses * self.radii**2

    # ------------------------------------------------------------------
    # Scalar contact-law helpers (magnitudes; directions handled by caller)
    # ------------------------------------------------------------------

    def _normal_magnitude(self, overlap: float, v_n: float, mass: float) -> float:
        """Normal contact-force magnitude, clamped ≥0 to avoid tension.

        Args:
            overlap: Penetration depth (positive when in contact).
            v_n: Normal component of relative velocity (approach negative).
            mass: Effective mass for the damping term.

        Returns:
            Non-negative normal force magnitude.
        """
        if self.contact_model == "dem-linear":
            f_n = self.k_n * overlap
        else:
            f_n = self.k_n * overlap**1.5
        f_damp = -self.damping * v_n * float(np.sqrt(self.k_n * mass))
        return max(f_n + f_damp, 0.0)

    def _tangential_magnitude(self, v_t: float, normal_force: float, mass: float) -> float:
        """Coulomb-limited tangential-force magnitude opposing slip.

        Args:
            v_t: Non-negative magnitude of the tangential slip velocity.
            normal_force: Current normal force magnitude.
            mass: Effective mass for the damping term.

        Returns:
            Tangential force component along the slip direction (already
            Coulomb-clamped to ±μ·normal_force).
        """
        if not self.rolling_friction or normal_force <= 0.0:
            return 0.0
        f_trial = -self.tangential_damping * float(np.sqrt(self.k_n * mass)) * v_t
        f_limit = self.sliding_friction * normal_force
        return float(np.clip(f_trial, -f_limit, f_limit))

    def _rolling_torque(
        self, omega_vec: np.ndarray, normal_force: float, radius: float, mass: float
    ) -> np.ndarray:
        """Coulomb-limited rolling-resistance torque vector opposing spin.

        Args:
            omega_vec: Angular-velocity vector of the particle.
            normal_force: Current normal force magnitude.
            radius: Particle radius.
            mass: Particle mass.

        Returns:
            Torque 3-vector opposing ``omega_vec`` (zero if no contact / disabled).
        """
        if not self.rolling_friction or normal_force <= 0.0:
            return np.zeros(3)
        omega_mag = float(np.linalg.norm(omega_vec))
        if omega_mag < 1e-12:
            return np.zeros(3)
        direction = omega_vec / omega_mag
        torque_trial = (
            self.rolling_damping * float(np.sqrt(self.k_n * mass)) * radius**2 * omega_mag
        )
        torque_limit = self.rolling_friction_coeff * normal_force * radius
        magnitude = min(torque_trial, torque_limit)
        return -magnitude * direction

    # ------------------------------------------------------------------
    # Contact assembly
    # ------------------------------------------------------------------

    def _add_contact(
        self,
        idx: int,
        normal: np.ndarray,
        overlap: float,
        other_vel: np.ndarray,
        other_omega: np.ndarray,
        other_radius: float,
        eff_mass: float,
        forces: np.ndarray,
        torques: np.ndarray,
        partner: int | None = None,
    ) -> None:
        """Apply normal + tangential + rolling loads of one contact.

        ``normal`` points from the contact surface toward particle ``idx`` (the
        direction it should be pushed).  Loads on ``idx`` are always applied.
        When ``partner`` is given (a sphere-sphere contact) the equal-and-opposite
        linear force, the partner's own tangential torque, and the partner's own
        rolling resistance are applied to ``partner`` too; ``other_vel`` /
        ``other_omega`` / ``other_radius`` must then describe that partner.

        Args:
            idx: Index of the particle receiving the load.
            normal: Unit contact normal pointing toward particle ``idx``.
            overlap: Penetration depth (positive in contact).
            other_vel: Velocity of the contacting body (0 for a static wall).
            other_omega: Angular velocity of the contacting body (0 for a wall).
            other_radius: Radius of the other body (0 for a wall/plane).
            eff_mass: Effective mass for damping terms.
            forces: (n,3) force accumulator (modified in place).
            torques: (n,3) torque accumulator (modified in place).
            partner: Index of the other particle for a sphere-sphere contact, or
                ``None`` for a contact against a static body (wall/cylinder).
        """
        r_i = self.radii[idx]
        # Relative velocity of i's surface w.r.t. the other body's surface.
        # Contact point on i is at -r_i * normal from its centre.
        contact_arm_i = -r_i * normal
        surf_vel_i = self.vel[idx] + np.cross(self.omega[idx], contact_arm_i)
        contact_arm_o = other_radius * normal
        surf_vel_o = other_vel + np.cross(other_omega, contact_arm_o)
        v_rel = surf_vel_i - surf_vel_o

        v_n = float(np.dot(v_rel, normal))
        f_n = self._normal_magnitude(overlap, v_n, eff_mass)
        forces[idx] += f_n * normal
        if partner is not None:
            forces[partner] -= f_n * normal  # Newton's 3rd law (normal part)

        # Tangential slip = relative surface velocity minus its normal part.
        v_t_vec = v_rel - v_n * normal
        v_t_mag = float(np.linalg.norm(v_t_vec))
        if v_t_mag > 1e-12 and f_n > 0.0:
            t_hat = v_t_vec / v_t_mag
            f_t_vec = self._tangential_magnitude(v_t_mag, f_n, eff_mass) * t_hat
            forces[idx] += f_t_vec
            torques[idx] += np.cross(contact_arm_i, f_t_vec)
            if partner is not None:
                forces[partner] -= f_t_vec  # Newton's 3rd law (tangential part)
                # Partner's contact arm points from its centre toward idx
                # (i.e. +other_radius * normal), and it feels -f_t_vec.
                torques[partner] += np.cross(contact_arm_o, -f_t_vec)

        torques[idx] += self._rolling_torque(self.omega[idx], f_n, r_i, self.masses[idx])
        if partner is not None:
            torques[partner] += self._rolling_torque(
                self.omega[partner], f_n, other_radius, self.masses[partner]
            )

    def compute_loads(self) -> tuple[np.ndarray, np.ndarray]:
        """Compute total force and torque on every particle.

        Returns:
            ``(forces, torques)`` each shaped (n, 3).
        """
        forces = np.zeros((self.n_p, 3))
        torques = np.zeros((self.n_p, 3))
        if self.n_p == 0:
            return forces, torques

        # Gravity (buoyancy-corrected).
        buoyancy_factor = 1.0 - 1.0 / self.density_ratio
        forces += (self.masses * self.g * buoyancy_factor)[:, None] * self.gravity_dir

        self._sphere_sphere_loads(forces, torques)
        self._wall_loads(forces, torques)
        self._cylinder_loads(forces, torques)
        return forces, torques

    def _sphere_sphere_loads(self, forces: np.ndarray, torques: np.ndarray) -> None:
        """Accumulate pairwise sphere-sphere contact loads (all-pairs)."""
        for i in range(self.n_p):
            for j in range(i + 1, self.n_p):
                dp = self.pos[j] - self.pos[i]
                dist = float(np.linalg.norm(dp))
                min_dist = self.radii[i] + self.radii[j]
                if dist <= 1e-10 or dist >= min_dist:
                    continue
                overlap = min_dist - dist
                # Normal pointing toward i (push i away from j).
                normal_i = -dp / dist
                eff_mass = (self.masses[i] + self.masses[j]) / 2.0
                self._add_contact(
                    idx=i,
                    normal=normal_i,
                    overlap=overlap,
                    other_vel=self.vel[j],
                    other_omega=self.omega[j],
                    other_radius=self.radii[j],
                    eff_mass=eff_mass,
                    forces=forces,
                    torques=torques,
                    partner=j,
                )

    def _wall_loads(self, forces: np.ndarray, torques: np.ndarray) -> None:
        """Accumulate sphere-wall contact loads for the configured walls."""
        extent = (self.nx, self.ny, self.nz)
        for wall in self.walls:
            axis, side, sign = _WALL_SPECS[wall]
            plane = 0.0 if side == "min" else float(extent[axis])
            normal = np.zeros(3)
            normal[axis] = sign  # points into the domain (toward the particle)
            for i in range(self.n_p):
                # Signed distance from particle centre to the wall plane.
                gap = (self.pos[i, axis] - plane) * sign
                overlap = self.radii[i] - gap
                if overlap <= 0.0:
                    continue
                self._add_contact(
                    idx=i,
                    normal=normal,
                    overlap=overlap,
                    other_vel=np.zeros(3),
                    other_omega=np.zeros(3),
                    other_radius=0.0,
                    eff_mass=self.masses[i],
                    forces=forces,
                    torques=torques,
                )

    def _cylinder_loads(self, forces: np.ndarray, torques: np.ndarray) -> None:
        """Accumulate sphere-(z-aligned finite cylinder) contact loads."""
        for cyl in self.cylinders:
            cx, cy, cr = cyl[0], cyl[1], cyl[2]
            z_lo = cyl[3] if len(cyl) > 3 else 0.0
            z_hi = cyl[4] if len(cyl) > 4 else float(self.nz)
            for i in range(self.n_p):
                if not (z_lo <= self.pos[i, 2] <= z_hi):
                    continue
                dx = self.pos[i, 0] - cx
                dy = self.pos[i, 1] - cy
                radial = float(np.hypot(dx, dy))
                min_dist = cr + self.radii[i]
                if radial <= 1e-10 or radial >= min_dist:
                    continue
                overlap = min_dist - radial
                # Normal points radially outward from the axis toward the sphere.
                normal = np.array([dx / radial, dy / radial, 0.0])
                self._add_contact(
                    idx=i,
                    normal=normal,
                    overlap=overlap,
                    other_vel=np.zeros(3),
                    other_omega=np.zeros(3),
                    other_radius=cr,
                    eff_mass=self.masses[i],
                    forces=forces,
                    torques=torques,
                )

    # ------------------------------------------------------------------
    # Time integration
    # ------------------------------------------------------------------

    def step(self, n_steps: int = 1) -> None:
        """Advance the DEM state by ``n_steps`` outer steps.

        Each outer step runs ``dem_substeps`` semi-implicit Euler sub-steps with
        ``dt_sub = 1 / dem_substeps`` so the effective outer ``dt`` is 1 lattice
        time unit, matching the 2D coupling convention.

        Args:
            n_steps: Number of outer time steps.
        """
        if self.n_p == 0:
            return
        dt_sub = 1.0 / self.dem_substeps
        for _ in range(n_steps):
            for _ in range(self.dem_substeps):
                forces, torques = self.compute_loads()
                self.vel += dt_sub * forces / self.masses[:, None]
                self.omega += dt_sub * torques / self.inertias[:, None]
                self.pos += dt_sub * self.vel
