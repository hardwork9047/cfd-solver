"""Fluid-independent 2D DEM contact law (issue #26).

The 2D :class:`particulate_flow.dem.solver.DEMSolver` historically computed the
contact-force magnitudes (normal Hertz, Coulomb tangential, Coulomb rolling) inline,
reading every material parameter off its coupled fluid solver (``self.sim``).  That
made the contact law impossible to exercise without constructing the whole fluid
solver.

``ContactLaw`` extracts those three pure functions into a small dataclass that owns
only the material parameters — no fluid-solver reference — so the contact physics can
be built and tested standalone.  ``DEMSolver`` delegates to it, so production results
are unchanged.  It is the 2D counterpart of the scalar-magnitude helpers on the 3D
:class:`particulate_flow.dem.contact3d.DEM3D` and shares the same clamp/Coulomb laws.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ContactLaw:
    """Material parameters and the pure 2D DEM contact-force laws.

    The three methods are the exact formulas previously inlined in ``DEMSolver``:
    a clamp-to-zero Hertz (or linear) normal force with velocity damping, a
    Coulomb-limited tangential force, and a Coulomb-limited rolling-resistance
    torque.  They depend only on the contact state passed in and the parameters
    held here — never on a fluid solver.

    Args:
        contact_model: ``"dem-hertz"`` (f_n ∝ overlap^1.5) or ``"dem-linear"``.
        k_n: Normal contact stiffness.
        damping: Normal viscous-damping coefficient.
        sliding_friction: Coulomb limit μ for the tangential force.
        tangential_damping: Viscous damping scale for tangential slip.
        rolling_friction: Master switch — when ``False`` both the tangential force
            and the rolling torque are disabled (matches the 2D solver semantics).
        rolling_friction_coeff: Rolling-resistance moment coefficient.
        rolling_damping: Viscous damping scale for angular velocity.
    """

    contact_model: str
    k_n: float
    damping: float
    sliding_friction: float
    tangential_damping: float
    rolling_friction: bool
    rolling_friction_coeff: float
    rolling_damping: float

    @classmethod
    def from_solver(cls, sim) -> ContactLaw:
        """Build a :class:`ContactLaw` from an object exposing the material params.

        Reads the seven scalar material attributes (``k_n, damping,
        sliding_friction, tangential_damping, rolling_friction,
        rolling_friction_coeff, rolling_damping``) off ``sim``.  The contact model
        is taken from ``sim.contact_model`` if present, else ``sim.particle_method``
        (which is where ``LBMDEMSolver`` stores it), else ``"dem-hertz"``.

        Note:
            ``DEMSolver`` does not use this constructor — it builds ``ContactLaw``
            directly so the contact model comes from its own constructor argument,
            which may differ from the coupled solver's.  ``from_solver`` is a
            convenience for callers that want the law a solver actually uses.

        Args:
            sim: An object carrying the material parameters (typically the LBM-DEM
                solver).

        Returns:
            A frozen ``ContactLaw`` carrying those parameter values.
        """
        contact_model = getattr(sim, "contact_model", None) or getattr(
            sim, "particle_method", "dem-hertz"
        )
        return cls(
            contact_model=contact_model,
            k_n=sim.k_n,
            damping=sim.damping,
            sliding_friction=sim.sliding_friction,
            tangential_damping=sim.tangential_damping,
            rolling_friction=sim.rolling_friction,
            rolling_friction_coeff=sim.rolling_friction_coeff,
            rolling_damping=sim.rolling_damping,
        )

    def normal_magnitude(self, overlap: float, v_n: float, mass: float) -> float:
        """Normal contact-force magnitude with damping, clamped to avoid tension.

        Args:
            overlap: Penetration depth (positive in contact).
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

    def tangential_magnitude(self, v_t: float, normal_force: float, mass: float) -> float:
        """Coulomb-limited tangential-force component opposing slip.

        Args:
            v_t: Tangential slip velocity along the contact tangent.
            normal_force: Current normal force magnitude.
            mass: Effective mass for the damping term.

        Returns:
            Tangential force clamped to ``±sliding_friction · normal_force``; zero
            when ``rolling_friction`` is disabled or there is no normal force.
        """
        if not self.rolling_friction or normal_force <= 0.0:
            return 0.0
        f_trial = -self.tangential_damping * float(np.sqrt(self.k_n * mass)) * v_t
        f_limit = self.sliding_friction * normal_force
        return float(np.clip(f_trial, -f_limit, f_limit))

    def rolling_torque(
        self, omega: float, normal_force: float, radius: float, mass: float
    ) -> float:
        """Coulomb-limited rolling-resistance torque opposing angular velocity.

        Args:
            omega: Particle angular velocity (scalar, 2D).
            normal_force: Current normal force magnitude.
            radius: Particle radius.
            mass: Particle mass.

        Returns:
            Torque clamped to ``±rolling_friction_coeff · normal_force · radius``;
            zero when ``rolling_friction`` is disabled or there is no normal force.
        """
        if not self.rolling_friction or normal_force <= 0.0:
            return 0.0
        torque_trial = -self.rolling_damping * float(np.sqrt(self.k_n * mass)) * radius**2 * omega
        torque_limit = self.rolling_friction_coeff * normal_force * radius
        return float(np.clip(torque_trial, -torque_limit, torque_limit))
