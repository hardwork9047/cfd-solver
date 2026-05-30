"""Tests for issue-26: 2D DEM contact law decoupled from the fluid solver.

Acceptance scenarios:
1. The contact law (normal Hertz, tangential friction, rolling resistance) runs
   and is correct WITHOUT constructing any fluid solver.
2. DEMSolver's contact-law methods return identical values to a ContactLaw with the
   same params (the delegation is behaviour-preserving — no numerical regression).
3. The 2D ContactLaw and the 3D DEM3D share the same contact-law concepts.
"""

from __future__ import annotations

import numpy as np
import pytest

from particulate_flow.dem.contact_law import ContactLaw


def _law(**overrides) -> ContactLaw:
    params = dict(
        contact_model="dem-hertz",
        k_n=50.0,
        damping=0.4,
        sliding_friction=0.5,
        tangential_damping=0.4,
        rolling_friction=True,
        rolling_friction_coeff=0.05,
        rolling_damping=0.2,
    )
    params.update(overrides)
    return ContactLaw(**params)


# ---------------------------------------------------------------------------
# Scenario 1: standalone contact law (no fluid solver)
# ---------------------------------------------------------------------------


class TestStandaloneContactLaw:
    def test_normal_hertz_positive_overlap(self):
        law = _law()
        # No approach velocity → pure Hertz f = k_n * overlap^1.5.
        f = law.normal_magnitude(overlap=0.5, v_n=0.0, mass=1.0)
        assert f == 50.0 * 0.5**1.5

    def test_normal_linear_model(self):
        law = _law(contact_model="dem-linear")
        f = law.normal_magnitude(overlap=0.5, v_n=0.0, mass=1.0)
        assert f == 50.0 * 0.5

    def test_normal_clamped_to_zero_under_strong_separation_damping(self):
        # Large separating velocity makes the damping term dominate → clamp at 0.
        law = _law()
        f = law.normal_magnitude(overlap=0.01, v_n=1e3, mass=1.0)
        assert f == 0.0

    def test_tangential_coulomb_bounded(self):
        law = _law(sliding_friction=0.5)
        f_n = 10.0
        # A huge slip velocity is clamped to ±μ·f_n.
        f_t = law.tangential_magnitude(v_t=1e6, normal_force=f_n, mass=1.0)
        assert abs(f_t) <= 0.5 * f_n + 1e-12
        assert f_t == -0.5 * f_n  # opposes +v_t at the limit

    def test_tangential_zero_without_contact_or_rolling(self):
        assert _law().tangential_magnitude(v_t=1.0, normal_force=0.0, mass=1.0) == 0.0
        assert (
            _law(rolling_friction=False).tangential_magnitude(v_t=1.0, normal_force=10.0, mass=1.0)
            == 0.0
        )

    def test_rolling_torque_coulomb_bounded(self):
        law = _law(rolling_friction_coeff=0.05)
        f_n, r = 10.0, 2.0
        t = law.rolling_torque(omega=1e6, normal_force=f_n, radius=r, mass=1.0)
        assert abs(t) <= 0.05 * f_n * r + 1e-12

    def test_rolling_zero_without_contact_or_rolling(self):
        assert _law().rolling_torque(omega=1.0, normal_force=0.0, radius=2.0, mass=1.0) == 0.0
        assert (
            _law(rolling_friction=False).rolling_torque(
                omega=1.0, normal_force=10.0, radius=2.0, mass=1.0
            )
            == 0.0
        )


# ---------------------------------------------------------------------------
# Scenario 2: DEMSolver delegates to ContactLaw (behaviour-preserving)
# ---------------------------------------------------------------------------


class TestDelegationEquivalence:
    def _solver_and_law(self, contact_model="dem-hertz"):
        from particulate_flow import LBMDEMSolver
        from particulate_flow.dem.solver import DEMSolver

        sim = LBMDEMSolver(
            nx=40,
            ny=20,
            Re=100.0,
            u_max=0.0,
            n_particles=2,
            particle_radius=2.0,
            density_ratio=2.0,
            gravity=0.0,
            k_n=70.0,
            damping=0.3,
            rolling_friction=True,
            sliding_friction=0.4,
            tangential_damping=0.35,
            rolling_friction_coeff=0.06,
            rolling_damping=0.25,
            particle_method=contact_model,
            seed=1,
        )
        dem = DEMSolver(sim, contact_model=contact_model)
        # Reference law built with the DEMSolver's own contact_model (matching how
        # DEMSolver builds its internal law).
        law = ContactLaw(
            contact_model=contact_model,
            k_n=sim.k_n,
            damping=sim.damping,
            sliding_friction=sim.sliding_friction,
            tangential_damping=sim.tangential_damping,
            rolling_friction=sim.rolling_friction,
            rolling_friction_coeff=sim.rolling_friction_coeff,
            rolling_damping=sim.rolling_damping,
        )
        return dem, law

    def test_from_solver_reads_params(self):
        _, law = self._solver_and_law()
        assert law.k_n == 70.0
        assert law.sliding_friction == 0.4
        assert law.rolling_friction_coeff == 0.06

    @pytest.mark.parametrize("model", ["dem-hertz", "dem-linear"])
    def test_normal_matches(self, model):
        dem, law = self._solver_and_law(model)
        for overlap, v_n, mass in [(0.5, 0.0, 1.0), (0.2, -0.1, 2.0), (0.01, 5.0, 0.5)]:
            assert dem.normal_contact_magnitude(overlap, v_n, mass) == law.normal_magnitude(
                overlap, v_n, mass
            )

    def test_tangential_matches(self):
        dem, law = self._solver_and_law()
        for v_t, f_n, mass in [(0.1, 5.0, 1.0), (10.0, 2.0, 0.5), (0.0, 0.0, 1.0)]:
            assert dem.tangential_force_magnitude(v_t, f_n, mass) == law.tangential_magnitude(
                v_t, f_n, mass
            )

    def test_rolling_matches(self):
        dem, law = self._solver_and_law()
        for omega, f_n, r, mass in [(0.1, 5.0, 2.0, 1.0), (100.0, 3.0, 1.5, 0.5)]:
            assert dem.rolling_resistance_torque(omega, f_n, r, mass) == law.rolling_torque(
                omega, f_n, r, mass
            )


# ---------------------------------------------------------------------------
# Scenario 3: concept parity with the 3D DEM3D contact law
# ---------------------------------------------------------------------------


class TestThreeDConceptParity:
    def test_normal_magnitude_matches_dem3d(self):
        from particulate_flow.dem.contact3d import DEM3D

        law = _law(k_n=50.0, damping=0.4, contact_model="dem-hertz")
        dem3d = DEM3D(
            pos=np.zeros((1, 3)),
            vel=np.zeros((1, 3)),
            radii=np.array([2.0]),
            nx=10,
            ny=10,
            nz=10,
            k_n=50.0,
            damping=0.4,
            contact_model="dem-hertz",
        )
        # Same clamp-to-0 Hertz normal law in both 2D and 3D.
        for overlap, v_n, mass in [(0.5, 0.0, 1.0), (0.2, -0.3, 2.0), (0.01, 1e3, 1.0)]:
            assert law.normal_magnitude(overlap, v_n, mass) == dem3d._normal_magnitude(
                overlap, v_n, mass
            )
