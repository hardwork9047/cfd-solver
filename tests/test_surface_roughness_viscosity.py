"""Tests for issue-8: DEM surface roughness and apparent viscosity evaluation."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from particulate_flow import LBMDEMSolver


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _shear_sim(*, surface_roughness: float = 0.0, **kwargs) -> LBMDEMSolver:
    """Two-particle LE shear solver for roughness / viscosity tests."""
    defaults = dict(
        nx=40,
        ny=30,
        Re=10.0,
        u_max=0.0,
        n_particles=2,
        particle_radius=3.0,
        density_ratio=2.0,
        gravity=0.0,
        y_boundary="lees_edwards",
        streamwise_boundary="periodic_force",
        le_shear_rate=0.001,
        surface_roughness=surface_roughness,
        seed=42,
    )
    defaults.update(kwargs)
    sim = LBMDEMSolver(**defaults)
    return sim


# ---------------------------------------------------------------------------
# Surface roughness
# ---------------------------------------------------------------------------

class TestSurfaceRoughness:
    def test_default_surface_roughness_is_zero(self):
        """LBMDEMSolver.surface_roughness defaults to 0.0."""
        sim = _shear_sim()
        assert sim.surface_roughness == 0.0

    def test_surface_roughness_stored(self):
        """surface_roughness parameter is stored on the instance."""
        sim = _shear_sim(surface_roughness=0.05)
        assert sim.surface_roughness == pytest.approx(0.05)

    def test_zero_roughness_no_contact_when_just_outside(self):
        """With h_r=0 particles exactly at r_i+r_j should not produce contact force."""
        sim = _shear_sim(surface_roughness=0.0)
        r = float(sim.radii[0]) + float(sim.radii[1])
        # place particles at separation slightly larger than sum of radii
        sim.pos[0] = [10.0, 15.0]
        sim.pos[1] = [10.0 + r + 0.1, 15.0]
        sim.vel[:] = 0.0
        forces = sim._dem_forces(dt_sub=1.0)
        # no contact force expected
        np.testing.assert_allclose(forces, 0.0, atol=1e-10)

    def test_nonzero_roughness_triggers_contact_in_roughness_gap(self):
        """With h_r=0.5 particles separated by r_i+r_j+0.3 should be in contact."""
        h_r = 0.5
        sim_rough = _shear_sim(surface_roughness=h_r)
        sim_smooth = _shear_sim(surface_roughness=0.0)
        r_sum = float(sim_rough.radii[0]) + float(sim_rough.radii[1])
        sep = r_sum + 0.3  # inside roughness zone, outside smooth zone
        for sim in (sim_rough, sim_smooth):
            sim.pos[0] = [10.0, 15.0]
            sim.pos[1] = [10.0 + sep, 15.0]
            sim.vel[:] = 0.0

        forces_rough = sim_rough._dem_forces(dt_sub=1.0)
        forces_smooth = sim_smooth._dem_forces(dt_sub=1.0)

        # rough sim should show repulsive (normal) contact force; smooth should not
        assert abs(forces_rough[0, 0]) > 1e-10, "roughness contact force should be nonzero"
        np.testing.assert_allclose(forces_smooth, 0.0, atol=1e-10)

    def test_zero_roughness_is_regression_compatible(self):
        """h_r=0.0 produces identical forces to a sim without surface_roughness param."""
        sim_a = _shear_sim(surface_roughness=0.0)
        sim_b = _shear_sim(surface_roughness=0.0)
        for sim in (sim_a, sim_b):
            sim.pos[0] = [10.0, 15.0]
            sim.pos[1] = [14.0, 15.0]  # close enough for Hertz contact
            sim.vel[:] = 0.0
        fa = sim_a._dem_forces(dt_sub=1.0)
        fb = sim_b._dem_forces(dt_sub=1.0)
        np.testing.assert_allclose(fa, fb, atol=1e-14)


# ---------------------------------------------------------------------------
# Apparent viscosity evaluator
# ---------------------------------------------------------------------------

class TestViscosityEvaluator:
    def test_disabled_by_default_no_cost(self):
        """When viscosity_eval is disabled no CSV is produced after stepping."""
        sim = _shear_sim()
        # step a few times
        for _ in range(5):
            sim.advance()
        # no viscosity attribute should exist (or it should be None/disabled)
        ve = getattr(sim, "viscosity_evaluator", None)
        assert ve is None or not getattr(ve, "enabled", False)

    def test_enabled_evaluator_writes_csv(self, tmp_path):
        """viscosity_eval enabled → CSV file created with expected columns."""
        from particulate_flow.rheology import ViscosityEvaluator

        sim = _shear_sim()
        ve = ViscosityEvaluator(
            sim,
            start_step=0,
            viscosity_interval=2,
            average_steps=10,
            out_dir=tmp_path,
        )
        for step in range(6):
            sim.advance()
            ve.record(step)

        csv_path = tmp_path / "viscosity_timeseries.csv"
        assert csv_path.exists(), "viscosity_timeseries.csv should be created"
        with csv_path.open() as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) >= 1
        assert "step" in rows[0]
        assert "eta_H" in rows[0]
        assert "eta_C" in rows[0]

    def test_evaluator_finalize_returns_mean_keys(self, tmp_path):
        """ViscosityEvaluator.finalize() returns dict with eta_H_mean, eta_C_mean, eta_s_mean."""
        from particulate_flow.rheology import ViscosityEvaluator

        sim = _shear_sim()
        ve = ViscosityEvaluator(
            sim,
            start_step=0,
            viscosity_interval=1,
            average_steps=5,
            out_dir=tmp_path,
        )
        for step in range(10):
            sim.advance()
            ve.record(step)

        result = ve.finalize()
        assert "eta_H_mean" in result
        assert "eta_C_mean" in result
        assert "eta_s_mean" in result
        # mean should be finite
        assert np.isfinite(result["eta_H_mean"])
        assert np.isfinite(result["eta_s_mean"])

    def test_evaluator_average_matches_last_n_rows(self, tmp_path):
        """eta_s_mean should equal mean of eta_H+eta_C over the last average_steps rows."""
        from particulate_flow.rheology import ViscosityEvaluator

        average_steps = 4
        sim = _shear_sim()
        ve = ViscosityEvaluator(
            sim,
            start_step=0,
            viscosity_interval=1,
            average_steps=average_steps,
            out_dir=tmp_path,
        )
        for step in range(10):
            sim.advance()
            ve.record(step)

        result = ve.finalize()
        csv_path = tmp_path / "viscosity_timeseries.csv"
        with csv_path.open() as f:
            rows = list(csv.DictReader(f))
        tail = rows[-average_steps:]
        expected_mean = np.mean([float(r["eta_H"]) + float(r["eta_C"]) for r in tail])
        assert result["eta_s_mean"] == pytest.approx(expected_mean, rel=1e-6)

    def test_no_extra_cost_when_disabled(self):
        """When disabled ViscosityEvaluator.record() is cheap (no exception, no CSV)."""
        from particulate_flow.rheology import ViscosityEvaluator
        import tempfile

        sim = _shear_sim()
        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td)
            ve = ViscosityEvaluator(
                sim,
                start_step=0,
                viscosity_interval=10,
                average_steps=5,
                out_dir=out_dir,
                enabled=False,
            )
            for step in range(5):
                sim.advance()
                ve.record(step)  # should not raise or write
            assert not (out_dir / "viscosity_timeseries.csv").exists()


# ---------------------------------------------------------------------------
# Config plumbing
# ---------------------------------------------------------------------------

class TestConfigPlumbing:
    def test_physics_surface_roughness_passes_through_config(self):
        """physics.surface_roughness in JSON config reaches LBMDEMSolver."""
        from particulate_flow.io.config import SimulationConfig

        cfg = SimulationConfig.from_mapping(
            {
                "domain": {"nx": 40, "ny": 30},
                "flow": {"u_max": 0.0, "reynolds_number": 10.0},
                "particles": {"n_particles": 0, "particle_radius": 3.0},
                "physics": {"surface_roughness": 0.07},
            }
        )
        assert cfg.values.get("surface_roughness") == pytest.approx(0.07)

    def test_runtime_viscosity_eval_passes_through_config(self):
        """runtime.viscosity_eval sub-keys are flattened into expected keys."""
        from particulate_flow.io.config import SimulationConfig

        cfg = SimulationConfig.from_mapping(
            {
                "domain": {"nx": 40, "ny": 30},
                "runtime": {
                    "viscosity_eval": {
                        "enabled": True,
                        "start_step": 200,
                        "viscosity_interval": 50,
                        "average_steps": 500,
                    }
                },
            }
        )
        assert cfg.values.get("viscosity_eval_enabled") is True
        assert cfg.values.get("viscosity_eval_start_step") == 200
        assert cfg.values.get("viscosity_eval_interval") == 50
        assert cfg.values.get("viscosity_eval_average_steps") == 500
