"""Tests for issue-15: 3D membrane-fouling config and runner plumbing.

Acceptance scenarios:
1. A 3D case extends fouling_supply_3d.json with minimal overrides.
2. run_lbm_dem.py --config <3d_case> completes and creates a result dir.
3. The result dir holds time_series.csv, summary.json, and ParaView VTK.
4. Adding a 3D cylinder geometry via extends enables obstacles.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from particulate_flow.io.config import SimulationConfig

_REPO_ROOT = Path(__file__).resolve().parents[1]
_TEMPLATES = _REPO_ROOT / "configs" / "lbm_dem" / "templates"
_GEOMETRIES = _REPO_ROOT / "configs" / "lbm_dem" / "geometries"
_CASES = _REPO_ROOT / "configs" / "lbm_dem" / "cases"


# ---------------------------------------------------------------------------
# Scenario 1 + 4: config templates exist and load
# ---------------------------------------------------------------------------


class TestConfigTemplates:
    def test_fouling_supply_3d_template_exists(self):
        assert (_TEMPLATES / "fouling_supply_3d.json").is_file()

    def test_template_declares_3d_pressure_flow(self):
        cfg = SimulationConfig.from_json(str(_TEMPLATES / "fouling_supply_3d.json"))
        v = cfg.values
        assert v.get("dimensions") == 3
        assert v.get("nz", 1) > 1
        assert v.get("streamwise_boundary") == "pressure"

    def test_three_d_cylinder_geometry_exists(self):
        # A z-aligned cylinder fragment that a case can extend.
        matches = list(_GEOMETRIES.glob("*3d*.json"))
        assert matches, "expected a 3D cylinder geometry fragment in geometries/"

    def test_minimal_case_extends_template(self):
        # A runnable case extends the 3D template with a tiny grid override.
        case = _CASES / "fouling_3d_smoke.json"
        assert case.is_file()
        cfg = SimulationConfig.from_json(str(case))
        assert cfg.values.get("dimensions") == 3


# ---------------------------------------------------------------------------
# Scenario 2 + 3: the 3D runner branch produces artifacts
# ---------------------------------------------------------------------------


class TestRunner3D:
    def _tiny_args(self, out_root: Path):
        import argparse

        return argparse.Namespace(
            dimensions=3,
            nx=24,
            ny=12,
            nz=12,
            reynolds_number=10.0,
            u_max=0.05,
            pressure_drop=3e-3,
            rho_out=1.0,
            streamwise_boundary="pressure",
            particle_radius=2.0,
            density_ratio=2.0,
            gravity=0.0,
            particle_source="left_inlet",
            particle_fluid_coupling="immersed_boundary",
            particle_volume_fraction=0.05,
            dem_substeps=4,
            total_steps=20,
            warmup_steps=5,
            snapshot_every=10,
            cylinder_spec=None,
            le_shear_rate=0.0,
            result_tag="test3d",
            output_root=str(out_root),
        )

    def test_run_3d_creates_artifacts(self, tmp_path):
        from particulate_flow.runner3d import run_3d

        out_dir = run_3d(self._tiny_args(tmp_path))
        out_dir = Path(out_dir)
        assert out_dir.is_dir()
        # Honors output_root: results land under tmp_path, not the repo tree.
        assert tmp_path in out_dir.parents
        assert (out_dir / "analysis" / "time_series.csv").is_file()
        assert (out_dir / "summary.json").is_file()
        assert (out_dir / "metadata.json").is_file()
        pv = out_dir / "paraview"
        assert pv.is_dir()
        assert list(pv.glob("*.vtk")), "expected at least one VTK file"
        assert list(pv.glob("*.pvd")), "expected a .pvd time series"

    def test_summary_reports_3d_metrics(self, tmp_path):
        from particulate_flow.runner3d import run_3d

        out_dir = Path(run_3d(self._tiny_args(tmp_path)))
        summary = json.loads((out_dir / "summary.json").read_text())
        assert summary.get("dimensions") == 3
        assert "final_particle_count" in summary
        assert np.isfinite(summary.get("max_speed", 0.0))


# ---------------------------------------------------------------------------
# Scenario 4: cylinder geometry enables a solid mask
# ---------------------------------------------------------------------------


class TestCylinderGeometry:
    def test_case_with_cylinder_geometry_has_solid(self, tmp_path):
        import argparse

        from particulate_flow.runner3d import build_3d_solver

        args = argparse.Namespace(
            dimensions=3,
            nx=32,
            ny=16,
            nz=12,
            reynolds_number=10.0,
            u_max=0.05,
            pressure_drop=3e-3,
            rho_out=1.0,
            streamwise_boundary="pressure",
            particle_radius=2.0,
            density_ratio=2.0,
            gravity=0.0,
            particle_source="none",
            particle_fluid_coupling="none",
            particle_volume_fraction=None,
            dem_substeps=4,
            total_steps=10,
            warmup_steps=0,
            snapshot_every=10,
            cylinder_spec=[(16.0, 8.0, 3.0)],
            le_shear_rate=0.0,
            result_tag=None,
            output_root=str(tmp_path),
        )
        sim = build_3d_solver(args)
        assert sim.solid.any()
