from pathlib import Path

from cfd_dem_lbm.cylinder_flow import (
    build_cylinder_flow_sim,
    compute_cylinder_flow_metrics,
    run_cylinder_flow,
)


def test_cylinder_flow_uses_fluid_only_fast_solver_path():
    sim = build_cylinder_flow_sim(
        nx=48,
        ny=24,
        reynolds_number=20.0,
        u_max=0.03,
        cylinders=[(18.0, 12.0, 3.0)],
        flow_condition="fixed_pressure",
        fluid_method="lbm-trt-guo",
    )
    assert sim.n_p == 0
    assert sim.fluid_method == "lbm-trt-guo"
    assert sim.cylinders == [(18.0, 12.0, 3.0)]
    assert sim.fixed_solid.any()

    sim.advance(3)
    metrics = compute_cylinder_flow_metrics(sim)
    assert metrics.step == 3
    assert metrics.max_speed >= 0.0
    assert metrics.reynolds_number >= 0.0


def test_cylinder_flow_runner_writes_outputs(tmp_path: Path):
    sim, metrics = run_cylinder_flow(
        output_dir=tmp_path,
        nx=36,
        ny=18,
        reynolds_number=10.0,
        u_max=0.02,
        cylinders=[(13.0, 9.0, 2.5)],
        reynolds_length=None,
        flow_condition="fixed_pressure",
        flow_control_gain=0.2,
        fluid_method="lbm-bgk-guo",
        steps=2,
        report_every=1,
    )
    assert sim.n_p == 0
    assert metrics[-1].step == 2
    assert (tmp_path / "metadata.json").exists()
    assert (tmp_path / "time_series.csv").exists()
    assert (tmp_path / "final_fields.npz").exists()
    assert (tmp_path / "paraview" / "final_fields.vtk").exists()
