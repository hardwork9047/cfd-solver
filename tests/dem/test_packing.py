import numpy as np

from particulate_flow.dem.packing import DEMPackingSimulation, write_metrics_csv, write_particles_vtk


def test_dem_packing_runs_and_keeps_particles_in_domain(tmp_path):
    sim = DEMPackingSimulation(
        nx=40,
        ny=60,
        n_particles=8,
        particle_radius=1.5,
        radius_variation=0.0,
        gravity=1e-3,
        dem_substeps=2,
        seed=7,
    )

    initial_y = sim.pos[:, 1].mean()
    sim.advance(5)
    metrics = sim.metrics()

    assert sim.pos[:, 1].mean() < initial_y
    assert metrics.n_particles == 8
    assert np.all(sim.pos[:, 0] >= sim.radii + 0.5)
    assert np.all(sim.pos[:, 0] <= sim.nx - 1.5 - sim.radii)
    assert np.all(sim.pos[:, 1] >= sim.radii + 0.5)

    write_metrics_csv(tmp_path / "metrics.csv", [metrics])
    write_particles_vtk(tmp_path / "particles.vtk", sim.snapshot())

    assert (tmp_path / "metrics.csv").exists()
    assert (tmp_path / "particles.vtk").exists()
