#!/usr/bin/env python3
"""Run a DEM-only gravity packing calculation."""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from particulate_flow.builder import build_dem_packing_solver
from particulate_flow.dem.packing import (
    write_metrics_csv,
    write_particles_vtk,
    write_pvd,
    write_summary,
)
from particulate_flow.geometry.pore import PoreGeometry
from particulate_flow.io.paths import program_results_dir


DEFAULTS = {
    "nx": 120,
    "ny": 220,
    "n_particles": 200,
    "particle_radius": 2.0,
    "radius_variation": 0.10,
    "density_ratio": 2.0,
    "gravity": 5e-3,
    "k_n": 120.0,
    "damping": 0.8,
    "linear_damping": 0.06,
    "dem_substeps": 5,
    "steps": 6000,
    "record_every": 250,
    "min_steps": 500,
    "settle_max_speed": 1e-2,
    "settle_kinetic_energy": 1.0,
    "seed": 42,
    "rolling_friction": True,
    "sliding_friction": 0.5,
    "tangential_damping": 0.5,
    "rolling_friction_coeff": 0.10,
    "rolling_damping": 0.35,
    "particle_attraction": False,
    "particle_repulsion": False,
    "attraction_strength": 1e-3,
    "repulsion_strength": 1e-3,
    "attraction_cutoff": 3.0,
    "repulsion_cutoff": 3.0,
    "attraction_min_gap": 0.05,
    "repulsion_min_gap": 0.05,
    "particle_method": "dem-hertz",
    "particle_search": "cell_list",
    "cylinders": [],
    "result_tag": "packing_200",
    "paraview_output": True,
}


def _load_config(path: str | None) -> dict:
    if not path:
        return {}
    config_path = Path(path)
    with config_path.open(encoding="utf-8") as handle:
        loaded = json.load(handle)
    return loaded


def _merge_config(args: argparse.Namespace) -> dict:
    config = DEFAULTS.copy()
    config.update(_load_config(args.config))
    for key, value in vars(args).items():
        if key == "config" or value is None:
            continue
        config[key] = value
    return config


def _cylinder_spec(value: list[list[float]] | list[tuple[float, float, float]]) -> list[tuple[float, float, float]]:
    return [tuple(float(v) for v in item) for item in value]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", help="JSON file with DEM packing settings")
    parser.add_argument("--nx", type=int)
    parser.add_argument("--ny", type=int)
    parser.add_argument("--n-particles", type=int, dest="n_particles")
    parser.add_argument("--particle-radius", type=float, dest="particle_radius")
    parser.add_argument("--radius-variation", type=float, dest="radius_variation")
    parser.add_argument("--density-ratio", type=float, dest="density_ratio")
    parser.add_argument("--gravity", type=float)
    parser.add_argument("--k-n", type=float, dest="k_n")
    parser.add_argument("--damping", type=float)
    parser.add_argument("--linear-damping", type=float, dest="linear_damping")
    parser.add_argument("--dem-substeps", type=int, dest="dem_substeps")
    parser.add_argument("--steps", type=int)
    parser.add_argument("--record-every", type=int, dest="record_every")
    parser.add_argument("--min-steps", type=int, dest="min_steps")
    parser.add_argument("--settle-max-speed", type=float, dest="settle_max_speed")
    parser.add_argument("--settle-kinetic-energy", type=float, dest="settle_kinetic_energy")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--rolling-friction", action=argparse.BooleanOptionalAction, dest="rolling_friction")
    parser.add_argument("--sliding-friction", type=float, dest="sliding_friction")
    parser.add_argument("--tangential-damping", type=float, dest="tangential_damping")
    parser.add_argument("--rolling-friction-coeff", type=float, dest="rolling_friction_coeff")
    parser.add_argument("--rolling-damping", type=float, dest="rolling_damping")
    parser.add_argument("--particle-attraction", action="store_true", default=None, dest="particle_attraction")
    parser.add_argument("--particle-repulsion", action="store_true", default=None, dest="particle_repulsion")
    parser.add_argument("--attraction-strength", type=float, dest="attraction_strength")
    parser.add_argument("--repulsion-strength", type=float, dest="repulsion_strength")
    parser.add_argument("--attraction-cutoff", type=float, dest="attraction_cutoff")
    parser.add_argument("--repulsion-cutoff", type=float, dest="repulsion_cutoff")
    parser.add_argument("--attraction-min-gap", type=float, dest="attraction_min_gap")
    parser.add_argument("--repulsion-min-gap", type=float, dest="repulsion_min_gap")
    parser.add_argument("--particle-method", choices=["dem-hertz", "dem-linear"], dest="particle_method")
    parser.add_argument("--particle-search", choices=["cell_list", "all_pairs"], dest="particle_search")
    parser.add_argument("--cylinder-spec", nargs=3, action="append", metavar=("X", "Y", "R"))
    parser.add_argument("--result-tag", dest="result_tag")
    parser.add_argument("--paraview-output", action=argparse.BooleanOptionalAction, dest="paraview_output")
    return parser.parse_args()


def _make_output_dir(config: dict) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = str(config["result_tag"]).strip() or "packing"
    return program_results_dir(__file__, f"{timestamp}_{tag}")


def main() -> None:
    args = parse_args()
    config = _merge_config(args)
    if args.cylinder_spec:
        config["cylinders"] = [[float(x), float(y), float(r)] for x, y, r in args.cylinder_spec]
    config["cylinders"] = _cylinder_spec(config.get("cylinders", []))

    out_dir = _make_output_dir(config)
    snapshot_dir = out_dir / "snapshots"
    paraview_dir = out_dir / "paraview"
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    if config["paraview_output"]:
        paraview_dir.mkdir(parents=True, exist_ok=True)
        PoreGeometry.from_cylinders(config["cylinders"]).write_vtk(paraview_dir / "cylinders.vtk")

    (out_dir / "config_resolved.json").write_text(
        json.dumps(config, indent=2, default=list),
        encoding="utf-8",
    )

    sim = build_dem_packing_solver(argparse.Namespace(**config))

    print("DEM gravity packing")
    print(f"  Particles : {config['n_particles']}")
    print(f"  Domain    : {config['nx']} x {config['ny']}")
    print(f"  Output    : {out_dir}")
    print(f"  Steps     : {config['steps']}  record_every={config['record_every']}")

    metrics_rows = [sim.metrics()]
    vtk_entries: list[tuple[int, Path]] = []
    settled = False
    t0 = time.perf_counter()

    def record(frame_index: int) -> None:
        snap = sim.snapshot()
        np.savez_compressed(
            snapshot_dir / f"snapshot_{frame_index:04d}_step_{sim.step_count:07d}.npz",
            **snap,
        )
        if config["paraview_output"]:
            vtk_path = paraview_dir / f"particles_{frame_index:04d}_step_{sim.step_count:07d}.vtk"
            write_particles_vtk(vtk_path, snap)
            vtk_entries.append((sim.step_count, vtk_path.relative_to(paraview_dir)))

    record(0)
    frame = 1
    for _ in range(config["steps"]):
        sim.advance(1)
        row = sim.metrics()
        if sim.step_count % config["record_every"] == 0:
            metrics_rows.append(row)
            record(frame)
            frame += 1
            print(
                f"step={row.step:6d} max_speed={row.max_speed:.4e} "
                f"KE={row.kinetic_energy:.4e} contacts={row.contact_count}"
            )

        if (
            sim.step_count >= config["min_steps"]
            and row.max_speed <= config["settle_max_speed"]
            and row.kinetic_energy <= config["settle_kinetic_energy"]
        ):
            settled = True
            metrics_rows.append(row)
            record(frame)
            break

    final = sim.metrics()
    if not metrics_rows or metrics_rows[-1].step != final.step:
        metrics_rows.append(final)
        record(frame)

    np.savez_compressed(
        out_dir / "particles_final.npz",
        pos=sim.pos,
        vel=sim.vel,
        radii=sim.radii,
        omega=sim.omega_p,
        force=sim.forces_p,
        torque=sim.torques_p,
    )
    write_metrics_csv(out_dir / "metrics.csv", metrics_rows)
    write_summary(out_dir / "packing_summary.json", config, final, settled)
    if config["paraview_output"]:
        write_particles_vtk(paraview_dir / "particles_final.vtk", sim.snapshot())
        write_pvd(paraview_dir / "particles_series.pvd", vtk_entries)

    elapsed = time.perf_counter() - t0
    print("\nDEM packing finished")
    print(f"  Settled       : {settled}")
    print(f"  Final step    : {sim.step_count}")
    print(f"  Final max vel : {final.max_speed:.4e}")
    print(f"  Contacts      : {final.contact_count}")
    print(f"  Bed height    : {final.bed_height:.3f}")
    print(f"  Packing frac. : {final.packing_fraction:.3f}")
    print(f"  Elapsed       : {elapsed:.2f} s")
    print(f"  Final NPZ     : {out_dir / 'particles_final.npz'}")
    if config["paraview_output"]:
        print(f"  ParaView      : {paraview_dir / 'particles_series.pvd'}")


if __name__ == "__main__":
    main()
