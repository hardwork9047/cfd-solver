#!/usr/bin/env python3
"""Benchmark Numba-accelerated LBM-DEM component kernels."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cfd_dem_lbm.result_paths import program_results_dir
from cfd_dem_lbm import FastLBMDEM


@dataclass
class ComponentRow:
    component: str
    backend: str
    active_backend: str
    repeat: int
    calls: int
    seconds: float
    calls_per_second: float


def _make_particles(sim: FastLBMDEM, n_particles: int) -> None:
    rng = np.random.default_rng(1234)
    sim._delete_particles(np.ones(sim.n_p, dtype=bool))
    sim.n_p = n_particles
    sim.radii = rng.uniform(1.6, 2.4, n_particles)
    sim.masses = 2.0 * np.pi * sim.radii**2
    sim.inertias = 0.5 * sim.masses * sim.radii**2
    sim.pos = np.column_stack(
        (
            rng.uniform(5.0, sim.nx - 5.0, n_particles),
            rng.uniform(4.0, sim.ny - 4.0, n_particles),
        )
    )
    sim.vel = rng.normal(0.0, 0.02, (n_particles, 2))
    sim.omega_p = rng.normal(0.0, 0.01, n_particles)
    sim.forces_p = np.zeros((n_particles, 2))
    sim.torques_p = np.zeros(n_particles)
    sim.ibm_forces_p = np.zeros((n_particles, 2))
    sim.ibm_torques_p = np.zeros(n_particles)


def _build_sim(args: argparse.Namespace, backend: str, coupling: str) -> FastLBMDEM:
    sim = FastLBMDEM(
        nx=args.nx,
        ny=args.ny,
        Re=args.reynolds_number,
        u_max=args.u_max,
        fluid_method="lbm-bgk-guo",
        fluid_accelerator="numpy",
        compute_accelerator=backend,
        n_particles=1,
        particle_radius=2.0,
        gravity=0.0,
        rolling_friction=True,
        particle_fluid_coupling=coupling,
        cylinders=[
            (0.35 * args.nx, 0.35 * args.ny, 0.08 * args.ny),
            (0.55 * args.nx, 0.65 * args.ny, 0.07 * args.ny),
            (0.72 * args.nx, 0.50 * args.ny, 0.06 * args.ny),
        ],
        seed=99,
    )
    _make_particles(sim, args.n_particles)
    return sim


def _time_calls(fn, calls: int) -> float:
    start = time.perf_counter()
    for _ in range(calls):
        fn()
    return time.perf_counter() - start


def _benchmark_component(
    args: argparse.Namespace,
    component: str,
    backend: str,
    repeat: int,
) -> ComponentRow:
    if component == "particle_pair_loads":
        sim = _build_sim(args, backend, "point_force")
        sim.particle_attraction = True
        sim.attraction_strength = 1e-3
        sim.attraction_cutoff = 3.0
        if backend == "numpy":
            sim.dem_solver.pair_kernel = None

        def fn() -> None:
            sim.dem_solver.compute_loads(1.0)

        calls = args.pair_calls
    elif component == "dem_boundary_loads":
        sim = _build_sim(args, backend, "point_force")

        def fn() -> None:
            forces = np.zeros((sim.n_p, 2))
            torques = np.zeros(sim.n_p)
            if backend == "numba" and sim.uses_numba_compute:
                sim.dem_solver.compute_loads(1.0)
            else:
                sim.dem_solver._apply_wall_loads(forces, torques)
                sim.dem_solver._apply_cylinder_loads(forces, torques)

        calls = args.dem_calls
    elif component == "solid_mask":
        sim = _build_sim(args, backend, "solid_boundary")

        def fn() -> None:
            sim._update_particle_solid_mask()

        calls = args.mask_calls
    elif component == "ibm_markers":
        sim = _build_sim(args, backend, "immersed_boundary")
        _, ux, uy = sim.get_fields()

        def fn() -> None:
            sim.Fx[:] = sim.F_drive
            sim.Fy[:] = 0.0
            sim._apply_immersed_boundary_forces(ux, uy)

        calls = args.ibm_calls
    else:
        raise ValueError(component)

    if backend == "numba" and sim.uses_numba_compute:
        fn()  # compile/warm cache
    seconds = _time_calls(fn, calls)
    return ComponentRow(
        component=component,
        backend=backend,
        active_backend="numba" if sim.uses_numba_compute and backend == "numba" else "numpy",
        repeat=repeat,
        calls=calls,
        seconds=seconds,
        calls_per_second=calls / max(seconds, 1e-12),
    )


def _write_outputs(output_dir: Path, rows: list[ComponentRow], args: argparse.Namespace) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "numba_component_benchmark.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)

    summary = {}
    for component in sorted({row.component for row in rows}):
        summary[component] = {}
        for backend in ("numpy", "numba"):
            group = [row for row in rows if row.component == component and row.backend == backend]
            if not group:
                continue
            cps = [row.calls_per_second for row in group]
            secs = [row.seconds for row in group]
            summary[component][backend] = {
                "active_backend": group[0].active_backend,
                "median_calls_per_second": statistics.median(cps),
                "median_seconds": statistics.median(secs),
            }
        if "numpy" in summary[component] and "numba" in summary[component]:
            summary[component]["numba_vs_numpy_speedup"] = (
                summary[component]["numba"]["median_calls_per_second"]
                / max(summary[component]["numpy"]["median_calls_per_second"], 1e-12)
            )

    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "command": sys.argv,
        "configuration": vars(args),
        "summary": summary,
    }
    (output_dir / "numba_component_benchmark.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    lines = [
        "# Numba Component Benchmark",
        "",
        f"Generated: {report['created_at']}",
        f"Particles: {args.n_particles}, grid: {args.nx} x {args.ny}",
        "",
        "| component | numpy calls/s | numba calls/s | speedup |",
        "|---|---:|---:|---:|",
    ]
    for component, item in summary.items():
        numpy_cps = item.get("numpy", {}).get("median_calls_per_second", 0.0)
        numba_cps = item.get("numba", {}).get("median_calls_per_second", 0.0)
        speedup = item.get("numba_vs_numpy_speedup", 0.0)
        lines.append(f"| {component} | {numpy_cps:.3f} | {numba_cps:.3f} | {speedup:.3f}x |")
    (output_dir / "numba_component_benchmark.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=180)
    parser.add_argument("--ny", type=int, default=70)
    parser.add_argument("--n-particles", type=int, default=300)
    parser.add_argument("--reynolds-number", type=float, default=50.0)
    parser.add_argument("--u-max", type=float, default=0.04)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--pair-calls", type=int, default=30)
    parser.add_argument("--dem-calls", type=int, default=200)
    parser.add_argument("--mask-calls", type=int, default=100)
    parser.add_argument("--ibm-calls", type=int, default=100)
    parser.add_argument("--result-tag", default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    timestamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    parts = [args.result_tag, timestamp] if args.result_tag else [timestamp]
    output_dir = args.out_dir or program_results_dir(__file__, *parts)
    rows: list[ComponentRow] = []
    for component in ("particle_pair_loads", "dem_boundary_loads", "solid_mask", "ibm_markers"):
        for backend in ("numpy", "numba"):
            for repeat in range(1, args.repeats + 1):
                row = _benchmark_component(args, component, backend, repeat)
                rows.append(row)
                print(
                    f"{component:>20} {backend:>5} repeat {repeat}: "
                    f"{row.calls_per_second:.2f} calls/s"
                )
    _write_outputs(output_dir, rows, args)
    print(f"\nSaved benchmark outputs: {output_dir}")


if __name__ == "__main__":
    main()
