#!/usr/bin/env python3
"""Benchmark LBM-DEM coupling modes across NumPy, Numba, and auto backends."""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import statistics
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from particulate_flow.result_paths import program_results_dir
from particulate_flow import FastLBMDEM


@dataclass
class CouplingBackendRow:
    coupling: str
    backend_config: str
    fluid_accelerator: str
    compute_accelerator: str
    active_fluid_backend: str
    active_compute_backend: str
    repeat: int
    timed_steps: int
    timed_seconds: float
    steps_per_second: float
    final_particles: int
    final_max_speed: float
    dynamic_solid_fraction: float


def _quiet_solver(**kwargs) -> FastLBMDEM:
    """Build the solver while suppressing the setup banner during benchmarks."""
    with contextlib.redirect_stdout(io.StringIO()):
        return FastLBMDEM(**kwargs)


def _make_particles(sim: FastLBMDEM, n_particles: int, seed: int) -> None:
    """Create a deterministic particle cloud shared by all benchmark cases."""
    rng = np.random.default_rng(seed)
    if sim.n_p:
        sim._delete_particles(np.ones(sim.n_p, dtype=bool))
    sim.n_p = n_particles
    sim.radii = rng.uniform(0.9 * sim.r_p, 1.1 * sim.r_p, n_particles)
    sim.masses = sim.density_ratio * np.pi * sim.radii**2
    sim.inertias = 0.5 * sim.masses * sim.radii**2
    spacing = max(2.8 * float(np.max(sim.radii)), 4.0)
    x_values = np.arange(7.0, sim.nx - 7.0, spacing)
    y_values = np.arange(5.0, sim.ny - 5.0, spacing)
    candidates = np.array([(x_pos, y_pos) for x_pos in x_values for y_pos in y_values])
    if len(candidates) < n_particles:
        raise ValueError(
            "Benchmark grid is too small for the requested non-overlapping particles. "
            "Increase --nx/--ny or reduce --n-particles."
        )
    order = rng.permutation(len(candidates))[:n_particles]
    jitter = rng.uniform(-0.12 * sim.r_p, 0.12 * sim.r_p, (n_particles, 2))
    sim.pos = candidates[order] + jitter
    sim.vel = rng.normal(0.0, 0.002, (n_particles, 2))
    sim.omega_p = rng.normal(0.0, 0.001, n_particles)
    sim.forces_p = np.zeros((n_particles, 2))
    sim.torques_p = np.zeros(n_particles)
    sim.ibm_forces_p = np.zeros((n_particles, 2))
    sim.ibm_torques_p = np.zeros(n_particles)
    sim._update_particle_solid_mask()


def _build_sim(
    args: argparse.Namespace,
    coupling: str,
    fluid_accelerator: str,
    compute_accelerator: str,
    seed: int,
) -> FastLBMDEM:
    sim = _quiet_solver(
        nx=args.nx,
        ny=args.ny,
        Re=args.reynolds_number,
        u_max=args.u_max,
        reynolds_length=2.0 * args.particle_radius,
        flow_control=args.flow_control,
        fluid_method=args.fluid_method,
        fluid_accelerator=fluid_accelerator,
        compute_accelerator=compute_accelerator,
        particle_fluid_coupling=coupling,
        ibm_marker_spacing=args.ibm_marker_spacing,
        n_particles=1,
        particle_radius=args.particle_radius,
        density_ratio=args.density_ratio,
        gravity=0.0,
        rolling_friction=True,
        cylinders=[(args.cylinder_x, args.cylinder_y, args.cylinder_radius)],
        seed=seed,
    )
    _make_particles(sim, args.n_particles, seed)
    return sim


def _run_case(
    args: argparse.Namespace,
    coupling: str,
    label: str,
    fluid_accelerator: str,
    compute_accelerator: str,
    repeat: int,
) -> CouplingBackendRow:
    seed = args.seed + 1000 * repeat
    sim = _build_sim(args, coupling, fluid_accelerator, compute_accelerator, seed)
    if args.warmup_steps:
        sim.advance(args.warmup_steps)
    start = time.perf_counter()
    sim.advance(args.timed_steps)
    timed_seconds = time.perf_counter() - start
    _, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)
    fluid_mask = ~sim.solid
    return CouplingBackendRow(
        coupling=coupling,
        backend_config=label,
        fluid_accelerator=fluid_accelerator,
        compute_accelerator=compute_accelerator,
        active_fluid_backend="numba" if sim.uses_numba_lbm else "numpy",
        active_compute_backend="numba" if sim.uses_numba_compute else "numpy",
        repeat=repeat,
        timed_steps=args.timed_steps,
        timed_seconds=timed_seconds,
        steps_per_second=args.timed_steps / max(timed_seconds, 1e-12),
        final_particles=sim.n_p,
        final_max_speed=float(speed[fluid_mask].max()) if np.any(fluid_mask) else 0.0,
        dynamic_solid_fraction=sim.dynamic_solid_fraction(),
    )


def _summarise(rows: list[CouplingBackendRow]) -> dict:
    summary: dict[str, dict] = {}
    for coupling in sorted({row.coupling for row in rows}):
        base_key = (coupling, "numpy/numpy")
        base_rows = [
            row.steps_per_second
            for row in rows
            if row.coupling == base_key[0] and row.backend_config == base_key[1]
        ]
        base_speed = statistics.median(base_rows) if base_rows else 0.0
        summary[coupling] = {}
        for label in sorted({row.backend_config for row in rows if row.coupling == coupling}):
            group = [row for row in rows if row.coupling == coupling and row.backend_config == label]
            speeds = [row.steps_per_second for row in group]
            summary[coupling][label] = {
                "active_fluid_backend": group[0].active_fluid_backend,
                "active_compute_backend": group[0].active_compute_backend,
                "median_steps_per_second": statistics.median(speeds),
                "min_steps_per_second": min(speeds),
                "max_steps_per_second": max(speeds),
                "speedup_vs_numpy_numpy": (
                    statistics.median(speeds) / base_speed if base_speed > 0.0 else 0.0
                ),
            }
    return summary


def _write_outputs(output_dir: Path, rows: list[CouplingBackendRow], args: argparse.Namespace) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "lbm_dem_coupling_backend_benchmark.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)

    summary = _summarise(rows)
    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "command": sys.argv,
        "configuration": vars(args),
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }
    (output_dir / "lbm_dem_coupling_backend_benchmark.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    lines = [
        "# LBM-DEM Coupling Backend Benchmark",
        "",
        f"Generated: {report['created_at']}",
        f"Grid: {args.nx} x {args.ny}, particles: {args.n_particles}, timed steps: {args.timed_steps}",
        "",
        "| coupling | config | active backend | median steps/s | range steps/s | vs numpy/numpy |",
        "|---|---|---|---:|---:|---:|",
    ]
    config_order = [label for label, _, _ in args.backend_configs]
    for coupling in args.couplings:
        for label in config_order:
            if label not in summary.get(coupling, {}):
                continue
            item = summary[coupling][label]
            lines.append(
                f"| {coupling} | {label} | "
                f"fluid={item['active_fluid_backend']}, compute={item['active_compute_backend']} | "
                f"{item['median_steps_per_second']:.3f} | "
                f"{item['min_steps_per_second']:.3f}-{item['max_steps_per_second']:.3f} | "
                f"{item['speedup_vs_numpy_numpy']:.3f}x |"
            )
    (output_dir / "lbm_dem_coupling_backend_benchmark.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def _parse_backend_config(value: str) -> tuple[str, str, str]:
    try:
        fluid, compute = value.split("/", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("backend config must be FLUID/COMPUTE") from exc
    allowed = {"numpy", "numba", "auto"}
    if fluid not in allowed or compute not in allowed:
        raise argparse.ArgumentTypeError("backend names must be numpy, numba, or auto")
    return value, fluid, compute


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=96)
    parser.add_argument("--ny", type=int, default=40)
    parser.add_argument("--n-particles", type=int, default=40)
    parser.add_argument("--particle-radius", type=float, default=1.5)
    parser.add_argument("--density-ratio", type=float, default=1.0e6)
    parser.add_argument("--reynolds-number", "--Re", dest="reynolds_number", type=float, default=40.0)
    parser.add_argument("--u-max", type=float, default=0.04)
    parser.add_argument("--flow-control", choices=["fixed_pressure", "target_max_velocity"], default="target_max_velocity")
    parser.add_argument("--fluid-method", choices=["lbm-bgk-guo", "lbm-trt-guo"], default="lbm-bgk-guo")
    parser.add_argument("--couplings", nargs="+", choices=["point_force", "immersed_boundary", "solid_boundary"], default=["point_force", "immersed_boundary", "solid_boundary"])
    parser.add_argument(
        "--backend-config",
        dest="backend_configs",
        action="append",
        type=_parse_backend_config,
        default=None,
        help="Backend pair FLUID/COMPUTE. Can be repeated. Default benchmarks common pairs.",
    )
    parser.add_argument("--timed-steps", type=int, default=250)
    parser.add_argument("--warmup-steps", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--ibm-marker-spacing", type=float, default=1.0)
    parser.add_argument("--cylinder-x", type=float, default=36.5)
    parser.add_argument("--cylinder-y", type=float, default=20.0)
    parser.add_argument("--cylinder-radius", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=20260509)
    parser.add_argument("--result-tag", default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()
    if args.backend_configs is None:
        args.backend_configs = [
            _parse_backend_config(item)
            for item in ("numpy/numpy", "numba/numba", "auto/auto", "auto/numba", "numba/auto")
        ]
    if args.nx <= 2 or args.ny <= 2:
        parser.error("--nx and --ny must be greater than 2")
    if args.n_particles <= 0:
        parser.error("--n-particles must be positive")
    if args.density_ratio <= 0.0:
        parser.error("--density-ratio must be positive")
    if args.timed_steps <= 0:
        parser.error("--timed-steps must be positive")
    if args.warmup_steps < 0:
        parser.error("--warmup-steps must be non-negative")
    if args.repeats <= 0:
        parser.error("--repeats must be positive")
    if args.ibm_marker_spacing <= 0.0:
        parser.error("--ibm-marker-spacing must be positive")
    return args


def main() -> None:
    args = parse_args()
    timestamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    parts = [args.result_tag, timestamp] if args.result_tag else [timestamp]
    output_dir = args.out_dir or program_results_dir(__file__, *parts)

    # Compile/warm the common Numba signatures before timing repeats.
    for coupling in args.couplings:
        for _, fluid_accelerator, compute_accelerator in args.backend_configs:
            sim = _build_sim(args, coupling, fluid_accelerator, compute_accelerator, args.seed)
            sim.advance(max(args.warmup_steps, 1))

    rows: list[CouplingBackendRow] = []
    for coupling in args.couplings:
        for label, fluid_accelerator, compute_accelerator in args.backend_configs:
            for repeat in range(1, args.repeats + 1):
                row = _run_case(
                    args,
                    coupling,
                    label,
                    fluid_accelerator,
                    compute_accelerator,
                    repeat,
                )
                rows.append(row)
                print(
                    f"{coupling:18s} {label:12s} repeat {repeat}: "
                    f"{row.steps_per_second:.2f} steps/s "
                    f"(fluid={row.active_fluid_backend}, compute={row.active_compute_backend})"
                )

    _write_outputs(output_dir, rows, args)
    print(f"\nSaved benchmark outputs: {output_dir}")


if __name__ == "__main__":
    main()
