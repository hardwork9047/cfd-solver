#!/usr/bin/env python3
"""Benchmark IBM marker-spacing convergence for a fixed particle in channel flow."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from particulate_flow.io.paths import program_results_dir
from particulate_flow import FastLBMDEM


@dataclass
class MarkerSpacingRow:
    marker_spacing: float
    marker_count: int
    markers_per_particle: float
    seconds: float
    steps_per_second: float
    force_x: float
    force_y: float
    force_magnitude: float
    mean_slip_speed: float
    max_slip_speed: float
    inlet_flux: float
    outlet_flux: float
    flux_imbalance: float
    pressure_drop: float
    relative_force_change_vs_finest: float = 0.0
    relative_pressure_drop_change_vs_finest: float = 0.0
    relative_outlet_flux_change_vs_finest: float = 0.0


def _positive_flux(ux: np.ndarray, solid: np.ndarray, x_idx: int) -> float:
    fluid = ~solid[x_idx, :]
    return float(np.sum(np.maximum(ux[x_idx, fluid], 0.0)))


def _boundary_mean(field: np.ndarray, solid: np.ndarray, x_idx: int) -> float:
    fluid = ~solid[x_idx, :]
    if not np.any(fluid):
        return float("nan")
    return float(np.mean(field[x_idx, fluid]))


def _marker_slip(sim: FastLBMDEM, ux: np.ndarray, uy: np.ndarray) -> tuple[float, float, int]:
    markers = sim.ibm_marker_points()
    marker_x = markers["x"]
    marker_y = markers["y"]
    if len(marker_x) == 0:
        return 0.0, 0.0, 0
    marker_owner = markers["owner"]
    marker_rx = markers["rx"]
    marker_ry = markers["ry"]
    uf_x, uf_y = sim._interp_velocity_many(marker_x, marker_y, ux=ux, uy=uy)
    owner_vel = sim.vel[marker_owner]
    owner_omega = sim.omega_p[marker_owner]
    ub_x = owner_vel[:, 0] - owner_omega * marker_ry
    ub_y = owner_vel[:, 1] + owner_omega * marker_rx
    slip = np.hypot(ub_x - uf_x, ub_y - uf_y)
    return float(np.mean(slip)), float(np.max(slip)), int(len(slip))


def _build_sim(args: argparse.Namespace, marker_spacing: float) -> FastLBMDEM:
    sim = FastLBMDEM(
        nx=args.nx,
        ny=args.ny,
        Re=args.reynolds_number,
        u_max=args.u_max,
        reynolds_length=float(args.ny),
        flow_control=args.flow_control,
        flow_control_gain=args.flow_control_gain,
        fluid_method=args.fluid_method,
        fluid_accelerator=args.fluid_accelerator,
        compute_accelerator=args.compute_accelerator,
        n_particles=1,
        particle_radius=args.particle_radius,
        radius_variation=0.0,
        density_ratio=args.density_ratio,
        gravity=0.0,
        particle_fluid_coupling="immersed_boundary",
        ibm_stiffness=args.ibm_stiffness,
        ibm_marker_spacing=marker_spacing,
        seed=22,
    )
    sim.pos[:] = np.array([[args.particle_x * args.nx, args.particle_y * args.ny]])
    sim.vel[:] = 0.0
    sim.omega_p[:] = 0.0
    return sim


def _run_one(args: argparse.Namespace, marker_spacing: float) -> MarkerSpacingRow:
    sim = _build_sim(args, marker_spacing)
    start = time.perf_counter()
    for _ in range(args.steps):
        sim.advance(1)
        if args.pin_particle:
            sim.pos[:] = np.array([[args.particle_x * args.nx, args.particle_y * args.ny]])
            sim.vel[:] = 0.0
            sim.omega_p[:] = 0.0
    seconds = time.perf_counter() - start
    rho, ux, uy = sim.get_fields()
    pressure = rho / 3.0
    mean_slip, max_slip, marker_count = _marker_slip(sim, ux, uy)
    inlet_flux = _positive_flux(ux, sim.solid, 1)
    outlet_flux = _positive_flux(ux, sim.solid, sim.nx - 2)
    pressure_drop = _boundary_mean(pressure, sim.solid, 1) - _boundary_mean(
        pressure, sim.solid, sim.nx - 2
    )
    flux_imbalance = abs(inlet_flux - outlet_flux) / max(abs(inlet_flux), abs(outlet_flux), 1e-12)
    force = sim.ibm_forces_p[0] if sim.n_p else np.zeros(2)
    return MarkerSpacingRow(
        marker_spacing=marker_spacing,
        marker_count=marker_count,
        markers_per_particle=float(marker_count / max(sim.n_p, 1)),
        seconds=seconds,
        steps_per_second=args.steps / max(seconds, 1e-12),
        force_x=float(force[0]),
        force_y=float(force[1]),
        force_magnitude=float(np.linalg.norm(force)),
        mean_slip_speed=mean_slip,
        max_slip_speed=max_slip,
        inlet_flux=inlet_flux,
        outlet_flux=outlet_flux,
        flux_imbalance=flux_imbalance,
        pressure_drop=pressure_drop,
    )


def _with_relative_changes(rows: list[MarkerSpacingRow]) -> list[MarkerSpacingRow]:
    if not rows:
        return rows
    finest = min(rows, key=lambda row: row.marker_spacing)
    for row in rows:
        row.relative_force_change_vs_finest = abs(row.force_magnitude - finest.force_magnitude) / max(
            abs(finest.force_magnitude), 1e-12
        )
        row.relative_pressure_drop_change_vs_finest = abs(row.pressure_drop - finest.pressure_drop) / max(
            abs(finest.pressure_drop), 1e-12
        )
        row.relative_outlet_flux_change_vs_finest = abs(row.outlet_flux - finest.outlet_flux) / max(
            abs(finest.outlet_flux), 1e-12
        )
    return rows


def _recommend_spacing(rows: list[MarkerSpacingRow], tolerance: float) -> MarkerSpacingRow:
    candidates = [
        row
        for row in sorted(rows, key=lambda item: item.marker_spacing, reverse=True)
        if row.relative_force_change_vs_finest <= tolerance
        and row.relative_pressure_drop_change_vs_finest <= tolerance
        and row.relative_outlet_flux_change_vs_finest <= tolerance
    ]
    return candidates[0] if candidates else min(rows, key=lambda row: row.marker_spacing)


def _write_outputs(output_dir: Path, rows: list[MarkerSpacingRow], args: argparse.Namespace) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "ibm_marker_spacing_convergence.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)

    recommended = _recommend_spacing(rows, args.tolerance)
    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "command": sys.argv,
        "configuration": vars(args),
        "recommended": asdict(recommended),
        "rows": [asdict(row) for row in rows],
    }
    (output_dir / "ibm_marker_spacing_convergence.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    lines = [
        "# IBM Marker Spacing Convergence",
        "",
        f"Generated: {report['created_at']}",
        f"Grid: {args.nx} x {args.ny}, particle radius: {args.particle_radius}",
        f"Reference: finest spacing = {min(row.marker_spacing for row in rows):.3g}",
        f"Recommended spacing at tolerance {args.tolerance:.1%}: {recommended.marker_spacing:.3g}",
        "",
        "| spacing | markers | force | rel force | pressure drop | rel dp | outlet flux | rel flux | mean slip | steps/s |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(rows, key=lambda item: item.marker_spacing, reverse=True):
        lines.append(
            f"| {row.marker_spacing:.3g} | {row.marker_count:d} | "
            f"{row.force_magnitude:.6g} | {row.relative_force_change_vs_finest:.3%} | "
            f"{row.pressure_drop:.6g} | {row.relative_pressure_drop_change_vs_finest:.3%} | "
            f"{row.outlet_flux:.6g} | {row.relative_outlet_flux_change_vs_finest:.3%} | "
            f"{row.mean_slip_speed:.6g} | {row.steps_per_second:.1f} |"
        )
    (output_dir / "ibm_marker_spacing_convergence.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=72)
    parser.add_argument("--ny", type=int, default=32)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--reynolds-number", type=float, default=40.0)
    parser.add_argument("--u-max", type=float, default=0.04)
    parser.add_argument("--particle-radius", type=float, default=3.0)
    parser.add_argument("--particle-x", type=float, default=0.45, help="Particle x as fraction of nx")
    parser.add_argument("--particle-y", type=float, default=0.50, help="Particle y as fraction of ny")
    parser.add_argument("--density-ratio", type=float, default=1e9)
    parser.add_argument("--ibm-stiffness", type=float, default=1.0)
    parser.add_argument("--flow-control", choices=["fixed_pressure", "target_max_velocity"], default="fixed_pressure")
    parser.add_argument("--flow-control-gain", type=float, default=0.2)
    parser.add_argument("--fluid-method", choices=["lbm-bgk-guo", "lbm-trt-guo"], default="lbm-bgk-guo")
    parser.add_argument("--fluid-accelerator", choices=["numpy", "numba", "auto"], default="numpy")
    parser.add_argument("--compute-accelerator", choices=["numpy", "numba", "auto"], default="auto")
    parser.add_argument("--pin-particle", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--tolerance", type=float, default=0.05)
    parser.add_argument(
        "--spacing",
        type=float,
        action="append",
        default=None,
        help="IBM marker spacing to test. Can be repeated.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    spacings = args.spacing or [2.0, 1.5, 1.0, 0.75, 0.5]
    rows = [_run_one(args, float(spacing)) for spacing in spacings]
    rows = _with_relative_changes(rows)
    tag = "spacing_" + "_".join(str(spacing).replace(".", "p") for spacing in spacings)
    output_dir = program_results_dir(__file__, tag, datetime.now().strftime("run_%Y%m%d_%H%M%S"))
    _write_outputs(output_dir, rows, args)
    recommended = _recommend_spacing(rows, args.tolerance)
    print(f"IBM marker spacing convergence outputs: {output_dir}")
    print(
        f"Recommended spacing: {recommended.marker_spacing:.3g} "
        f"({recommended.marker_count} markers, tolerance={args.tolerance:.1%})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
