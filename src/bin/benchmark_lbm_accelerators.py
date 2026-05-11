#!/usr/bin/env python3
"""Benchmark NumPy and Numba LBM backends on the shared LBM-DEM solver path."""

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

from particulate_flow.result_paths import program_results_dir
from particulate_flow import FastLBMDEM


@dataclass
class BenchmarkRow:
    accelerator: str
    repeat: int
    active_backend: str
    compile_warmup_seconds: float
    timed_seconds: float
    timed_steps: int
    steps_per_second: float
    final_max_speed: float
    final_mean_density: float


def _fluid_only(sim: FastLBMDEM) -> FastLBMDEM:
    if sim.n_p:
        sim._delete_particles(np.ones(sim.n_p, dtype=bool))
        sim._invalidate_macroscopic_cache()
    sim.total_particles_requested = 0
    sim.generated_particles = 0
    return sim


def _build_sim(args: argparse.Namespace, accelerator: str) -> FastLBMDEM:
    cylinder = (args.cylinder_x, args.cylinder_y, args.cylinder_radius)
    return _fluid_only(
        FastLBMDEM(
            nx=args.nx,
            ny=args.ny,
            Re=args.reynolds_number,
            u_max=args.u_max,
            reynolds_length=2.0 * args.cylinder_radius,
            flow_control="fixed_pressure",
            fluid_method=args.fluid_method,
            fluid_accelerator=accelerator,
            n_particles=1,
            particle_radius=1.0,
            gravity=0.0,
            cylinders=[cylinder],
            seed=100 + args.nx + args.ny,
        )
    )


def _run_one(args: argparse.Namespace, accelerator: str, repeat: int) -> BenchmarkRow:
    sim = _build_sim(args, accelerator)

    t0 = time.perf_counter()
    if args.compile_warmup_steps:
        sim.advance(args.compile_warmup_steps)
    compile_warmup_seconds = time.perf_counter() - t0

    t1 = time.perf_counter()
    sim.advance(args.timed_steps)
    timed_seconds = time.perf_counter() - t1

    rho, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)
    return BenchmarkRow(
        accelerator=accelerator,
        repeat=repeat,
        active_backend="numba" if sim.uses_numba_lbm else "numpy",
        compile_warmup_seconds=compile_warmup_seconds,
        timed_seconds=timed_seconds,
        timed_steps=args.timed_steps,
        steps_per_second=args.timed_steps / max(timed_seconds, 1e-12),
        final_max_speed=float(speed[~sim.solid].max()),
        final_mean_density=float(np.mean(rho[~sim.solid])),
    )


def _write_outputs(output_dir: Path, rows: list[BenchmarkRow], args: argparse.Namespace) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = list(asdict(rows[0]).keys())
    csv_path = output_dir / "lbm_accelerator_benchmark.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)

    grouped: dict[str, list[BenchmarkRow]] = {}
    for row in rows:
        grouped.setdefault(row.accelerator, []).append(row)

    summary = {}
    for name, group in grouped.items():
        speeds = [row.steps_per_second for row in group]
        times = [row.timed_seconds for row in group]
        summary[name] = {
            "active_backend": group[0].active_backend,
            "mean_steps_per_second": statistics.fmean(speeds),
            "median_steps_per_second": statistics.median(speeds),
            "mean_timed_seconds": statistics.fmean(times),
            "median_timed_seconds": statistics.median(times),
            "mean_compile_warmup_seconds": statistics.fmean(
                row.compile_warmup_seconds for row in group
            ),
        }
    numpy_speed = summary.get("numpy", {}).get("median_steps_per_second")
    numba_speed = summary.get("numba", {}).get("median_steps_per_second")
    if numpy_speed and numba_speed:
        summary["numba_vs_numpy_speedup"] = numba_speed / numpy_speed

    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "command": sys.argv,
        "configuration": vars(args),
        "summary": summary,
    }
    (output_dir / "lbm_accelerator_benchmark.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    lines = [
        "# LBM Accelerator Benchmark",
        "",
        f"Generated: {report['created_at']}",
        f"Grid: {args.nx} x {args.ny}, timed steps: {args.timed_steps}",
        "",
        "| accelerator | active backend | median steps/s | median timed s | warmup s |",
        "|---|---:|---:|---:|---:|",
    ]
    for name in ("numpy", "numba"):
        if name not in summary:
            continue
        item = summary[name]
        lines.append(
            f"| {name} | {item['active_backend']} | "
            f"{item['median_steps_per_second']:.3f} | "
            f"{item['median_timed_seconds']:.6f} | "
            f"{item['mean_compile_warmup_seconds']:.6f} |"
        )
    if "numba_vs_numpy_speedup" in summary:
        lines.extend(
            [
                "",
                f"Numba/NumPy median speed ratio: {summary['numba_vs_numpy_speedup']:.3f}x",
            ]
        )
    (output_dir / "lbm_accelerator_benchmark.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=160)
    parser.add_argument("--ny", type=int, default=64)
    parser.add_argument("--reynolds-number", "--Re", dest="reynolds_number", type=float, default=50.0)
    parser.add_argument("--u-max", type=float, default=0.04)
    parser.add_argument("--fluid-method", choices=["lbm-bgk-guo", "lbm-trt-guo"], default="lbm-bgk-guo")
    parser.add_argument("--timed-steps", type=int, default=200)
    parser.add_argument("--compile-warmup-steps", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cylinder-x", type=float, default=56.0)
    parser.add_argument("--cylinder-y", type=float, default=32.0)
    parser.add_argument("--cylinder-radius", type=float, default=6.0)
    parser.add_argument("--result-tag", default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()
    if args.nx <= 2 or args.ny <= 2:
        parser.error("--nx and --ny must be greater than 2")
    if args.timed_steps <= 0:
        parser.error("--timed-steps must be positive")
    if args.compile_warmup_steps < 0:
        parser.error("--compile-warmup-steps must be non-negative")
    if args.repeats <= 0:
        parser.error("--repeats must be positive")
    return args


def main() -> None:
    args = parse_args()
    timestamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    parts = [args.result_tag, timestamp] if args.result_tag else [timestamp]
    output_dir = args.out_dir or program_results_dir(__file__, *parts)

    rows: list[BenchmarkRow] = []
    for accelerator in ("numpy", "numba"):
        for repeat in range(1, args.repeats + 1):
            row = _run_one(args, accelerator, repeat)
            rows.append(row)
            print(
                f"{accelerator:>5} repeat {repeat}: "
                f"backend={row.active_backend}, {row.steps_per_second:.2f} steps/s, "
                f"timed={row.timed_seconds:.3f}s, warmup={row.compile_warmup_seconds:.3f}s"
            )

    _write_outputs(output_dir, rows, args)
    print(f"\nSaved benchmark outputs: {output_dir}")


if __name__ == "__main__":
    main()
