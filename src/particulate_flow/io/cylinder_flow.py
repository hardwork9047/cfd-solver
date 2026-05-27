"""Cylinder-flow runner using the production LBM solver path.

This module intentionally builds a fluid-only cylinder calculation from
``FastLBMDEM``.  That keeps standalone cylinder-flow runs on the same LBM
collision, forcing, solid-boundary, and control implementation used by
``run_lbm_dem.py`` and the fluid-verification suite.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..fast_solver import FastLBMDEM
from ..lbm_dem import COMPUTE_ACCELERATORS, CS2, FLUID_ACCELERATORS, FLUID_METHODS
from .paths import program_results_dir


@dataclass
class CylinderFlowMetrics:
    step: int
    max_speed: float
    mean_inlet_flux: float
    mean_outlet_flux: float
    relative_flux_imbalance: float
    pressure_drop: float
    reynolds_number: float


def remove_dem_particles(sim: FastLBMDEM) -> FastLBMDEM:
    """Convert a coupled solver instance into a fluid-only calculation."""
    if sim.n_p:
        sim._delete_particles(np.ones(sim.n_p, dtype=bool))
        sim._invalidate_macroscopic_cache()
    sim.total_particles_requested = 0
    sim.generated_particles = 0
    return sim


def default_cylinder(nx: int, ny: int) -> tuple[float, float, float]:
    """Return a conservative single-cylinder default in lattice units."""
    return (0.35 * float(nx), 0.5 * float(ny), 0.10 * float(ny))


def build_cylinder_flow_sim(
    *,
    nx: int,
    ny: int,
    reynolds_number: float,
    u_max: float,
    cylinders: list[tuple[float, float, float]],
    reynolds_length: float | None = None,
    flow_condition: str = "fixed_pressure",
    flow_control_gain: float = 0.2,
    fluid_method: str = "lbm-bgk-guo",
    fluid_accelerator: str = "numpy",
    compute_accelerator: str = "auto",
) -> FastLBMDEM:
    """Create a fluid-only cylinder-flow simulation with the shared solver."""
    if not cylinders:
        cylinders = [default_cylinder(nx, ny)]
    if reynolds_length is None:
        reynolds_length = 2.0 * max(radius for _, _, radius in cylinders)

    sim = FastLBMDEM(
        nx=nx,
        ny=ny,
        Re=reynolds_number,
        u_max=u_max,
        reynolds_length=reynolds_length,
        flow_control=flow_condition,
        flow_control_gain=flow_control_gain,
        fluid_method=fluid_method,
        fluid_accelerator=fluid_accelerator,
        compute_accelerator=compute_accelerator,
        n_particles=1,
        particle_radius=1.0,
        gravity=0.0,
        cylinders=cylinders,
        seed=123,
    )
    return remove_dem_particles(sim)


def _section_flux(ux: np.ndarray, solid: np.ndarray, x_idx: int) -> float:
    fluid = ~solid[x_idx, :]
    return float(np.sum(ux[x_idx, fluid]))


def compute_cylinder_flow_metrics(sim: FastLBMDEM) -> CylinderFlowMetrics:
    """Compute scalar diagnostics for a cylinder-flow state."""
    rho, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)
    pressure = CS2 * rho
    q_in = _section_flux(ux, sim.solid, 1)
    q_out = _section_flux(ux, sim.solid, sim.nx - 2)
    imbalance = abs(q_in - q_out) / max(abs(q_in), abs(q_out), 1e-12)
    inlet_fluid = ~sim.solid[1, :]
    outlet_fluid = ~sim.solid[sim.nx - 2, :]
    pressure_drop = float(np.mean(pressure[1, inlet_fluid]) - np.mean(pressure[sim.nx - 2, outlet_fluid]))
    max_speed = float(speed[~sim.solid].max()) if np.any(~sim.solid) else 0.0
    re = max_speed * sim.reynolds_length / max(sim.nu, 1e-12)
    return CylinderFlowMetrics(
        step=sim.step_count,
        max_speed=max_speed,
        mean_inlet_flux=q_in,
        mean_outlet_flux=q_out,
        relative_flux_imbalance=imbalance,
        pressure_drop=pressure_drop,
        reynolds_number=float(re),
    )


def write_structured_vtk(path: Path, sim: FastLBMDEM) -> None:
    """Write fluid fields as legacy VTK structured-points data."""
    rho, ux, uy = sim.get_fields()
    pressure = CS2 * rho
    speed = np.sqrt(ux**2 + uy**2)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write("particulate_flow cylinder flow\n")
        handle.write("ASCII\n")
        handle.write("DATASET STRUCTURED_POINTS\n")
        handle.write(f"DIMENSIONS {sim.nx} {sim.ny} 1\n")
        handle.write("ORIGIN 0 0 0\n")
        handle.write("SPACING 1 1 1\n")
        handle.write(f"POINT_DATA {sim.nx * sim.ny}\n")

        for name, data in (
            ("pressure", pressure),
            ("speed", speed),
            ("density", rho),
            ("solid", sim.solid.astype(float)),
        ):
            handle.write(f"SCALARS {name} float 1\n")
            handle.write("LOOKUP_TABLE default\n")
            for y in range(sim.ny):
                for x in range(sim.nx):
                    handle.write(f"{float(data[x, y]):.9e}\n")

        handle.write("VECTORS velocity float\n")
        for y in range(sim.ny):
            for x in range(sim.nx):
                handle.write(f"{float(ux[x, y]):.9e} {float(uy[x, y]):.9e} 0.0\n")


def write_plots(output_dir: Path, sim: FastLBMDEM) -> None:
    """Save velocity and pressure figures for quick inspection."""
    rho, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)
    pressure = CS2 * rho

    for name, field, cmap, label in (
        ("speed", speed, "viridis", "|u|"),
        ("pressure", pressure, "coolwarm", "pressure"),
    ):
        fig, ax = plt.subplots(figsize=(8, 3.5))
        image = ax.imshow(field.T, origin="lower", cmap=cmap, aspect="auto")
        ax.contour(sim.solid.T, levels=[0.5], colors="black", linewidths=0.8)
        ax.set_xlabel("x [lattice]")
        ax.set_ylabel("y [lattice]")
        ax.set_title(f"Cylinder flow {label}, step={sim.step_count}")
        fig.colorbar(image, ax=ax, label=label)
        fig.tight_layout()
        fig.savefig(output_dir / f"{name}.png", dpi=180)
        plt.close(fig)


def run_cylinder_flow(
    *,
    output_dir: Path,
    nx: int,
    ny: int,
    reynolds_number: float,
    u_max: float,
    cylinders: list[tuple[float, float, float]],
    reynolds_length: float | None,
    flow_condition: str,
    flow_control_gain: float,
    fluid_method: str,
    fluid_accelerator: str,
    compute_accelerator: str = "auto",
    steps: int,
    report_every: int,
) -> tuple[FastLBMDEM, list[CylinderFlowMetrics]]:
    """Run a cylinder-flow calculation and write standard outputs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    sim = build_cylinder_flow_sim(
        nx=nx,
        ny=ny,
        reynolds_number=reynolds_number,
        u_max=u_max,
        cylinders=cylinders,
        reynolds_length=reynolds_length,
        flow_condition=flow_condition,
        flow_control_gain=flow_control_gain,
        fluid_method=fluid_method,
        fluid_accelerator=fluid_accelerator,
        compute_accelerator=compute_accelerator,
    )

    metrics: list[CylinderFlowMetrics] = [compute_cylinder_flow_metrics(sim)]
    report = max(report_every, 1)
    while sim.step_count < steps:
        chunk = min(report, steps - sim.step_count)
        sim.advance(chunk)
        metrics.append(compute_cylinder_flow_metrics(sim))
        latest = metrics[-1]
        print(
            f"step {latest.step:>8}  max|u|={latest.max_speed:.6g}  "
            f"Re={latest.reynolds_number:.4g}  flux_imb={latest.relative_flux_imbalance:.3g}"
        )

    config = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "solver_path": "particulate_flow.FastLBMDEM",
        "calculation": "fluid_only_cylinder_flow",
        "nx": nx,
        "ny": ny,
        "reynolds_number_input": reynolds_number,
        "u_max": u_max,
        "reynolds_length": sim.reynolds_length,
        "flow_condition": flow_condition,
        "flow_control_gain": flow_control_gain,
        "fluid_method": fluid_method,
        "fluid_accelerator": fluid_accelerator,
        "uses_numba_lbm": sim.uses_numba_lbm,
        "compute_accelerator": compute_accelerator,
        "uses_numba_compute": sim.uses_numba_compute,
        "cylinders": [{"x": x, "y": y, "radius": r} for x, y, r in sim.cylinders],
        "steps": steps,
    }
    (output_dir / "metadata.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    with (output_dir / "time_series.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(metrics[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in metrics)

    rho, ux, uy = sim.get_fields()
    np.savez_compressed(
        output_dir / "final_fields.npz",
        rho=rho,
        ux=ux,
        uy=uy,
        solid=sim.solid,
        pressure=CS2 * rho,
        cylinders=np.array(sim.cylinders, dtype=float),
    )
    write_structured_vtk(output_dir / "paraview" / "final_fields.vtk", sim)
    write_plots(output_dir, sim)
    return sim, metrics


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fluid-only cylinder-flow calculation using particulate_flow.FastLBMDEM"
    )
    parser.add_argument("--nx", type=int, default=180, help="Grid width in lattice nodes")
    parser.add_argument("--ny", type=int, default=70, help="Grid height in lattice nodes")
    parser.add_argument("--Re", "--reynolds-number", dest="reynolds_number", type=float, default=100.0)
    parser.add_argument("--u-max", type=float, default=0.05)
    parser.add_argument(
        "--reynolds-length",
        type=float,
        default=None,
        help="Reference length for Re. Default: diameter of largest cylinder.",
    )
    parser.add_argument(
        "--flow-condition",
        choices=["fixed_pressure", "target_max_velocity", "fixed-pressure", "target-max-velocity"],
        default="fixed_pressure",
    )
    parser.add_argument("--flow-control-gain", type=float, default=0.2)
    parser.add_argument("--fluid-method", choices=FLUID_METHODS, default="lbm-bgk-guo")
    parser.add_argument("--fluid-accelerator", choices=FLUID_ACCELERATORS, default="numpy")
    parser.add_argument("--compute-accelerator", choices=COMPUTE_ACCELERATORS, default="auto")
    parser.add_argument(
        "--cylinder",
        "--cylinder-spec",
        dest="cylinders",
        type=float,
        nargs=3,
        action="append",
        metavar=("X", "Y", "R"),
        default=None,
        help="Add fixed cylinder as x y radius. Can be repeated.",
    )
    parser.add_argument("--steps", type=int, default=10000)
    parser.add_argument("--report-every", type=int, default=1000)
    parser.add_argument("--result-tag", default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args(argv)
    if args.nx <= 2 or args.ny <= 2:
        parser.error("--nx and --ny must be greater than 2")
    if args.reynolds_number <= 0.0:
        parser.error("--Re must be positive")
    if args.u_max < 0.0:
        parser.error("--u-max must be non-negative")
    if args.reynolds_length is not None and args.reynolds_length <= 0.0:
        parser.error("--reynolds-length must be positive")
    if args.steps < 0:
        parser.error("--steps must be non-negative")
    if args.report_every <= 0:
        parser.error("--report-every must be positive")
    args.flow_condition = args.flow_condition.replace("-", "_")
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    timestamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    parts = [args.result_tag, timestamp] if args.result_tag else [timestamp]
    output_dir = args.out_dir or program_results_dir(__file__, *parts)
    cylinders = args.cylinders if args.cylinders is not None else []
    sim, metrics = run_cylinder_flow(
        output_dir=output_dir,
        nx=args.nx,
        ny=args.ny,
        reynolds_number=args.reynolds_number,
        u_max=args.u_max,
        cylinders=[tuple(item) for item in cylinders],
        reynolds_length=args.reynolds_length,
        flow_condition=args.flow_condition,
        flow_control_gain=args.flow_control_gain,
        fluid_method=args.fluid_method,
        fluid_accelerator=args.fluid_accelerator,
        compute_accelerator=args.compute_accelerator,
        steps=args.steps,
        report_every=args.report_every,
    )
    latest = metrics[-1]
    print("")
    print(f"Saved: {output_dir}")
    print(
        "Solver path: particulate_flow.FastLBMDEM, "
        f"fluid_method={sim.fluid_method}, accelerator={sim.fluid_accelerator}, "
        f"compute={sim.compute_accelerator}"
    )
    print(
        f"Final: step={latest.step}, max|u|={latest.max_speed:.6g}, "
        f"Re={latest.reynolds_number:.6g}, pressure_drop={latest.pressure_drop:.6g}"
    )


if __name__ == "__main__":
    main()
