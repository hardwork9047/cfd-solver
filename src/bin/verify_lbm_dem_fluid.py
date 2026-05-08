#!/usr/bin/env python3
"""Fluid-only verification checks for the LBM-DEM solver.

The checks intentionally isolate the LBM part from DEM particles:

1. Plane Poiseuille profile against a parabolic reference.
2. Mass-flux consistency between two streamwise sections.
3. Target maximum-velocity control under a fixed cylinder geometry.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cfd.result_paths import program_results_dir
from cfd_dem_lbm import LBMDEMSolver


@dataclass
class VerificationResult:
    name: str
    passed: bool
    metric: float
    tolerance: float
    details: str


def _fluid_only(sim: LBMDEMSolver) -> LBMDEMSolver:
    """Remove initialized DEM particles so only the fluid solver advances."""
    if sim.n_p:
        sim._delete_particles(np.ones(sim.n_p, dtype=bool))
    sim.total_particles_requested = 0
    sim.generated_particles = 0
    return sim


def _positive_flux(ux: np.ndarray, solid: np.ndarray, x_idx: int) -> float:
    fluid = ~solid[x_idx, :]
    return float(np.sum(np.maximum(ux[x_idx, fluid], 0.0)))


def _poiseuille_reference(ny: int, u_max: float) -> tuple[np.ndarray, np.ndarray]:
    y = np.arange(ny, dtype=float)
    channel_height = max(ny - 2.0, 1.0)
    eta = np.clip((y - 0.5) / channel_height, 0.0, 1.0)
    profile = 4.0 * u_max * eta * (1.0 - eta)
    profile[0] = 0.0
    profile[-1] = 0.0
    return y, profile


def verify_plane_poiseuille(
    nx: int,
    ny: int,
    steps: int,
    re: float,
    u_max: float,
    tolerance: float,
    output_dir: Path,
) -> VerificationResult:
    sim = _fluid_only(
        LBMDEMSolver(
            nx=nx,
            ny=ny,
            Re=re,
            u_max=u_max,
            reynolds_length=float(ny),
            flow_control="fixed_pressure",
            n_particles=1,
            particle_radius=2.0,
            gravity=0.0,
            seed=11,
        )
    )
    sim.advance(steps)
    _, ux, _ = sim.get_fields()
    measured = np.mean(ux[2:-2, :], axis=0)
    y, reference = _poiseuille_reference(ny, u_max)
    interior = slice(1, ny - 1)
    denom = max(float(np.dot(reference[interior], reference[interior])), 1e-12)
    fitted_scale = float(np.dot(measured[interior], reference[interior]) / denom)
    fitted_reference = fitted_scale * reference
    rel_l2 = float(
        np.linalg.norm(measured[interior] - fitted_reference[interior])
        / max(float(np.linalg.norm(fitted_reference[interior])), 1e-12)
    )

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(fitted_reference, y, label="best-fit parabolic reference")
    ax.plot(measured, y, "--", label="LBM")
    ax.set_xlabel("ux")
    ax.set_ylabel("y")
    ax.set_title(f"Plane Poiseuille verification, rel L2={rel_l2:.3g}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "plane_poiseuille_profile.png", dpi=180)
    plt.close(fig)

    return VerificationResult(
        name="plane_poiseuille_profile",
        passed=rel_l2 <= tolerance,
        metric=rel_l2,
        tolerance=tolerance,
        details=(
            "Relative L2 shape error between cross-section averaged ux and "
            f"best-fit parabolic profile. fitted_scale={fitted_scale:.6g}."
        ),
    )


def verify_flux_balance(
    nx: int,
    ny: int,
    steps: int,
    re: float,
    u_max: float,
    tolerance: float,
) -> VerificationResult:
    sim = _fluid_only(
        LBMDEMSolver(
            nx=nx,
            ny=ny,
            Re=re,
            u_max=u_max,
            reynolds_length=float(ny),
            flow_control="fixed_pressure",
            n_particles=1,
            particle_radius=2.0,
            gravity=0.0,
            cylinders=[(0.45 * nx, 0.5 * ny, 0.11 * ny)],
            seed=12,
        )
    )
    sim.advance(steps)
    _, ux, _ = sim.get_fields()
    q_left = _positive_flux(ux, sim.solid, 1)
    q_right = _positive_flux(ux, sim.solid, nx - 2)
    imbalance = abs(q_left - q_right) / max(abs(q_left), abs(q_right), 1e-12)
    return VerificationResult(
        name="streamwise_flux_balance",
        passed=imbalance <= tolerance,
        metric=imbalance,
        tolerance=tolerance,
        details="Relative positive-flux imbalance between near-inlet and near-outlet sections.",
    )


def verify_target_max_velocity(
    nx: int,
    ny: int,
    steps: int,
    re: float,
    u_max: float,
    tolerance: float,
) -> VerificationResult:
    sim = _fluid_only(
        LBMDEMSolver(
            nx=nx,
            ny=ny,
            Re=re,
            u_max=u_max,
            reynolds_length=4.0,
            flow_control="target_max_velocity",
            flow_control_gain=0.2,
            n_particles=1,
            particle_radius=2.0,
            gravity=0.0,
            cylinders=[(0.45 * nx, 0.5 * ny, 0.10 * ny)],
            seed=13,
        )
    )
    sim.advance(steps)
    _, ux, uy = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2)
    observed = float(speed[~sim.solid].max())
    rel_error = abs(observed - u_max) / max(u_max, 1e-12)
    return VerificationResult(
        name="target_max_velocity_control",
        passed=rel_error <= tolerance,
        metric=rel_error,
        tolerance=tolerance,
        details=f"Observed max speed={observed:.6g}, target={u_max:.6g}.",
    )


def write_outputs(results: list[VerificationResult], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "fluid_verification_results.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(result) for result in results)

    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "passed": all(result.passed for result in results),
        "results": [asdict(result) for result in results],
    }
    (output_dir / "fluid_verification_results.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    lines = [
        "# LBM-DEM Fluid Verification",
        "",
        f"Generated: {report['created_at']}",
        f"Overall: {'PASS' if report['passed'] else 'FAIL'}",
        "",
    ]
    for result in results:
        status = "PASS" if result.passed else "FAIL"
        lines.append(
            f"- {status} `{result.name}`: metric={result.metric:.6g}, "
            f"tolerance={result.tolerance:.6g}. {result.details}"
        )
    (output_dir / "fluid_verification_report.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify LBM fluid accuracy used by LBM-DEM")
    parser.add_argument("--nx", type=int, default=72)
    parser.add_argument("--ny", type=int, default=32)
    parser.add_argument("--steps", type=int, default=3000)
    parser.add_argument("--reynolds-number", "--Re", dest="re", type=float, default=40.0)
    parser.add_argument("--u-max", type=float, default=0.04)
    parser.add_argument("--poiseuille-tol", type=float, default=0.20)
    parser.add_argument("--flux-tol", type=float, default=0.05)
    parser.add_argument("--control-tol", type=float, default=0.20)
    args = parser.parse_args()

    output_dir = program_results_dir(__file__, datetime.now().strftime("run_%Y%m%d_%H%M%S"))
    results = [
        verify_plane_poiseuille(
            args.nx,
            args.ny,
            args.steps,
            args.re,
            args.u_max,
            args.poiseuille_tol,
            output_dir,
        ),
        verify_flux_balance(args.nx, args.ny, args.steps, args.re, args.u_max, args.flux_tol),
        verify_target_max_velocity(args.nx, args.ny, args.steps, args.re, args.u_max, args.control_tol),
    ]
    write_outputs(results, output_dir)
    for result in results:
        status = "PASS" if result.passed else "FAIL"
        print(f"{status} {result.name}: metric={result.metric:.6g}, tol={result.tolerance:.6g}")
    print(f"Verification outputs: {output_dir}")
    return 0 if all(result.passed for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
