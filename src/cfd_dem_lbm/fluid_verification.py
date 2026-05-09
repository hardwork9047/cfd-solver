"""Fluid-only verification routines for the production LBM-DEM solver path."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .fast_solver import FastLBMDEM


@dataclass
class VerificationResult:
    name: str
    passed: bool
    metric: float
    tolerance: float
    details: str


def fluid_only(sim: FastLBMDEM) -> FastLBMDEM:
    """Remove initialized DEM particles so only the fluid solver advances."""
    if sim.n_p:
        sim._delete_particles(np.ones(sim.n_p, dtype=bool))
        sim._invalidate_macroscopic_cache()
    sim.total_particles_requested = 0
    sim.generated_particles = 0
    return sim


def positive_flux(ux: np.ndarray, solid: np.ndarray, x_idx: int) -> float:
    """Return positive x-direction flux through one vertical section."""
    fluid = ~solid[x_idx, :]
    return float(np.sum(np.maximum(ux[x_idx, fluid], 0.0)))


def poiseuille_reference(ny: int, u_max: float) -> tuple[np.ndarray, np.ndarray]:
    """Return the parabolic channel-flow shape used by the production setup."""
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
    """Check that the no-cylinder fluid path recovers a parabolic profile."""
    sim = fluid_only(
        FastLBMDEM(
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
    y, reference = poiseuille_reference(ny, u_max)
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
    """Check streamwise flux consistency through a cylinder-containing domain."""
    sim = fluid_only(
        FastLBMDEM(
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
    q_left = positive_flux(ux, sim.solid, 1)
    q_right = positive_flux(ux, sim.solid, nx - 2)
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
    """Check the production max-speed control loop under a cylinder geometry."""
    sim = fluid_only(
        FastLBMDEM(
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


def verify_particle_solid_boundary_reduces_flux(
    nx: int,
    ny: int,
    steps: int,
    re: float,
    u_max: float,
    min_reduction: float,
) -> VerificationResult:
    """Check that dynamic particle solids reduce through-channel flux."""
    common = dict(
        nx=nx,
        ny=ny,
        Re=re,
        u_max=u_max,
        reynolds_length=float(ny),
        flow_control="fixed_pressure",
        n_particles=1,
        particle_radius=max(2.0, 0.18 * ny),
        gravity=0.0,
        seed=14,
    )
    open_sim = fluid_only(FastLBMDEM(**common))
    blocked_sim = FastLBMDEM(
        **common,
        particle_fluid_coupling="solid_boundary",
    )
    blocked_sim.pos[:] = np.array([[0.5 * nx, 0.5 * ny]])
    blocked_sim.vel[:] = 0.0
    blocked_sim._update_particle_solid_mask()

    open_sim.advance(steps)
    blocked_sim.advance(steps)
    _, ux_open, _ = open_sim.get_fields()
    _, ux_blocked, _ = blocked_sim.get_fields()
    section = nx // 2
    q_open = positive_flux(ux_open, open_sim.solid, section)
    q_blocked = positive_flux(ux_blocked, blocked_sim.solid, section)
    reduction = 1.0 - q_blocked / max(q_open, 1e-12)
    return VerificationResult(
        name="particle_solid_boundary_flux_reduction",
        passed=reduction >= min_reduction,
        metric=reduction,
        tolerance=min_reduction,
        details=(
            "Flux reduction from overlaying one DEM particle as a moving LBM solid. "
            f"section={section}, open_flux={q_open:.6g}, blocked_flux={q_blocked:.6g}."
        ),
    )


def verify_immersed_boundary_particle_load(
    nx: int,
    ny: int,
    steps: int,
    re: float,
    u_max: float,
    min_force: float,
) -> VerificationResult:
    """Check that IBM produces a hydrodynamic load on a fixed particle."""
    sim = FastLBMDEM(
        nx=nx,
        ny=ny,
        Re=re,
        u_max=u_max,
        reynolds_length=float(ny),
        flow_control="fixed_pressure",
        n_particles=1,
        particle_radius=max(2.0, 0.14 * ny),
        gravity=0.0,
        particle_fluid_coupling="immersed_boundary",
        ibm_stiffness=1.0,
        ibm_marker_spacing=1.0,
        seed=15,
    )
    sim.pos[:] = np.array([[0.5 * nx, 0.5 * ny]])
    sim.vel[:] = 0.0
    sim.omega_p[:] = 0.0
    sim.advance(steps)
    force_mag = float(np.linalg.norm(sim.ibm_forces_p[0])) if sim.n_p else 0.0
    return VerificationResult(
        name="immersed_boundary_particle_load",
        passed=force_mag >= min_force,
        metric=force_mag,
        tolerance=min_force,
        details=(
            "Magnitude of the particle reaction force from direct-forcing IBM "
            "for a stationary particle in pressure-driven channel flow."
        ),
    )


def run_fluid_verification(
    *,
    nx: int,
    ny: int,
    steps: int,
    re: float,
    u_max: float,
    poiseuille_tol: float,
    flux_tol: float,
    control_tol: float,
    particle_solid_min_reduction: float,
    ibm_min_force: float,
    output_dir: Path,
) -> list[VerificationResult]:
    """Run the standard fluid-only verification suite."""
    output_dir.mkdir(parents=True, exist_ok=True)
    return [
        verify_plane_poiseuille(nx, ny, steps, re, u_max, poiseuille_tol, output_dir),
        verify_flux_balance(nx, ny, steps, re, u_max, flux_tol),
        verify_target_max_velocity(nx, ny, steps, re, u_max, control_tol),
        verify_particle_solid_boundary_reduces_flux(
            nx,
            ny,
            max(steps // 3, 1),
            re,
            u_max,
            particle_solid_min_reduction,
        ),
        verify_immersed_boundary_particle_load(
            nx,
            ny,
            max(steps // 3, 1),
            re,
            u_max,
            ibm_min_force,
        ),
    ]


def write_verification_outputs(results: list[VerificationResult], output_dir: Path) -> None:
    """Write CSV, JSON, and Markdown verification summaries."""
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "fluid_verification_results.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(results[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(result) for result in results)

    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "passed": all(result.passed for result in results),
        "solver_path": "cfd_dem_lbm.FastLBMDEM",
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
        "Solver path: `cfd_dem_lbm.FastLBMDEM`",
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
