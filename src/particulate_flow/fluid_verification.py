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
            reynolds_length=float(ny),
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


def verify_near_clogged_solid_boundary_suppresses_flux(
    nx: int,
    ny: int,
    steps: int,
    re: float,
    u_max: float,
    max_leakage_ratio: float,
) -> VerificationResult:
    """Check flux leakage when a solid-boundary particle spans the channel height."""
    common = dict(
        nx=nx,
        ny=ny,
        Re=re,
        u_max=u_max,
        reynolds_length=float(ny),
        flow_control="fixed_pressure",
        n_particles=1,
        particle_radius=2.0,
        density_ratio=1e9,
        gravity=0.0,
        seed=18,
    )
    open_sim = fluid_only(FastLBMDEM(**common))
    clogged_sim = FastLBMDEM(
        **common,
        particle_fluid_coupling="solid_boundary",
    )
    radius = 0.65 * ny
    clogged_sim.radii[:] = radius
    clogged_sim.masses[:] = 2.0 * np.pi * radius**2 * clogged_sim.density_ratio
    clogged_sim.inertias[:] = 0.5 * clogged_sim.masses * radius**2
    clogged_sim.pos[:] = np.array([[0.5 * nx, 0.5 * ny]])
    clogged_sim.vel[:] = 0.0
    clogged_sim.omega_p[:] = 0.0
    clogged_sim._update_particle_solid_mask()

    open_sim.advance(steps)
    clogged_sim.advance(steps)
    rho_open, ux_open, _ = open_sim.get_fields()
    rho_clogged, ux_clogged, _ = clogged_sim.get_fields()
    outlet = nx - 2
    inlet = 1
    q_open = positive_flux(ux_open, open_sim.solid, outlet)
    q_clogged = positive_flux(ux_clogged, clogged_sim.solid, outlet)
    leakage_ratio = q_clogged / max(q_open, 1e-12)
    q_in_clogged = positive_flux(ux_clogged, clogged_sim.solid, inlet)
    flux_imbalance = abs(q_in_clogged - q_clogged) / max(q_in_clogged, q_clogged, 1e-12)
    p_open = rho_open / 3.0
    p_clogged = rho_clogged / 3.0
    open_dp = float(np.mean(p_open[inlet, ~open_sim.solid[inlet, :]]) - np.mean(p_open[outlet, ~open_sim.solid[outlet, :]]))
    clogged_dp = float(
        np.mean(p_clogged[inlet, ~clogged_sim.solid[inlet, :]])
        - np.mean(p_clogged[outlet, ~clogged_sim.solid[outlet, :]])
    )
    return VerificationResult(
        name="near_clogged_solid_boundary_leakage",
        passed=leakage_ratio <= max_leakage_ratio,
        metric=leakage_ratio,
        tolerance=max_leakage_ratio,
        details=(
            "A channel-spanning DEM particle represented as solid_boundary should "
            "strongly suppress outlet flux. "
            f"open_outlet_flux={q_open:.6g}, clogged_outlet_flux={q_clogged:.6g}, "
            f"clogged_inlet_flux={q_in_clogged:.6g}, flux_imbalance={flux_imbalance:.6g}, "
            f"open_pressure_drop={open_dp:.6g}, clogged_pressure_drop={clogged_dp:.6g}, "
            f"dynamic_solid_fraction={clogged_sim.dynamic_solid_fraction():.6g}."
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


def verify_porous_resistance_reduces_flux(
    nx: int,
    ny: int,
    steps: int,
    re: float,
    u_max: float,
    min_reduction: float,
) -> VerificationResult:
    """Check that particle-occupancy porous resistance lowers local flux."""
    common = dict(
        nx=nx,
        ny=ny,
        Re=re,
        u_max=u_max,
        reynolds_length=float(ny),
        flow_control="fixed_pressure",
        n_particles=1,
        particle_radius=max(2.0, 0.16 * ny),
        density_ratio=1e6,
        gravity=0.0,
        seed=16,
    )
    open_sim = fluid_only(FastLBMDEM(**common))
    porous_sim = FastLBMDEM(
        **common,
        porous_resistance=True,
        porous_resistance_coeff=0.08,
    )
    porous_sim.pos[:] = np.array([[0.5 * nx, 0.5 * ny]])
    porous_sim.vel[:] = 0.0

    open_sim.advance(steps)
    porous_sim.advance(steps)
    _, ux_open, _ = open_sim.get_fields()
    _, ux_porous, _ = porous_sim.get_fields()
    section = nx // 2
    q_open = positive_flux(ux_open, open_sim.solid, section)
    q_porous = positive_flux(ux_porous, porous_sim.solid, section)
    reduction = 1.0 - q_porous / max(q_open, 1e-12)
    return VerificationResult(
        name="porous_resistance_flux_reduction",
        passed=reduction >= min_reduction,
        metric=reduction,
        tolerance=min_reduction,
        details=(
            "Flux reduction from applying Brinkman-like drag in particle-occupied cells. "
            f"section={section}, open_flux={q_open:.6g}, porous_flux={q_porous:.6g}."
        ),
    )


def verify_cylinder_surface_force(
    nx: int,
    ny: int,
    re: float,
    u_max: float,
    min_force: float,
) -> VerificationResult:
    """Check that particle-cylinder Hamaker attraction produces a surface load."""
    sim = FastLBMDEM(
        nx=nx,
        ny=ny,
        Re=re,
        u_max=u_max,
        reynolds_length=float(ny),
        flow_control="fixed_pressure",
        n_particles=1,
        particle_radius=max(2.0, 0.08 * ny),
        gravity=0.0,
        particle_attraction=True,
        attraction_strength=1e-2,
        attraction_cutoff=4.0,
        attraction_min_gap=0.05,
        cylinders=[(0.45 * nx, 0.5 * ny, 0.10 * ny)],
        seed=17,
    )
    cx, cy, cr = sim.cylinders[0]
    gap = 1.0
    sim.pos[:] = np.array([[cx + cr + sim.radii[0] + gap, cy]])
    sim.vel[:] = 0.0
    forces, _ = sim._dem_loads(1.0 / sim.dem_substeps)
    force_x = float(forces[0, 0])
    return VerificationResult(
        name="particle_cylinder_hamaker_attraction",
        passed=force_x <= -min_force,
        metric=abs(force_x),
        tolerance=min_force,
        details=(
            "A particle just outside the cylinder attraction cutoff should feel a "
            f"force toward the cylinder. force_x={force_x:.6g}."
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
    porous_min_reduction: float = 0.01,
    cylinder_surface_min_force: float = 1e-6,
    near_clogged_max_leakage: float = 0.02,
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
        verify_near_clogged_solid_boundary_suppresses_flux(
            nx,
            ny,
            max(steps // 3, 1),
            re,
            u_max,
            near_clogged_max_leakage,
        ),
        verify_immersed_boundary_particle_load(
            nx,
            ny,
            max(steps // 3, 1),
            re,
            u_max,
            ibm_min_force,
        ),
        verify_porous_resistance_reduces_flux(
            nx,
            ny,
            max(steps // 3, 1),
            re,
            u_max,
            porous_min_reduction,
        ),
        verify_cylinder_surface_force(
            nx,
            ny,
            re,
            u_max,
            cylinder_surface_min_force,
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
        "solver_path": "particulate_flow.FastLBMDEM",
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
        "Solver path: `particulate_flow.FastLBMDEM`",
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
