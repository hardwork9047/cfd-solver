"""Solver builder: assemble FastLBMDEM or DEMPackingSimulation from parsed args.

Centralises the solver construction logic that previously lived inline in the
runner scripts.  Runners become thin entry-points that call one builder function.
"""

from __future__ import annotations

import argparse

import numpy as np

from particulate_flow.dem.packing import DEMPackingSimulation
from particulate_flow.fast_solver import FastLBMDEM
from particulate_flow.geometry.pore import PoreGeometry
from particulate_flow.lbm.constants import FLUID_ACCELERATORS, FLUID_METHODS

# ---------------------------------------------------------------------------
# Flow-control mapping (shared between runner and builder)
# ---------------------------------------------------------------------------

FLOW_CONTROL_MAP: dict[str, str] = {
    "fixed-pressure": "fixed_pressure",
    "constant-pressure": "constant_pressure",
    "target-max-velocity": "target_max_velocity",
    "constant-flux": "constant_flux",
}


def _poiseuille_flow_rate(ny: int, u_max: float) -> float:
    """Approximate left-boundary flow rate for a developed channel profile."""
    y_nodes = np.arange(1, ny - 1)
    channel_height = max(ny - 2.0, 1.0)
    eta = np.clip((y_nodes - 0.5) / channel_height, 0.0, 1.0)
    return float(np.sum(4.0 * u_max * eta * (1.0 - eta)))


def _build_geometry(args: argparse.Namespace, nx: int) -> PoreGeometry:
    """Construct PoreGeometry from parsed cylinder args."""
    cylinders: list[tuple[float, float, float]] = []
    if getattr(args, "cylinder", False):
        cyl_x = args.cyl_x if getattr(args, "cyl_x", None) is not None else nx / 4
        cyl_y = args.cyl_y if getattr(args, "cyl_y", None) is not None else args.ny / 2
        cyl_r = args.cyl_r
        cylinders.append((cyl_x, cyl_y, cyl_r))
    if getattr(args, "cylinder_spec", None) is not None:
        cylinders.extend(
            (float(x), float(y), float(r)) for x, y, r in args.cylinder_spec
        )
    return PoreGeometry.from_cylinders(cylinders)


def build_lbm_dem_solver(args: argparse.Namespace) -> FastLBMDEM:
    """Assemble a FastLBMDEM from argparse Namespace values.

    Extracts geometry construction and particle-count derivation from the
    runner, so the runner only needs to call this function and handle I/O.

    Args:
        args: Parsed namespace from the LBM-DEM runner's ArgumentParser,
              after config defaults have been applied.

    Returns:
        Fully initialised FastLBMDEM ready for time-stepping.
    """
    nx: int = args.nx
    ny: int = args.ny
    re: float = args.reynolds_number
    u_max: float = args.u_max
    radius: float = args.particle_radius
    reynolds_length: float = 2.0 * radius
    radius_variation: float = args.radius_variation
    density_ratio: float = args.density_ratio
    gravity: float = args.gravity

    fluid_method: str = getattr(args, "fluid_method", "lbm-trt-guo")
    if fluid_method not in FLUID_METHODS:
        raise ValueError(f"Unknown fluid_method {fluid_method!r}. Valid: {sorted(FLUID_METHODS)}")
    fluid_accelerator: str = getattr(args, "fluid_accelerator", "auto")
    if fluid_accelerator not in FLUID_ACCELERATORS:
        raise ValueError(
            f"Unknown fluid_accelerator {fluid_accelerator!r}. Valid: {sorted(FLUID_ACCELERATORS)}"
        )

    flow_condition = (
        args.flow_condition
        if args.flow_condition is not None
        else ("target-max-velocity" if args.flow_control else "fixed-pressure")
    )
    flow_control = FLOW_CONTROL_MAP[flow_condition]

    geometry = _build_geometry(args, nx)

    # Derive particle count from volume fraction if given.
    particle_volume_fraction = getattr(args, "particle_volume_fraction", None)
    if particle_volume_fraction is not None and particle_volume_fraction > 1.0:
        particle_volume_fraction /= 100.0

    expected_particle_area = np.pi * radius**2 * (1.0 + radius_variation**2 / 3.0)
    total_steps: int = args.total_steps
    warmup_steps: int = getattr(args, "warmup_steps", 0)
    particle_source: str = getattr(args, "particle_source", "domain").replace("-", "_")
    y_boundary: str = getattr(args, "y_boundary", "wall")

    if particle_volume_fraction is None:
        default_n = 40
        n_particles: int = args.n_particles if args.n_particles is not None else default_n
        if particle_source == "left_inlet":
            expected_flow_area = _poiseuille_flow_rate(ny, u_max) * (warmup_steps + total_steps)
            source_volume_fraction: float | None = min(
                0.95,
                n_particles * expected_particle_area / max(expected_flow_area, 1.0),
            )
        else:
            source_volume_fraction = None
    else:
        source_volume_fraction = particle_volume_fraction
        if particle_source == "left_inlet":
            expected_flow_area = _poiseuille_flow_rate(ny, u_max) * (warmup_steps + total_steps)
            n_particles = max(
                1,
                int(np.ceil(1.25 * particle_volume_fraction * expected_flow_area  # noqa: E501
                    / expected_particle_area)),
            )
        else:
            n_particles = max(
                1,
                int(round(
                    particle_volume_fraction
                    * geometry.water_area(nx, ny, wall_y=(y_boundary == "wall"))
                    / expected_particle_area
                )),
            )

    return FastLBMDEM(
        nx=nx, ny=ny,
        Re=re, u_max=u_max,
        reynolds_length=reynolds_length,
        flow_control=flow_control,
        flow_control_gain=getattr(args, "flow_control_gain", 1.0),
        y_boundary=y_boundary,
        streamwise_boundary=getattr(args, "streamwise_boundary", "periodic").replace("-", "_"),
        pressure_drop=getattr(args, "pressure_drop", None),
        rho_out=getattr(args, "rho_out", 1.0),
        fluid_method=fluid_method,
        fluid_accelerator=fluid_accelerator,
        compute_accelerator=getattr(args, "compute_accelerator", "auto"),
        particle_method=getattr(args, "particle_method", "dem-hertz"),
        particle_search=getattr(args, "particle_search", "cell_list"),
        particle_fluid_coupling=getattr(args, "particle_fluid_coupling", "immersed_boundary"),
        ibm_stiffness=getattr(args, "ibm_stiffness", 0.5),
        ibm_marker_spacing=getattr(args, "ibm_marker_spacing", 2.0),
        n_particles=n_particles,
        particle_radius=radius,
        radius_variation=radius_variation,
        density_ratio=density_ratio,
        gravity=gravity,
        dem_substeps=getattr(args, "dem_substeps", 5),
        seed=42,
        rolling_friction=getattr(args, "rolling_friction", True),
        sliding_friction=getattr(args, "sliding_friction", 0.5),
        tangential_damping=getattr(args, "tangential_damping", 0.5),
        rolling_friction_coeff=getattr(args, "rolling_friction_coeff", 0.1),
        rolling_damping=getattr(args, "rolling_damping", 0.35),
        particle_attraction=getattr(args, "particle_attraction", False),
        particle_repulsion=getattr(args, "particle_repulsion", False),
        attraction_strength=getattr(args, "attraction_strength", 1e-3),
        repulsion_strength=getattr(args, "repulsion_strength", 1e-3),
        attraction_cutoff=getattr(args, "attraction_cutoff", 3.0),
        repulsion_cutoff=getattr(args, "repulsion_cutoff", 3.0),
        attraction_min_gap=getattr(args, "attraction_min_gap", 0.05),
        repulsion_min_gap=getattr(args, "repulsion_min_gap", 0.05),
        porous_resistance=getattr(args, "porous_resistance", False),
        porous_resistance_coeff=getattr(args, "porous_resistance_coeff", 0.0),
        geometry=geometry,
        particle_source=particle_source,
        source_volume_fraction=source_volume_fraction,
        le_shear_rate=getattr(args, "le_shear_rate", 0.0),
        le_shear_axis=getattr(args, "le_shear_axis", 0),
        le_boundary_axis=getattr(args, "le_boundary_axis", 1),
        le_interpolation_order=getattr(args, "le_interpolation_order", 3),
        surface_roughness=getattr(args, "surface_roughness", 0.0),
    )


def build_dem_packing_solver(args: argparse.Namespace) -> DEMPackingSimulation:
    """Assemble a DEMPackingSimulation from argparse Namespace values.

    Args:
        args: Parsed namespace after config defaults have been applied.

    Returns:
        Fully initialised DEMPackingSimulation ready for time-stepping.
    """
    cylinders: list[tuple[float, float, float]] = []
    raw_cylinders = getattr(args, "cylinders", None) or []
    for item in raw_cylinders:
        if isinstance(item, (list, tuple)):
            cylinders.append((float(item[0]), float(item[1]), float(item[2])))
        else:
            cylinders.append(item)

    return DEMPackingSimulation(
        nx=args.nx,
        ny=args.ny,
        n_particles=args.n_particles,
        particle_radius=args.particle_radius,
        radius_variation=getattr(args, "radius_variation", 0.1),
        density_ratio=getattr(args, "density_ratio", 2.0),
        gravity=args.gravity,
        k_n=getattr(args, "k_n", 120.0),
        damping=getattr(args, "damping", 0.8),
        linear_damping=getattr(args, "linear_damping", 0.06),
        dem_substeps=getattr(args, "dem_substeps", 5),
        seed=getattr(args, "seed", 42),
        rolling_friction=getattr(args, "rolling_friction", True),
        sliding_friction=getattr(args, "sliding_friction", 0.5),
        tangential_damping=getattr(args, "tangential_damping", 0.5),
        rolling_friction_coeff=getattr(args, "rolling_friction_coeff", 0.1),
        rolling_damping=getattr(args, "rolling_damping", 0.35),
        particle_attraction=getattr(args, "particle_attraction", False),
        particle_repulsion=getattr(args, "particle_repulsion", False),
        attraction_strength=getattr(args, "attraction_strength", 1e-3),
        repulsion_strength=getattr(args, "repulsion_strength", 1e-3),
        attraction_cutoff=getattr(args, "attraction_cutoff", 3.0),
        repulsion_cutoff=getattr(args, "repulsion_cutoff", 3.0),
        attraction_min_gap=getattr(args, "attraction_min_gap", 0.05),
        repulsion_min_gap=getattr(args, "repulsion_min_gap", 0.05),
        cylinders=cylinders,
        particle_method=getattr(args, "particle_method", "dem-hertz"),
        particle_search=getattr(args, "particle_search", "cell_list"),
    )
