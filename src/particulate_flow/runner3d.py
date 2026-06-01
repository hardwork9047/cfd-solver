"""3D LBM-DEM execution path (issue #15).

The main 2D runner (``run_lbm_dem.py``) is a monolithic pipeline built around the
2D ``FastLBMDEM`` (matplotlib animation, 2D VTK, fouling-specific analysis).  This
module provides a streamlined 3D path that ``run_lbm_dem.py`` dispatches to when
``--dimensions 3``: build the 3D solver, advance in snapshot chunks, record
dimension-agnostic scalar time-series, and write ParaView VTK for the 3D fields and
particles.  It reuses the 3D solver from :mod:`particulate_flow.lbm3d` and the DEM
contact from :mod:`particulate_flow.dem.contact3d`; the 2D path is untouched.
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

import numpy as np

from particulate_flow.io.paths import program_results_dir
from particulate_flow.lbm3d import LBMDEMSolver3D


def build_3d_solver(args: argparse.Namespace) -> LBMDEMSolver3D:
    """Construct an :class:`LBMDEMSolver3D` from parsed runner args.

    Mirrors the field names the 2D builder reads, but targets the 3D solver:
    pressure-driven flow, optional ``left_inlet`` particle source, immersed-boundary
    coupling, and z-aligned cylinder obstacles from ``cylinder_spec``.

    Args:
        args: Parsed namespace with at least ``nx, ny, nz, reynolds_number, u_max,
            streamwise_boundary, pressure_drop, rho_out`` and the particle/source
            fields.

    Returns:
        A ready-to-advance 3D solver.
    """
    cylinders = []
    for spec in getattr(args, "cylinder_spec", None) or []:
        cylinders.append((float(spec[0]), float(spec[1]), float(spec[2])))

    source = getattr(args, "particle_source", "none")
    source = "left_inlet" if source in ("left-inlet", "left_inlet") else "none"
    coupling = getattr(args, "particle_fluid_coupling", "none")
    coupling = "immersed_boundary" if coupling == "immersed_boundary" else "none"

    # The 3D solver only knows "periodic" / "pressure"; map the 2D-style
    # "periodic-force" alias to "periodic" so a config copied from a 2D case
    # does not crash deep inside the solver.
    streamwise = getattr(args, "streamwise_boundary", "pressure")
    if streamwise == "periodic-force":
        streamwise = "periodic"

    phi = getattr(args, "particle_volume_fraction", None)
    if phi is not None and phi > 1.0:
        phi /= 100.0

    return LBMDEMSolver3D(
        nx=args.nx,
        ny=args.ny,
        nz=getattr(args, "nz", 1),
        Re=getattr(args, "reynolds_number", 10.0),
        u_max=getattr(args, "u_max", 0.05),
        streamwise_boundary=streamwise,
        pressure_drop=getattr(args, "pressure_drop", None) or 0.0,
        rho_out=getattr(args, "rho_out", 1.0),
        cylinders=cylinders,
        particle_source=source,
        particle_fluid_coupling=coupling,
        ibm_stiffness=getattr(args, "ibm_stiffness", 1.0),
        ibm_marker_spacing=getattr(args, "ibm_marker_spacing", 1.0),
        particle_radius=getattr(args, "particle_radius", 2.0),
        density_ratio=getattr(args, "density_ratio", 2.0),
        gravity=getattr(args, "gravity", 0.0),
        dem_substeps=getattr(args, "dem_substeps", 4),
        sliding_friction=getattr(args, "sliding_friction", 0.5),
        rolling_friction_coeff=getattr(args, "rolling_friction_coeff", 0.05),
        particle_attraction=getattr(args, "particle_attraction", False),
        particle_repulsion=getattr(args, "particle_repulsion", False),
        attraction_strength=getattr(args, "attraction_strength", 1e-3),
        repulsion_strength=getattr(args, "repulsion_strength", 1e-3),
        attraction_cutoff=getattr(args, "attraction_cutoff", 3.0),
        repulsion_cutoff=getattr(args, "repulsion_cutoff", 3.0),
        attraction_min_gap=getattr(args, "attraction_min_gap", 0.05),
        repulsion_min_gap=getattr(args, "repulsion_min_gap", 0.05),
        source_volume_fraction=phi if source == "left_inlet" else None,
    )


def _scalar_row(sim: LBMDEMSolver3D, step: int) -> dict[str, float | int]:
    """Return a dimension-agnostic time-series row for the current state."""
    rho, ux, uy, uz = sim.get_fields()
    speed = np.sqrt(ux**2 + uy**2 + uz**2)
    n_p = sim.dem.n_p if sim.dem is not None else 0
    if n_p:
        pspeed = np.sqrt((sim.dem.vel**2).sum(axis=1))
        p_mean = float(pspeed.mean())
        p_max = float(pspeed.max())
    else:
        p_mean = p_max = 0.0
    return {
        "step": int(step),
        "u_max": float(speed.max()),
        "rho_mean": float(rho.mean()),
        "particle_count": int(n_p),
        "particle_speed_mean": p_mean,
        "particle_speed_max": p_max,
        "generated_particles": int(getattr(sim, "generated_particles", 0)),
        "removed_particles": int(getattr(sim, "removed_particles", 0)),
    }


def _vtk_float(value: float) -> str:
    """Format a float for VTK ASCII output."""
    return f"{float(value):.6g}"


def _write_fluid_vtk_3d(path: Path, sim: LBMDEMSolver3D) -> None:
    """Write the 3D fluid field as a STRUCTURED_POINTS VTK file.

    Records density and velocity-magnitude scalars plus the velocity vector and a
    solid mask over the (nx, ny, nz) lattice.

    Args:
        path: Output ``.vtk`` path.
        sim: The 3D solver to snapshot.
    """
    rho, ux, uy, uz = sim.get_fields()
    nx, ny, nz = sim.nx, sim.ny, sim.nz
    n = nx * ny * nz
    # VTK STRUCTURED_POINTS iterates x fastest, then y, then z.
    order = (2, 1, 0)  # transpose (nx,ny,nz) -> iterate z slowest
    rho_f = np.transpose(rho, order).ravel()
    ux_f = np.transpose(ux, order).ravel()
    uy_f = np.transpose(uy, order).ravel()
    uz_f = np.transpose(uz, order).ravel()
    solid_f = np.transpose(sim.solid, order).ravel().astype(int)
    speed_f = np.sqrt(ux_f**2 + uy_f**2 + uz_f**2)

    lines = [
        "# vtk DataFile Version 3.0",
        "LBMDEMSolver3D fluid field",
        "ASCII",
        "DATASET STRUCTURED_POINTS",
        f"DIMENSIONS {nx} {ny} {nz}",
        "ORIGIN 0 0 0",
        "SPACING 1 1 1",
        f"POINT_DATA {n}",
        "SCALARS density float 1",
        "LOOKUP_TABLE default",
    ]
    lines.extend(_vtk_float(v) for v in rho_f)
    lines.append("SCALARS speed float 1")
    lines.append("LOOKUP_TABLE default")
    lines.extend(_vtk_float(v) for v in speed_f)
    lines.append("SCALARS solid int 1")
    lines.append("LOOKUP_TABLE default")
    lines.extend(str(int(v)) for v in solid_f)
    lines.append("VECTORS velocity float")
    lines.extend(
        f"{_vtk_float(a)} {_vtk_float(b)} {_vtk_float(c)}" for a, b, c in zip(ux_f, uy_f, uz_f)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_particles_vtk_3d(path: Path, sim: LBMDEMSolver3D) -> None:
    """Write the 3D particles as a POLYDATA point cloud with radius and velocity.

    Args:
        path: Output ``.vtk`` path.
        sim: The 3D solver whose particles to write.
    """
    n_p = sim.dem.n_p if sim.dem is not None else 0
    lines = [
        "# vtk DataFile Version 3.0",
        "LBMDEMSolver3D particles",
        "ASCII",
        "DATASET POLYDATA",
        f"POINTS {n_p} float",
    ]
    if n_p:
        for p in sim.dem.pos:
            lines.append(f"{_vtk_float(p[0])} {_vtk_float(p[1])} {_vtk_float(p[2])}")
        lines.append(f"VERTICES {n_p} {2 * n_p}")
        lines.extend(f"1 {i}" for i in range(n_p))
        lines.append(f"POINT_DATA {n_p}")
        lines.append("SCALARS radius float 1")
        lines.append("LOOKUP_TABLE default")
        lines.extend(_vtk_float(r) for r in sim.dem.radii)
        lines.append("VECTORS velocity float")
        lines.extend(
            f"{_vtk_float(v[0])} {_vtk_float(v[1])} {_vtk_float(v[2])}" for v in sim.dem.vel
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_pvd(path: Path, entries: list[tuple[int, str]]) -> None:
    """Write a ParaView .pvd time-series index over (step, filename) entries."""
    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        "  <Collection>",
    ]
    for step, name in entries:
        lines.append(f'    <DataSet timestep="{step}" file="{name}"/>')
    lines.append("  </Collection>")
    lines.append("</VTKFile>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_3d(args: argparse.Namespace) -> Path:
    """Run a 3D LBM-DEM simulation end-to-end and write standard artifacts.

    Builds the solver, optionally warms up, advances in ``snapshot_every`` chunks to
    ``total_steps``, records a scalar time-series, and writes ``analysis/
    time_series.csv`` + ``.npz``, ``summary.json``, ``metadata.json``,
    ``run_status.json``, and ParaView VTK (fluid + particles + a ``.pvd`` index) to a
    timestamped result directory.

    Args:
        args: Parsed runner namespace (``dimensions == 3``).  An optional
            ``output_root`` attribute redirects the result directory (used by tests);
            otherwise results land under the standard ``results/run_lbm_dem/`` tree.

    Returns:
        The result directory path.
    """
    sim = build_3d_solver(args)

    tag_parts = ["dim3"]
    if getattr(args, "result_tag", None):
        tag_parts.append(args.result_tag)
    tag_parts.append(datetime.now().strftime("run_%Y%m%d_%H%M%S_%f"))
    output_root = getattr(args, "output_root", None)
    if output_root:
        out_dir = Path(output_root).joinpath(*tag_parts)
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        out_dir = program_results_dir("run_lbm_dem.py", *tag_parts)
    paraview_dir = out_dir / "paraview"
    paraview_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir = out_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    warmup = int(getattr(args, "warmup_steps", 0))
    total = int(getattr(args, "total_steps", 0))
    every = max(1, int(getattr(args, "snapshot_every", 100)))

    print("=" * 60)
    print("  LBM-DEM 3D シミュレーション")
    print("=" * 60)
    print(f"  Grid     : {sim.nx} × {sim.ny} × {sim.nz}")
    print(f"  Boundary : x={sim.streamwise_boundary}, y=periodic, z=periodic")
    print(f"  Cylinders: {len(sim.cylinders)} fixed")
    print(f"  Source   : {sim.particle_source}")

    if warmup:
        print(f"\n[1/2] ウォームアップ {warmup} ステップ ...")
        sim.advance(warmup)

    rows: list[dict[str, float | int]] = []
    fluid_entries: list[tuple[int, str]] = []
    particle_entries: list[tuple[int, str]] = []
    print(f"\n[2/2] メインシミュレーション {total} ステップ ...")
    done = 0
    while done < total:
        chunk = min(every, total - done)
        sim.advance(chunk)
        done += chunk
        rows.append(_scalar_row(sim, sim.step_count))
        fluid_name = f"fluid_{sim.step_count:06d}.vtk"
        particle_name = f"particles_{sim.step_count:06d}.vtk"
        _write_fluid_vtk_3d(paraview_dir / fluid_name, sim)
        _write_particles_vtk_3d(paraview_dir / particle_name, sim)
        fluid_entries.append((sim.step_count, fluid_name))
        particle_entries.append((sim.step_count, particle_name))

    _write_pvd(paraview_dir / "fluid_series.pvd", fluid_entries)
    _write_pvd(paraview_dir / "particles_series.pvd", particle_entries)

    # Time-series CSV + NPZ.
    csv_path = analysis_dir / "time_series.csv"
    if rows:
        fieldnames = list(rows[0].keys())
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        arrays = {k: np.array([r[k] for r in rows]) for k in fieldnames}
        np.savez(analysis_dir / "time_series.npz", **arrays)
    else:
        csv_path.write_text("", encoding="utf-8")

    # Summary.
    final = rows[-1] if rows else {}
    summary = {
        "dimensions": 3,
        "grid": [sim.nx, sim.ny, sim.nz],
        "streamwise_boundary": sim.streamwise_boundary,
        "total_steps": total,
        "warmup_steps": warmup,
        "n_cylinders": len(sim.cylinders),
        "particle_source": sim.particle_source,
        "final_particle_count": int(final.get("particle_count", 0)),
        "generated_particles": int(getattr(sim, "generated_particles", 0)),
        "removed_particles": int(getattr(sim, "removed_particles", 0)),
        "max_speed": float(final.get("u_max", 0.0)),
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Metadata (reproducibility): grid, physics, and the resolved solver config.
    metadata = {
        "dimensions": 3,
        "grid": [sim.nx, sim.ny, sim.nz],
        "reynolds_number": float(getattr(args, "reynolds_number", 0.0)),
        "u_max": float(getattr(args, "u_max", 0.0)),
        "pressure_drop": float(getattr(args, "pressure_drop", 0.0) or 0.0),
        "streamwise_boundary": sim.streamwise_boundary,
        "particle_source": sim.particle_source,
        "particle_radius": float(getattr(args, "particle_radius", 0.0)),
        "density_ratio": float(getattr(args, "density_ratio", 0.0)),
        "n_cylinders": len(sim.cylinders),
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    (out_dir / "run_status.json").write_text(
        json.dumps({"status": "completed", "dimensions": 3}, indent=2),
        encoding="utf-8",
    )

    print("\n完了!")
    print(f"  結果ディレクトリ: {out_dir}")
    print(f"  時系列CSV: {csv_path}")
    print(f"  ParaView : {paraview_dir / 'fluid_series.pvd'}")
    return out_dir
