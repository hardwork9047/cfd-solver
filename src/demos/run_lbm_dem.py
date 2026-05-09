"""
LBM-DEM 連成シミュレーション 実行スクリプト
==============================================
- 流体: D2Q9 LBM (Poiseuille チャネル流)
- 粒子: DEM (Hertz接触 + Stokes抗力)
- 動画: MP4 (ffmpeg) を出力
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import numpy as np

# パッケージを src/ から読み込む
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cfd_dem_lbm import FastLBMDEM, LBMDEMSolver
from cfd.result_paths import program_results_dir

# ---------------------------------------------------------------------------
# コマンドライン引数
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="LBM-DEM coupled simulation runner")
parser.add_argument("--nx", type=int, default=180,
                    help="Grid width in lattice nodes (default: 180)")
parser.add_argument("--ny", type=int, default=70,
                    help="Grid height in lattice nodes (default: 70)")
parser.add_argument("--reynolds-number", "--Re", dest="reynolds_number",
                    type=float, default=100.0,
                    help="Particle-diameter Reynolds number based on max flow speed "
                         "(default: 100)")
parser.add_argument("--u-max", type=float, default=0.05,
                    help="Target maximum fluid speed in lattice units (default: 0.05)")
parser.add_argument("--flow-control", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="Adjust drive force to keep max fluid speed near --u-max "
                         "(default: enabled)")
parser.add_argument("--flow-control-gain", type=float, default=0.2,
                    help="Relaxation gain for max-speed flow control (default: 0.2)")
parser.add_argument("--fluid-method",
                    choices=["lbm-bgk-guo", "lbm-trt-guo"],
                    default="lbm-bgk-guo",
                    help="LBM collision/forcing method (default: lbm-bgk-guo)")
parser.add_argument("--particle-method",
                    choices=["dem-hertz", "dem-linear"],
                    default="dem-hertz",
                    help="DEM normal contact model (default: dem-hertz)")
parser.add_argument("--particle-fluid-coupling",
                    choices=["point_force", "solid_boundary", "immersed_boundary"],
                    default="point_force",
                    help="Fluid-particle coupling mode: point force, dynamic "
                         "solid-boundary mask, or immersed boundary (default: point_force)")
parser.add_argument("--ibm-stiffness", type=float, default=1.0,
                    help="Direct-forcing immersed-boundary stiffness (default: 1.0)")
parser.add_argument("--ibm-marker-spacing", type=float, default=1.0,
                    help="Approximate immersed-boundary marker spacing in lattice units "
                         "(default: 1.0)")
parser.add_argument("--cylinder", action="store_true", help="Add a fixed solid cylinder")
parser.add_argument("--cyl-x", type=float, default=None,
                    help="Cylinder center x [lattice] (default: NX/4)")
parser.add_argument("--cyl-y", type=float, default=None,
                    help="Cylinder center y [lattice] (default: NY/2)")
parser.add_argument("--cyl-r", type=float, default=6.0,
                    help="Cylinder radius [lattice] (default: 6.0)")
parser.add_argument("--cylinder-spec", type=float, nargs=3, action="append",
                    metavar=("X", "Y", "R"), default=None,
                    help="Add a fixed cylinder as x y radius. Can be repeated "
                         "to define membrane-pore geometry")
parser.add_argument("--n-particles", type=int, default=None,
                    help="Particle count. Ignored when --particle-volume-fraction is set")
parser.add_argument("--particle-radius", type=float, default=3.0,
                    help="Mean particle radius in lattice nodes (default: 3.0)")
parser.add_argument("--radius-variation", type=float, default=0.15,
                    help="Uniform particle radius variation fraction (default: 0.15)")
parser.add_argument("--particle-volume-fraction", "--particle-area-fraction",
                    type=float, default=None,
                    help="Target 2-D particle area fraction relative to water area. "
                         "Use 0.05 or 5 for 5%%")
parser.add_argument("--particle-source", choices=["initial", "left-inlet"],
                    default="initial",
                    help="Particle supply mode: initial distribution or uniform left inlet")
parser.add_argument("--total-steps", type=int, default=20000,
                    help="Main simulation steps (default: 20000)")
parser.add_argument("--snapshot-every", type=int, default=40,
                    help="Steps between snapshots (default: 40)")
parser.add_argument("--warmup-steps", type=int, default=1000,
                    help="Warmup steps before recording (default: 1000)")
parser.add_argument("--no-video", action="store_true",
                    help="Skip MP4 generation")
parser.add_argument("--paraview-output", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="Export ParaView-readable VTK files (default: enabled)")
parser.add_argument("--paraview-every", type=int, default=0,
                    help="Export every Nth recorded snapshot to ParaView VTK. "
                         "0 exports only the final snapshot (default: 0)")
parser.add_argument("--snapshot-storage", choices=["disk", "memory"], default="disk",
                    help="Store recorded snapshots on disk or in memory (default: disk)")
parser.add_argument("--result-tag", default=None,
                    help="Extra result directory tag to avoid overwriting parameter sweeps")
parser.add_argument("--rolling-friction", action=argparse.BooleanOptionalAction,
                    default=False,
                    help="Enable angular motion, tangential friction, and rolling resistance")
parser.add_argument("--sliding-friction", type=float, default=0.5,
                    help="Coulomb limit for tangential contact force (default: 0.5)")
parser.add_argument("--tangential-damping", type=float, default=0.4,
                    help="Tangential slip damping scale (default: 0.4)")
parser.add_argument("--rolling-friction-coeff", type=float, default=0.05,
                    help="Rolling-resistance moment coefficient (default: 0.05)")
parser.add_argument("--rolling-damping", type=float, default=0.2,
                    help="Angular velocity damping scale for rolling resistance (default: 0.2)")
parser.add_argument("--particle-attraction", action="store_true",
                    help="Enable Hamaker-like particle-particle attraction")
parser.add_argument("--particle-repulsion", action="store_true",
                    help="Enable Hamaker-like particle-particle repulsion")
parser.add_argument("--attraction-strength", type=float, default=1e-3,
                    help="Dimensionless Hamaker-like attraction strength (default: 1e-3)")
parser.add_argument("--repulsion-strength", type=float, default=1e-3,
                    help="Dimensionless Hamaker-like repulsion strength (default: 1e-3)")
parser.add_argument("--attraction-cutoff", type=float, default=3.0,
                    help="Surface gap cutoff for attraction [lattice] (default: 3.0)")
parser.add_argument("--repulsion-cutoff", type=float, default=3.0,
                    help="Surface gap cutoff for repulsion [lattice] (default: 3.0)")
parser.add_argument("--attraction-min-gap", type=float, default=0.05,
                    help="Minimum surface gap for attraction regularisation (default: 0.05)")
parser.add_argument("--repulsion-min-gap", type=float, default=0.05,
                    help="Minimum surface gap for repulsion regularisation (default: 0.05)")
args = parser.parse_args()
if args.particle_attraction and args.particle_repulsion:
    parser.error("--particle-attraction and --particle-repulsion are mutually exclusive")
if args.n_particles is not None and args.n_particles <= 0:
    parser.error("--n-particles must be positive")
if args.nx <= 2:
    parser.error("--nx must be greater than 2")
if args.ny <= 2:
    parser.error("--ny must be greater than 2")
if args.reynolds_number <= 0.0:
    parser.error("--reynolds-number must be positive")
if args.u_max < 0.0:
    parser.error("--u-max must be non-negative")
if not (0.0 < args.flow_control_gain <= 1.0):
    parser.error("--flow-control-gain must be in the range (0, 1]")
if args.ibm_stiffness <= 0.0:
    parser.error("--ibm-stiffness must be positive")
if args.ibm_marker_spacing <= 0.0:
    parser.error("--ibm-marker-spacing must be positive")
if args.particle_radius <= 0.0:
    parser.error("--particle-radius must be positive")
if args.radius_variation < 0.0:
    parser.error("--radius-variation must be non-negative")
if args.particle_volume_fraction is not None and args.particle_volume_fraction <= 0.0:
    parser.error("--particle-volume-fraction must be positive")
if args.total_steps <= 0:
    parser.error("--total-steps must be positive")
if args.snapshot_every <= 0:
    parser.error("--snapshot-every must be positive")
if args.warmup_steps < 0:
    parser.error("--warmup-steps must be non-negative")
if args.paraview_every < 0:
    parser.error("--paraview-every must be non-negative")
if args.result_tag is not None and any(ch in args.result_tag for ch in "/\\"):
    parser.error("--result-tag must be a single directory name")
if args.cylinder_spec is not None:
    for spec in args.cylinder_spec:
        if spec[2] <= 0.0:
            parser.error("--cylinder-spec radius must be positive")
if args.cyl_r <= 0.0:
    parser.error("--cyl-r must be positive")

# ---------------------------------------------------------------------------
# シミュレーション設定
# ---------------------------------------------------------------------------

NX              = args.nx    # 格子幅
NY              = args.ny     # 格子高さ
RE              = args.reynolds_number  # 粒子径基準 Reynolds 数
U_MAX           = args.u_max   # 目標最大流速 (格子単位)
DEFAULT_N_PARTICLES = 40 # 粒子数 (分率未指定時)
RADIUS          = args.particle_radius    # 粒子平均半径 (格子単位)
REYNOLDS_LENGTH = 2.0 * RADIUS  # 代表長さ: 粒子径
RADIUS_VARIATION = args.radius_variation  # 粒子径バリエーション
DENSITY_RATIO   = 2.0    # ρ_p / ρ_f
GRAVITY         = 3e-5   # 重力加速度 (格子単位)
TOTAL_STEPS  = args.total_steps  # 総 LBM ステップ数
SNAPSHOT_EVERY = args.snapshot_every   # 何ステップごとにスナップショット取得
FPS          = 24     # 動画フレームレート

# 円柱設定:
#   --cylinder は従来互換の単円柱
#   --cylinder-spec X Y R は複数円柱を明示配置
CYLINDERS = []
if args.cylinder:
    CYL_X = args.cyl_x if args.cyl_x is not None else NX / 4
    CYL_Y = args.cyl_y if args.cyl_y is not None else NY / 2
    CYL_R = args.cyl_r
    CYLINDERS.append((CYL_X, CYL_Y, CYL_R))
if args.cylinder_spec is not None:
    CYLINDERS.extend((float(x), float(y), float(r)) for x, y, r in args.cylinder_spec)
CYLINDER = CYLINDERS[0] if CYLINDERS else None

def _normalise_fraction(value: float | None) -> float | None:
    """Accept 0.05 or 5 as a 5% area fraction."""
    if value is None:
        return None
    if value > 1.0:
        value /= 100.0
    if value <= 0.0 or value >= 1.0:
        parser.error("--particle-volume-fraction must be between 0 and 1, or 0 and 100%")
    return value


def _water_area(nx: int, ny: int, cylinders: list[tuple[float, float, float]]) -> float:
    """Approximate available water area in 2-D lattice units."""
    area = float(nx * (ny - 2))
    for _, _, radius in cylinders:
        area -= float(np.pi * radius ** 2)
    return max(area, 1.0)


def _poiseuille_flow_rate(ny: int, u_max: float) -> float:
    """Approximate left-boundary flow rate for the developed channel profile."""
    y_nodes = np.arange(1, ny - 1)
    channel_height = max(ny - 2.0, 1.0)
    eta = np.clip((y_nodes - 0.5) / channel_height, 0.0, 1.0)
    return float(np.sum(4.0 * u_max * eta * (1.0 - eta)))


PARTICLE_VOLUME_FRACTION = _normalise_fraction(args.particle_volume_fraction)
expected_particle_area = np.pi * RADIUS**2 * (1.0 + RADIUS_VARIATION**2 / 3.0)
if PARTICLE_VOLUME_FRACTION is None:
    N_PARTICLES = args.n_particles if args.n_particles is not None else DEFAULT_N_PARTICLES
    if args.particle_source == "left-inlet":
        expected_flow_area = _poiseuille_flow_rate(NY, U_MAX) * (args.warmup_steps + TOTAL_STEPS)
        SOURCE_VOLUME_FRACTION = min(
            0.95,
            N_PARTICLES * expected_particle_area / max(expected_flow_area, 1.0),
        )
    else:
        SOURCE_VOLUME_FRACTION = None
else:
    SOURCE_VOLUME_FRACTION = PARTICLE_VOLUME_FRACTION
    if args.particle_source == "left-inlet":
        expected_flow_area = _poiseuille_flow_rate(NY, U_MAX) * (args.warmup_steps + TOTAL_STEPS)
        N_PARTICLES = max(
            1,
            int(np.ceil(1.25 * PARTICLE_VOLUME_FRACTION * expected_flow_area / expected_particle_area)),
        )
    else:
        N_PARTICLES = max(
            1,
            int(round(PARTICLE_VOLUME_FRACTION * _water_area(NX, NY, CYLINDERS) / expected_particle_area)),
        )

if len(CYLINDERS) > 1:
    GEOMETRY_MODE = "multi_cylinder"
elif len(CYLINDERS) == 1:
    GEOMETRY_MODE = "cylinder"
else:
    GEOMETRY_MODE = "channel"
if args.particle_attraction:
    SURFACE_FORCE_MODE = "with_attraction"
elif args.particle_repulsion:
    SURFACE_FORCE_MODE = "with_repulsion"
else:
    SURFACE_FORCE_MODE = "no_surface_force"
ROLLING_MODE = "rolling" if args.rolling_friction else "free_roll"
SOURCE_MODE = "left_inlet" if args.particle_source == "left-inlet" else "initial"
out_parts = [f"{GEOMETRY_MODE}_{SURFACE_FORCE_MODE}_{ROLLING_MODE}_{SOURCE_MODE}"]
if args.result_tag:
    out_parts.append(args.result_tag)
if PARTICLE_VOLUME_FRACTION is not None:
    pct = int(round(PARTICLE_VOLUME_FRACTION * 100))
    out_parts.append(f"phi_{pct:02d}pct")
run_timestamp = datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
out_parts.append(run_timestamp)
OUT_DIR = program_results_dir(__file__, *out_parts)
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _git_value(args_: list[str]) -> str | None:
    """Return a git command result, or None when unavailable."""
    try:
        completed = subprocess.run(
            ["git", *args_],
            cwd=Path(__file__).resolve().parents[2],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()


def _write_metadata(path: Path, sim: LBMDEMSolver | None = None) -> None:
    """Write reproducibility metadata next to numerical outputs."""
    git_status = _git_value(["status", "--short"])
    metadata = {
        "schema_version": 1,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "program": str(Path(__file__).resolve()),
        "command": sys.argv,
        "output_dir": str(OUT_DIR),
        "python": {
            "version": sys.version,
            "executable": sys.executable,
            "platform": platform.platform(),
        },
        "git": {
            "commit": _git_value(["rev-parse", "HEAD"]),
            "branch": _git_value(["branch", "--show-current"]),
            "is_dirty": bool(git_status),
            "status_short": git_status or "",
        },
        "arguments": vars(args),
        "configuration": {
            "nx": NX,
            "ny": NY,
            "reynolds_number": RE,
            "target_u_max": U_MAX,
            "reynolds_length": REYNOLDS_LENGTH,
            "fluid_method": args.fluid_method,
            "particle_method": args.particle_method,
            "particle_fluid_coupling": args.particle_fluid_coupling,
            "ibm_stiffness": args.ibm_stiffness,
            "ibm_marker_spacing": args.ibm_marker_spacing,
            "particle_radius": RADIUS,
            "radius_variation": RADIUS_VARIATION,
            "density_ratio": DENSITY_RATIO,
            "gravity": GRAVITY,
            "total_steps": TOTAL_STEPS,
            "snapshot_every": SNAPSHOT_EVERY,
            "warmup_steps": args.warmup_steps,
            "particle_volume_fraction": PARTICLE_VOLUME_FRACTION,
            "source_volume_fraction": SOURCE_VOLUME_FRACTION,
            "n_particles_requested": N_PARTICLES,
            "cylinders": [
                {"x": float(cx), "y": float(cy), "radius": float(cr)}
                for cx, cy, cr in CYLINDERS
            ],
        },
    }
    if sim is not None:
        metadata["solver"] = {
            "fluid_method": sim.fluid_method,
            "particle_method": sim.particle_method,
            "particle_fluid_coupling": sim.particle_fluid_coupling,
            "ibm_stiffness": sim.ibm_stiffness,
            "ibm_marker_spacing": sim.ibm_marker_spacing,
            "nu": sim.nu,
            "tau": sim.tau,
            "omega": sim.omega,
            "initial_drive_force": sim.initial_F_drive,
            "current_drive_force": sim.F_drive,
            "fluid_area": sim.fluid_area,
            "solid_area_fraction": float(np.count_nonzero(sim.solid) / sim.solid.size),
        }
    path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")


def _vtk_float(value: float) -> str:
    """Format a float compactly for ASCII VTK."""
    return f"{float(value):.9g}"


def _write_fluid_vtk(path: Path, snap: dict, solid: np.ndarray) -> None:
    """Write 2-D fluid fields as a legacy VTK structured-points file."""
    nx, ny = snap["speed"].shape
    n_points = nx * ny
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write(f"LBM-DEM fluid field step {snap['step']}\n")
        handle.write("ASCII\n")
        handle.write("DATASET STRUCTURED_POINTS\n")
        handle.write(f"DIMENSIONS {nx} {ny} 1\n")
        handle.write("ORIGIN 0 0 0\n")
        handle.write("SPACING 1 1 1\n")
        handle.write(f"POINT_DATA {n_points}\n")
        handle.write("VECTORS velocity float\n")
        for y_idx in range(ny):
            for x_idx in range(nx):
                handle.write(
                    f"{_vtk_float(snap['ux'][x_idx, y_idx])} "
                    f"{_vtk_float(snap['uy'][x_idx, y_idx])} 0\n"
                )
        handle.write("SCALARS speed float 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for y_idx in range(ny):
            for x_idx in range(nx):
                handle.write(f"{_vtk_float(snap['speed'][x_idx, y_idx])}\n")
        handle.write("SCALARS pressure float 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for y_idx in range(ny):
            for x_idx in range(nx):
                handle.write(f"{_vtk_float(snap['pressure'][x_idx, y_idx])}\n")
        handle.write("SCALARS pressure_gauge float 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for y_idx in range(ny):
            for x_idx in range(nx):
                handle.write(f"{_vtk_float(snap['pressure_gauge'][x_idx, y_idx])}\n")
        handle.write("SCALARS solid int 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for y_idx in range(ny):
            for x_idx in range(nx):
                handle.write(f"{int(solid[x_idx, y_idx])}\n")


def _write_particles_vtk(path: Path, snap: dict) -> None:
    """Write DEM particles as VTK polydata points with radius and force scalars."""
    positions = snap["pos"]
    velocities = snap["vel"]
    radii = snap["radii"]
    forces = snap["total_force"]
    n_particles = len(positions)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write(f"LBM-DEM particles step {snap['step']}\n")
        handle.write("ASCII\n")
        handle.write("DATASET POLYDATA\n")
        handle.write(f"POINTS {n_particles} float\n")
        for pos in positions:
            handle.write(f"{_vtk_float(pos[0])} {_vtk_float(pos[1])} 0\n")
        handle.write(f"VERTICES {n_particles} {2 * n_particles}\n")
        for idx in range(n_particles):
            handle.write(f"1 {idx}\n")
        handle.write(f"POINT_DATA {n_particles}\n")
        handle.write("SCALARS radius float 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for radius in radii:
            handle.write(f"{_vtk_float(radius)}\n")
        handle.write("VECTORS velocity float\n")
        for vel in velocities:
            handle.write(f"{_vtk_float(vel[0])} {_vtk_float(vel[1])} 0\n")
        handle.write("SCALARS force_magnitude float 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for force in forces:
            handle.write(f"{_vtk_float(force)}\n")


def _write_cylinders_vtk(
    path: Path,
    cylinders: list[tuple[float, float, float]],
    n_segments: int = 96,
) -> None:
    """Write fixed cylinders as polygonal discs for ParaView."""
    points: list[tuple[float, float, float]] = []
    polygons: list[list[int]] = []
    for cx, cy, radius in cylinders:
        start = len(points)
        polygon = []
        for seg in range(n_segments):
            theta = 2.0 * np.pi * seg / n_segments
            points.append((cx + radius * np.cos(theta), cy + radius * np.sin(theta), 0.0))
            polygon.append(start + seg)
        polygons.append(polygon)

    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write("LBM-DEM fixed cylinders\n")
        handle.write("ASCII\n")
        handle.write("DATASET POLYDATA\n")
        handle.write(f"POINTS {len(points)} float\n")
        for x_pos, y_pos, z_pos in points:
            handle.write(f"{_vtk_float(x_pos)} {_vtk_float(y_pos)} {_vtk_float(z_pos)}\n")
        handle.write(f"POLYGONS {len(polygons)} {len(polygons) * (n_segments + 1)}\n")
        for polygon in polygons:
            handle.write(f"{len(polygon)} {' '.join(str(idx) for idx in polygon)}\n")
        handle.write(f"CELL_DATA {len(polygons)}\n")
        handle.write("SCALARS radius float 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for _, _, radius in cylinders:
            handle.write(f"{_vtk_float(radius)}\n")


def _write_pvd(path: Path, entries: list[tuple[int, Path]]) -> None:
    """Write a ParaView collection file for a VTK time series."""
    with path.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0"?>\n')
        handle.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        handle.write("  <Collection>\n")
        for step, vtk_path in entries:
            handle.write(
                f'    <DataSet timestep="{step}" group="" part="0" '
                f'file="{vtk_path.as_posix()}"/>\n'
            )
        handle.write("  </Collection>\n")
        handle.write("</VTKFile>\n")


def _paraview_snapshot_indices(n_snapshots: int, every: int) -> list[int]:
    """Return snapshot indices to export. 0 means final snapshot only."""
    if n_snapshots <= 0:
        return []
    if every == 0:
        return [n_snapshots - 1]
    indices = list(range(0, n_snapshots, every))
    if indices[-1] != n_snapshots - 1:
        indices.append(n_snapshots - 1)
    return indices


def _export_paraview_data(
    out_dir: Path,
    snapshots: list[dict],
    solid: np.ndarray,
    cylinders: list[tuple[float, float, float]],
    every: int,
) -> Path:
    """Export fluid, particle, and cylinder data for ParaView."""
    paraview_dir = out_dir / "paraview"
    paraview_dir.mkdir(parents=True, exist_ok=True)

    _write_cylinders_vtk(paraview_dir / "cylinders.vtk", cylinders)
    fluid_entries: list[tuple[int, Path]] = []
    particle_entries: list[tuple[int, Path]] = []
    for frame_idx in _paraview_snapshot_indices(len(snapshots), every):
        snap = snapshots[frame_idx]
        step = int(snap["step"])
        fluid_name = f"fluid_{frame_idx:04d}_step_{step:07d}.vtk"
        particle_name = f"particles_{frame_idx:04d}_step_{step:07d}.vtk"
        _write_fluid_vtk(paraview_dir / fluid_name, snap, solid)
        _write_particles_vtk(paraview_dir / particle_name, snap)
        fluid_entries.append((step, Path(fluid_name)))
        particle_entries.append((step, Path(particle_name)))

    _write_pvd(paraview_dir / "fluid_series.pvd", fluid_entries)
    _write_pvd(paraview_dir / "particles_series.pvd", particle_entries)
    return paraview_dir


def _write_snapshot_npz(path: Path, snap: dict) -> None:
    """Persist one snapshot so long runs do not keep every field in RAM."""
    np.savez(
        path,
        step=np.array(snap["step"], dtype=np.int64),
        ux=snap["ux"],
        uy=snap["uy"],
        speed=snap["speed"],
        pressure=snap["pressure"],
        pressure_gauge=snap["pressure_gauge"],
        pos=snap["pos"],
        vel=snap["vel"],
        radii=snap["radii"],
        masses=snap["masses"],
        total_force=snap["total_force"],
        removed_particles=np.array(snap["removed_particles"], dtype=np.int64),
    )


def _load_snapshot_npz(path: Path) -> dict:
    """Load a snapshot written by ``_write_snapshot_npz``."""
    with np.load(path) as data:
        return {
            "step": int(data["step"]),
            "ux": data["ux"],
            "uy": data["uy"],
            "speed": data["speed"],
            "pressure": data["pressure"],
            "pressure_gauge": data["pressure_gauge"],
            "pos": data["pos"],
            "vel": data["vel"],
            "radii": data["radii"],
            "masses": data["masses"],
            "total_force": data["total_force"],
            "removed_particles": int(data["removed_particles"]),
        }


def _boundary_flux(ux: np.ndarray, solid: np.ndarray, x_idx: int) -> float:
    """Positive x-direction volumetric flux through one vertical boundary."""
    fluid = ~solid[x_idx, :]
    return float(np.sum(np.maximum(ux[x_idx, fluid], 0.0)))


def _boundary_mean(field: np.ndarray, solid: np.ndarray, x_idx: int) -> float:
    """Mean scalar value on fluid nodes of one vertical boundary."""
    fluid = ~solid[x_idx, :]
    if not bool(np.any(fluid)):
        return float("nan")
    return float(np.mean(field[x_idx, fluid]))


def _contact_counts(sim: FastLBMDEM) -> dict[str, int]:
    """Count active DEM contact points for post-run analysis."""
    particle_particle = 0
    for i, j in sim._particle_pair_candidates():
        dist = float(np.hypot(sim.pos[j, 0] - sim.pos[i, 0], sim.pos[j, 1] - sim.pos[i, 1]))
        if dist < sim.radii[i] + sim.radii[j]:
            particle_particle += 1

    wall = 0
    if sim.n_p:
        wall_bot = sim.radii + 0.5
        wall_top = sim.ny - 1.5 - sim.radii
        wall = int(np.count_nonzero(sim.pos[:, 1] < wall_bot))
        wall += int(np.count_nonzero(sim.pos[:, 1] > wall_top))

    cylinder = 0
    cylinder_contacts_by_id: dict[str, int] = {}
    for idx, (cx, cy, cr) in enumerate(sim.cylinders, start=1):
        if sim.n_p == 0:
            cylinder_contacts_by_id[f"cylinder_{idx:02d}_contacts"] = 0
            continue
        dist = np.hypot(sim.pos[:, 0] - cx, sim.pos[:, 1] - cy)
        contacts = int(np.count_nonzero(dist < cr + sim.radii))
        cylinder_contacts_by_id[f"cylinder_{idx:02d}_contacts"] = contacts
        cylinder += contacts

    return {
        "particle_particle_contacts": particle_particle,
        "wall_contacts": wall,
        "cylinder_contacts": cylinder,
        "total_contacts": particle_particle + wall + cylinder,
        **cylinder_contacts_by_id,
    }


def _write_analysis_outputs(out_dir: Path, rows: list[dict[str, float | int]]) -> tuple[Path, Path]:
    """Write scalar time-series data as CSV and NPZ for later analysis."""
    analysis_dir = out_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    csv_path = analysis_dir / "time_series.csv"
    npz_path = analysis_dir / "time_series.npz"
    if not rows:
        csv_path.write_text("", encoding="utf-8")
        np.savez(npz_path)
        return csv_path, npz_path

    fieldnames = list(rows[0].keys())
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    arrays = {name: np.array([row[name] for row in rows]) for name in fieldnames}
    np.savez(npz_path, **arrays)
    return csv_path, npz_path

# ---------------------------------------------------------------------------
# 初期化
# ---------------------------------------------------------------------------

print("=" * 60)
print("  LBM-DEM 連成シミュレーション")
print("=" * 60)
print(
    f"Reynolds definition: Re_p = U_max * d_p / nu = {RE:.1f}, "
    f"d_p={REYNOLDS_LENGTH:.3f}, U_max={U_MAX:.4f}"
)
print(
    "Flow control: "
    + (
        f"target max velocity (gain={args.flow_control_gain:.2f})"
        if args.flow_control
        else "fixed pressure/body force"
    )
)
print(f"Fluid method: {args.fluid_method}")
print(f"Particle method: {args.particle_method}")
print(f"Particle-fluid coupling: {args.particle_fluid_coupling}")
if PARTICLE_VOLUME_FRACTION is not None:
    if args.particle_source == "left-inlet":
        print(
            f"Source particle area fraction: {PARTICLE_VOLUME_FRACTION:.1%} "
            f"-> up to {N_PARTICLES} queued particles"
        )
    else:
        print(
            f"Target particle area fraction: {PARTICLE_VOLUME_FRACTION:.1%} "
            f"-> {N_PARTICLES} particles"
        )
elif args.particle_source == "left-inlet":
    print(
        f"Derived source particle area fraction: {SOURCE_VOLUME_FRACTION:.3%} "
        f"for {N_PARTICLES} queued particles"
    )
print(f"Particle source: {args.particle_source}")
if CYLINDERS:
    print("Fixed cylinders:")
    for idx, (cx, cy, cr) in enumerate(CYLINDERS, start=1):
        print(f"  #{idx}: x={cx:.1f}, y={cy:.1f}, r={cr:.1f}")

sim = FastLBMDEM(
    nx=NX, ny=NY,
    Re=RE, u_max=U_MAX,
    reynolds_length=REYNOLDS_LENGTH,
    flow_control="target_max_velocity" if args.flow_control else "fixed_pressure",
    flow_control_gain=args.flow_control_gain,
    fluid_method=args.fluid_method,
    particle_method=args.particle_method,
    particle_fluid_coupling=args.particle_fluid_coupling,
    ibm_stiffness=args.ibm_stiffness,
    ibm_marker_spacing=args.ibm_marker_spacing,
    n_particles=N_PARTICLES,
    particle_radius=RADIUS,
    radius_variation=RADIUS_VARIATION,
    density_ratio=DENSITY_RATIO,
    gravity=GRAVITY,
    dem_substeps=4,
    seed=42,
    rolling_friction=args.rolling_friction,
    sliding_friction=args.sliding_friction,
    tangential_damping=args.tangential_damping,
    rolling_friction_coeff=args.rolling_friction_coeff,
    rolling_damping=args.rolling_damping,
    particle_attraction=args.particle_attraction,
    particle_repulsion=args.particle_repulsion,
    attraction_strength=args.attraction_strength,
    repulsion_strength=args.repulsion_strength,
    attraction_cutoff=args.attraction_cutoff,
    repulsion_cutoff=args.repulsion_cutoff,
    attraction_min_gap=args.attraction_min_gap,
    repulsion_min_gap=args.repulsion_min_gap,
    cylinders=CYLINDERS,
    particle_source=args.particle_source.replace("-", "_"),
    source_volume_fraction=SOURCE_VOLUME_FRACTION,
)
metadata_path = OUT_DIR / "metadata.json"
_write_metadata(metadata_path, sim)
print(f"Metadata saved: {metadata_path}")

# ---------------------------------------------------------------------------
# ウォームアップ (流れを発達させる)
# ---------------------------------------------------------------------------

WARMUP = args.warmup_steps
print(f"\n[1/3] ウォームアップ {WARMUP} ステップ ...")
t0 = time.perf_counter()
sim.advance(WARMUP)
print(f"      完了 ({time.perf_counter()-t0:.1f} s)")

# ---------------------------------------------------------------------------
# スナップショットを取りながらシミュレーション
# ---------------------------------------------------------------------------

n_frames = TOTAL_STEPS // SNAPSHOT_EVERY
if n_frames == 0:
    n_frames = 1
    SNAPSHOT_EVERY = TOTAL_STEPS
print(f"\n[2/3] メインシミュレーション {TOTAL_STEPS} ステップ "
      f"(スナップショット {n_frames} 枚) ...")

snapshots = [] if args.snapshot_storage == "memory" else None
snapshot_paths: list[Path] = []
snapshot_dir = OUT_DIR / "snapshots"
if args.snapshot_storage == "disk":
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    print(f"      snapshots: {snapshot_dir} に逐次保存")

paraview_dir = None
fluid_entries: list[tuple[int, Path]] = []
particle_entries: list[tuple[int, Path]] = []
if args.paraview_output:
    paraview_dir = OUT_DIR / "paraview"
    paraview_dir.mkdir(parents=True, exist_ok=True)
    _write_cylinders_vtk(paraview_dir / "cylinders.vtk", CYLINDERS)

last = None
speed_global_max = 0.0
force_global_min = np.inf
force_global_max = -np.inf
max_particles_in_frames = 0
analysis_rows: list[dict[str, float | int]] = []
reference_outlet_flux: float | None = None
reference_pressure_drop: float | None = None

t1 = time.perf_counter()
for frame_idx in range(n_frames):
    sim.advance(SNAPSHOT_EVERY)

    rho, ux, uy = sim.get_fields()

    # 各粒子に加わる合力の大きさ (重力 + 抗力 + 接触力)
    total_force_mags = np.linalg.norm(sim.forces_p, axis=1)
    pressure = rho / 3.0
    pressure_gauge = pressure - float(np.mean(pressure[~sim.solid]))
    speed = np.sqrt(ux**2 + uy**2)
    fluid_mask = ~sim.solid
    inlet_flux = _boundary_flux(ux, sim.solid, 0)
    outlet_flux = _boundary_flux(ux, sim.solid, NX - 1)
    inlet_pressure = _boundary_mean(pressure, sim.solid, 0)
    outlet_pressure = _boundary_mean(pressure, sim.solid, NX - 1)
    pressure_drop = inlet_pressure - outlet_pressure
    if reference_outlet_flux is None and outlet_flux > 1e-12:
        reference_outlet_flux = outlet_flux
    if reference_pressure_drop is None and np.isfinite(pressure_drop):
        reference_pressure_drop = pressure_drop
    contacts = _contact_counts(sim)
    generated_for_ratio = max(sim.generated_particles, 1)
    pass_ratio = sim.removed_particles / generated_for_ratio
    retained_ratio = max(sim.generated_particles - sim.removed_particles, 0) / generated_for_ratio
    normalized_permeate_flux = (
        outlet_flux / reference_outlet_flux
        if reference_outlet_flux is not None and reference_outlet_flux > 1e-12
        else 0.0
    )
    pressure_drop_ratio = (
        pressure_drop / reference_pressure_drop
        if reference_pressure_drop is not None and abs(reference_pressure_drop) > 1e-12
        else 1.0
    )
    inlet_outlet_flux_ratio = outlet_flux / inlet_flux if inlet_flux > 1e-12 else 0.0
    dynamic_solid_fraction = sim.dynamic_solid_fraction()
    total_solid_fraction = float(np.count_nonzero(sim.solid) / sim.solid.size)

    snap = {
        "step": sim.step_count,
        "ux":   ux.copy(),
        "uy":   uy.copy(),
        "speed": speed.copy(),
        "pressure": pressure.copy(),
        "pressure_gauge": pressure_gauge.copy(),
        "pos":  sim.pos.copy(),
        "vel":  sim.vel.copy(),
        "radii": sim.radii.copy(),
        "masses": sim.masses.copy(),
        "total_force": total_force_mags.copy(),
        "removed_particles": sim.removed_particles,
    }
    analysis_rows.append({
        "frame": frame_idx + 1,
        "step": sim.step_count,
        "time_lattice": float(sim.step_count),
        "active_particles": sim.n_p,
        "generated_particles": sim.generated_particles,
        "pending_particles": int(len(sim._pending_radii)),
        "passed_particles": sim.removed_particles,
        "inlet_flux": inlet_flux,
        "outlet_flux": outlet_flux,
        "permeate_flux": outlet_flux,
        "reynolds_number": RE,
        "reynolds_length": REYNOLDS_LENGTH,
        "target_u_max": U_MAX,
        "drive_force": sim.F_drive,
        "mean_ux": float(np.mean(ux[fluid_mask])),
        "max_speed": float(speed[fluid_mask].max()),
        "mean_pressure": float(np.mean(pressure[fluid_mask])),
        "inlet_pressure": inlet_pressure,
        "outlet_pressure": outlet_pressure,
        "pressure_drop": pressure_drop,
        "pressure_drop_ratio": pressure_drop_ratio,
        "mean_pressure_gauge": float(np.mean(pressure_gauge[fluid_mask])),
        "particle_area_fraction": float(np.sum(np.pi * sim.radii**2) / sim.fluid_area),
        "retained_particle_ratio": retained_ratio,
        "passed_particle_ratio": pass_ratio,
        "normalized_permeate_flux": normalized_permeate_flux,
        "inlet_outlet_flux_ratio": inlet_outlet_flux_ratio,
        "dynamic_particle_solid_fraction": dynamic_solid_fraction,
        "total_solid_fraction": total_solid_fraction,
        "fouling_resistance_index": pressure_drop_ratio / max(normalized_permeate_flux, 1e-12),
        "particle_ke": 0.5 * float(np.sum(sim.masses * np.sum(sim.vel**2, axis=1))),
        "mean_force": float(total_force_mags.mean()) if len(total_force_mags) else 0.0,
        **contacts,
    })
    last = snap
    speed_global_max = max(speed_global_max, float(snap["speed"].max()))
    if len(total_force_mags):
        force_global_min = min(force_global_min, float(total_force_mags.min()))
        force_global_max = max(force_global_max, float(total_force_mags.max()))
    max_particles_in_frames = max(max_particles_in_frames, len(snap["pos"]))

    if snapshots is not None:
        snapshots.append(snap)
    else:
        snapshot_path = snapshot_dir / f"snapshot_{frame_idx:04d}_step_{sim.step_count:07d}.npz"
        _write_snapshot_npz(snapshot_path, snap)
        snapshot_paths.append(snapshot_path)

    if paraview_dir is not None and args.paraview_every > 0 and frame_idx % args.paraview_every == 0:
        step = int(snap["step"])
        fluid_name = f"fluid_{frame_idx:04d}_step_{step:07d}.vtk"
        particle_name = f"particles_{frame_idx:04d}_step_{step:07d}.vtk"
        _write_fluid_vtk(paraview_dir / fluid_name, snap, sim.solid)
        _write_particles_vtk(paraview_dir / particle_name, snap)
        fluid_entries.append((step, Path(fluid_name)))
        particle_entries.append((step, Path(particle_name)))

    elapsed = time.perf_counter() - t1
    frac = (frame_idx + 1) / n_frames
    eta = elapsed / frac - elapsed if frac > 0 else 0
    speed_max = float(np.sqrt(ux**2 + uy**2).max())
    p_ke = 0.5 * float(np.sum(sim.masses * np.sum(sim.vel**2, axis=1)))
    f_mean = float(total_force_mags.mean()) if len(total_force_mags) else 0.0
    print(f"  frame {frame_idx+1:>3}/{n_frames}  step={sim.step_count:>6,}"
          f"  |u|_max={speed_max:.5f}  KE_p={p_ke:.3e}  |F|_mean={f_mean:.3e}  ETA {eta:.0f}s")

total_time = time.perf_counter() - t1
print(f"      完了 ({total_time:.1f} s, {TOTAL_STEPS/total_time:.0f} steps/s)")

analysis_csv, analysis_npz = _write_analysis_outputs(OUT_DIR, analysis_rows)
print(f"\n後解析データ保存: {analysis_csv}")
print(f"後解析NPZ保存  : {analysis_npz}")

if last is None:
    raise RuntimeError("No snapshots were recorded")
if not np.isfinite(force_global_min):
    force_global_min = 0.0
    force_global_max = 1.0
if force_global_max <= force_global_min:
    force_global_max = force_global_min + 1e-12

if paraview_dir is not None:
    should_export_final = (
        args.paraview_every == 0
        or not fluid_entries
        or fluid_entries[-1][0] != int(last["step"])
    )
    if should_export_final:
        frame_idx = n_frames - 1
        step = int(last["step"])
        fluid_name = f"fluid_{frame_idx:04d}_step_{step:07d}.vtk"
        particle_name = f"particles_{frame_idx:04d}_step_{step:07d}.vtk"
        _write_fluid_vtk(paraview_dir / fluid_name, last, sim.solid)
        _write_particles_vtk(paraview_dir / particle_name, last)
        fluid_entries.append((step, Path(fluid_name)))
        particle_entries.append((step, Path(particle_name)))
    _write_pvd(paraview_dir / "fluid_series.pvd", fluid_entries)
    _write_pvd(paraview_dir / "particles_series.pvd", particle_entries)
    if args.paraview_every == 0:
        print(f"\nParaViewデータ保存: {paraview_dir} (final snapshot only)")
    else:
        print(f"\nParaViewデータ保存: {paraview_dir} (every {args.paraview_every} snapshots)")

# ---------------------------------------------------------------------------
# 統計出力
# ---------------------------------------------------------------------------

print("\n--- 最終統計 ---")
p_vel_mag = np.linalg.norm(last["vel"], axis=1)
print(f"  粒子数          : {sim.n_p}")
print(f"  左入口投入済み  : {sim.generated_particles}/{sim.total_particles_requested}")
print(f"  左入口待機粒子数: {len(sim._pending_radii)}")
print(f"  右出口削除粒子数: {sim.removed_particles}")
print(f"  流体最大速度    : {last['speed'].max():.5f} (格子単位)")
if len(p_vel_mag):
    print(f"  粒子平均速度    : {p_vel_mag.mean():.4e} (格子単位)")
    print(f"  粒子最大速度    : {p_vel_mag.max():.4e} (格子単位)")
    print(f"  粒子KE          : {0.5*np.sum(last['masses']*np.sum(last['vel']**2, axis=1)):.4e}")
    print(f"  粒子面積分率    : {np.sum(np.pi*last['radii']**2)/sim.fluid_area:.3f}")
    print(f"  粒子Y重心       : {last['pos'][:,1].mean():.2f} / {NY} (格子単位)")
    print(f"  粒子Y重心 (正規): {last['pos'][:,1].mean()/NY:.3f}")
else:
    if len(sim._pending_radii):
        print("  粒子統計        : 現時点では有効粒子なし (入口投入待ちあり)")
    else:
        print("  粒子統計        : 全粒子が右出口から流出済み")
if args.particle_source == "left-inlet":
    print(f"  最終入口流量    : {sim.last_inlet_flow_rate:.4e} (格子面積/step)")
    print(f"  投入済み粒子面積: {sim.injected_particle_area:.4e}")
    print(f"  未投入面積予算  : {sim.inlet_particle_area_budget:.4e}")

# 静止画を保存
fig_stat, axes = plt.subplots(1, 3, figsize=(16, 5))
fig_stat.suptitle(
    f"LBM-DEM  Re={RE:.0f}  step={sim.step_count:,}  {sim.n_p} particles  ({NX}×{NY})",
    fontsize=12,
)
snap = last
x = np.arange(NX)
y = np.arange(NY)

def _add_cylinder_patch(ax, color="cyan", alpha=0.6, zorder=4):
    """静止画用: 円柱パッチを追加する。"""
    for cx, cy, cr in CYLINDERS:
        ax.add_patch(plt.Circle(
            (cx, cy), cr,
            color=color, alpha=alpha, zorder=zorder,
        ))

ax = axes[0]
im = ax.imshow(snap["speed"].T, origin="lower", cmap="inferno",
               extent=[0, NX, 0, NY], aspect="auto")
for i in range(len(snap["pos"])):
    c = plt.Circle((snap["pos"][i,0], snap["pos"][i,1]), snap["radii"][i],
                   color="white", lw=0.8, fill=False)
    ax.add_patch(c)
_add_cylinder_patch(ax)
ax.set_title("速度大きさ |u|")
ax.set_xlabel("x [格子]"); ax.set_ylabel("y [格子]")
fig_stat.colorbar(im, ax=ax, shrink=0.7)

ax = axes[1]
lw_arr = 1.5 * snap["speed"].T / (snap["speed"].T.max() + 1e-12)
ax.streamplot(x, y, snap["ux"].T, snap["uy"].T,
              color=snap["speed"].T, cmap="cool",
              linewidth=lw_arr, density=1.2, arrowsize=0.8)
for i in range(len(snap["pos"])):
    c = plt.Circle((snap["pos"][i,0], snap["pos"][i,1]), snap["radii"][i],
                   color="white", lw=0.8, fill=False)
    ax.add_patch(c)
_add_cylinder_patch(ax)
ax.set_xlim(0, NX); ax.set_ylim(0, NY)
ax.set_title("流線"); ax.set_xlabel("x [格子]")

ax = axes[2]
if len(snap["pos"]):
    sc = ax.scatter(snap["pos"][:,0], snap["pos"][:,1],
                    c=snap["total_force"], cmap="plasma", s=(snap["radii"]*4)**2,
                    edgecolors="k", lw=0.5, zorder=5)
else:
    sc = ax.scatter([], [], c=[], cmap="plasma", edgecolors="k", lw=0.5, zorder=5)
ax.imshow(snap["speed"].T, origin="lower", cmap="Blues",
          extent=[0, NX, 0, NY], aspect="auto", alpha=0.5)
_add_cylinder_patch(ax)
fig_stat.colorbar(sc, ax=ax, label="合力 |F_total| [格子単位]", shrink=0.7)
ax.set_xlim(0, NX); ax.set_ylim(0, NY)
ax.set_title("粒子位置 (合力でカラー)"); ax.set_xlabel("x [格子]")

plt.tight_layout()
static_path = OUT_DIR / "lbm_dem_final.png"
fig_stat.savefig(static_path, dpi=150, bbox_inches="tight")
print(f"\n静止画保存: {static_path}")
plt.close(fig_stat)

# ---------------------------------------------------------------------------
# 動画作成
# ---------------------------------------------------------------------------

if args.no_video:
    print("\n[3/3] 動画作成をスキップ (--no-video)")
    _write_metadata(metadata_path, sim)
    print(f"\n出力ファイル:")
    print(f"  静止画: {static_path}")
    print(f"  時系列CSV: {analysis_csv}")
    print(f"  時系列NPZ: {analysis_npz}")
    if paraview_dir is not None:
        print(f"  ParaView: {paraview_dir}")
        print(f"    流体時系列: {paraview_dir / 'fluid_series.pvd'}")
        print(f"    粒子時系列: {paraview_dir / 'particles_series.pvd'}")
        print(f"    固定円柱  : {paraview_dir / 'cylinders.vtk'}")
    print(f"\n完了!")
    sys.exit(0)

print(f"\n[3/3] 動画作成 ({n_frames} フレーム, {FPS} fps) ...")

fig_anim, ax_anim = plt.subplots(figsize=(10, 4))
fig_anim.patch.set_facecolor("#0a0a0a")
ax_anim.set_facecolor("#0a0a0a")

def _snapshot_at(frame_idx: int) -> dict:
    if snapshots is not None:
        return snapshots[frame_idx]
    return _load_snapshot_npz(snapshot_paths[frame_idx])


snap0 = _snapshot_at(0)

# 合力の全フレームにわたるグローバルmin/max（カラースケール固定）
force_cmap = plt.cm.plasma
force_norm = plt.Normalize(vmin=force_global_min, vmax=force_global_max)

im_fluid = ax_anim.imshow(
    snap0["speed"].T,
    origin="lower", cmap="inferno",
    extent=[0, NX, 0, NY], aspect="auto",
    vmin=0, vmax=speed_global_max,
    animated=True,
)
cbar_fluid = fig_anim.colorbar(im_fluid, ax=ax_anim, shrink=0.75, pad=0.01)
cbar_fluid.set_label("|u| [格子単位]", color="white")
cbar_fluid.ax.yaxis.set_tick_params(color="white")
plt.setp(cbar_fluid.ax.yaxis.get_ticklabels(), color="white")

# 合力カラーバー（ScalarMappable で追加）
sm_force = plt.cm.ScalarMappable(cmap=force_cmap, norm=force_norm)
sm_force.set_array([])
cbar_force = fig_anim.colorbar(sm_force, ax=ax_anim, shrink=0.75, pad=0.12)
cbar_force.set_label("|F_total| [格子単位]", color="white")
cbar_force.ax.yaxis.set_tick_params(color="white")
plt.setp(cbar_force.ax.yaxis.get_ticklabels(), color="white")

# 固定円柱パッチ (アニメーション中は動かないので animated=False)
for cx, cy, cr in CYLINDERS:
    ax_anim.add_patch(mpatches.Circle(
        (cx, cy), cr,
        linewidth=1.5, edgecolor="cyan", facecolor="cyan", alpha=0.55, zorder=4,
    ))

# 粒子の円パッチ (合力に応じた初期色、per-particle半径でサイズ)
particle_circles = []
for i in range(max_particles_in_frames):
    visible = i < len(snap0["pos"])
    center = (snap0["pos"][i, 0], snap0["pos"][i, 1]) if visible else (-100.0, -100.0)
    radius = snap0["radii"][i] if visible else RADIUS
    force_value = snap0["total_force"][i] if visible else force_global_min
    c = mpatches.Circle(
        center,
        radius,
        linewidth=1.2,
        edgecolor="white",
        facecolor=force_cmap(force_norm(force_value)),
        alpha=0.85 if visible else 0.0,
        animated=True,
    )
    ax_anim.add_patch(c)
    particle_circles.append(c)

title_txt = ax_anim.set_title(
    f"LBM-DEM  step={snap0['step']:,}", color="white", fontsize=11
)
ax_anim.set_xlabel("x [格子]", color="white")
ax_anim.set_ylabel("y [格子]", color="white")
ax_anim.tick_params(colors="white")
for spine in ax_anim.spines.values():
    spine.set_edgecolor("white")

plt.tight_layout()


def update(frame_idx: int):
    snap = _snapshot_at(frame_idx)
    im_fluid.set_data(snap["speed"].T)
    for i, c in enumerate(particle_circles):
        if i < len(snap["pos"]):
            c.center = (snap["pos"][i, 0], snap["pos"][i, 1])
            c.radius = snap["radii"][i]
            # 合力の大きさに応じて色を更新
            rgba = force_cmap(force_norm(snap["total_force"][i]))
            c.set_facecolor(rgba)
            c.set_alpha(0.85)
        else:
            c.center = (-100.0, -100.0)
            c.set_alpha(0.0)
    p_ke = 0.5 * float(np.sum(snap["masses"] * np.sum(snap["vel"] ** 2, axis=1)))
    force_mean = float(snap["total_force"].mean()) if len(snap["total_force"]) else 0.0
    title_txt.set_text(
        f"LBM-DEM  Re={RE:.0f}  step={snap['step']:,}"
        f"  KE_p={p_ke:.3e}  |F_total|_mean={force_mean:.3e}"
    )
    artists = [im_fluid, title_txt] + particle_circles
    return artists


ani = animation.FuncAnimation(
    fig_anim, update,
    frames=n_frames,
    interval=1000 / FPS,
    blit=True,
)

video_path = OUT_DIR / "lbm_dem_simulation.mp4"
writer = animation.FFMpegWriter(fps=FPS, bitrate=2000,
                                 extra_args=["-pix_fmt", "yuv420p"])
ani.save(str(video_path), writer=writer, dpi=150)
plt.close(fig_anim)
_write_metadata(metadata_path, sim)

print(f"動画保存: {video_path}")
print(f"\n出力ファイル:")
print(f"  静止画: {static_path}")
print(f"  動  画: {video_path}")
print(f"  時系列CSV: {analysis_csv}")
print(f"  時系列NPZ: {analysis_npz}")
if paraview_dir is not None:
    print(f"  ParaView: {paraview_dir}")
    print(f"    流体時系列: {paraview_dir / 'fluid_series.pvd'}")
    print(f"    粒子時系列: {paraview_dir / 'particles_series.pvd'}")
    print(f"    固定円柱  : {paraview_dir / 'cylinders.vtk'}")
print(f"\n完了!")
