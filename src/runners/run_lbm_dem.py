"""
LBM-DEM 連成シミュレーション 実行スクリプト
==============================================
- 流体: D2Q9 LBM
- 標準境界: 上下周期境界、左圧力入口、右圧力出口
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
import traceback
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import numpy as np

# パッケージを src/ から読み込む
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from particulate_flow import FastLBMDEM, LBMDEMSolver, PoreGeometry
from particulate_flow.builder import build_lbm_dem_solver
from particulate_flow.rheology import ViscosityEvaluator
from particulate_flow.io.paths import program_results_dir
from particulate_flow.io.config import SimulationConfig

# ---------------------------------------------------------------------------
# コマンドライン引数
# ---------------------------------------------------------------------------

bootstrap_parser = argparse.ArgumentParser(add_help=False)
bootstrap_parser.add_argument(
    "--config",
    type=str,
    default=None,
    help="JSON file with run defaults. Explicit CLI options override config values.",
)
bootstrap_args, _ = bootstrap_parser.parse_known_args()
CONFIG = (
    SimulationConfig.from_json(bootstrap_args.config)
    if bootstrap_args.config is not None
    else SimulationConfig()
)

parser = argparse.ArgumentParser(
    description="LBM-DEM coupled simulation runner",
    parents=[bootstrap_parser],
)
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
parser.add_argument("--flow-condition",
                    choices=[
                        "fixed-pressure",
                        "target-max-velocity",
                        "constant-pressure",
                        "constant-flux",
                    ],
                    default=None,
                    help="Explicit flow boundary/control condition. Overrides "
                         "--flow-control when set")
parser.add_argument("--flow-control-gain", type=float, default=0.2,
                    help="Relaxation gain for max-speed flow control (default: 0.2)")
parser.add_argument("--y-boundary",
                    choices=["wall", "periodic"],
                    default="periodic",
                    help="Boundary condition in the transverse y direction "
                         "(default: periodic for membrane unit-cell calculations)")
parser.add_argument("--streamwise-boundary",
                    choices=["periodic-force", "pressure"],
                    default="pressure",
                    help="Streamwise x boundary: pressure inlet/outlet, or periodic with body force "
                         "(default: pressure)")
parser.add_argument("--pressure-drop", type=float, default=1e-4,
                    help="Pressure drop p_in-p_out for --streamwise-boundary pressure")
parser.add_argument("--rho-out", type=float, default=1.0,
                    help="Outlet density for pressure inlet/outlet boundary (default: 1.0)")
parser.add_argument("--fluid-method",
                    choices=["lbm-bgk-guo", "lbm-trt-guo"],
                    default="lbm-bgk-guo",
                    help="LBM collision/forcing method (default: lbm-bgk-guo)")
parser.add_argument("--fluid-accelerator",
                    choices=["numpy", "numba", "auto"],
                    default="numpy",
                    help="LBM execution backend. 'numba' uses the compiled step "
                         "when available; default keeps the stable NumPy path")
parser.add_argument("--compute-accelerator",
                    choices=["numpy", "numba", "auto"],
                    default="auto",
                    help="Backend for DEM boundary loads, solid masks, and IBM "
                         "bookkeeping (default: auto)")
parser.add_argument("--particle-method",
                    choices=["dem-hertz", "dem-linear"],
                    default="dem-hertz",
                    help="DEM normal contact model (default: dem-hertz)")
parser.add_argument("--particle-search",
                    choices=["cell_list", "all_pairs"],
                    default="cell_list",
                    help="Particle-pair neighbour search (default: cell_list)")
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
parser.add_argument("--show-ibm-markers", action=argparse.BooleanOptionalAction,
                    default=False,
                    help="Export and overlay IBM Lagrangian marker points when using "
                         "--particle-fluid-coupling immersed_boundary")
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
parser.add_argument("--density-ratio", type=float, default=2.0,
                    help="Particle/fluid density ratio (default: 2.0)")
parser.add_argument("--gravity", type=float, default=3e-5,
                    help="Gravity acceleration in lattice units (default: 3e-5)")
parser.add_argument("--dem-substeps", type=int, default=4,
                    help="DEM substeps per LBM step (default: 4)")
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
parser.add_argument("--snapshot-storage", choices=["disk", "memory", "none"], default="disk",
                    help="Store recorded snapshots on disk, in memory, or not at all "
                         "(default: disk)")
parser.add_argument("--output-profile",
                    choices=["full", "analysis", "minimal"],
                    default="full",
                    help="Output workload preset: full keeps snapshots/video defaults, "
                         "analysis keeps CSV/NPZ/final ParaView and skips video, "
                         "minimal keeps only analysis/final fields (default: full)")
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
parser.add_argument("--surface-roughness", type=float, default=0.0,
                    help="DEM surface roughness h_r; extends contact threshold by h_r [lattice] (default: 0.0)")
parser.add_argument("--viscosity-eval-enabled", action=argparse.BooleanOptionalAction,
                    default=False,
                    help="Enable apparent viscosity evaluation (requires LE shear BC)")
parser.add_argument("--viscosity-eval-start-step", type=int, default=0,
                    help="Step at which to start viscosity accumulation (default: 0)")
parser.add_argument("--viscosity-eval-interval", type=int, default=100,
                    help="Write viscosity CSV row every N steps (default: 100)")
parser.add_argument("--viscosity-eval-average-steps", type=int, default=1000,
                    help="Number of trailing rows for time-averaged η_s (default: 1000)")
parser.add_argument("--porous-resistance", action=argparse.BooleanOptionalAction,
                    default=False,
                    help="Apply Brinkman-like drag in particle-occupied cake cells")
parser.add_argument("--porous-resistance-coeff", type=float, default=0.0,
                    help="Porous drag coefficient for cake-layer resistance (default: 0)")
parser.add_argument("--max-stable-speed", type=float, default=1.0,
                    help="Abort and mark run failed if max fluid or particle speed exceeds this value")
parser.add_argument("--max-stable-pressure", type=float, default=10.0,
                    help="Abort and mark run failed if absolute pressure exceeds this value")
parser.set_defaults(**CONFIG.argparse_defaults())
args = parser.parse_args()
if args.particle_attraction and args.particle_repulsion:
    parser.error("--particle-attraction and --particle-repulsion are mutually exclusive")
if args.n_particles is not None and args.n_particles < 0:
    parser.error("--n-particles must be non-negative")
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
if args.pressure_drop < 0.0:
    parser.error("--pressure-drop must be non-negative")
if args.rho_out <= 0.0:
    parser.error("--rho-out must be positive")
if args.ibm_stiffness <= 0.0:
    parser.error("--ibm-stiffness must be positive")
if args.ibm_marker_spacing <= 0.0:
    parser.error("--ibm-marker-spacing must be positive")
if args.particle_radius <= 0.0:
    parser.error("--particle-radius must be positive")
if args.radius_variation < 0.0:
    parser.error("--radius-variation must be non-negative")
if args.density_ratio <= 0.0:
    parser.error("--density-ratio must be positive")
if args.dem_substeps <= 0:
    parser.error("--dem-substeps must be positive")
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
if args.porous_resistance_coeff < 0.0:
    parser.error("--porous-resistance-coeff must be non-negative")
if args.max_stable_speed <= 0.0:
    parser.error("--max-stable-speed must be positive")
if args.max_stable_pressure <= 0.0:
    parser.error("--max-stable-pressure must be positive")

if args.output_profile == "analysis":
    args.no_video = True
    args.snapshot_storage = "none"
    if args.paraview_every != 0:
        args.paraview_every = 0
elif args.output_profile == "minimal":
    args.no_video = True
    args.paraview_output = False
    args.snapshot_storage = "none"
if args.snapshot_storage == "none" and not args.no_video:
    parser.error("--snapshot-storage none requires --no-video or an output profile that skips video")
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
FLOW_CONDITION  = (
    args.flow_condition
    if args.flow_condition is not None
    else ("target-max-velocity" if args.flow_control else "fixed-pressure")
)
FLOW_CONTROL_MAP = {
    "fixed-pressure": "fixed_pressure",
    "constant-pressure": "constant_pressure",
    "target-max-velocity": "target_max_velocity",
    "constant-flux": "constant_flux",
}
DEFAULT_N_PARTICLES = 40 # 粒子数 (分率未指定時)
RADIUS          = args.particle_radius    # 粒子平均半径 (格子単位)
REYNOLDS_LENGTH = 2.0 * RADIUS  # 代表長さ: 粒子径
RADIUS_VARIATION = args.radius_variation  # 粒子径バリエーション
DENSITY_RATIO   = args.density_ratio    # ρ_p / ρ_f
GRAVITY         = args.gravity   # 重力加速度 (格子単位)
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
GEOMETRY = PoreGeometry.from_cylinders(CYLINDERS)
CYLINDERS = GEOMETRY.as_tuples()
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
            int(
                round(
                    PARTICLE_VOLUME_FRACTION
                    * GEOMETRY.water_area(NX, NY, wall_y=args.y_boundary == "wall")
                    / expected_particle_area
                )
            ),
        )

GEOMETRY_MODE = GEOMETRY.mode_name()
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
RUN_STARTED_AT = datetime.now()


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


def _write_run_status(
    status: str,
    *,
    error: str | None = None,
    traceback_text: str | None = None,
    extra: dict | None = None,
) -> Path:
    """Write a small status file so sweeps can be resumed or audited."""
    now = datetime.now()
    payload = {
        "schema_version": 1,
        "status": status,
        "created_at": RUN_STARTED_AT.isoformat(timespec="seconds"),
        "updated_at": now.isoformat(timespec="seconds"),
        "elapsed_seconds": (now - RUN_STARTED_AT).total_seconds(),
        "program": str(Path(__file__).resolve()),
        "command": sys.argv,
        "output_dir": str(OUT_DIR),
        "git": {
            "commit": _git_value(["rev-parse", "HEAD"]),
            "branch": _git_value(["branch", "--show-current"]),
        },
        "configuration": {
            "nx": NX,
            "ny": NY,
            "reynolds_number": RE,
            "target_u_max": U_MAX,
            "total_steps": TOTAL_STEPS,
            "snapshot_every": SNAPSHOT_EVERY,
            "y_boundary": args.y_boundary,
            "streamwise_boundary": args.streamwise_boundary,
            "pressure_drop": args.pressure_drop,
            "rho_out": args.rho_out,
            "particle_fluid_coupling": args.particle_fluid_coupling,
            "fluid_accelerator": args.fluid_accelerator,
            "compute_accelerator": args.compute_accelerator,
            "particle_volume_fraction": PARTICLE_VOLUME_FRACTION,
            "rolling_friction": args.rolling_friction,
            "particle_attraction": args.particle_attraction,
            "particle_repulsion": args.particle_repulsion,
            "cylinders": [
                {"x": float(cx), "y": float(cy), "radius": float(cr)}
                for cx, cy, cr in CYLINDERS
            ],
        },
    }
    if error:
        payload["error"] = error
    if traceback_text:
        payload["traceback"] = traceback_text
    if extra:
        payload["extra"] = extra
    path = OUT_DIR / "run_status.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path


def _record_uncaught_exception(exc_type, exc_value, exc_tb) -> None:
    """Persist failed run state before Python prints the original traceback."""
    try:
        existing_status = OUT_DIR / "run_status.json"
        if existing_status.exists():
            try:
                payload = json.loads(existing_status.read_text(encoding="utf-8"))
                if str(payload.get("status", "")).startswith("failed"):
                    return
            except (OSError, json.JSONDecodeError):
                pass
        _write_run_status(
            "failed",
            error=f"{exc_type.__name__}: {exc_value}",
            traceback_text="".join(traceback.format_exception(exc_type, exc_value, exc_tb)),
        )
    finally:
        sys.__excepthook__(exc_type, exc_value, exc_tb)


sys.excepthook = _record_uncaught_exception


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
            "flow_condition": FLOW_CONDITION,
            "y_boundary": args.y_boundary,
            "streamwise_boundary": args.streamwise_boundary,
            "pressure_drop": args.pressure_drop,
            "rho_out": args.rho_out,
            "fluid_method": args.fluid_method,
            "fluid_accelerator": args.fluid_accelerator,
            "compute_accelerator": args.compute_accelerator,
            "particle_method": args.particle_method,
            "particle_search": args.particle_search,
            "particle_fluid_coupling": args.particle_fluid_coupling,
            "ibm_stiffness": args.ibm_stiffness,
            "ibm_marker_spacing": args.ibm_marker_spacing,
            "show_ibm_markers": args.show_ibm_markers,
            "particle_radius": RADIUS,
            "radius_variation": RADIUS_VARIATION,
            "density_ratio": DENSITY_RATIO,
            "gravity": GRAVITY,
            "total_steps": TOTAL_STEPS,
            "snapshot_every": SNAPSHOT_EVERY,
            "warmup_steps": args.warmup_steps,
            "paraview_output": args.paraview_output,
            "paraview_every": args.paraview_every,
            "snapshot_storage": args.snapshot_storage,
            "output_profile": args.output_profile,
            "particle_volume_fraction": PARTICLE_VOLUME_FRACTION,
            "source_volume_fraction": SOURCE_VOLUME_FRACTION,
            "n_particles_requested": N_PARTICLES,
            "porous_resistance": args.porous_resistance,
            "porous_resistance_coeff": args.porous_resistance_coeff,
            "max_stable_speed": args.max_stable_speed,
            "max_stable_pressure": args.max_stable_pressure,
            "cylinders": [
                {"x": float(cx), "y": float(cy), "radius": float(cr)}
                for cx, cy, cr in CYLINDERS
            ],
        },
    }
    if sim is not None:
        metadata["solver"] = {
            "fluid_method": sim.fluid_method,
            "fluid_accelerator": sim.fluid_accelerator,
            "uses_numba_lbm": sim.uses_numba_lbm,
            "compute_accelerator": sim.compute_accelerator,
            "uses_numba_compute": sim.uses_numba_compute,
            "particle_method": sim.particle_method,
            "particle_search": sim.particle_search,
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
            "dimensionless_groups": sim.dimensionless_groups(),
        }
    path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")


def _write_standard_run_artifacts(out_dir: Path, sim: LBMDEMSolver | None = None) -> None:
    """Write small, stable files that identify how this run was produced."""
    CONFIG.write_effective(out_dir / "config.json", vars(args))
    git_commit = _git_value(["rev-parse", "HEAD"]) or "unknown"
    (out_dir / "git_commit.txt").write_text(git_commit + "\n", encoding="utf-8")
    environment = {
        "schema_version": 1,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "python": {
            "version": sys.version,
            "executable": sys.executable,
            "platform": platform.platform(),
        },
        "git": {
            "commit": git_commit,
            "branch": _git_value(["branch", "--show-current"]),
            "status_short": _git_value(["status", "--short"]) or "",
        },
        "accelerators": {
            "fluid_accelerator": args.fluid_accelerator,
            "compute_accelerator": args.compute_accelerator,
            "uses_numba_lbm": bool(sim.uses_numba_lbm) if sim is not None else None,
            "uses_numba_compute": bool(sim.uses_numba_compute) if sim is not None else None,
        },
    }
    (out_dir / "environment.json").write_text(
        json.dumps(environment, indent=2, sort_keys=True),
        encoding="utf-8",
    )


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


def _write_ibm_markers_vtk(path: Path, snap: dict) -> None:
    """Write IBM Lagrangian boundary markers as VTK polydata points."""
    marker_x = snap.get("ibm_marker_x", np.empty(0))
    marker_y = snap.get("ibm_marker_y", np.empty(0))
    marker_owner = snap.get("ibm_marker_owner", np.empty(0, dtype=np.int64))
    marker_ds = snap.get("ibm_marker_ds", np.empty(0))
    n_markers = len(marker_x)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write(f"LBM-DEM IBM markers step {snap['step']}\n")
        handle.write("ASCII\n")
        handle.write("DATASET POLYDATA\n")
        handle.write(f"POINTS {n_markers} float\n")
        for x_pos, y_pos in zip(marker_x, marker_y):
            handle.write(f"{_vtk_float(x_pos)} {_vtk_float(y_pos)} 0\n")
        handle.write(f"VERTICES {n_markers} {2 * n_markers}\n")
        for idx in range(n_markers):
            handle.write(f"1 {idx}\n")
        handle.write(f"POINT_DATA {n_markers}\n")
        handle.write("SCALARS particle_id int 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for owner in marker_owner:
            handle.write(f"{int(owner)}\n")
        handle.write("SCALARS marker_spacing float 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for ds in marker_ds:
            handle.write(f"{_vtk_float(ds)}\n")


def _write_cylinders_vtk(
    path: Path,
    cylinders: list[tuple[float, float, float]],
    n_segments: int = 96,
) -> None:
    """Write fixed cylinders as polygonal discs for ParaView."""
    PoreGeometry.from_cylinders(cylinders).write_vtk(path, n_segments=n_segments)


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
        ibm_marker_x=snap.get("ibm_marker_x", np.empty(0)),
        ibm_marker_y=snap.get("ibm_marker_y", np.empty(0)),
        ibm_marker_owner=snap.get("ibm_marker_owner", np.empty(0, dtype=np.int64)),
        ibm_marker_ds=snap.get("ibm_marker_ds", np.empty(0)),
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
            "ibm_marker_x": data["ibm_marker_x"] if "ibm_marker_x" in data else np.empty(0),
            "ibm_marker_y": data["ibm_marker_y"] if "ibm_marker_y" in data else np.empty(0),
            "ibm_marker_owner": data["ibm_marker_owner"] if "ibm_marker_owner" in data else np.empty(0, dtype=np.int64),
            "ibm_marker_ds": data["ibm_marker_ds"] if "ibm_marker_ds" in data else np.empty(0),
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


def _section_mean(field: np.ndarray, solid: np.ndarray, x_idx: int) -> float:
    """Mean scalar value on a clamped vertical section."""
    x_idx = int(np.clip(x_idx, 0, field.shape[0] - 1))
    return _boundary_mean(field, solid, x_idx)


def _pressure_sections(cylinders: list[tuple[float, float, float]]) -> dict[str, int]:
    """Return standard pressure sampling sections for global and pore-local loss."""
    return PoreGeometry.from_cylinders(cylinders).pressure_sections(NX)


def _instability_reason(
    *,
    rho: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    pressure: np.ndarray,
    speed: np.ndarray,
    sim: FastLBMDEM,
) -> str | None:
    """Return a reason string when fields or particle states have diverged."""
    arrays = {
        "rho": rho,
        "ux": ux,
        "uy": uy,
        "pressure": pressure,
        "particle_position": sim.pos,
        "particle_velocity": sim.vel,
        "particle_force": sim.forces_p,
    }
    for name, array in arrays.items():
        if not np.all(np.isfinite(array)):
            return f"non-finite values detected in {name}"
    fluid_mask = ~sim.solid
    max_fluid_speed = float(np.max(speed[fluid_mask])) if np.any(fluid_mask) else 0.0
    if max_fluid_speed > args.max_stable_speed:
        return (
            f"fluid speed exceeded --max-stable-speed "
            f"({max_fluid_speed:.6g} > {args.max_stable_speed:.6g})"
        )
    max_pressure = float(np.max(np.abs(pressure[fluid_mask]))) if np.any(fluid_mask) else 0.0
    if max_pressure > args.max_stable_pressure:
        return (
            f"pressure exceeded --max-stable-pressure "
            f"({max_pressure:.6g} > {args.max_stable_pressure:.6g})"
        )
    if sim.n_p:
        particle_speed = np.linalg.norm(sim.vel, axis=1)
        max_particle_speed = float(np.max(particle_speed))
        if max_particle_speed > args.max_stable_speed:
            return (
                f"particle speed exceeded --max-stable-speed "
                f"({max_particle_speed:.6g} > {args.max_stable_speed:.6g})"
            )
    return None


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


def _write_fouling_summary(
    out_dir: Path,
    rows: list[dict[str, float | int]],
    sim: FastLBMDEM,
    analysis_csv: Path,
    analysis_npz: Path,
    *,
    apparent_viscosity: dict | None = None,
) -> tuple[Path, Path]:
    """Write compact final metrics for comparing large condition sweeps."""
    summary_dir = out_dir / "analysis"
    summary_dir.mkdir(parents=True, exist_ok=True)
    json_path = summary_dir / "summary.json"
    md_path = summary_dir / "summary.md"
    if not rows:
        summary = {"schema_version": 1, "status": "no_rows"}
    else:
        final = rows[-1]
        max_pressure_drop = max(float(row["pressure_drop"]) for row in rows)
        max_local_pressure_drop = max(float(row["local_pressure_drop"]) for row in rows)
        min_normalized_flux = min(float(row["normalized_permeate_flux"]) for row in rows)
        max_contacts = max(int(row["total_contacts"]) for row in rows)
        max_cylinder_contacts = max(int(row["cylinder_contacts"]) for row in rows)
        summary = {
            "schema_version": 1,
            "status": "completed",
            "analysis_csv": str(analysis_csv),
            "analysis_npz": str(analysis_npz),
            "final": final,
            "extrema": {
                "max_pressure_drop": max_pressure_drop,
                "max_local_pressure_drop": max_local_pressure_drop,
                "min_normalized_permeate_flux": min_normalized_flux,
                "max_total_contacts": max_contacts,
                "max_cylinder_contacts": max_cylinder_contacts,
            },
            "fouling_indicators": {
                "final_retained_particle_ratio": final["retained_particle_ratio"],
                "final_passed_particle_ratio": final["passed_particle_ratio"],
                "final_normalized_permeate_flux": final["normalized_permeate_flux"],
                "final_pressure_drop_ratio": final["pressure_drop_ratio"],
                "final_local_pressure_drop_ratio": final["local_pressure_drop_ratio"],
                "final_fouling_resistance_index": final["fouling_resistance_index"],
                "final_cylinder_contacts": final["cylinder_contacts"],
                "final_effective_inlet_particle_fraction": final["effective_inlet_particle_fraction"],
                "final_inlet_particle_fraction_error": final["inlet_particle_fraction_error"],
                "max_total_contacts": max_contacts,
            },
            "solver": {
                "fluid_accelerator": sim.fluid_accelerator,
                "uses_numba_lbm": sim.uses_numba_lbm,
                "compute_accelerator": sim.compute_accelerator,
                "uses_numba_compute": sim.uses_numba_compute,
                "particle_fluid_coupling": sim.particle_fluid_coupling,
            },
        }
    if apparent_viscosity is not None:
        summary["apparent_viscosity"] = apparent_viscosity
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    lines = ["# LBM-DEM Fouling Summary", ""]
    if rows:
        indicators = summary["fouling_indicators"]
        lines.extend(
            [
                f"- Final normalized permeate flux: {indicators['final_normalized_permeate_flux']:.6g}",
                f"- Final pressure drop ratio: {indicators['final_pressure_drop_ratio']:.6g}",
                f"- Final local pressure drop ratio: {indicators['final_local_pressure_drop_ratio']:.6g}",
                f"- Fouling resistance index: {indicators['final_fouling_resistance_index']:.6g}",
                f"- Retained particle ratio: {indicators['final_retained_particle_ratio']:.6g}",
                f"- Passed particle ratio: {indicators['final_passed_particle_ratio']:.6g}",
                f"- Effective inlet particle fraction: {indicators['final_effective_inlet_particle_fraction']:.6g}",
                f"- Inlet particle fraction error: {indicators['final_inlet_particle_fraction_error']:.6g}",
                f"- Final cylinder contacts: {indicators['final_cylinder_contacts']}",
                f"- Max total contacts: {indicators['max_total_contacts']}",
                "",
                f"Time-series CSV: `{analysis_csv}`",
                f"Time-series NPZ: `{analysis_npz}`",
            ]
        )
    else:
        lines.append("No time-series rows were recorded.")
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, md_path

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
        f"target max velocity / constant flux proxy (gain={args.flow_control_gain:.2f})"
        if FLOW_CONDITION in {"target-max-velocity", "constant-flux"}
        else "fixed pressure/body force"
    )
)
print(f"Fluid method: {args.fluid_method}")
print(f"Boundary condition: y={args.y_boundary}, x={args.streamwise_boundary}")
print(f"Fluid accelerator: {args.fluid_accelerator}")
print(f"Compute accelerator: {args.compute_accelerator}")
print(f"Particle method: {args.particle_method}")
print(f"Particle search: {args.particle_search}")
print(f"Particle-fluid coupling: {args.particle_fluid_coupling}")
print(
    "Porous resistance: "
    + (f"on (coeff={args.porous_resistance_coeff:g})" if args.porous_resistance else "off")
)
print(f"Output profile: {args.output_profile}")
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

sim = build_lbm_dem_solver(args)
viscosity_evaluator = ViscosityEvaluator(
    sim,
    start_step=getattr(args, "viscosity_eval_start_step", 0),
    viscosity_interval=getattr(args, "viscosity_eval_interval", 100),
    average_steps=getattr(args, "viscosity_eval_average_steps", 1000),
    out_dir=OUT_DIR,
    enabled=getattr(args, "viscosity_eval_enabled", False),
)
# Sync derived values with what the builder actually computed to avoid divergence.
N_PARTICLES = sim.n_p
SOURCE_VOLUME_FRACTION = sim.source_volume_fraction
metadata_path = OUT_DIR / "metadata.json"
_write_metadata(metadata_path, sim)
_write_standard_run_artifacts(OUT_DIR, sim)
_write_run_status("running", extra={"metadata": str(metadata_path)})
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
elif args.snapshot_storage == "none":
    print("      snapshots: 中間NPZ保存なし")

paraview_dir = None
fluid_entries: list[tuple[int, Path]] = []
particle_entries: list[tuple[int, Path]] = []
marker_entries: list[tuple[int, Path]] = []
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
reference_local_pressure_drop: float | None = None
pressure_sections = _pressure_sections(CYLINDERS)

t1 = time.perf_counter()
for frame_idx in range(n_frames):
    sim.advance(SNAPSHOT_EVERY)
    viscosity_evaluator.record(sim.step_count)

    rho, ux, uy = sim.get_fields()

    # 各粒子に加わる合力の大きさ (重力 + 抗力 + 接触力)
    total_force_mags = np.linalg.norm(sim.forces_p, axis=1)
    pressure = rho / 3.0
    pressure_gauge = pressure - float(np.mean(pressure[~sim.solid]))
    speed = np.sqrt(ux**2 + uy**2)
    fluid_mask = ~sim.solid
    instability = _instability_reason(
        rho=rho,
        ux=ux,
        uy=uy,
        pressure=pressure,
        speed=speed,
        sim=sim,
    )
    if instability is not None:
        _write_run_status(
            "failed_unstable",
            error=instability,
            extra={"frame": frame_idx + 1, "step": sim.step_count},
        )
        raise RuntimeError(instability)
    inlet_flux = _boundary_flux(ux, sim.solid, pressure_sections["inlet_x"])
    outlet_flux = _boundary_flux(ux, sim.solid, pressure_sections["outlet_x"])
    inlet_pressure = _section_mean(pressure, sim.solid, pressure_sections["inlet_x"])
    outlet_pressure = _section_mean(pressure, sim.solid, pressure_sections["outlet_x"])
    pore_upstream_pressure = _section_mean(
        pressure,
        sim.solid,
        pressure_sections["pore_upstream_x"],
    )
    pore_downstream_pressure = _section_mean(
        pressure,
        sim.solid,
        pressure_sections["pore_downstream_x"],
    )
    pressure_drop = inlet_pressure - outlet_pressure
    local_pressure_drop = pore_upstream_pressure - pore_downstream_pressure
    if reference_outlet_flux is None and outlet_flux > 1e-12:
        reference_outlet_flux = outlet_flux
    if reference_pressure_drop is None and np.isfinite(pressure_drop):
        reference_pressure_drop = pressure_drop
    if reference_local_pressure_drop is None and np.isfinite(local_pressure_drop):
        reference_local_pressure_drop = local_pressure_drop
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
    local_pressure_drop_ratio = (
        local_pressure_drop / reference_local_pressure_drop
        if reference_local_pressure_drop is not None and abs(reference_local_pressure_drop) > 1e-12
        else 1.0
    )
    inlet_outlet_flux_ratio = outlet_flux / inlet_flux if inlet_flux > 1e-12 else 0.0
    dynamic_solid_fraction = sim.dynamic_solid_fraction()
    porous_fraction = sim.porous_resistance_fraction()
    total_solid_fraction = float(np.count_nonzero(sim.solid) / sim.solid.size)
    dimensionless = sim.dimensionless_groups(float(speed[fluid_mask].max()))
    effective_inlet_particle_fraction = (
        sim.injected_particle_area / sim.cumulative_inlet_flow_area
        if sim.cumulative_inlet_flow_area > 1e-12
        else 0.0
    )
    target_source_fraction = SOURCE_VOLUME_FRACTION if SOURCE_VOLUME_FRACTION is not None else 0.0

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
    if args.show_ibm_markers:
        markers = sim.ibm_marker_points()
        snap["ibm_marker_x"] = markers["x"].copy()
        snap["ibm_marker_y"] = markers["y"].copy()
        snap["ibm_marker_owner"] = markers["owner"].copy()
        snap["ibm_marker_ds"] = markers["ds"].copy()
    else:
        snap["ibm_marker_x"] = np.empty(0)
        snap["ibm_marker_y"] = np.empty(0)
        snap["ibm_marker_owner"] = np.empty(0, dtype=np.int64)
        snap["ibm_marker_ds"] = np.empty(0)
    analysis_rows.append({
        "frame": frame_idx + 1,
        "step": sim.step_count,
        "time_lattice": float(sim.step_count),
        "active_particles": sim.n_p,
        "generated_particles": sim.generated_particles,
        "pending_particles": int(len(sim._pending_radii)),
        "passed_particles": sim.removed_particles,
        "inlet_x": pressure_sections["inlet_x"],
        "outlet_x": pressure_sections["outlet_x"],
        "pore_upstream_x": pressure_sections["pore_upstream_x"],
        "pore_downstream_x": pressure_sections["pore_downstream_x"],
        "inlet_flux": inlet_flux,
        "outlet_flux": outlet_flux,
        "permeate_flux": outlet_flux,
        "cumulative_inlet_flow_area": sim.cumulative_inlet_flow_area,
        "target_source_particle_fraction": target_source_fraction,
        "effective_inlet_particle_fraction": effective_inlet_particle_fraction,
        "injected_particle_area": sim.injected_particle_area,
        "inlet_particle_area_budget": sim.inlet_particle_area_budget,
        "inlet_particle_fraction_error": effective_inlet_particle_fraction - target_source_fraction,
        "reynolds_number": RE,
        "observed_reynolds_number": dimensionless["observed_reynolds_number"],
        "particle_reynolds_number": dimensionless["particle_reynolds_number"],
        "stokes_number_estimate": dimensionless["stokes_number_estimate"],
        "particle_to_length_ratio": dimensionless["particle_to_length_ratio"],
        "brinkman_resistance_number": dimensionless["brinkman_resistance_number"],
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
        "pore_upstream_pressure": pore_upstream_pressure,
        "pore_downstream_pressure": pore_downstream_pressure,
        "local_pressure_drop": local_pressure_drop,
        "local_pressure_drop_ratio": local_pressure_drop_ratio,
        "mean_pressure_gauge": float(np.mean(pressure_gauge[fluid_mask])),
        "particle_area_fraction": float(np.sum(np.pi * sim.radii**2) / sim.fluid_area),
        "retained_particle_ratio": retained_ratio,
        "passed_particle_ratio": pass_ratio,
        "normalized_permeate_flux": normalized_permeate_flux,
        "inlet_outlet_flux_ratio": inlet_outlet_flux_ratio,
        "dynamic_particle_solid_fraction": dynamic_solid_fraction,
        "porous_resistance_fraction": porous_fraction,
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
    elif args.snapshot_storage == "disk":
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
        if args.show_ibm_markers:
            marker_name = f"ibm_markers_{frame_idx:04d}_step_{step:07d}.vtk"
            _write_ibm_markers_vtk(paraview_dir / marker_name, snap)
            marker_entries.append((step, Path(marker_name)))

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
viscosity_result = viscosity_evaluator.finalize()
summary_json, summary_md = _write_fouling_summary(
    OUT_DIR,
    analysis_rows,
    sim,
    analysis_csv,
    analysis_npz,
    apparent_viscosity=viscosity_result if getattr(args, "viscosity_eval_enabled", False) else None,
)
print(f"\n後解析データ保存: {analysis_csv}")
print(f"後解析NPZ保存  : {analysis_npz}")
print(f"サマリ保存      : {summary_json}")

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
        if args.show_ibm_markers:
            marker_name = f"ibm_markers_{frame_idx:04d}_step_{step:07d}.vtk"
            _write_ibm_markers_vtk(paraview_dir / marker_name, last)
            marker_entries.append((step, Path(marker_name)))
    _write_pvd(paraview_dir / "fluid_series.pvd", fluid_entries)
    _write_pvd(paraview_dir / "particles_series.pvd", particle_entries)
    if args.show_ibm_markers:
        _write_pvd(paraview_dir / "ibm_markers_series.pvd", marker_entries)
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
    print(f"  累積入口流量    : {sim.cumulative_inlet_flow_area:.4e} (格子面積)")
    print(f"  投入済み粒子面積: {sim.injected_particle_area:.4e}")
    print(f"  未投入面積予算  : {sim.inlet_particle_area_budget:.4e}")
    if sim.cumulative_inlet_flow_area > 1e-12:
        print(
            f"  実効投入分率    : "
            f"{sim.injected_particle_area / sim.cumulative_inlet_flow_area:.4%}"
        )

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
if args.show_ibm_markers and len(snap["ibm_marker_x"]):
    ax.scatter(
        snap["ibm_marker_x"], snap["ibm_marker_y"],
        s=10, c="#00ff88", marker=".", zorder=7, label="IBM markers",
    )
_add_cylinder_patch(ax)
ax.set_title("Velocity magnitude |u|")
ax.set_xlabel("x [lattice units]"); ax.set_ylabel("y [lattice units]")
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
if args.show_ibm_markers and len(snap["ibm_marker_x"]):
    ax.scatter(
        snap["ibm_marker_x"], snap["ibm_marker_y"],
        s=10, c="#00ff88", marker=".", zorder=7,
    )
_add_cylinder_patch(ax)
ax.set_xlim(0, NX); ax.set_ylim(0, NY)
ax.set_title("Streamlines"); ax.set_xlabel("x [lattice units]")

ax = axes[2]
if len(snap["pos"]):
    sc = ax.scatter(snap["pos"][:,0], snap["pos"][:,1],
                    c=snap["total_force"], cmap="plasma", s=(snap["radii"]*4)**2,
                    edgecolors="k", lw=0.5, zorder=5)
else:
    sc = ax.scatter([], [], c=[], cmap="plasma", edgecolors="k", lw=0.5, zorder=5)
if args.show_ibm_markers and len(snap["ibm_marker_x"]):
    ax.scatter(
        snap["ibm_marker_x"], snap["ibm_marker_y"],
        s=12, c="#00ff88", marker=".", zorder=7,
    )
ax.imshow(snap["speed"].T, origin="lower", cmap="Blues",
          extent=[0, NX, 0, NY], aspect="auto", alpha=0.5)
_add_cylinder_patch(ax)
fig_stat.colorbar(sc, ax=ax, label="Total force |F_total| [lattice units]", shrink=0.7)
ax.set_xlim(0, NX); ax.set_ylim(0, NY)
ax.set_title("Particle positions colored by total force")
ax.set_xlabel("x [lattice units]")

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
    _write_run_status(
        "completed",
        extra={
            "metadata": str(metadata_path),
            "static_image": str(static_path),
            "analysis_csv": str(analysis_csv),
            "analysis_npz": str(analysis_npz),
            "summary_json": str(summary_json),
            "summary_md": str(summary_md),
            "paraview_dir": str(paraview_dir) if paraview_dir is not None else None,
            "video": None,
        },
    )
    print(f"\n出力ファイル:")
    print(f"  静止画: {static_path}")
    print(f"  時系列CSV: {analysis_csv}")
    print(f"  時系列NPZ: {analysis_npz}")
    if paraview_dir is not None:
        print(f"  ParaView: {paraview_dir}")
        print(f"    流体時系列: {paraview_dir / 'fluid_series.pvd'}")
        print(f"    粒子時系列: {paraview_dir / 'particles_series.pvd'}")
        if args.show_ibm_markers:
            print(f"    IBMマーカー: {paraview_dir / 'ibm_markers_series.pvd'}")
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
    if not snapshot_paths:
        raise RuntimeError("Snapshots were not stored; rerun without --snapshot-storage none")
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
cbar_fluid.set_label("|u| [lattice units]", color="white")
cbar_fluid.ax.yaxis.set_tick_params(color="white")
plt.setp(cbar_fluid.ax.yaxis.get_ticklabels(), color="white")

# 合力カラーバー（ScalarMappable で追加）
sm_force = plt.cm.ScalarMappable(cmap=force_cmap, norm=force_norm)
sm_force.set_array([])
cbar_force = fig_anim.colorbar(sm_force, ax=ax_anim, shrink=0.75, pad=0.12)
cbar_force.set_label("|F_total| [lattice units]", color="white")
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

if args.show_ibm_markers:
    marker_offsets = (
        np.column_stack([snap0["ibm_marker_x"], snap0["ibm_marker_y"]])
        if len(snap0["ibm_marker_x"])
        else np.empty((0, 2))
    )
    ibm_marker_scatter = ax_anim.scatter(
        marker_offsets[:, 0] if len(marker_offsets) else [],
        marker_offsets[:, 1] if len(marker_offsets) else [],
        s=12,
        c="#00ff88",
        marker=".",
        alpha=0.95,
        zorder=7,
        animated=True,
    )
else:
    ibm_marker_scatter = None

title_txt = ax_anim.set_title(
    f"LBM-DEM  step={snap0['step']:,}", color="white", fontsize=11
)
ax_anim.set_xlabel("x [lattice units]", color="white")
ax_anim.set_ylabel("y [lattice units]", color="white")
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
    if ibm_marker_scatter is not None:
        marker_offsets = (
            np.column_stack([snap["ibm_marker_x"], snap["ibm_marker_y"]])
            if len(snap["ibm_marker_x"])
            else np.empty((0, 2))
        )
        ibm_marker_scatter.set_offsets(marker_offsets)
    p_ke = 0.5 * float(np.sum(snap["masses"] * np.sum(snap["vel"] ** 2, axis=1)))
    force_mean = float(snap["total_force"].mean()) if len(snap["total_force"]) else 0.0
    title_txt.set_text(
        f"LBM-DEM  Re={RE:.0f}  step={snap['step']:,}"
        f"  KE_p={p_ke:.3e}  |F_total|_mean={force_mean:.3e}"
    )
    artists = [im_fluid, title_txt] + particle_circles
    if ibm_marker_scatter is not None:
        artists.append(ibm_marker_scatter)
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
_write_run_status(
    "completed",
    extra={
        "metadata": str(metadata_path),
        "static_image": str(static_path),
        "analysis_csv": str(analysis_csv),
        "analysis_npz": str(analysis_npz),
        "summary_json": str(summary_json),
        "summary_md": str(summary_md),
        "paraview_dir": str(paraview_dir) if paraview_dir is not None else None,
        "video": str(video_path),
    },
)

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
    if args.show_ibm_markers:
        print(f"    IBMマーカー: {paraview_dir / 'ibm_markers_series.pvd'}")
    print(f"    固定円柱  : {paraview_dir / 'cylinders.vtk'}")
print(f"\n完了!")
