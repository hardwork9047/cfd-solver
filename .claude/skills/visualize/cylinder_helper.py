#!/usr/bin/env python3
"""
Analytic cylinder mesh helper for ParaView (pvpython) visualization.

The LBM-DEM 3D solver voxelizes membrane pores as cylinders on a coarse grid
(e.g. radius=3.0 LU on a 48x24x24 domain). Extracting those voxels from the
`solid` VTK field produces a blocky, square-looking surface. This module
instead builds *analytic* smooth cylinders directly from the geometry JSON
(center x, y + radius), matching the solver's `_add_cylinder` convention:
z-aligned cylinders spanning the full domain depth.

Use `add_analytic_cylinders(view, geometry_json, nz, ...)` inside a pvpython
visualization script in place of the threshold/extract-surface approach.
"""

import json
from pathlib import Path

try:
    from paraview.simple import Cylinder, Transform, Show
except ImportError:  # allow import outside pvpython for unit checks
    Cylinder = Transform = Show = None


def load_cylinder_specs(geometry_json: str | Path) -> list[dict]:
    """
    Load cylinder specs (x, y, radius) from a geometry fragment JSON.

    Args:
        geometry_json: Path to e.g. configs/lbm_dem/geometries/four_cylinder_3d.json

    Returns:
        List of dicts, each with keys 'x', 'y', 'radius'.

    The solver treats each cylinder as z-aligned (axis parallel to z), centered
    at (x, y), spanning the full domain depth. Only x, y, radius are read here;
    the z extent comes from the domain nz passed to add_analytic_cylinders.
    """
    with open(geometry_json) as f:
        data = json.load(f)
    return list(data["geometry"]["cylinders"])


def add_analytic_cylinders(
    view,
    geometry_json: str | Path,
    nz: int,
    opacity: float = 0.55,
    color: tuple[float, float, float] = (0.25, 0.55, 0.95),
    resolution: int = 48,
):
    """
    Build smooth analytic z-aligned cylinders and show them in the view.

    This replaces the voxel-based Threshold(solid)->ExtractSurface approach,
    which looks square on coarse grids. Each cylinder is a ParaView Cylinder
    source (radius from geometry, height = nz) rotated so its axis is the
    z-axis and translated to (x, y, nz/2), matching the solver convention in
    geometry3d.py `_add_cylinder` (z-aligned, full depth).

    Args:
        view: The render view to show the cylinders in.
        geometry_json: Path to the geometry fragment JSON (e.g. four_cylinder_3d.json).
        nz: Domain depth in lattice units (cylinder height / z-span).
        opacity: Surface opacity (0-1). 0.55 lets particles behind show through.
        color: RGB tuple (0-1) for the cylinder surface.
        resolution: Number of facets around the cylinder. 48 looks smooth.

    Returns:
        List of display proxies (one per cylinder), so the caller can further
        tweak them if desired.
    """
    if Cylinder is None:
        raise RuntimeError("paraview.simple not available — run under pvpython")

    specs = load_cylinder_specs(geometry_json)
    displays = []
    for spec in specs:
        cyl = Cylinder()
        cyl.Radius = float(spec["radius"])
        cyl.Height = float(nz)
        cyl.Resolution = resolution
        # ParaView's Cylinder axis is +Y by default. The solver's cylinders are
        # z-aligned, so rotate -90 deg about X to point the axis along +Z, then
        # translate to the cell center (x, y, nz/2).
        tf = Transform(Input=cyl)
        tf.Transform = "Transform"
        tf.Transform.Rotate = [-90.0, 0.0, 0.0]
        tf.Transform.Translate = [float(spec["x"]), float(spec["y"]), float(nz) / 2.0]

        disp = Show(tf, view)
        disp.Representation = "Surface"
        disp.AmbientColor = list(color)
        disp.DiffuseColor = list(color)
        disp.Opacity = opacity
        disp.Specular = 0.4
        disp.SpecularPower = 30
        displays.append(disp)

    return displays


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python cylinder_helper.py <geometry_json>")
        sys.exit(1)
    specs = load_cylinder_specs(sys.argv[1])
    print(f"Loaded {len(specs)} cylinders from {sys.argv[1]}:")
    for s in specs:
        print(f"  center=({s['x']}, {s['y']}), radius={s['radius']}")
