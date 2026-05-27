"""Geometry definitions for pore-scale LBM-DEM simulations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


def _vtk_float(value: float) -> str:
    return f"{float(value):.9g}"


@dataclass(frozen=True)
class Cylinder:
    """A fixed circular solid obstacle in lattice coordinates."""

    x: float
    y: float
    radius: float

    def __post_init__(self) -> None:
        if self.radius <= 0.0:
            raise ValueError("cylinder radius must be positive")

    @classmethod
    def from_tuple(cls, value: Iterable[float]) -> "Cylinder":
        items = tuple(value)
        if len(items) != 3:
            raise ValueError("each cylinder must be (x, y, radius)")
        return cls(float(items[0]), float(items[1]), float(items[2]))

    def as_tuple(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.radius)


@dataclass(frozen=True)
class PoreGeometry:
    """Pore geometry shared by solvers, runners, and result writers."""

    cylinders: tuple[Cylinder, ...] = ()

    @classmethod
    def from_cylinders(
        cls,
        cylinders: Iterable[Cylinder | Iterable[float]] | None = None,
    ) -> "PoreGeometry":
        if cylinders is None:
            return cls()
        normalised = []
        for item in cylinders:
            if isinstance(item, Cylinder):
                normalised.append(item)
            else:
                normalised.append(Cylinder.from_tuple(item))
        return cls(tuple(normalised))

    @classmethod
    def from_legacy_inputs(
        cls,
        *,
        cylinder: Iterable[float] | None = None,
        cylinders: Iterable[Cylinder | Iterable[float]] | None = None,
    ) -> "PoreGeometry":
        items: list[Cylinder | Iterable[float]] = []
        if cylinder is not None:
            items.append(cylinder)
        if cylinders is not None:
            items.extend(cylinders)
        return cls.from_cylinders(items)

    def as_tuples(self) -> list[tuple[float, float, float]]:
        return [cylinder.as_tuple() for cylinder in self.cylinders]

    def is_empty(self) -> bool:
        return len(self.cylinders) == 0

    def mode_name(self) -> str:
        if len(self.cylinders) > 1:
            return "multi_cylinder"
        if len(self.cylinders) == 1:
            return "cylinder"
        return "channel"

    def water_area(self, nx: int, ny: int, *, wall_y: bool = False) -> float:
        """Approximate available 2-D water area in lattice units."""
        height = ny - 2 if wall_y else ny
        area = float(nx * max(height, 1))
        for cylinder in self.cylinders:
            area -= float(np.pi * cylinder.radius**2)
        return max(area, 1.0)

    def cylinder_solid_mask(self, nx: int, ny: int, *, y_boundary: str) -> np.ndarray:
        """Return solid lattice cells occupied by fixed cylinders."""
        solid = np.zeros((nx, ny), dtype=bool)
        if not self.cylinders:
            return solid
        ix = np.arange(nx)
        iy = np.arange(ny)
        xx, yy = np.meshgrid(ix, iy, indexing="ij")
        for cylinder in self.cylinders:
            dy = yy - cylinder.y
            if y_boundary == "periodic":
                dy = dy - ny * np.rint(dy / ny)
            solid |= (xx - cylinder.x) ** 2 + dy**2 <= cylinder.radius**2
        return solid

    def pressure_sections(self, nx: int) -> dict[str, int]:
        """Return global and pore-local pressure sampling sections."""
        if self.cylinders:
            upstream_edge = min(cylinder.x - cylinder.radius for cylinder in self.cylinders)
            downstream_edge = max(cylinder.x + cylinder.radius for cylinder in self.cylinders)
            margin = max(2.0, 0.05 * nx)
            upstream = int(np.floor(upstream_edge - margin))
            downstream = int(np.ceil(downstream_edge + margin))
        else:
            upstream = int(0.25 * nx)
            downstream = int(0.75 * nx)
        return {
            "inlet_x": 1,
            "outlet_x": nx - 2,
            "pore_upstream_x": int(np.clip(upstream, 1, nx - 2)),
            "pore_downstream_x": int(np.clip(downstream, 1, nx - 2)),
        }

    def write_vtk(self, path: str | Path, *, n_segments: int = 96) -> None:
        """Write fixed cylinders as polygonal discs for ParaView."""
        target = Path(path)
        points: list[tuple[float, float, float]] = []
        polygons: list[list[int]] = []
        for cylinder in self.cylinders:
            start = len(points)
            polygon = []
            for seg in range(n_segments):
                theta = 2.0 * np.pi * seg / n_segments
                points.append(
                    (
                        cylinder.x + cylinder.radius * np.cos(theta),
                        cylinder.y + cylinder.radius * np.sin(theta),
                        0.0,
                    )
                )
                polygon.append(start + seg)
            polygons.append(polygon)

        with target.open("w", encoding="utf-8") as handle:
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
            for cylinder in self.cylinders:
                handle.write(f"{_vtk_float(cylinder.radius)}\n")
