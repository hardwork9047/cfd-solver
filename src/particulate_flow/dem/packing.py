"""Gravity packing calculations using the shared DEM contact solver."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .solver import DEMSolver
from ..geometry.pore import PoreGeometry


def _vtk_float(value: float) -> str:
    """Format a float compactly for ASCII VTK."""
    return f"{float(value):.9g}"


@dataclass
class PackingMetrics:
    """Compact time-series metrics for a DEM packing calculation."""

    step: int
    n_particles: int
    kinetic_energy: float
    mean_speed: float
    max_speed: float
    contact_count: int
    bed_height: float
    packing_fraction: float


class DEMPackingSimulation:
    """DEM-only gravity settling and packing simulation.

    The class intentionally exposes the attribute names expected by
    :class:`particulate_flow.dem_solver.DEMSolver`, so packing calculations and
    coupled LBM-DEM calculations share the same particle contact, rolling
    friction, and cylinder surface-force implementation.
    """

    def __init__(
        self,
        *,
        nx: int = 120,
        ny: int = 220,
        n_particles: int = 200,
        particle_radius: float = 2.0,
        radius_variation: float = 0.10,
        density_ratio: float = 2.0,
        gravity: float = 5e-4,
        k_n: float = 120.0,
        damping: float = 0.8,
        linear_damping: float = 0.03,
        dem_substeps: int = 5,
        seed: int = 42,
        rolling_friction: bool = True,
        sliding_friction: float = 0.5,
        tangential_damping: float = 0.5,
        rolling_friction_coeff: float = 0.10,
        rolling_damping: float = 0.35,
        particle_attraction: bool = False,
        particle_repulsion: bool = False,
        attraction_strength: float = 1e-3,
        repulsion_strength: float = 1e-3,
        attraction_cutoff: float = 3.0,
        repulsion_cutoff: float = 3.0,
        attraction_min_gap: float = 0.05,
        repulsion_min_gap: float = 0.05,
        cylinders: Iterable[tuple[float, float, float]] | None = None,
        particle_method: str = "dem-hertz",
        particle_search: str = "cell_list",
    ):
        if particle_attraction and particle_repulsion:
            raise ValueError("particle_attraction and particle_repulsion are mutually exclusive")
        if particle_search not in {"cell_list", "all_pairs"}:
            raise ValueError("particle_search must be 'cell_list' or 'all_pairs'")
        if n_particles < 0:
            raise ValueError("n_particles must be non-negative")
        if nx <= 0 or ny <= 0:
            raise ValueError("nx and ny must be positive")
        if particle_radius <= 0.0:
            raise ValueError("particle_radius must be positive")

        self.nx = int(nx)
        self.ny = int(ny)
        self.n_p = int(n_particles)
        self.r_p = float(particle_radius)
        self.radius_variation = float(radius_variation)
        self.density_ratio = float(density_ratio)
        self.g = float(gravity)
        self.k_n = float(k_n)
        self.damping = float(damping)
        self.linear_damping = float(linear_damping)
        self.dem_substeps = int(dem_substeps)
        self.y_boundary = "wall"
        self.particle_fluid_coupling = "point_force"
        self.uses_numba_compute = False
        self.particle_search = particle_search
        self.rolling_friction = bool(rolling_friction)
        self.sliding_friction = float(sliding_friction)
        self.tangential_damping = float(tangential_damping)
        self.rolling_friction_coeff = float(rolling_friction_coeff)
        self.rolling_damping = float(rolling_damping)
        self.particle_attraction = bool(particle_attraction)
        self.particle_repulsion = bool(particle_repulsion)
        self.attraction_strength = float(attraction_strength)
        self.repulsion_strength = float(repulsion_strength)
        self.attraction_cutoff = float(attraction_cutoff)
        self.repulsion_cutoff = float(repulsion_cutoff)
        self.attraction_min_gap = float(attraction_min_gap)
        self.repulsion_min_gap = float(repulsion_min_gap)
        self.geometry = PoreGeometry.from_cylinders(cylinders or [])
        self.cylinders = self.geometry.as_tuples()
        self.step_count = 0

        rng = np.random.default_rng(seed)
        if radius_variation > 0.0:
            self.radii = particle_radius * (
                1.0 + rng.uniform(-radius_variation, radius_variation, self.n_p)
            )
        else:
            self.radii = np.full(self.n_p, particle_radius)
        self.masses = self.density_ratio * np.pi * self.radii**2
        self.inertias = 0.5 * self.masses * self.radii**2
        self.pos = np.zeros((self.n_p, 2))
        self.vel = np.zeros((self.n_p, 2))
        self.omega_p = np.zeros(self.n_p)
        self.forces_p = np.zeros((self.n_p, 2))
        self.torques_p = np.zeros(self.n_p)
        self.ibm_forces_p = np.zeros((self.n_p, 2))
        self.ibm_torques_p = np.zeros(self.n_p)

        self._init_loose_cloud(rng)
        self.dem_solver = DEMSolver(self, contact_model=particle_method)

    def _periodic_y_delta(self, dy: float) -> float:
        """Compatibility hook for the shared DEM solver."""
        return float(dy)

    def _particle_drag_forces(self) -> np.ndarray:
        """Return numerical settling drag, not fluid drag."""
        return -self.linear_damping * self.masses[:, np.newaxis] * self.vel

    def _init_loose_cloud(self, rng: np.random.Generator) -> None:
        """Place particles in a loose, non-overlapping cloud above the bottom."""
        if self.n_p == 0:
            return

        order = np.argsort(-self.radii)
        placed_positions = np.zeros_like(self.pos)
        placed_radii = np.zeros_like(self.radii)
        placed = 0
        max_r = float(np.max(self.radii))
        y_min = max(self.ny * 0.35, max_r + 1.0)

        for idx in order:
            radius = float(self.radii[idx])
            found = False
            for _ in range(5000):
                x = rng.uniform(radius + 1.0, self.nx - radius - 1.0)
                y = rng.uniform(y_min + radius, self.ny - radius - 1.0)
                if self._overlaps_fixed_solid(x, y, radius, clearance=1.05):
                    continue
                if placed:
                    d = np.hypot(
                        placed_positions[:placed, 0] - x,
                        placed_positions[:placed, 1] - y,
                    )
                    if np.any(d < 1.05 * (placed_radii[:placed] + radius)):
                        continue
                found = True
                break
            if not found:
                x, y = self._fallback_lattice_position(placed, radius, max_r)
            placed_positions[placed] = [x, y]
            placed_radii[placed] = radius
            self.pos[idx] = [x, y]
            placed += 1

    def _fallback_lattice_position(self, index: int, radius: float, max_r: float) -> tuple[float, float]:
        """Deterministic fallback placement for crowded initial clouds."""
        dx = 2.2 * max_r
        dy = np.sqrt(3.0) * max_r
        columns = max(1, int((self.nx - 2.0 * max_r - 2.0) // dx))
        row = index // columns
        col = index % columns
        x = max_r + 1.0 + col * dx + (0.5 * dx if row % 2 else 0.0)
        x = min(max(x, radius + 1.0), self.nx - radius - 1.0)
        y = self.ny - radius - 1.0 - row * dy
        y = min(max(y, radius + 1.0), self.ny - radius - 1.0)
        return x, y

    def _overlaps_fixed_solid(self, x: float, y: float, radius: float, clearance: float) -> bool:
        """Return whether a candidate particle overlaps walls or cylinders."""
        if x < radius + 0.5 or x > self.nx - 1.5 - radius:
            return True
        if y < radius + 0.5 or y > self.ny - 1.5 - radius:
            return True
        for cx, cy, cr in self.cylinders:
            if np.hypot(x - cx, y - cy) < cr + clearance * radius:
                return True
        return False

    def _neighbor_search_radius(self) -> float:
        interaction_cutoff = 0.0
        if self.particle_attraction:
            interaction_cutoff = max(interaction_cutoff, self.attraction_cutoff)
        if self.particle_repulsion:
            interaction_cutoff = max(interaction_cutoff, self.repulsion_cutoff)
        return 2.0 * float(np.max(self.radii)) + interaction_cutoff if self.n_p else 1.0

    def _particle_pair_candidates(self):
        """Yield nearby particle pairs."""
        if self.particle_search == "all_pairs" or self.n_p < 2:
            for i in range(self.n_p):
                for j in range(i + 1, self.n_p):
                    yield i, j
            return

        cell_size = max(self._neighbor_search_radius(), 1e-6)
        bins: dict[tuple[int, int], list[int]] = {}
        for idx, (x_pos, y_pos) in enumerate(self.pos):
            key = (int(x_pos // cell_size), int(y_pos // cell_size))
            bins.setdefault(key, []).append(idx)

        seen: set[tuple[int, int]] = set()
        for key, indices in bins.items():
            kx, ky = key
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for i in indices:
                        for j in bins.get((kx + dx, ky + dy), []):
                            if j <= i:
                                continue
                            pair = (i, j)
                            if pair not in seen:
                                seen.add(pair)
                                yield pair

    def compute_loads(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        """Compute DEM loads, including left/right wall contact."""
        forces, torques = self.dem_solver.compute_loads(dt)
        self._apply_side_wall_loads(forces, torques)
        return forces, torques

    def _apply_side_wall_loads(self, forces: np.ndarray, torques: np.ndarray) -> None:
        """Apply vertical side-wall contact loads for packing boxes."""
        for i in range(self.n_p):
            left = self.radii[i] + 0.5
            right = self.nx - 1.5 - self.radii[i]
            if self.pos[i, 0] < left:
                overlap = left - self.pos[i, 0]
                v_n = -self.vel[i, 0]
                f_mag = self.dem_solver.normal_contact_magnitude(overlap, v_n, self.masses[i])
                forces[i, 0] += f_mag
                self.dem_solver._apply_single_body_tangential_load(
                    i, np.array([1.0, 0.0]), f_mag, forces, torques
                )
            if self.pos[i, 0] > right:
                overlap = self.pos[i, 0] - right
                v_n = self.vel[i, 0]
                f_mag = self.dem_solver.normal_contact_magnitude(overlap, v_n, self.masses[i])
                forces[i, 0] -= f_mag
                self.dem_solver._apply_single_body_tangential_load(
                    i, np.array([-1.0, 0.0]), f_mag, forces, torques
                )

    def advance(self, steps: int = 1) -> None:
        """Advance the DEM packing simulation."""
        dt = 1.0 / max(self.dem_substeps, 1)
        for _ in range(steps):
            for _ in range(max(self.dem_substeps, 1)):
                self._substep(dt)
            self.step_count += 1

    def _substep(self, dt: float) -> None:
        if self.n_p == 0:
            return
        forces, torques = self.compute_loads(dt)
        acc = forces / self.masses[:, np.newaxis]
        alpha = torques / self.inertias
        self.vel += 0.5 * dt * acc
        self.omega_p += 0.5 * dt * alpha
        self.pos += dt * self.vel
        self._clamp_boundaries()

        forces_new, torques_new = self.compute_loads(dt)
        acc_new = forces_new / self.masses[:, np.newaxis]
        alpha_new = torques_new / self.inertias
        self.vel += 0.5 * dt * acc_new
        self.omega_p += 0.5 * dt * alpha_new
        self.forces_p = forces_new
        self.torques_p = torques_new
        self._clamp_boundaries()

    def _clamp_boundaries(self) -> None:
        """Keep particles inside walls and outside fixed cylinders."""
        for i in range(self.n_p):
            left = self.radii[i] + 0.5
            right = self.nx - 1.5 - self.radii[i]
            bottom = self.radii[i] + 0.5
            top = self.ny - 1.5 - self.radii[i]
            if self.pos[i, 0] < left:
                self.pos[i, 0] = left
                if self.vel[i, 0] < 0.0:
                    self.vel[i, 0] *= -0.15
            if self.pos[i, 0] > right:
                self.pos[i, 0] = right
                if self.vel[i, 0] > 0.0:
                    self.vel[i, 0] *= -0.15
            if self.pos[i, 1] < bottom:
                self.pos[i, 1] = bottom
                if self.vel[i, 1] < 0.0:
                    self.vel[i, 1] *= -0.15
            if self.pos[i, 1] > top:
                self.pos[i, 1] = top
                if self.vel[i, 1] > 0.0:
                    self.vel[i, 1] *= -0.15

            for cx, cy, cr in self.cylinders:
                dx = self.pos[i, 0] - cx
                dy = self.pos[i, 1] - cy
                dist = float(np.hypot(dx, dy))
                min_dist = cr + self.radii[i]
                if 1e-10 < dist < min_dist:
                    normal = np.array([dx / dist, dy / dist])
                    self.pos[i] = np.array([cx, cy]) + min_dist * normal
                    v_n = float(np.dot(self.vel[i], normal))
                    if v_n < 0.0:
                        self.vel[i] -= 1.15 * v_n * normal

    def metrics(self) -> PackingMetrics:
        """Return current packing metrics."""
        speeds = np.linalg.norm(self.vel, axis=1) if self.n_p else np.empty(0)
        kinetic = 0.5 * float(np.sum(self.masses * speeds**2))
        top = float(np.max(self.pos[:, 1] + self.radii)) if self.n_p else 0.0
        bottom = float(np.min(self.pos[:, 1] - self.radii)) if self.n_p else 0.0
        bed_height = max(top - bottom, 0.0)
        particle_area = float(np.sum(np.pi * self.radii**2))
        packing_area = max(self.nx * max(bed_height, 1e-12), 1e-12)
        return PackingMetrics(
            step=self.step_count,
            n_particles=self.n_p,
            kinetic_energy=kinetic,
            mean_speed=float(np.mean(speeds)) if self.n_p else 0.0,
            max_speed=float(np.max(speeds)) if self.n_p else 0.0,
            contact_count=self.contact_count(),
            bed_height=bed_height,
            packing_fraction=particle_area / packing_area,
        )

    def contact_count(self) -> int:
        """Count particle-particle and particle-wall/cylinder contacts."""
        count = 0
        for i, j in self._particle_pair_candidates():
            if np.linalg.norm(self.pos[j] - self.pos[i]) < self.radii[i] + self.radii[j]:
                count += 1
        for i in range(self.n_p):
            if self.pos[i, 0] <= self.radii[i] + 0.501:
                count += 1
            if self.pos[i, 0] >= self.nx - 1.501 - self.radii[i]:
                count += 1
            if self.pos[i, 1] <= self.radii[i] + 0.501:
                count += 1
            if self.pos[i, 1] >= self.ny - 1.501 - self.radii[i]:
                count += 1
            for cx, cy, cr in self.cylinders:
                if np.hypot(self.pos[i, 0] - cx, self.pos[i, 1] - cy) <= cr + self.radii[i] + 1e-3:
                    count += 1
        return count

    def snapshot(self) -> dict[str, np.ndarray | int]:
        """Return a serializable snapshot dictionary."""
        return {
            "step": self.step_count,
            "pos": self.pos.copy(),
            "vel": self.vel.copy(),
            "radii": self.radii.copy(),
            "omega": self.omega_p.copy(),
            "force": self.forces_p.copy(),
            "torque": self.torques_p.copy(),
        }


def write_particles_vtk(path: Path, snap: dict[str, np.ndarray | int]) -> None:
    """Write DEM particle centers as legacy VTK polydata."""
    positions = np.asarray(snap["pos"])
    velocities = np.asarray(snap["vel"])
    radii = np.asarray(snap["radii"])
    forces = np.asarray(snap["force"])
    n_particles = len(positions)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write(f"DEM packing particles step {snap['step']}\n")
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
        for force in np.linalg.norm(forces, axis=1):
            handle.write(f"{_vtk_float(force)}\n")


def write_metrics_csv(path: Path, rows: list[PackingMetrics]) -> None:
    """Write packing metrics to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def write_pvd(path: Path, entries: list[tuple[int, Path]]) -> None:
    """Write a ParaView collection file for particle VTK snapshots."""
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


def write_summary(path: Path, config: dict, metrics: PackingMetrics, settled: bool) -> None:
    """Write a compact JSON summary for later reuse."""
    summary = {
        "config": config,
        "final_metrics": asdict(metrics),
        "settled": settled,
    }
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
