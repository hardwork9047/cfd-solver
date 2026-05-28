"""Configuration helpers for LBM-DEM production runs.

The command-line runner is still the main execution entrypoint.  This module
keeps file-based run definitions small and explicit so sweeps can be reproduced
without hiding options in shell scripts.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

META_KEYS = {
    "schema_version",
    "name",
    "description",
    "case_type",
    "notes",
}
SECTION_KEYS = {
    "domain",
    "flow",
    "numerics",
    "particles",
    "physics",
    "runtime",
    "outputs",
    "stability",
    # Phase 2 additions — solver separates physical model choice; accelerator separates backend choice.
    "solver",
    "accelerator",
}


def _normalise_key(key: str) -> str:
    """Map JSON-friendly keys to argparse destination names."""
    return key.strip().replace("-", "_")


def _normalise_value(value: Any) -> Any:
    """Recursively convert tuples and paths into JSON-serialisable values."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return [_normalise_value(item) for item in value]
    if isinstance(value, list):
        return [_normalise_value(item) for item in value]
    if isinstance(value, dict):
        return {_normalise_key(str(key)): _normalise_value(item) for key, item in value.items()}
    return value


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge two dictionaries without mutating either input."""
    merged = dict(base)
    for key, value in override.items():
        if (
            key in merged
            and isinstance(merged[key], dict)
            and isinstance(value, dict)
        ):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _load_json_object(path: Path) -> dict[str, Any]:
    """Load a JSON object from disk."""
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("simulation config JSON must contain an object at the top level")
    return payload


def _resolve_relative_path(value: str | Path, *, base_dir: Path) -> Path:
    """Resolve config references relative to the including file."""
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return path


def _load_with_extends(path: Path, seen: set[Path] | None = None) -> dict[str, Any]:
    """Load a config and recursively merge any ``extends`` parents."""
    resolved = path.resolve()
    seen = set() if seen is None else set(seen)
    if resolved in seen:
        raise ValueError(f"cyclic config extends detected at {path}")
    seen.add(resolved)

    payload = _load_json_object(path)
    parents = payload.pop("extends", None)
    if parents is None:
        return payload
    if isinstance(parents, (str, Path)):
        parent_list = [parents]
    elif isinstance(parents, list):
        parent_list = parents
    else:
        raise ValueError("extends must be a string or list of strings")

    merged: dict[str, Any] = {}
    for parent in parent_list:
        parent_path = _resolve_relative_path(parent, base_dir=path.parent)
        merged = _deep_merge(merged, _load_with_extends(parent_path, seen))
    return _deep_merge(merged, payload)


# Keys in the ``accelerator`` section use short names; map to argparse destinations.
_ACCELERATOR_KEY_MAP = {
    "fluid": "fluid_accelerator",
    "compute": "compute_accelerator",
}

# Keys in the ``solver.lees_edwards`` sub-section mapped to solver __init__ params.
_LE_KEY_MAP = {
    "shear_rate": "le_shear_rate",
    "shear_axis": "le_shear_axis",
    "boundary_axis": "le_boundary_axis",
    "interpolation_order": "le_interpolation_order",
}


def _expand_solver_section(solver: dict[str, Any]) -> dict[str, Any]:
    """Expand ``solver.lees_edwards`` sub-dict into flat ``le_*`` keys.

    When ``lees_edwards.enabled`` is ``true``, also injects ``y_boundary``
    and ``streamwise_boundary`` defaults (``"lees_edwards"`` and
    ``"periodic_force"`` respectively) unless those keys are already present.
    """
    result = dict(solver)
    le_config = result.pop("lees_edwards", None)
    if le_config is None:
        return result
    if not isinstance(le_config, dict):
        raise ValueError("solver.lees_edwards must be an object when provided")
    if le_config.get("enabled", False):
        result.setdefault("y_boundary", "lees_edwards")
        result.setdefault("streamwise_boundary", "periodic_force")
    for k, v in le_config.items():
        if k == "enabled":
            continue
        dest = _LE_KEY_MAP.get(k, f"le_{k}")
        result[dest] = v
    return result


_VISCOSITY_EVAL_KEY_MAP = {
    "enabled": "viscosity_eval_enabled",
    "start_step": "viscosity_eval_start_step",
    "viscosity_interval": "viscosity_eval_interval",
    "average_steps": "viscosity_eval_average_steps",
}


def _expand_runtime_section(runtime: dict[str, Any]) -> dict[str, Any]:
    """Expand ``runtime.viscosity_eval`` sub-dict into flat ``viscosity_eval_*`` keys."""
    result = dict(runtime)
    ve_config = result.pop("viscosity_eval", None)
    if ve_config is None:
        return result
    if not isinstance(ve_config, dict):
        raise ValueError("runtime.viscosity_eval must be an object when provided")
    for k, v in ve_config.items():
        dest = _VISCOSITY_EVAL_KEY_MAP.get(k, f"viscosity_eval_{k}")
        result[dest] = v
    return result


def _flatten_sections(values: dict[str, Any]) -> dict[str, Any]:
    """Flatten research-friendly config sections into runner argument keys."""
    payload: dict[str, Any] = {}
    for key, value in values.items():
        if key in META_KEYS:
            continue
        if key in SECTION_KEYS:
            if not isinstance(value, dict):
                raise ValueError(f"{key} must be an object when provided")
            if key == "accelerator":
                # Remap short keys (fluid, compute) to argparse destinations.
                value = {_ACCELERATOR_KEY_MAP.get(k, k): v for k, v in value.items()}
            elif key == "solver":
                value = _expand_solver_section(value)
            elif key == "runtime":
                value = _expand_runtime_section(value)
            payload = _deep_merge(payload, value)
        else:
            payload[key] = value
    return payload


def _normalise_geometry_cylinders(values: dict[str, Any]) -> dict[str, Any]:
    """Allow structured geometry blocks while keeping argparse compatibility."""
    payload = dict(values)
    geometry = payload.pop("geometry", None)
    if geometry is None or "cylinder_spec" in payload:
        return payload
    if not isinstance(geometry, dict):
        raise ValueError("geometry must be an object when provided")
    cylinders = geometry.get("cylinders")
    if cylinders is None:
        return payload
    cylinder_spec = []
    for item in cylinders:
        if isinstance(item, dict):
            cylinder_spec.append([item["x"], item["y"], item["radius"]])
        else:
            cylinder_spec.append(item)
    payload["cylinder_spec"] = cylinder_spec
    return payload


def _prepare_values(mapping: dict[str, Any]) -> dict[str, Any]:
    """Normalize, flatten, and adapt a user-facing config mapping."""
    values = _normalise_value(mapping)
    values = _flatten_sections(values)
    return _normalise_geometry_cylinders(values)


@dataclass(frozen=True)
class SimulationConfig:
    """A file-backed simulation configuration.

    Unknown fields are preserved so newer runners can add options without
    changing the config loader.  Keys may use either argparse destination names
    such as ``streamwise_boundary`` or CLI-style names such as
    ``streamwise-boundary``.
    """

    values: dict[str, Any] = field(default_factory=dict)
    source_path: str | None = None

    @classmethod
    def from_mapping(
        cls,
        mapping: dict[str, Any],
        *,
        source_path: str | None = None,
    ) -> "SimulationConfig":
        values = _prepare_values(mapping)
        return cls(values=values, source_path=source_path)

    @classmethod
    def from_json(cls, path: str | Path) -> "SimulationConfig":
        config_path = Path(path)
        payload = _load_with_extends(config_path)
        return cls.from_mapping(payload, source_path=str(config_path))

    def argparse_defaults(self) -> dict[str, Any]:
        """Return values suitable for ``ArgumentParser.set_defaults``."""
        return dict(self.values)

    def to_dict(self) -> dict[str, Any]:
        """Return a stable JSON-serialisable representation."""
        payload = dict(self.values)
        if self.source_path is not None:
            payload["_source_path"] = self.source_path
        return payload

    def write_effective(self, path: str | Path, effective_arguments: dict[str, Any]) -> None:
        """Write the resolved run config, including CLI overrides."""
        target = Path(path)
        payload = {
            "schema_version": 1,
            "source": self.to_dict() if self.values else None,
            "effective_arguments": _normalise_value(effective_arguments),
        }
        target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
