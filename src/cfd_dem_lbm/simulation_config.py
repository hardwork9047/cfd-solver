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
        return cls(values=_normalise_value(mapping), source_path=source_path)

    @classmethod
    def from_json(cls, path: str | Path) -> "SimulationConfig":
        config_path = Path(path)
        with config_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            raise ValueError("simulation config JSON must contain an object at the top level")
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
