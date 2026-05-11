"""Helpers for keeping generated results under ``src/results``."""

from __future__ import annotations

from pathlib import Path

RESULTS_ROOT = Path(__file__).resolve().parents[1] / "results"


def program_results_dir(program_file: str | Path, *parts: str) -> Path:
    """Return ``src/results/<program-name>/`` and create it if needed."""
    program_name = Path(program_file).stem
    output_dir = RESULTS_ROOT / program_name
    if parts:
        output_dir = output_dir.joinpath(*parts)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir
