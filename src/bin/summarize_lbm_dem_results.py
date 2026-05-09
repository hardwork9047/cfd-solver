#!/usr/bin/env python3
"""Summarize LBM-DEM run outputs into CSV and Markdown reports."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
RESULT_ROOT = REPO_ROOT / "src" / "results" / "run_lbm_dem"


def _safe_float(value: Any, default: float = math.nan) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _read_timeseries(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return [{key: _safe_float(value) for key, value in row.items()} for row in rows]


def _final(rows: list[dict[str, float]], key: str) -> float:
    for row in reversed(rows):
        value = row.get(key, math.nan)
        if math.isfinite(value):
            return value
    return math.nan


def _mean(rows: list[dict[str, float]], key: str) -> float:
    values = [row.get(key, math.nan) for row in rows]
    values = [value for value in values if math.isfinite(value)]
    return float(sum(values) / len(values)) if values else math.nan


def _slope(rows: list[dict[str, float]], key: str) -> float:
    values = [(row.get("step", idx), row.get(key, math.nan)) for idx, row in enumerate(rows)]
    values = [(float(step), float(value)) for step, value in values if math.isfinite(value)]
    if len(values) < 2:
        return math.nan
    x0 = values[0][0]
    xs = [step - x0 for step, _ in values]
    ys = [value for _, value in values]
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    denom = sum((x - x_mean) ** 2 for x in xs)
    if denom <= 1e-30:
        return math.nan
    return sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / denom


def _metadata_for(timeseries: Path) -> dict[str, Any]:
    metadata_path = timeseries.parents[1] / "metadata.json"
    if not metadata_path.exists():
        return {}
    try:
        return json.loads(metadata_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _run_record(timeseries: Path) -> dict[str, Any] | None:
    rows = _read_timeseries(timeseries)
    if not rows:
        return None
    metadata = _metadata_for(timeseries)
    config = metadata.get("configuration", {})
    solver = metadata.get("solver", {})
    args = metadata.get("arguments", {})
    final_generated = _final(rows, "generated_particles")
    final_passed = _final(rows, "passed_particles")
    passed_ratio = (
        final_passed / final_generated
        if math.isfinite(final_generated) and final_generated > 0
        else math.nan
    )
    return {
        "run_dir": str(timeseries.parents[1]),
        "created_at": metadata.get("created_at", ""),
        "result_tag": args.get("result_tag", ""),
        "fluid_method": solver.get("fluid_method", config.get("fluid_method", "")),
        "particle_method": solver.get("particle_method", config.get("particle_method", "")),
        "particle_fluid_coupling": solver.get(
            "particle_fluid_coupling",
            config.get("particle_fluid_coupling", ""),
        ),
        "flow_condition": config.get("flow_condition", args.get("flow_condition", "")),
        "reynolds_number": config.get("reynolds_number", math.nan),
        "particle_volume_fraction": config.get("particle_volume_fraction", math.nan),
        "source_volume_fraction": config.get("source_volume_fraction", math.nan),
        "n_cylinders": len(config.get("cylinders", [])),
        "final_step": _final(rows, "step"),
        "final_active_particles": _final(rows, "active_particles"),
        "final_generated_particles": final_generated,
        "final_passed_particles": final_passed,
        "final_adhered_particles": _final(rows, "adhered_particles"),
        "final_adhesion_events": _final(rows, "adhesion_events"),
        "final_detachment_events": _final(rows, "detachment_events"),
        "passed_particle_ratio": passed_ratio,
        "final_normalized_permeate_flux": _final(rows, "normalized_permeate_flux"),
        "mean_normalized_permeate_flux": _mean(rows, "normalized_permeate_flux"),
        "final_pressure_drop": _final(rows, "pressure_drop"),
        "mean_pressure_drop": _mean(rows, "pressure_drop"),
        "pressure_drop_slope": _slope(rows, "pressure_drop"),
        "final_fouling_resistance_index": _final(rows, "fouling_resistance_index"),
        "mean_fouling_resistance_index": _mean(rows, "fouling_resistance_index"),
        "final_total_contacts": _final(rows, "total_contacts"),
        "mean_total_contacts": _mean(rows, "total_contacts"),
        "contact_growth_rate": _slope(rows, "total_contacts"),
        "final_dynamic_particle_solid_fraction": _final(rows, "dynamic_particle_solid_fraction"),
        "final_porous_resistance_fraction": _final(rows, "porous_resistance_fraction"),
        "mean_observed_reynolds_number": _mean(rows, "observed_reynolds_number"),
        "mean_particle_reynolds_number": _mean(rows, "particle_reynolds_number"),
        "mean_stokes_number_estimate": _mean(rows, "stokes_number_estimate"),
        "mean_brinkman_resistance_number": _mean(rows, "brinkman_resistance_number"),
    }


def collect_records(root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for timeseries in sorted(root.glob("**/analysis/time_series.csv")):
        record = _run_record(timeseries)
        if record is not None:
            records.append(record)
    return records


def _write_csv(records: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(records[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def _write_markdown(records: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    high_risk = sorted(
        records,
        key=lambda row: _safe_float(row.get("mean_fouling_resistance_index")),
        reverse=True,
    )[:10]
    best_flux = sorted(
        records,
        key=lambda row: _safe_float(row.get("final_normalized_permeate_flux"), -math.inf),
        reverse=True,
    )[:10]
    lines = [
        "# LBM-DEM Result Summary",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Runs analysed: {len(records)}",
        "",
        "## Highest Fouling Resistance",
        "",
    ]
    for record in high_risk:
        lines.append(
            "- "
            f"`{record['result_tag'] or Path(record['run_dir']).name}`: "
            f"FRI={_safe_float(record.get('mean_fouling_resistance_index')):.3g}, "
            f"J/J0={_safe_float(record.get('final_normalized_permeate_flux')):.3g}, "
            f"coupling={record.get('particle_fluid_coupling')}, "
            f"flow={record.get('flow_condition')}"
        )
    lines.extend(["", "## Highest Final Permeate Flux", ""])
    for record in best_flux:
        lines.append(
            "- "
            f"`{record['result_tag'] or Path(record['run_dir']).name}`: "
            f"J/J0={_safe_float(record.get('final_normalized_permeate_flux')):.3g}, "
            f"FRI={_safe_float(record.get('mean_fouling_resistance_index')):.3g}, "
            f"coupling={record.get('particle_fluid_coupling')}, "
            f"flow={record.get('flow_condition')}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize LBM-DEM result directories")
    parser.add_argument("--root", type=Path, default=RESULT_ROOT)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=RESULT_ROOT / "summary",
        help="Directory for summary CSV and Markdown",
    )
    args = parser.parse_args()

    records = collect_records(args.root)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "lbm_dem_run_summary.csv"
    md_path = args.out_dir / "lbm_dem_run_summary.md"
    _write_csv(records, csv_path)
    _write_markdown(records, md_path)
    print(f"Runs analysed: {len(records)}")
    print(f"CSV: {csv_path}")
    print(f"Report: {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
