#!/usr/bin/env python3
"""Aggregate LBM-DEM design sweeps and rank likely fouling drivers."""

from __future__ import annotations

import csv
import math
import re
import shlex
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
RESULT_ROOT = REPO_ROOT / "src" / "results" / "run_lbm_dem"
SWEEP_ROOTS = [
    RESULT_ROOT / "pore_geometry_16_parallel",
    RESULT_ROOT / "design_sweep_100_parallel",
]
OUTPUT_DIR = RESULT_ROOT / "design_sweep_analysis"


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _safe_float(value: str | float | int | None, default: float = math.nan) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _command_tokens(log_file: Path) -> list[str]:
    if not log_file.exists():
        return []
    for line in log_file.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("Command:"):
            return shlex.split(line.removeprefix("Command:").strip())
    return []


def _arg_value(tokens: list[str], option: str, default: str | None = None) -> str | None:
    for idx, token in enumerate(tokens):
        if token == option and idx + 1 < len(tokens):
            return tokens[idx + 1]
        if token.startswith(option + "="):
            return token.split("=", 1)[1]
    return default


def _has_arg(tokens: list[str], option: str) -> bool:
    return option in tokens or any(token.startswith(option + "=") for token in tokens)


def _cylinders_from_tokens(tokens: list[str]) -> list[tuple[float, float, float]]:
    cylinders: list[tuple[float, float, float]] = []
    idx = 0
    while idx < len(tokens):
        if tokens[idx] == "--cylinder-spec" and idx + 3 < len(tokens):
            cylinders.append((float(tokens[idx + 1]), float(tokens[idx + 2]), float(tokens[idx + 3])))
            idx += 4
        else:
            idx += 1
    return cylinders


def _find_time_series(tag: str) -> Path | None:
    candidates = list(RESULT_ROOT.glob(f"**/{tag}/**/analysis/time_series.csv"))
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def _read_timeseries(path: Path) -> list[dict[str, float]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    parsed: list[dict[str, float]] = []
    for row in rows:
        parsed.append({key: _safe_float(value) for key, value in row.items()})
    return parsed


def _slope(rows: list[dict[str, float]], key: str) -> float:
    values = np.array([row.get(key, math.nan) for row in rows], dtype=float)
    steps = np.array([row.get("step", idx) for idx, row in enumerate(rows)], dtype=float)
    mask = np.isfinite(values) & np.isfinite(steps)
    if np.count_nonzero(mask) < 2:
        return math.nan
    x = steps[mask] - steps[mask][0]
    if np.allclose(x, 0.0):
        return math.nan
    return float(np.polyfit(x, values[mask], 1)[0])


def _mean(rows: list[dict[str, float]], key: str) -> float:
    values = np.array([row.get(key, math.nan) for row in rows], dtype=float)
    values = values[np.isfinite(values)]
    return float(values.mean()) if values.size else math.nan


def _final(rows: list[dict[str, float]], key: str) -> float:
    for row in reversed(rows):
        value = row.get(key, math.nan)
        if math.isfinite(value):
            return float(value)
    return math.nan


def _geometry_features(cylinders: list[tuple[float, float, float]], nx: float, ny: float) -> dict[str, float]:
    count = len(cylinders)
    radii = np.array([radius for _, _, radius in cylinders], dtype=float)
    solid_area = float(np.sum(np.pi * radii**2)) if count else 0.0
    channel_area = max(nx * max(ny - 2.0, 1.0), 1.0)
    min_gap = math.nan
    if count >= 2:
        gaps = []
        for i, (x0, y0, r0) in enumerate(cylinders):
            for x1, y1, r1 in cylinders[i + 1 :]:
                gaps.append(math.hypot(x1 - x0, y1 - y0) - r0 - r1)
        min_gap = float(min(gaps))
    return {
        "cylinder_count": float(count),
        "mean_cylinder_radius": float(radii.mean()) if count else 0.0,
        "solid_area_fraction": solid_area / channel_area,
        "min_cylinder_gap": min_gap,
    }


def _normalise(values: list[float], invert: bool = False) -> list[float]:
    arr = np.array(values, dtype=float)
    mask = np.isfinite(arr)
    out = np.zeros_like(arr, dtype=float)
    if np.count_nonzero(mask) == 0:
        return out.tolist()
    min_v = float(arr[mask].min())
    max_v = float(arr[mask].max())
    if math.isclose(min_v, max_v):
        out[mask] = 0.5
    else:
        out[mask] = (arr[mask] - min_v) / (max_v - min_v)
    if invert:
        out[mask] = 1.0 - out[mask]
    return out.tolist()


def _rank(values: list[float]) -> np.ndarray:
    arr = np.array(values, dtype=float)
    order = np.argsort(arr)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(len(arr), dtype=float)
    return ranks


def _spearman(x_values: list[float], y_values: list[float]) -> float:
    x = np.array(x_values, dtype=float)
    y = np.array(y_values, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 3:
        return math.nan
    xr = _rank(x[mask].tolist())
    yr = _rank(y[mask].tolist())
    if np.std(xr) == 0.0 or np.std(yr) == 0.0:
        return math.nan
    return float(np.corrcoef(xr, yr)[0, 1])


def _collect_records() -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    records: list[dict[str, object]] = []
    missing: list[dict[str, str]] = []

    for sweep_root in SWEEP_ROOTS:
        for row in _read_tsv(sweep_root / "status.tsv"):
            tag = row.get("condition", "")
            if not tag:
                continue
            status = row.get("status", "")
            if status != "ok":
                missing.append({"condition": tag, "reason": f"status={status}", "sweep": str(sweep_root)})
                continue

            log_file = REPO_ROOT / row.get("log_file", "")
            tokens = _command_tokens(log_file)
            time_series = _find_time_series(tag)
            if time_series is None:
                missing.append({"condition": tag, "reason": "time_series.csv not found", "sweep": str(sweep_root)})
                continue

            rows = _read_timeseries(time_series)
            if not rows:
                missing.append({"condition": tag, "reason": "empty time_series.csv", "sweep": str(sweep_root)})
                continue

            nx = _safe_float(row.get("nx"), _safe_float(_arg_value(tokens, "--nx"), 540.0))
            ny = _safe_float(row.get("ny"), _safe_float(_arg_value(tokens, "--ny"), 210.0))
            cylinders = _cylinders_from_tokens(tokens)
            geometry = _geometry_features(cylinders, nx, ny)

            reynolds_number = _safe_float(
                row.get("reynolds_number"),
                _safe_float(_arg_value(tokens, "--reynolds-number") or _arg_value(tokens, "--Re"), 100.0),
            )
            phi = _safe_float(
                row.get("phi"),
                _safe_float(_arg_value(tokens, "--particle-volume-fraction"), 0.40),
            )
            attraction_strength = _safe_float(
                row.get("attraction_strength"),
                _safe_float(_arg_value(tokens, "--attraction-strength"), 0.0),
            )
            rolling_coeff = _safe_float(
                row.get("rolling_friction_coeff"),
                _safe_float(_arg_value(tokens, "--rolling-friction-coeff"), 0.0),
            )
            surface_force = row.get("surface_force") or ("attraction" if _has_arg(tokens, "--particle-attraction") else "none")
            rolling = row.get("rolling_friction") or ("rolling" if _has_arg(tokens, "--rolling-friction") else "free")

            final_generated = _final(rows, "generated_particles")
            final_passed = _final(rows, "passed_particles")
            passed_ratio = final_passed / final_generated if final_generated and final_generated > 0 else math.nan
            flux_values = [record.get("permeate_flux", math.nan) for record in rows]
            flux_norm = [value for value in flux_values if math.isfinite(value)]
            flux_decline = (flux_norm[0] - flux_norm[-1]) / abs(flux_norm[0]) if len(flux_norm) >= 2 and abs(flux_norm[0]) > 1e-12 else math.nan

            record: dict[str, object] = {
                "condition": tag,
                "sweep": sweep_root.name,
                "time_series": str(time_series),
                "reynolds_number": reynolds_number,
                "phi": phi,
                "surface_force": surface_force,
                "attraction_strength": attraction_strength,
                "rolling_friction": rolling,
                "rolling_friction_coeff": rolling_coeff,
                "final_generated_particles": final_generated,
                "final_passed_particles": final_passed,
                "passed_particle_ratio": passed_ratio,
                "mean_pressure_drop": _mean(rows, "pressure_drop"),
                "pressure_drop_slope": _slope(rows, "pressure_drop"),
                "mean_permeate_flux": _mean(rows, "permeate_flux"),
                "flux_decline_ratio": flux_decline,
                "final_total_contacts": _final(rows, "total_contacts"),
                "mean_total_contacts": _mean(rows, "total_contacts"),
                "contact_growth_rate": _slope(rows, "total_contacts"),
                "final_particle_area_fraction": _final(rows, "particle_area_fraction"),
                "mean_particle_area_fraction": _mean(rows, "particle_area_fraction"),
                "mean_max_speed": _mean(rows, "max_speed"),
                "final_max_speed": _final(rows, "max_speed"),
                **geometry,
            }
            records.append(record)
    return records, missing


def _add_scores(records: list[dict[str, object]]) -> None:
    components = {
        "contact_growth_rate": False,
        "mean_total_contacts": False,
        "pressure_drop_slope": False,
        "mean_pressure_drop": False,
        "final_particle_area_fraction": False,
        "mean_permeate_flux": True,
        "passed_particle_ratio": True,
    }
    normalised = {
        key: _normalise([float(record.get(key, math.nan)) for record in records], invert=invert)
        for key, invert in components.items()
    }
    for idx, record in enumerate(records):
        values = [normalised[key][idx] for key in components]
        record["fouling_risk_score"] = float(np.mean(values)) if values else math.nan


def _write_records_csv(records: list[dict[str, object]], path: Path) -> None:
    if not records:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(records[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def _group_mean(records: list[dict[str, object]], key: str, metric: str) -> list[tuple[str, float, int]]:
    groups: dict[str, list[float]] = defaultdict(list)
    for record in records:
        value = _safe_float(record.get(metric))
        if math.isfinite(value):
            groups[str(record.get(key))].append(value)
    return sorted(
        ((group, float(np.mean(values)), len(values)) for group, values in groups.items()),
        key=lambda item: item[1],
        reverse=True,
    )


def _make_plots(records: list[dict[str, object]], output_dir: Path) -> list[Path]:
    plot_paths: list[Path] = []
    if not records:
        return plot_paths

    numeric_features = [
        "reynolds_number",
        "phi",
        "attraction_strength",
        "rolling_friction_coeff",
        "cylinder_count",
        "mean_cylinder_radius",
        "solid_area_fraction",
        "min_cylinder_gap",
    ]
    correlations = [
        (feature, _spearman([_safe_float(r.get(feature)) for r in records], [_safe_float(r.get("fouling_risk_score")) for r in records]))
        for feature in numeric_features
    ]
    correlations = [(feature, value) for feature, value in correlations if math.isfinite(value)]
    if correlations:
        labels, values = zip(*sorted(correlations, key=lambda item: abs(item[1]), reverse=True))
        fig, ax = plt.subplots(figsize=(9, 5))
        colors = ["#b23a48" if value > 0 else "#2f6f9f" for value in values]
        ax.barh(labels, values, color=colors)
        ax.axvline(0.0, color="black", linewidth=0.8)
        ax.set_xlabel("Spearman correlation with fouling risk score")
        ax.set_title("Likely Drivers Of Fouling Risk")
        ax.invert_yaxis()
        fig.tight_layout()
        path = output_dir / "01_factor_correlations.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        plot_paths.append(path)

    fig, ax = plt.subplots(figsize=(7, 5))
    for re_value in sorted({str(r.get("reynolds_number")) for r in records}):
        subset = [r for r in records if str(r.get("reynolds_number")) == re_value]
        ax.scatter(
            [_safe_float(r.get("mean_pressure_drop")) for r in subset],
            [_safe_float(r.get("mean_total_contacts")) for r in subset],
            label=f"Re={re_value}",
            alpha=0.75,
        )
    ax.set_xlabel("Mean pressure drop")
    ax.set_ylabel("Mean total contacts")
    ax.set_title("Pressure Drop vs Contact Accumulation")
    ax.legend()
    fig.tight_layout()
    path = output_dir / "02_pressure_drop_vs_contacts.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    plot_paths.append(path)

    re_groups = _group_mean(records, "reynolds_number", "fouling_risk_score")
    if re_groups:
        fig, ax = plt.subplots(figsize=(6, 4))
        labels = [item[0] for item in re_groups]
        values = [item[1] for item in re_groups]
        ax.bar(labels, values, color="#6a8d73")
        ax.set_xlabel("Reynolds number")
        ax.set_ylabel("Mean fouling risk score")
        ax.set_title("Fouling Severity By Re")
        fig.tight_layout()
        path = output_dir / "03_reynolds_severity.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        plot_paths.append(path)

    return plot_paths


def _write_markdown(records: list[dict[str, object]], missing: list[dict[str, str]], plots: list[Path], path: Path) -> None:
    completed = len(records)
    top_risk = sorted(records, key=lambda row: _safe_float(row.get("fouling_risk_score")), reverse=True)[:10]
    low_risk = sorted(records, key=lambda row: _safe_float(row.get("fouling_risk_score")))[:10]
    numeric_features = [
        "reynolds_number",
        "phi",
        "attraction_strength",
        "rolling_friction_coeff",
        "cylinder_count",
        "mean_cylinder_radius",
        "solid_area_fraction",
        "min_cylinder_gap",
    ]
    correlations = sorted(
        [
            (
                feature,
                _spearman(
                    [_safe_float(record.get(feature)) for record in records],
                    [_safe_float(record.get("fouling_risk_score")) for record in records],
                ),
            )
            for feature in numeric_features
        ],
        key=lambda item: abs(item[1]) if math.isfinite(item[1]) else -1.0,
        reverse=True,
    )

    lines = [
        "# LBM-DEM Design Sweep Analysis",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        f"Completed conditions analysed: {completed}",
        f"Missing or failed conditions: {len(missing)}",
        "",
        "## Fouling Risk Score",
        "",
        "The score is a normalized composite of contact growth, contact count, pressure-drop growth, mean pressure drop, deposited particle fraction, low permeate flux, and low passed-particle ratio. It is a screening score, not a replacement for visual inspection or final LBM-DEM validation.",
        "",
        "## Strongest Numeric Associations",
        "",
    ]
    for feature, corr in correlations:
        if math.isfinite(corr):
            lines.append(f"- `{feature}`: Spearman correlation {corr:+.3f}")
    lines.extend(["", "## Highest Fouling-Risk Conditions", ""])
    for record in top_risk:
        lines.append(
            "- "
            f"`{record['condition']}`: score={_safe_float(record.get('fouling_risk_score')):.3f}, "
            f"Re={_safe_float(record.get('reynolds_number')):.0f}, "
            f"phi={_safe_float(record.get('phi')):.2f}, "
            f"cylinders={_safe_float(record.get('cylinder_count')):.0f}, "
            f"radius={_safe_float(record.get('mean_cylinder_radius')):.1f}, "
            f"attr={_safe_float(record.get('attraction_strength')):.4g}, "
            f"mu_r={_safe_float(record.get('rolling_friction_coeff')):.3f}"
        )
    lines.extend(["", "## Lowest Fouling-Risk Conditions", ""])
    for record in low_risk:
        lines.append(
            "- "
            f"`{record['condition']}`: score={_safe_float(record.get('fouling_risk_score')):.3f}, "
            f"Re={_safe_float(record.get('reynolds_number')):.0f}, "
            f"phi={_safe_float(record.get('phi')):.2f}, "
            f"cylinders={_safe_float(record.get('cylinder_count')):.0f}, "
            f"radius={_safe_float(record.get('mean_cylinder_radius')):.1f}, "
            f"attr={_safe_float(record.get('attraction_strength')):.4g}, "
            f"mu_r={_safe_float(record.get('rolling_friction_coeff')):.3f}"
        )
    lines.extend(["", "## Re Severity", ""])
    for group, mean_value, count in _group_mean(records, "reynolds_number", "fouling_risk_score"):
        lines.append(f"- Re={group}: mean fouling risk={mean_value:.3f} (n={count})")
    lines.extend(["", "## Categorical Group Means", ""])
    for key in ["surface_force", "rolling_friction", "sweep"]:
        lines.append(f"### {key}")
        for group, mean_value, count in _group_mean(records, key, "fouling_risk_score"):
            lines.append(f"- {group}: mean fouling risk={mean_value:.3f} (n={count})")
        lines.append("")
    if plots:
        lines.extend(["## Plots", ""])
        for plot in plots:
            lines.append(f"- `{plot.relative_to(REPO_ROOT)}`")
    if missing:
        lines.extend(["", "## Missing Or Failed Conditions", ""])
        for item in missing[:200]:
            lines.append(f"- `{item['condition']}`: {item['reason']}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    records, missing = _collect_records()
    _add_scores(records)
    records.sort(key=lambda row: str(row.get("condition")))

    csv_path = OUTPUT_DIR / "design_sweep_metrics.csv"
    missing_path = OUTPUT_DIR / "missing_or_failed_conditions.tsv"
    report_path = OUTPUT_DIR / "analysis_report.md"

    _write_records_csv(records, csv_path)
    with missing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["condition", "reason", "sweep"], delimiter="\t")
        writer.writeheader()
        writer.writerows(missing)
    plots = _make_plots(records, OUTPUT_DIR)
    _write_markdown(records, missing, plots, report_path)

    print(f"Analysed conditions: {len(records)}")
    print(f"Metrics: {csv_path}")
    print(f"Report:  {report_path}")
    if missing:
        print(f"Missing/failed: {missing_path}")


if __name__ == "__main__":
    main()
