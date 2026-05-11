#!/usr/bin/env python3
"""Plot fouling metrics from LBM-DEM analysis/time_series.csv files."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cfd_dem_lbm.result_paths import program_results_dir


METRICS = [
    ("normalized_permeate_flux", "Normalized permeate flux"),
    ("pressure_drop_ratio", "Global pressure drop ratio"),
    ("local_pressure_drop_ratio", "Local pore pressure drop ratio"),
    ("fouling_resistance_index", "Fouling resistance index"),
    ("retained_particle_ratio", "Retained particle ratio"),
    ("cylinder_contacts", "Cylinder contacts"),
    ("effective_inlet_particle_fraction", "Effective inlet particle fraction"),
]


def _safe_float(value: str | None, default: float = np.nan) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _find_time_series(roots: list[Path]) -> list[Path]:
    paths: list[Path] = []
    for root in roots:
        if root.is_file() and root.name == "time_series.csv":
            paths.append(root)
        elif root.is_dir():
            paths.extend(root.glob("**/analysis/time_series.csv"))
    return sorted(set(paths))


def _run_label(path: Path) -> str:
    run_dir = path.parent.parent
    parts = run_dir.parts
    if len(parts) >= 3:
        return "/".join(parts[-3:])
    return run_dir.name


def _summary_path(time_series: Path) -> Path:
    return time_series.parent / "summary.json"


def _load_summary(time_series: Path) -> dict:
    path = _summary_path(time_series)
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_dataset(output_dir: Path, series_paths: list[Path]) -> Path:
    rows = []
    for path in series_paths:
        data_rows = _read_rows(path)
        if not data_rows:
            continue
        final = data_rows[-1]
        summary = _load_summary(path)
        run_status = path.parent.parent / "run_status.json"
        rows.append(
            {
                "label": _run_label(path),
                "time_series": str(path),
                "summary_json": str(_summary_path(path)) if _summary_path(path).exists() else "",
                "run_status": str(run_status) if run_status.exists() else "",
                "final_step": final.get("step", ""),
                "final_normalized_permeate_flux": final.get("normalized_permeate_flux", ""),
                "final_pressure_drop_ratio": final.get("pressure_drop_ratio", ""),
                "final_local_pressure_drop_ratio": final.get("local_pressure_drop_ratio", ""),
                "final_fouling_resistance_index": final.get("fouling_resistance_index", ""),
                "final_retained_particle_ratio": final.get("retained_particle_ratio", ""),
                "final_cylinder_contacts": final.get("cylinder_contacts", ""),
                "final_effective_inlet_particle_fraction": final.get(
                    "effective_inlet_particle_fraction",
                    "",
                ),
                "summary_status": summary.get("status", ""),
            }
        )
    dataset_path = output_dir / "fouling_metrics_dataset.csv"
    if not rows:
        dataset_path.write_text("", encoding="utf-8")
        return dataset_path
    with dataset_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return dataset_path


def _plot_time_series(output_dir: Path, series_paths: list[Path], max_runs: int) -> list[Path]:
    selected = series_paths[:max_runs]
    output_paths: list[Path] = []
    if not selected:
        return output_paths

    for metric, ylabel in METRICS:
        fig, ax = plt.subplots(figsize=(8, 4.8))
        plotted = False
        for path in selected:
            rows = _read_rows(path)
            if not rows or metric not in rows[0]:
                continue
            steps = np.array([_safe_float(row.get("step")) for row in rows])
            values = np.array([_safe_float(row.get(metric)) for row in rows])
            mask = np.isfinite(steps) & np.isfinite(values)
            if not np.any(mask):
                continue
            ax.plot(steps[mask], values[mask], lw=1.7, label=_run_label(path))
            plotted = True
        if not plotted:
            plt.close(fig)
            continue
        ax.set_xlabel("Step")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.grid(True, alpha=0.3)
        if len(selected) <= 8:
            ax.legend(fontsize=8)
        fig.tight_layout()
        out_path = output_dir / f"{metric}.png"
        fig.savefig(out_path, dpi=180)
        plt.close(fig)
        output_paths.append(out_path)
    return output_paths


def _plot_final_scatter(output_dir: Path, dataset_path: Path) -> Path | None:
    rows = _read_rows(dataset_path)
    if not rows:
        return None
    x = np.array([_safe_float(row.get("final_retained_particle_ratio")) for row in rows])
    y = np.array([_safe_float(row.get("final_pressure_drop_ratio")) for row in rows])
    c = np.array([_safe_float(row.get("final_normalized_permeate_flux")) for row in rows])
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
    if not np.any(mask):
        return None
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(x[mask], y[mask], c=c[mask], cmap="viridis", edgecolor="black", linewidth=0.4)
    ax.set_xlabel("Retained particle ratio")
    ax.set_ylabel("Global pressure drop ratio")
    ax.set_title("Final fouling state")
    ax.grid(True, alpha=0.3)
    fig.colorbar(sc, ax=ax, label="Normalized permeate flux")
    fig.tight_layout()
    out_path = output_dir / "final_pressure_drop_vs_retention.png"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "roots",
        nargs="*",
        type=Path,
        default=[Path("src/results/run_lbm_dem")],
        help="Result roots or analysis/time_series.csv files to plot",
    )
    parser.add_argument("--max-runs", type=int, default=12)
    parser.add_argument("--result-tag", default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()
    if args.max_runs <= 0:
        parser.error("--max-runs must be positive")
    return args


def main() -> int:
    args = parse_args()
    timestamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    parts = [args.result_tag, timestamp] if args.result_tag else [timestamp]
    output_dir = args.out_dir or program_results_dir(__file__, *parts)
    output_dir.mkdir(parents=True, exist_ok=True)

    series_paths = _find_time_series(args.roots)
    dataset_path = _write_dataset(output_dir, series_paths)
    plot_paths = _plot_time_series(output_dir, series_paths, args.max_runs)
    scatter_path = _plot_final_scatter(output_dir, dataset_path)

    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "roots": [str(root) for root in args.roots],
        "n_time_series": len(series_paths),
        "dataset": str(dataset_path),
        "plots": [str(path) for path in plot_paths],
        "final_scatter": str(scatter_path) if scatter_path is not None else None,
    }
    (output_dir / "fouling_metrics_plot_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(f"Found time-series files: {len(series_paths)}")
    print(f"Dataset: {dataset_path}")
    for path in plot_paths:
        print(f"Plot: {path}")
    if scatter_path is not None:
        print(f"Plot: {scatter_path}")
    print(f"Report: {output_dir / 'fouling_metrics_plot_report.json'}")
    return 0 if series_paths else 1


if __name__ == "__main__":
    raise SystemExit(main())
