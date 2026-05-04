#!/usr/bin/env python3
"""Plot summary figures for the 168-condition LBM-DEM sweep."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[2]
BASE_SWEEP = REPO_ROOT / "src/results/run_lbm_dem/sweep_56_conditions_3x"
ROLLING_SWEEP = (
    REPO_ROOT / "src/results/run_lbm_dem/sweep_56_conditions_rolling_friction_2x_3x"
)
OUTPUT_DIR = REPO_ROOT / "src/results/run_lbm_dem/sweep_168_analysis"

FRAME_RE = re.compile(
    r"frame\s+(\d+).*?\|u\|_max=([0-9.eE+-]+).*?"
    r"KE_p=([0-9.eE+-]+).*?\|F\|_mean=([0-9.eE+-]+)"
)


def parse_tag(name: str, path: Path) -> dict[str, object] | None:
    if name.startswith("t3x_muR2x_"):
        rolling_multiplier = 2
        tag = name.removeprefix("t3x_muR2x_")
    elif name.startswith("t3x_muR3x_"):
        rolling_multiplier = 3
        tag = name.removeprefix("t3x_muR3x_")
    elif name.startswith("t3x_") and BASE_SWEEP in path.parents:
        rolling_multiplier = 1
        tag = name.removeprefix("t3x_")
    else:
        return None

    phi_match = re.search(r"_(05|10|20|40)pct$", tag)
    if phi_match is None:
        return None
    phi = int(phi_match.group(1))
    rolling = "_rolling_" in tag and "_free_roll_" not in tag

    if tag.startswith("surface_none_") or tag.startswith("none_"):
        force_type = "none"
        force_multiplier = 0.0
        force_label = "none"
    elif tag.startswith("attraction_0p5x_"):
        force_type = "attraction"
        force_multiplier = 0.5
        force_label = "attraction 0.5x"
    elif tag.startswith("attraction_2x_"):
        force_type = "attraction"
        force_multiplier = 2.0
        force_label = "attraction 2x"
    elif tag.startswith("attraction_"):
        force_type = "attraction"
        force_multiplier = 1.0
        force_label = "attraction 1x"
    elif tag.startswith("repulsion_0p5x_"):
        force_type = "repulsion"
        force_multiplier = 0.5
        force_label = "repulsion 0.5x"
    elif tag.startswith("repulsion_2x_"):
        force_type = "repulsion"
        force_multiplier = 2.0
        force_label = "repulsion 2x"
    elif tag.startswith("repulsion_"):
        force_type = "repulsion"
        force_multiplier = 1.0
        force_label = "repulsion 1x"
    else:
        return None

    return {
        "tag": tag,
        "phi": phi,
        "rolling": rolling,
        "rolling_multiplier": rolling_multiplier,
        "force_type": force_type,
        "force_multiplier": force_multiplier,
        "force_label": force_label,
    }


def parse_log(path: Path) -> dict[str, object] | None:
    meta = parse_tag(path.stem, path)
    if meta is None:
        return None

    text = path.read_text(errors="replace")
    frames = [
        (int(match.group(1)), float(match.group(2)), float(match.group(3)), float(match.group(4)))
        for match in FRAME_RE.finditer(text)
    ]
    if not frames:
        return None

    last_n = min(50, len(frames))
    last = frames[-last_n:]
    outlet_match = re.search(r"右出口削除粒子数:\s*(\d+)", text)
    queued_match = re.search(r"up to\s+(\d+)\s+queued particles", text)

    row = {
        "path": str(path.relative_to(REPO_ROOT)),
        "name": path.stem,
        "n_frames": len(frames),
        "u_last50": sum(frame[1] for frame in last) / last_n,
        "ke_last50": sum(frame[2] for frame in last) / last_n,
        "force_last50": sum(frame[3] for frame in last) / last_n,
        "u_final": frames[-1][1],
        "outlet_removed": int(outlet_match.group(1)) if outlet_match else None,
        "queued_particles": int(queued_match.group(1)) if queued_match else None,
    }
    row.update(meta)
    if row["outlet_removed"] is not None and row["queued_particles"]:
        row["pass_ratio"] = row["outlet_removed"] / row["queued_particles"]
    else:
        row["pass_ratio"] = None
    return row


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def collect_rows() -> list[dict[str, object]]:
    rows = []
    for root in (BASE_SWEEP, ROLLING_SWEEP):
        for path in sorted(root.rglob("*.log")):
            row = parse_log(path)
            if row is not None:
                rows.append(row)
    return rows


def grouped(rows: list[dict[str, object]], keys: tuple[str, ...]):
    bins = defaultdict(list)
    for row in rows:
        bins[tuple(row[key] for key in keys)].append(row)
    return bins


def save_csv(rows: list[dict[str, object]]) -> None:
    fields = [
        "name",
        "phi",
        "force_label",
        "force_type",
        "force_multiplier",
        "rolling",
        "rolling_multiplier",
        "u_last50",
        "ke_last50",
        "force_last50",
        "outlet_removed",
        "queued_particles",
        "pass_ratio",
        "path",
    ]
    with (OUTPUT_DIR / "summary_168_conditions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in fields})


def plot_phi_force(rows: list[dict[str, object]]) -> None:
    order = [
        "none",
        "attraction 0.5x",
        "attraction 1x",
        "attraction 2x",
        "repulsion 0.5x",
        "repulsion 1x",
        "repulsion 2x",
    ]
    colors = {
        "none": "#4b5563",
        "attraction 0.5x": "#f59e0b",
        "attraction 1x": "#dc2626",
        "attraction 2x": "#7f1d1d",
        "repulsion 0.5x": "#38bdf8",
        "repulsion 1x": "#2563eb",
        "repulsion 2x": "#1e3a8a",
    }
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    force_bins = grouped(rows, ("phi", "force_label"))
    phis = [5, 10, 20, 40]
    for label in order:
        y = []
        for phi in phis:
            bucket = force_bins[(phi, label)]
            y.append(mean([row["u_last50"] for row in bucket]))
        ax.plot(phis, y, marker="o", linewidth=2, color=colors[label], label=label)

    ax.set_title("Flow reduction appears only when attraction and high loading coincide")
    ax.set_xlabel("Particle volume fraction in feed (%)")
    ax.set_ylabel("Mean |u|max over last 50 frames")
    ax.set_xticks(phis)
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "01_flow_vs_fraction_by_interaction.png", dpi=180)
    plt.close(fig)


def plot_attraction_rolling(rows: list[dict[str, object]]) -> None:
    subset = [
        row
        for row in rows
        if row["phi"] == 40 and row["force_type"] == "attraction" and row["rolling"]
    ]
    strengths = [0.5, 1.0, 2.0]
    mu_values = [1, 2, 3]
    width = 0.24
    x = list(range(len(strengths)))
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    colors = {1: "#f97316", 2: "#ef4444", 3: "#7f1d1d"}
    for i, mu in enumerate(mu_values):
        y = []
        for strength in strengths:
            bucket = [
                row
                for row in subset
                if row["force_multiplier"] == strength and row["rolling_multiplier"] == mu
            ]
            y.append(bucket[0]["u_last50"])
        positions = [value + (i - 1) * width for value in x]
        ax.bar(positions, y, width=width, color=colors[mu], label=f"rolling friction {mu}x")
    ax.set_title("At 40%, rolling friction stabilizes attractive deposits")
    ax.set_xlabel("Attraction strength")
    ax.set_ylabel("Mean |u|max over last 50 frames")
    ax.set_xticks(x)
    ax.set_xticklabels(["0.5x", "1x", "2x"])
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "02_40pct_attraction_rolling_friction.png", dpi=180)
    plt.close(fig)


def plot_pass_ratio(rows: list[dict[str, object]]) -> None:
    subset = [row for row in rows if row["phi"] == 40 and row["pass_ratio"] is not None]
    labels = ["none", "attraction 0.5x", "attraction 1x", "attraction 2x", "repulsion 0.5x", "repulsion 1x", "repulsion 2x"]
    series = [
        ("free rolling", False, None, "#94a3b8"),
        ("rolling friction 1x", True, 1, "#14b8a6"),
        ("rolling friction 2x", True, 2, "#0f766e"),
        ("rolling friction 3x", True, 3, "#134e4a"),
    ]
    x = list(range(len(labels)))
    width = 0.2
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    for index, (series_label, rolling, rolling_multiplier, color) in enumerate(series):
        values = []
        for label in labels:
            bucket = [
                row["pass_ratio"]
                for row in subset
                if row["force_label"] == label
                and row["rolling"] == rolling
                and (rolling_multiplier is None or row["rolling_multiplier"] == rolling_multiplier)
            ]
            values.append(mean(bucket))
        offset = (index - 1.5) * width
        ax.bar([value + offset for value in x], values, width, label=series_label, color=color)
    ax.set_title("At 40%, attraction traps particles while repulsion preserves throughput")
    ax.set_xlabel("Particle-particle interaction")
    ax.set_ylabel("Right-outlet pass ratio")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=32, ha="right")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "03_40pct_pass_ratio_by_interaction.png", dpi=180)
    plt.close(fig)


def plot_heatmap(rows: list[dict[str, object]]) -> None:
    labels = [
        "repulsion 2x",
        "repulsion 1x",
        "repulsion 0.5x",
        "none",
        "attraction 0.5x",
        "attraction 1x",
        "attraction 2x",
    ]
    phis = [5, 10, 20, 40]
    data = []
    for label in labels:
        line = []
        for phi in phis:
            bucket = [
                row["u_last50"]
                for row in rows
                if row["phi"] == phi and row["force_label"] == label and row["rolling"]
            ]
            line.append(mean(bucket))
        data.append(line)

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    vmin = min(map(min, data))
    vmax = max(map(max, data))
    image = ax.imshow(data, cmap="viridis", aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_title("Rolling-friction cases: mean of 1x/2x/3x")
    ax.set_xlabel("Particle volume fraction (%)")
    ax.set_ylabel("Interaction")
    ax.set_xticks(range(len(phis)))
    ax.set_xticklabels([str(phi) for phi in phis])
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    for y, line in enumerate(data):
        for x, value in enumerate(line):
            normalized = (value - vmin) / (vmax - vmin)
            text_color = "white" if normalized < 0.45 else "black"
            ax.text(x, y, f"{value:.4f}", ha="center", va="center", fontsize=8, color=text_color)
    fig.colorbar(image, ax=ax, label="Mean |u|max over last 50 frames")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "04_low_flow_heatmap_rolling_cases.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = collect_rows()
    if len(rows) != 168:
        raise SystemExit(f"Expected 168 parsed logs, got {len(rows)}")
    save_csv(rows)
    plot_phi_force(rows)
    plot_attraction_rolling(rows)
    plot_pass_ratio(rows)
    plot_heatmap(rows)
    print(f"Wrote 4 figures and summary CSV to {OUTPUT_DIR.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
