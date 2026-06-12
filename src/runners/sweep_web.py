#!/usr/bin/env python3
"""Web UI backend for monitoring and controlling parameter sweeps.

A thin FastAPI control+view layer over ``sweep_cli.py``. The "Run" action
launches the existing ``./bin/sweep run`` command as a background subprocess
(the CLI already does sequential execution, STARTED/FAILED markers, light-result
harvest, and git push); the web layer only tracks the PID, tails the log, counts
VTK snapshot files for a live progress %, and can kill the process group to stop.

Control is LOCAL (this machine launches/stops runs); the progress *view* is
multi-machine because status is read from the git-synced ``sweep_results/`` dir.

Launch via ``./bin/sweep-web`` (needs ``poetry install --with web``).
"""

from __future__ import annotations

import json
import math
import os
import signal
import subprocess
import time
from datetime import datetime
from pathlib import Path

from runners.sweep_cli import (  # noqa: F401  (re-exported constants used in tests)
    MANIFESTS_DIR,
    REPO_ROOT,
    RESULTS_ROOT,
    SHARED_ROOT,
    case_name_for,
    case_status,
    find_completed_run_dir,
    load_manifest,
    resolve_cases,
)

WEB_DIR = REPO_ROOT / ".sweep_web"
LOG_DIR = WEB_DIR / "logs"
REGISTRY = WEB_DIR / "runs.json"
STATIC_DIR = Path(__file__).resolve().parent / "sweep_web_static"

# Mirrors the heavy artifacts the runner writes one-per-snapshot; counting these
# in the active run dir gives a live progress estimate without touching the solver.
_FLUID_GLOB = "paraview/fluid_*.vtk"

# mtime-keyed cache for expected_snapshots (configs don't change mid-run).
_expected_cache: dict[Path, tuple[float, int]] = {}

# run_id -> Popen for processes launched in THIS web-app process. Polled in
# _reconcile() so finished children are reaped (otherwise they linger as
# zombies and os.kill(pid, 0) keeps reporting them alive).
_PROCS: dict[str, "subprocess.Popen"] = {}


# ---------------------------------------------------------------------------
# Progress (non-invasive: derived from on-disk snapshot files)
# ---------------------------------------------------------------------------


def active_run_dir(case_name: str, results_root: Path | None = None) -> Path | None:
    """Newest run dir for a case that is NOT yet completed (running or crashed).

    Returns None when the newest run is already completed (nothing active) or
    when no run dir exists. ``results_root`` defaults to the module-level
    RESULTS_ROOT, resolved at call time so tests can monkeypatch it.
    """
    if results_root is None:
        results_root = RESULTS_ROOT
    candidates = []
    for parent in (results_root / "dim3" / case_name, results_root / case_name):
        if parent.is_dir():
            candidates.extend(parent.glob("run_*"))
    if not candidates:
        return None
    newest = max(candidates, key=lambda p: p.name)
    status_file = newest / "run_status.json"
    if status_file.is_file():
        try:
            if json.loads(status_file.read_text()).get("status") == "completed":
                return None
        except json.JSONDecodeError:
            pass
    return newest


def expected_snapshots(config_path: Path) -> int:
    """Expected snapshot count = ceil(total_steps / snapshot_every).

    Resolves the ``extends`` chain via SimulationConfig so inherited runtime
    sections (variants that don't set runtime inline) are read correctly.
    """
    config_path = Path(config_path)
    try:
        mtime = config_path.stat().st_mtime
    except OSError:
        mtime = 0.0
    cached = _expected_cache.get(config_path)
    if cached and cached[0] == mtime:
        return cached[1]

    # Imported lazily so the pure-progress helpers don't pull the solver stack
    # until a config actually needs resolving.
    from particulate_flow.io.config import SimulationConfig

    values = SimulationConfig.from_json(config_path).values
    total = int(values.get("total_steps", 0))
    every = max(1, int(values.get("snapshot_every", 1)))
    expected = math.ceil(total / every) if total > 0 else 0
    _expected_cache[config_path] = (mtime, expected)
    return expected


def live_progress(case_name: str, config_path: Path) -> dict:
    """Progress for a (possibly) running case, from snapshot file count."""
    expected = expected_snapshots(config_path)
    run_dir = active_run_dir(case_name)
    done = len(list(run_dir.glob(_FLUID_GLOB))) if run_dir else 0
    percent = round(100 * done / expected, 1) if expected else 0.0
    return {
        "snapshots_done": done,
        "snapshots_expected": expected,
        "percent": percent,
        "run_dir": _display_path(run_dir) if run_dir else None,
    }


def _display_path(path: Path) -> str:
    """Path relative to the repo when possible, else absolute (test-safe)."""
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


# ---------------------------------------------------------------------------
# Per-case settings table (compare conditions of a sweep)
# ---------------------------------------------------------------------------

# Meta keys excluded from the settings comparison (identity, not parameters).
_META_KEYS = {"name", "description", "result_tag", "source_path", "schema_version"}


def _jsonable(value):
    """Coerce a flattened-config value into something JSON-serializable."""
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    return str(value)


def sweep_config_table(name: str) -> dict:
    """Resolve every case config and report values + which keys vary.

    Returns: {name, cases:[{alias, case_name, values:{...}, error?}],
              varying_keys:[...], all_keys:[...]}.
    The ``extends`` chain is flattened via SimulationConfig so inherited
    parameters (e.g. geometry fragments) are included.
    """
    from particulate_flow.io.config import SimulationConfig

    manifest = load_manifest(name)
    cases = []
    for entry in resolve_cases(manifest, ["all"]):
        case = {"alias": entry["alias"], "case_name": entry["case_name"]}
        try:
            values = SimulationConfig.from_json(entry["config"]).values
            case["values"] = {
                k: _jsonable(v) for k, v in values.items() if k not in _META_KEYS
            }
        except Exception as exc:  # noqa: BLE001 - surface config errors in the UI
            case["values"] = {}
            case["error"] = str(exc)
        cases.append(case)

    all_keys = sorted({k for c in cases for k in c["values"]})
    varying = []
    for key in all_keys:
        seen = {
            json.dumps(c["values"].get(key), sort_keys=True, default=str) for c in cases
        }
        if len(seen) > 1:
            varying.append(key)

    return {
        "name": manifest["name"],
        "cases": cases,
        "varying_keys": varying,
        "all_keys": all_keys,
    }


# ---------------------------------------------------------------------------
# Process registry
# ---------------------------------------------------------------------------


def _load_registry() -> dict:
    if not REGISTRY.is_file():
        return {}
    try:
        return json.loads(REGISTRY.read_text())
    except (json.JSONDecodeError, OSError):
        return {}


def _save_registry(reg: dict) -> None:
    WEB_DIR.mkdir(parents=True, exist_ok=True)
    tmp = REGISTRY.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(reg, indent=2) + "\n")
    os.replace(tmp, REGISTRY)


def _pid_alive(pid: int) -> bool:
    """True if the PID is a live (non-zombie) process.

    A zombie still answers ``os.kill(pid, 0)`` until reaped, so we additionally
    reject processes in state ``Z`` (defunct children of a still-running parent
    that we couldn't poll, e.g. after a web-app restart — init reaps these).
    """
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    stat = Path(f"/proc/{pid}/stat")
    if stat.exists():
        try:
            state = stat.read_text().rsplit(")", 1)[1].split()[0]
            if state == "Z":
                return False
        except (OSError, IndexError):
            pass
    return True


def _reconcile() -> dict:
    """Mark registry entries whose process has exited as finished; persist if changed."""
    # Reap any children we launched in this process (clears zombies, gives the
    # authoritative exit code) before consulting the registry.
    for run_id, proc in list(_PROCS.items()):
        if proc.poll() is not None:
            _PROCS.pop(run_id, None)

    reg = _load_registry()
    changed = False
    for entry in reg.values():
        if entry.get("status") == "running" and not _pid_alive(entry.get("pid", -1)):
            entry["status"] = "finished"
            changed = True
    if changed:
        _save_registry(reg)
    return reg


def _make_run_id(sweep: str, case_name: str) -> str:
    return f"{sweep}__{case_name}__{int(time.time() * 1000)}"


# ---------------------------------------------------------------------------
# Results parsing
# ---------------------------------------------------------------------------

# summary.json keys surfaced as headline metrics in the UI (others still shown).
_METRIC_KEYS = [
    "final_particle_count",
    "generated_particles",
    "removed_particles",
    "max_speed",
    "total_steps",
    "n_cylinders",
]


def _parse_time_series(csv_path: Path, max_points: int = 300) -> dict:
    """Parse time_series.csv into column arrays, downsampled to max_points."""
    lines = csv_path.read_text().strip().splitlines()
    if not lines:
        return {"columns": []}
    columns = lines[0].split(",")
    rows = [ln.split(",") for ln in lines[1:]]
    if len(rows) > max_points:
        step = math.ceil(len(rows) / max_points)
        rows = rows[::step]
    data: dict[str, list] = {col: [] for col in columns}
    for row in rows:
        for col, val in zip(columns, row):
            try:
                data[col].append(float(val))
            except ValueError:
                data[col].append(val)
    return {"columns": columns, **data}


def result_payload(sweep: str, case_name: str) -> dict:
    """Result dir + parsed summary/time-series, preferring git-synced harvest.

    Falls back to the local completed run dir when not yet harvested.
    """
    harvested = SHARED_ROOT / sweep / case_name
    source = harvested if (harvested / "summary.json").is_file() else None
    if source is None:
        local = find_completed_run_dir(case_name)
        source = local
    if source is None:
        return {"error": f"no completed results for {case_name}"}

    summary = None
    summary_file = source / "summary.json"
    if summary_file.is_file():
        summary = json.loads(summary_file.read_text())
    metrics = {k: summary[k] for k in _METRIC_KEYS if summary and k in summary}

    time_series = {"columns": []}
    ts_file = source / "analysis" / "time_series.csv"
    if ts_file.is_file():
        time_series = _parse_time_series(ts_file)

    return {
        "result_dir": str(source.resolve()),
        "summary": summary,
        "metrics": metrics,
        "time_series": time_series,
    }


# ---------------------------------------------------------------------------
# Run / stop
# ---------------------------------------------------------------------------


def _build_run_command(sweep: str, cases: list[str], options: dict) -> list[str]:
    cmd = ["./bin/sweep", "run", sweep, *cases]
    if options.get("force"):
        cmd.append("--force")
    if options.get("no_push"):
        cmd.append("--no-push")
    if options.get("no_install"):
        cmd.append("--no-install")
    return cmd


def start_run(sweep: str, cases: list[str], options: dict) -> dict:
    """Validate, launch ``./bin/sweep run`` as a process group, register it.

    Raises ValueError on validation failure (unknown case / overlap).
    """
    manifest = load_manifest(sweep)
    resolved = resolve_cases(manifest, cases)  # raises SystemExit on unknown
    resolved_names = {e["case_name"] for e in resolved}

    reg = _reconcile()
    for entry in reg.values():
        if entry.get("status") == "running" and resolved_names.intersection(
            entry.get("case_names", [])
        ):
            raise ValueError(
                f"a run is already active for: "
                f"{sorted(resolved_names.intersection(entry.get('case_names', [])))}"
            )

    primary = resolved[0]["case_name"]
    run_id = _make_run_id(sweep, primary)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"{run_id}.log"
    cmd = _build_run_command(sweep, cases, options)

    log_fh = open(log_path, "w")
    proc = subprocess.Popen(
        cmd,
        cwd=REPO_ROOT,
        stdout=log_fh,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )

    reg[run_id] = {
        "run_id": run_id,
        "sweep": sweep,
        "cases": cases,
        "case_names": sorted(resolved_names),
        "pid": proc.pid,
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "log_path": str(log_path),
        "status": "running",
        "cmd": " ".join(cmd),
    }
    _save_registry(reg)
    _PROCS[run_id] = proc
    return reg[run_id]


def stop_run(run_id: str) -> bool:
    """SIGTERM (then SIGKILL) the process group for run_id. Returns True if signalled."""
    reg = _load_registry()
    entry = reg.get(run_id)
    if not entry:
        return False
    pid = entry.get("pid", -1)
    if pid > 0 and _pid_alive(pid):
        try:
            pgid = os.getpgid(pid)
            os.killpg(pgid, signal.SIGTERM)
            for _ in range(20):  # up to ~2s grace
                if not _pid_alive(pid):
                    break
                time.sleep(0.1)
            if _pid_alive(pid):
                os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    entry["status"] = "stopped"
    _save_registry(reg)
    return True


def tail_log(run_id: str, tail: int = 200) -> dict:
    reg = _load_registry()
    entry = reg.get(run_id)
    if not entry:
        return {"error": f"unknown run_id {run_id}"}
    tail = max(1, min(int(tail), 2000))
    log_path = Path(entry["log_path"])
    lines: list[str] = []
    if log_path.is_file():
        lines = log_path.read_text(errors="replace").splitlines()[-tail:]
    return {"run_id": run_id, "lines": lines, "status": entry.get("status")}


# ---------------------------------------------------------------------------
# FastAPI app
# ---------------------------------------------------------------------------


def create_app():
    from fastapi import FastAPI
    from fastapi.responses import JSONResponse
    from fastapi.staticfiles import StaticFiles

    app = FastAPI(title="Sweep Web UI")

    def _running_case_index() -> dict[str, dict]:
        """Map case_name -> registry entry for currently-running local runs."""
        reg = _reconcile()
        index = {}
        for entry in reg.values():
            if entry.get("status") == "running":
                for name in entry.get("case_names", []):
                    index[name] = entry
        return index

    @app.get("/api/sweeps")
    def list_sweeps():
        out = []
        for path in sorted(MANIFESTS_DIR.glob("*.json")):
            manifest = load_manifest(path.stem)
            out.append(
                {
                    "name": manifest["name"],
                    "description": manifest.get("description", ""),
                    "case_count": len(manifest["cases"]),
                }
            )
        return out

    @app.get("/api/sweeps/{name}")
    def sweep_detail(name: str):
        try:
            manifest = load_manifest(name)
        except SystemExit as exc:
            return JSONResponse({"error": str(exc)}, status_code=404)
        running = _running_case_index()
        cases = []
        for entry in resolve_cases(manifest, ["all"]):
            case_name = entry["case_name"]
            st = case_status(manifest["name"], case_name)
            run_entry = running.get(case_name)
            local_running = run_entry is not None
            progress = None
            if local_running or st["status"] == "in-progress":
                progress = live_progress(case_name, entry["config"])
            cases.append(
                {
                    "alias": entry["alias"],
                    "case_name": case_name,
                    "config": str(entry["config"].relative_to(REPO_ROOT)),
                    "status": st["status"],
                    "machine": st["machine"],
                    "when": st["when"],
                    "local_running": local_running,
                    "run_id": run_entry["run_id"] if run_entry else None,
                    "progress": progress,
                }
            )
        return {"name": manifest["name"], "cases": cases}

    @app.get("/api/sweeps/{name}/config")
    def sweep_config(name: str):
        try:
            return sweep_config_table(name)
        except SystemExit as exc:
            return JSONResponse({"error": str(exc)}, status_code=404)

    @app.post("/api/run")
    def post_run(body: dict):
        sweep = body.get("sweep")
        cases = body.get("cases") or []
        options = body.get("options") or {}
        if not sweep or not cases:
            return JSONResponse({"error": "sweep and cases are required"}, status_code=400)
        try:
            entry = start_run(sweep, cases, options)
        except SystemExit as exc:
            return JSONResponse({"error": str(exc)}, status_code=400)
        except ValueError as exc:
            return JSONResponse({"error": str(exc)}, status_code=409)
        return {
            "run_id": entry["run_id"],
            "pid": entry["pid"],
            "log_path": entry["log_path"],
            "cmd": entry["cmd"],
        }

    @app.post("/api/stop")
    def post_stop(body: dict):
        run_id = body.get("run_id")
        if not run_id:
            return JSONResponse({"error": "run_id is required"}, status_code=400)
        if not stop_run(run_id):
            return JSONResponse({"error": f"unknown run_id {run_id}"}, status_code=404)
        return {"stopped": True, "run_id": run_id}

    @app.get("/api/runs")
    def list_runs():
        reg = _reconcile()
        return list(reg.values())

    @app.get("/api/logs/{run_id}")
    def get_logs(run_id: str, tail: int = 200):
        result = tail_log(run_id, tail)
        if "error" in result:
            return JSONResponse(result, status_code=404)
        return result

    @app.get("/api/results/{sweep}/{case}")
    def get_results(sweep: str, case: str):
        payload = result_payload(sweep, case)
        if "error" in payload:
            return JSONResponse(payload, status_code=404)
        return payload

    @app.get("/api/results/{sweep}/{case}/open")
    def open_results(sweep: str, case: str):
        payload = result_payload(sweep, case)
        if "error" in payload:
            return JSONResponse(payload, status_code=404)
        path = payload["result_dir"]
        opened = False
        try:
            subprocess.Popen(["xdg-open", path])
            opened = True
        except Exception:
            opened = False
        return {"path": path, "opened": opened}

    # Mount static LAST so /api/* routes are not shadowed.
    if STATIC_DIR.is_dir():
        app.mount("/", StaticFiles(directory=STATIC_DIR, html=True), name="static")

    return app


# uvicorn target: ``runners.sweep_web:app``
app = create_app()
