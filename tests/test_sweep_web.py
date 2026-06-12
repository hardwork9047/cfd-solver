"""Tests for the sweep web backend (src/runners/sweep_web.py).

No real simulations are launched: subprocess.Popen and os.killpg are mocked.
"""

import json

import pytest

pytest.importorskip("fastapi")

from fastapi.testclient import TestClient  # noqa: E402

from runners import sweep_web  # noqa: E402

CONFIGS = sweep_web.REPO_ROOT / "configs" / "lbm_dem" / "cases"


def _write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


@pytest.fixture
def client():
    return TestClient(sweep_web.app)


@pytest.fixture(autouse=True)
def _isolate_registry(tmp_path, monkeypatch):
    """Point the registry/log dir at tmp so tests never touch real .sweep_web/."""
    monkeypatch.setattr(sweep_web, "WEB_DIR", tmp_path / "web")
    monkeypatch.setattr(sweep_web, "LOG_DIR", tmp_path / "web" / "logs")
    monkeypatch.setattr(sweep_web, "REGISTRY", tmp_path / "web" / "runs.json")
    sweep_web._expected_cache.clear()
    sweep_web._PROCS.clear()


# ---------------------------------------------------------------------------
# expected_snapshots / live_progress / active_run_dir
# ---------------------------------------------------------------------------


def test_expected_snapshots_inline():
    assert sweep_web.expected_snapshots(CONFIGS / "fouling_3d_smoke.json") == 2


def test_expected_snapshots_inherited():
    # a1 inherits runtime (125000/1250) from b0 via extends -> proves resolution.
    assert (
        sweep_web.expected_snapshots(CONFIGS / "fouling_3d_sweep_20260609_a1_attr_0015.json")
        == 100
    )


def _make_run(results_root, case, run_name, *, completed=False, n_vtk=0):
    run_dir = results_root / "dim3" / case / run_name
    (run_dir / "paraview").mkdir(parents=True, exist_ok=True)
    for i in range(n_vtk):
        (run_dir / "paraview" / f"fluid_{i:06d}.vtk").write_text("x")
    if completed:
        _write_json(run_dir / "run_status.json", {"status": "completed"})
    return run_dir


def test_live_progress_counts_fake_vtks(tmp_path, monkeypatch):
    monkeypatch.setattr(sweep_web, "RESULTS_ROOT", tmp_path)
    _make_run(tmp_path, "smoke3d", "run_20260610_000000_000000", n_vtk=1)
    prog = sweep_web.live_progress("smoke3d", CONFIGS / "fouling_3d_smoke.json")
    assert prog["snapshots_done"] == 1
    assert prog["snapshots_expected"] == 2
    assert prog["percent"] == 50.0
    assert prog["run_dir"] is not None


def test_live_progress_no_active_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(sweep_web, "RESULTS_ROOT", tmp_path)
    prog = sweep_web.live_progress("smoke3d", CONFIGS / "fouling_3d_smoke.json")
    assert prog["snapshots_done"] == 0
    assert prog["run_dir"] is None
    assert prog["percent"] == 0.0


def test_active_run_dir_ignores_completed(tmp_path, monkeypatch):
    monkeypatch.setattr(sweep_web, "RESULTS_ROOT", tmp_path)
    _make_run(tmp_path, "c", "run_20260610_000000_000000", completed=True, n_vtk=2)
    assert sweep_web.active_run_dir("c", results_root=tmp_path) is None
    running = _make_run(tmp_path, "c", "run_20260611_000000_000000", n_vtk=1)
    assert sweep_web.active_run_dir("c", results_root=tmp_path) == running


# ---------------------------------------------------------------------------
# API: sweeps / detail
# ---------------------------------------------------------------------------


def test_api_sweeps_lists_smoke(client):
    names = {s["name"]: s for s in client.get("/api/sweeps").json()}
    assert "smoke_test" in names
    assert names["smoke_test"]["case_count"] == 1


def test_api_sweep_detail_status(client):
    data = client.get("/api/sweeps/smoke_test").json()
    assert data["name"] == "smoke_test"
    case = data["cases"][0]
    assert case["case_name"] == "smoke3d"
    assert case["status"] in {"done", "failed", "in-progress", "pending"}


def test_api_sweep_detail_unknown_404(client):
    assert client.get("/api/sweeps/does_not_exist").status_code == 404


# ---------------------------------------------------------------------------
# API: config / conditions table
# ---------------------------------------------------------------------------


def test_sweep_config_table_detects_varying_keys():
    table = sweep_web.sweep_config_table("fouling_3d_20260609")
    assert len(table["cases"]) == 12
    # attraction_strength and particle_volume_fraction differ across cases.
    assert "attraction_strength" in table["varying_keys"]
    assert "particle_volume_fraction" in table["varying_keys"]
    # result_tag is identity, not a parameter -> excluded.
    assert "result_tag" not in table["varying_keys"]
    assert "result_tag" not in table["all_keys"]
    # varying is a subset of all.
    assert set(table["varying_keys"]).issubset(table["all_keys"])


def test_api_sweep_config_single_case_has_no_varying():
    table = sweep_web.sweep_config_table("smoke_test")
    assert len(table["cases"]) == 1
    # One case -> nothing varies.
    assert table["varying_keys"] == []
    assert table["all_keys"]  # but it still reports its parameters


def test_api_sweep_config_endpoint(client):
    data = client.get("/api/sweeps/fouling_3d_20260609/config").json()
    assert data["name"] == "fouling_3d_20260609"
    assert data["cases"][0]["values"]  # flattened params present


def test_api_sweep_config_unknown_404(client):
    assert client.get("/api/sweeps/does_not_exist/config").status_code == 404


# ---------------------------------------------------------------------------
# API: results parsing
# ---------------------------------------------------------------------------


def test_api_results_parses_summary_and_timeseries(tmp_path, monkeypatch, client):
    monkeypatch.setattr(sweep_web, "SHARED_ROOT", tmp_path)
    case_dir = tmp_path / "smoke_test" / "smoke3d"
    _write_json(case_dir / "summary.json", {"final_particle_count": 42, "max_speed": 0.01})
    (case_dir / "analysis").mkdir(parents=True)
    (case_dir / "analysis" / "time_series.csv").write_text(
        "step,particle_count\n100,5\n200,9\n"
    )
    data = client.get("/api/results/smoke_test/smoke3d").json()
    assert data["metrics"]["final_particle_count"] == 42
    assert data["time_series"]["columns"] == ["step", "particle_count"]
    assert data["time_series"]["particle_count"] == [5.0, 9.0]


# ---------------------------------------------------------------------------
# API: run / stop (mocked subprocess)
# ---------------------------------------------------------------------------


class _FakePopen:
    def __init__(self, *a, **k):
        self.pid = 999999

    def poll(self):
        return None  # pretend still running (no real process to reap)


def test_api_run_validates_unknown_case(monkeypatch, client):
    spy = {"called": False}

    def _spy(*a, **k):
        spy["called"] = True
        return _FakePopen()

    monkeypatch.setattr(sweep_web.subprocess, "Popen", _spy)
    resp = client.post("/api/run", json={"sweep": "smoke_test", "cases": ["bogus"]})
    assert resp.status_code == 400
    assert spy["called"] is False


def test_api_run_launch_is_mocked(monkeypatch, client):
    monkeypatch.setattr(sweep_web.subprocess, "Popen", _FakePopen)
    resp = client.post(
        "/api/run",
        json={"sweep": "smoke_test", "cases": ["smoke"], "options": {"no_push": True}},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["pid"] == 999999
    runs = client.get("/api/runs").json()
    assert any(r["run_id"] == body["run_id"] for r in runs)


def test_stop_killpg_called(monkeypatch, client):
    # Seed a registry entry with a fake live pid.
    sweep_web._save_registry(
        {
            "rid": {
                "run_id": "rid",
                "sweep": "smoke_test",
                "cases": ["smoke"],
                "case_names": ["smoke3d"],
                "pid": 4242,
                "status": "running",
                "log_path": ".sweep_web/logs/rid.log",
            }
        }
    )
    killed = {}
    monkeypatch.setattr(sweep_web, "_pid_alive", lambda pid: pid == 4242 and not killed)
    monkeypatch.setattr(sweep_web.os, "getpgid", lambda pid: pid)

    def _killpg(pgid, sig):
        killed["pgid"] = pgid

    monkeypatch.setattr(sweep_web.os, "killpg", _killpg)
    resp = client.post("/api/stop", json={"run_id": "rid"})
    assert resp.status_code == 200
    assert killed["pgid"] == 4242


def test_stop_unknown_run_404(client):
    assert client.post("/api/stop", json={"run_id": "nope"}).status_code == 404


def test_reconcile_marks_dead_pid(monkeypatch):
    sweep_web._save_registry(
        {"rid": {"run_id": "rid", "pid": 4242, "status": "running"}}
    )
    monkeypatch.setattr(sweep_web, "_pid_alive", lambda pid: False)
    reg = sweep_web._reconcile()
    assert reg["rid"]["status"] == "finished"
