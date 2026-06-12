"""Tests for the multi-machine sweep CLI (src/runners/sweep_cli.py)."""

import json

import pytest

from runners import sweep_cli


def _write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


# ---------------------------------------------------------------------------
# Manifest loading
# ---------------------------------------------------------------------------


def test_load_manifest_reads_cases(tmp_path):
    _write_json(
        tmp_path / "demo.json",
        {"description": "d", "cases": [{"alias": "a", "config": "c.json"}]},
    )
    manifest = sweep_cli.load_manifest("demo", manifests_dir=tmp_path)
    assert manifest["name"] == "demo"
    assert manifest["cases"][0]["alias"] == "a"


def test_load_manifest_missing_file_exits(tmp_path):
    with pytest.raises(SystemExit, match="not found"):
        sweep_cli.load_manifest("nope", manifests_dir=tmp_path)


def test_load_manifest_rejects_empty_cases(tmp_path):
    _write_json(tmp_path / "demo.json", {"cases": []})
    with pytest.raises(SystemExit, match="no 'cases'"):
        sweep_cli.load_manifest("demo", manifests_dir=tmp_path)


def test_load_manifest_rejects_duplicate_alias(tmp_path):
    _write_json(
        tmp_path / "demo.json",
        {"cases": [{"alias": "x", "config": "a.json"}, {"alias": "x", "config": "b.json"}]},
    )
    with pytest.raises(SystemExit, match="duplicate alias"):
        sweep_cli.load_manifest("demo", manifests_dir=tmp_path)


# ---------------------------------------------------------------------------
# Case-name resolution
# ---------------------------------------------------------------------------


def test_case_name_prefers_result_tag(tmp_path):
    cfg = tmp_path / "case.json"
    _write_json(cfg, {"outputs": {"result_tag": "tag_name"}, "name": "other"})
    assert sweep_cli.case_name_for(cfg) == "tag_name"


def test_case_name_falls_back_to_name_then_stem(tmp_path):
    cfg = tmp_path / "case.json"
    _write_json(cfg, {"name": "config_name"})
    assert sweep_cli.case_name_for(cfg) == "config_name"
    _write_json(cfg, {})
    assert sweep_cli.case_name_for(cfg) == "case"


# ---------------------------------------------------------------------------
# Selector resolution
# ---------------------------------------------------------------------------


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """Fake repo root with two case configs."""
    monkeypatch.setattr(sweep_cli, "REPO_ROOT", tmp_path)
    _write_json(tmp_path / "configs" / "one.json", {"outputs": {"result_tag": "case_one"}})
    _write_json(tmp_path / "configs" / "two.json", {"outputs": {"result_tag": "case_two"}})
    return tmp_path


@pytest.fixture
def manifest(repo):
    return {
        "name": "demo",
        "cases": [
            {"alias": "c1", "config": "configs/one.json"},
            {"alias": "c2", "config": "configs/two.json"},
        ],
    }


def test_resolve_all(manifest):
    cases = sweep_cli.resolve_cases(manifest, ["all"])
    assert [c["case_name"] for c in cases] == ["case_one", "case_two"]


def test_resolve_by_alias_and_full_name(manifest):
    cases = sweep_cli.resolve_cases(manifest, ["c2", "case_one"])
    assert [c["case_name"] for c in cases] == ["case_two", "case_one"]


def test_resolve_deduplicates(manifest):
    cases = sweep_cli.resolve_cases(manifest, ["c1", "case_one"])
    assert len(cases) == 1


def test_resolve_unknown_selector_exits(manifest):
    with pytest.raises(SystemExit, match="unknown case"):
        sweep_cli.resolve_cases(manifest, ["zz"])


# ---------------------------------------------------------------------------
# Status classification
# ---------------------------------------------------------------------------


def test_status_pending_when_nothing_harvested(tmp_path):
    st = sweep_cli.case_status("s", "c", shared_root=tmp_path)
    assert st["status"] == "pending"


def test_status_done_with_machine_from_harvest(tmp_path):
    case_dir = tmp_path / "s" / "c"
    _write_json(case_dir / "run_status.json", {"status": "completed"})
    _write_json(
        case_dir / "harvest.json",
        {"machine": "pc-a", "harvested_at": "2026-06-10T01:00:00"},
    )
    st = sweep_cli.case_status("s", "c", shared_root=tmp_path)
    assert st == {"status": "done", "machine": "pc-a", "when": "2026-06-10T01:00:00"}


def test_status_incomplete_run_status_is_not_done(tmp_path):
    case_dir = tmp_path / "s" / "c"
    _write_json(case_dir / "run_status.json", {"status": "running"})
    assert sweep_cli.case_status("s", "c", shared_root=tmp_path)["status"] == "pending"


def test_status_failed_beats_in_progress(tmp_path):
    case_dir = tmp_path / "s" / "c"
    _write_json(case_dir / "FAILED.json", {"machine": "pc-b", "failed_at": "t"})
    _write_json(case_dir / "STARTED.json", {"machine": "pc-b", "started_at": "t"})
    assert sweep_cli.case_status("s", "c", shared_root=tmp_path)["status"] == "failed"


def test_status_in_progress(tmp_path):
    case_dir = tmp_path / "s" / "c"
    _write_json(case_dir / "STARTED.json", {"machine": "pc-b", "started_at": "t"})
    st = sweep_cli.case_status("s", "c", shared_root=tmp_path)
    assert st["status"] == "in-progress"
    assert st["machine"] == "pc-b"


# ---------------------------------------------------------------------------
# Harvest run-dir selection
# ---------------------------------------------------------------------------


def _make_run(results_root, case, run_name, status="completed", dim3=True):
    base = results_root / "dim3" if dim3 else results_root
    run_dir = base / case / run_name
    _write_json(run_dir / "run_status.json", {"status": status})
    return run_dir


def test_find_run_dir_picks_newest_completed(tmp_path):
    _make_run(tmp_path, "c", "run_20260601_000000_000000")
    newest = _make_run(tmp_path, "c", "run_20260609_000000_000000")
    _make_run(tmp_path, "c", "run_20260610_000000_000000", status="running")
    assert sweep_cli.find_completed_run_dir("c", results_root=tmp_path) == newest


def test_find_run_dir_respects_since(tmp_path):
    import os

    old = _make_run(tmp_path, "c", "run_20260601_000000_000000")
    os.utime(old / "run_status.json", (1000, 1000))
    assert sweep_cli.find_completed_run_dir("c", since=2000, results_root=tmp_path) is None


def test_find_run_dir_searches_non_dim3_layout(tmp_path):
    run_dir = _make_run(tmp_path, "c", "run_20260601_000000_000000", dim3=False)
    assert sweep_cli.find_completed_run_dir("c", results_root=tmp_path) == run_dir


def test_find_run_dir_none_when_no_results(tmp_path):
    assert sweep_cli.find_completed_run_dir("c", results_root=tmp_path) is None


# ---------------------------------------------------------------------------
# Harvest copy
# ---------------------------------------------------------------------------


def test_harvest_copies_light_files_and_clears_markers(tmp_path, monkeypatch):
    monkeypatch.setattr(sweep_cli, "REPO_ROOT", tmp_path)
    run_dir = tmp_path / "results" / "dim3" / "c" / "run_x"
    _write_json(run_dir / "run_status.json", {"status": "completed"})
    _write_json(run_dir / "summary.json", {"metric": 1})
    _write_json(run_dir / "analysis" / "time_series.csv", {})
    (run_dir / "paraview").mkdir()
    (run_dir / "paraview" / "fluid_000001.vtk").write_text("huge")

    shared = tmp_path / "sweep_results"
    case_dir = shared / "s" / "c"
    _write_json(case_dir / "STARTED.json", {"machine": "m"})

    sweep_cli.harvest_case("s", "c", run_dir, shared_root=shared)

    assert (case_dir / "run_status.json").is_file()
    assert (case_dir / "summary.json").is_file()
    assert (case_dir / "analysis" / "time_series.csv").is_file()
    assert not (case_dir / "paraview").exists()
    assert not (case_dir / "STARTED.json").exists()
    harvest = json.loads((case_dir / "harvest.json").read_text())
    assert harvest["source_run_dir"] == "results/dim3/c/run_x"
