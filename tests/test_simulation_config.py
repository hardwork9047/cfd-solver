"""Tests for file-backed simulation run configuration."""

import json

from cfd_dem_lbm.simulation_config import SimulationConfig


def test_simulation_config_normalises_cli_style_keys(tmp_path):
    """Config files may use CLI-style dashed keys or argparse-style keys."""
    path = tmp_path / "case.json"
    path.write_text(
        json.dumps(
            {
                "streamwise-boundary": "pressure",
                "y_boundary": "periodic",
                "cylinder-spec": [[10, 5, 2]],
            }
        ),
        encoding="utf-8",
    )

    config = SimulationConfig.from_json(path)

    assert config.argparse_defaults()["streamwise_boundary"] == "pressure"
    assert config.argparse_defaults()["y_boundary"] == "periodic"
    assert config.argparse_defaults()["cylinder_spec"] == [[10, 5, 2]]


def test_simulation_config_writes_effective_arguments(tmp_path):
    """Each run can persist source config and final effective arguments."""
    config = SimulationConfig.from_mapping({"nx": 40}, source_path="case.json")
    out = tmp_path / "config.json"

    config.write_effective(out, {"nx": 80, "ny": 20})
    payload = json.loads(out.read_text(encoding="utf-8"))

    assert payload["schema_version"] == 1
    assert payload["source"]["nx"] == 40
    assert payload["source"]["_source_path"] == "case.json"
    assert payload["effective_arguments"]["nx"] == 80
    assert payload["effective_arguments"]["ny"] == 20
