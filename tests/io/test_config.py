"""Tests for file-backed simulation run configuration."""

import json

from particulate_flow.io.config import SimulationConfig


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


def test_simulation_config_accepts_structured_geometry_block(tmp_path):
    """Research configs can define pore geometry without CLI-shaped keys."""
    path = tmp_path / "case.json"
    path.write_text(
        json.dumps(
            {
                "geometry": {
                    "cylinders": [
                        {"x": 34, "y": 13, "radius": 4.0},
                        {"x": 54, "y": 25, "radius": 4.5},
                    ]
                }
            }
        ),
        encoding="utf-8",
    )

    config = SimulationConfig.from_json(path)

    assert config.argparse_defaults()["cylinder_spec"] == [
        [34, 13, 4.0],
        [54, 25, 4.5],
    ]


def test_simulation_config_flattens_research_sections():
    """Configs can group settings by purpose while still feeding argparse."""
    config = SimulationConfig.from_mapping(
        {
            "case_type": "fouling_supply",
            "domain": {"nx": 100, "ny": 50},
            "flow": {"streamwise-boundary": "pressure"},
            "particles": {"n_particles": 0, "particle_radius": 1.5},
            "outputs": {"no_video": True},
        }
    )

    defaults = config.argparse_defaults()

    assert defaults["nx"] == 100
    assert defaults["ny"] == 50
    assert defaults["streamwise_boundary"] == "pressure"
    assert defaults["n_particles"] == 0
    assert defaults["particle_radius"] == 1.5
    assert defaults["no_video"] is True
    assert "case_type" not in defaults


def test_simulation_config_extends_parent_configs(tmp_path):
    """Case files can extend reusable templates."""
    parent = tmp_path / "template.json"
    parent.write_text(
        json.dumps(
            {
                "domain": {"nx": 80, "ny": 40},
                "outputs": {"no_video": True, "paraview_output": False},
            }
        ),
        encoding="utf-8",
    )
    child = tmp_path / "case.json"
    child.write_text(
        json.dumps(
            {
                "extends": "template.json",
                "domain": {"nx": 120},
                "outputs": {"paraview_output": True},
            }
        ),
        encoding="utf-8",
    )

    defaults = SimulationConfig.from_json(child).argparse_defaults()

    assert defaults["nx"] == 120
    assert defaults["ny"] == 40
    assert defaults["no_video"] is True
    assert defaults["paraview_output"] is True
