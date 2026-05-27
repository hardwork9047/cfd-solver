"""Tests for pore geometry definitions."""

import numpy as np

from particulate_flow import PoreGeometry


def test_pore_geometry_builds_periodic_cylinder_mask_across_y_boundary():
    """Cylinders near y=0 should wrap their solid mask through the top side."""
    geometry = PoreGeometry.from_cylinders([(10.0, 0.5, 2.0)])

    mask = geometry.cylinder_solid_mask(24, 16, y_boundary="periodic")

    assert mask[10, 0]
    assert mask[10, -1]


def test_pore_geometry_pressure_sections_bracket_cylinder_group():
    """Pressure probe sections should surround the fixed pore geometry."""
    geometry = PoreGeometry.from_cylinders([(30.0, 10.0, 4.0), (60.0, 10.0, 5.0)])

    sections = geometry.pressure_sections(100)

    assert sections["inlet_x"] == 1
    assert sections["outlet_x"] == 98
    assert sections["pore_upstream_x"] < 30
    assert sections["pore_downstream_x"] > 60


def test_pore_geometry_water_area_removes_cylinder_area():
    """Water-area estimates should subtract fixed cylinder areas."""
    geometry = PoreGeometry.from_cylinders([(10.0, 10.0, 2.0)])

    area = geometry.water_area(20, 10, wall_y=False)

    assert np.isclose(area, 200.0 - np.pi * 4.0)
