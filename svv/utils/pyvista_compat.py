"""Compatibility helpers for PyVista API changes."""

from __future__ import annotations


def compute_cell_quality(mesh, quality_measure="scaled_jacobian", **kwargs):
    """Return a mesh with a legacy-compatible ``CellQuality`` cell array."""
    if hasattr(mesh, "compute_cell_quality"):
        quality = mesh.compute_cell_quality(quality_measure=quality_measure, **kwargs)
    else:
        quality = mesh.cell_quality(quality_measure=quality_measure, **kwargs)

    if "CellQuality" not in quality.cell_data and quality_measure in quality.cell_data:
        quality.cell_data["CellQuality"] = quality.cell_data[quality_measure]
    return quality
