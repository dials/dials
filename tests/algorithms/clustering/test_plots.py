from __future__ import annotations

import random

from cctbx import sgtbx
from scitbx.array_family import flex

from dials.algorithms.clustering import plots
from dials.algorithms.clustering.unit_cell import cluster_unit_cells


def test_plot_uc_histograms():
    sgi = sgtbx.space_group_info("P1")
    unit_cells = [
        sgi.any_compatible_unit_cell(volume=random.uniform(990, 1010)).parameters()
        for i in range(10)
    ]
    uc_params = [flex.double(a) for a in zip(*unit_cells)]
    d = plots.plot_uc_histograms(uc_params)
    assert set(d) == {"uc_scatter", "uc_hist"}
    assert len(d["uc_scatter"]["data"]) == 3
    assert len(d["uc_hist"]["data"]) == 3


def test_scipy_dendrogram_to_plotly_json():
    # generate some random unit cells
    sgi = sgtbx.space_group_info("P1")
    crystal_symmetries = [
        sgi.any_compatible_crystal_symmetry(volume=random.uniform(990, 1010))
        for i in range(10)
    ]
    lattice_ids = list(range(len(crystal_symmetries)))
    result = cluster_unit_cells(crystal_symmetries, lattice_ids)
    d = plots.scipy_dendrogram_to_plotly_json(
        result.dendrogram, title="Unit cell clustering"
    )
    assert set(d) == {"layout", "data"}
