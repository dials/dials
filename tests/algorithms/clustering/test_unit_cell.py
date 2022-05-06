from __future__ import annotations

import random

from cctbx import sgtbx

from dials.algorithms.clustering.unit_cell import cluster_unit_cells


def test_unit_cell():
    # generate some random unit cells
    sgi = sgtbx.space_group_info("P1")
    crystal_symmetries = [
        sgi.any_compatible_crystal_symmetry(volume=random.uniform(990, 1010))
        for i in range(10)
    ]
    lattice_ids = list(range(len(crystal_symmetries)))
    result = cluster_unit_cells(crystal_symmetries, lattice_ids)
    assert "cluster_1" in str(result)
    assert len(result.clusters) == 1
    assert "dcoord" in result.dendrogram.keys()
