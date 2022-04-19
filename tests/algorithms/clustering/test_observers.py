from __future__ import annotations

import random
from unittest import mock

from cctbx import sgtbx
from dxtbx.model import Crystal, Experiment, ExperimentList
from scitbx import matrix

from dials.algorithms.clustering import observers
from dials.algorithms.clustering.unit_cell import cluster_unit_cells


def test_UnitCellAnalysisObserver():
    # generate some random unit cells
    sgi = sgtbx.space_group_info("P1")
    unit_cells = [
        sgi.any_compatible_unit_cell(volume=random.uniform(990, 1010))
        for i in range(10)
    ]

    # generate experiment list
    experiments = ExperimentList()
    U = matrix.identity(3)
    for uc in unit_cells:
        B = matrix.sqr(uc.fractionalization_matrix()).transpose()
        direct_matrix = (U * B).inverse()
        experiments.append(
            Experiment(
                crystal=Crystal(
                    direct_matrix[:3],
                    direct_matrix[3:6],
                    direct_matrix[6:9],
                    space_group=sgi.group(),
                )
            )
        )

    # generate dendrogram
    crystal_symmetries = [expt.crystal.get_crystal_symmetry() for expt in experiments]
    result = cluster_unit_cells(
        crystal_symmetries, lattice_ids=list(experiments.identifiers())
    )

    # setup script
    script = mock.Mock()
    script._experiments = experiments
    script.unit_cell_dendrogram = result.dendrogram

    # test the observer
    observer = observers.UnitCellAnalysisObserver()
    observer.update(script)
    assert set(observer.data) == {"experiments", "dendrogram"}
    d = observer.make_plots()
    assert "unit_cell_graphs" in d
