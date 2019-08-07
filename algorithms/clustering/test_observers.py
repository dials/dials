from __future__ import absolute_import, division, print_function
import mock
import random

from cctbx import sgtbx
from scitbx import matrix
from dxtbx.model import Crystal, Experiment, ExperimentList
from dials.algorithms.clustering.unit_cell import UnitCellCluster
from dials.algorithms.clustering import observers


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
    lattice_ids = experiments.identifiers()
    ucs = UnitCellCluster.from_crystal_symmetries(
        crystal_symmetries, lattice_ids=lattice_ids
    )
    _, dendrogram, _ = ucs.ab_cluster(write_file_lists=False, doplot=False)

    # setup script
    script = mock.Mock()
    script._experiments = experiments
    script.unit_cell_dendrogram = dendrogram

    # test the observer
    observer = observers.UnitCellAnalysisObserver()
    observer.update(script)
    assert set(observer.data) == {"experiments", "dendrogram"}
    d = observer.make_plots()
    assert "unit_cell_graphs" in d
