from __future__ import absolute_import, division, print_function

import os

import pytest


def test_dials_cluster_unit_cell_command_line(dials_regression, run_in_tmpdir):
    pytest.importorskip("scipy")
    pytest.importorskip("xfel")

    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "multi_narrow_wedges"
    )
    import glob

    experiments = glob.glob(os.path.join(data_dir, "data/sweep_*/experiments.json"))

    import procrunner

    result = procrunner.run(
        command=["dials.cluster_unit_cell", "plot.show=False"] + experiments,
        print_stdout=False,
    )
    assert not result.returncode

    # print result
    assert os.path.exists("cluster_unit_cell.png")

    from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory
    from dials.command_line import cluster_unit_cell

    experiments = ExperimentList(
        [
            ExperimentListFactory.from_json_file(expt, check_format=False)[0]
            for expt in experiments
        ]
    )
    from cctbx import crystal

    crystal_symmetries = [
        crystal.symmetry(
            unit_cell=expt.crystal.get_unit_cell(),
            space_group=expt.crystal.get_space_group(),
        )
        for expt in experiments
    ]

    params = cluster_unit_cell.phil_scope.extract()
    params.plot.show = False
    params.plot.name = None
    clusters = cluster_unit_cell.do_cluster_analysis(crystal_symmetries, params)
    assert len(clusters) == 1
    cluster = clusters[0]
    assert len(cluster.members) == 40
    assert cluster.medians == pytest.approx(
        [
            90.9430182020995,
            90.9430182020995,
            90.9430182020995,
            109.47122063449069,
            109.47122063449069,
            109.47122063449069,
        ],
        abs=1e-6,
    )
    assert cluster.stdevs == pytest.approx(
        [0.09509739126548639, 0.09509739126548526, 0.0950973912654865, 0, 0, 0],
        abs=1e-6,
    )
