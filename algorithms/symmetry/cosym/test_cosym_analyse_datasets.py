from __future__ import absolute_import, division, print_function

import libtbx
import pytest

from cctbx import sgtbx
from dials.algorithms.symmetry.cosym._generate_test_data import generate_test_data
from dials.algorithms.symmetry.cosym import phil_scope
from dials.algorithms.symmetry.cosym import CosymAnalysis


@pytest.mark.parametrize(
    (
        "space_group",
        "unit_cell",
        "dimensions",
        "sample_size",
        "use_known_space_group",
        "use_known_lattice_group",
        "best_monoclinic_beta",
    ),
    [
        ("P2", None, None, 10, False, False, True),
        ("P3", None, None, 20, False, False, True),
        ("I23", None, libtbx.Auto, 10, False, False, True),
        ("I23", None, libtbx.Auto, 10, True, False, True),
        ("P422", (79, 79, 37, 90, 90, 90), None, 20, True, True, True),
        ("P321", (59.39, 59.39, 28.35, 90, 90, 120), None, 5, False, False, True),
        (
            "C121",
            (112.67, 52.85, 44.47, 90.00, 102.97, 90.00),
            None,
            5,
            False,
            False,
            False,
        ),
        (
            "C121",
            (112.67, 52.85, 44.47, 90.00, 102.97, 90.00),
            None,
            5,
            True,
            True,
            False,
        ),
        (
            "I121",
            (44.47, 52.85, 111.46, 90.00, 99.91, 90.00),
            None,
            5,
            False,
            False,
            True,
        ),
        (
            "I121",
            (44.47, 52.85, 111.46, 90.00, 99.91, 90.00),
            None,
            5,
            True,
            True,
            True,
        ),
    ],
)
def test_cosym(
    space_group,
    unit_cell,
    dimensions,
    sample_size,
    use_known_space_group,
    use_known_lattice_group,
    best_monoclinic_beta,
    run_in_tmpdir,
):
    import matplotlib

    matplotlib.use("Agg")

    datasets, expected_reindexing_ops = generate_test_data(
        space_group=sgtbx.space_group_info(symbol=space_group).group(),
        unit_cell=unit_cell,
        unit_cell_volume=10000,
        d_min=1.5,
        map_to_p1=True,
        sample_size=sample_size,
        seed=1,
    )
    expected_space_group = sgtbx.space_group_info(symbol=space_group).group()

    params = phil_scope.extract()
    params.cluster.n_clusters = len(expected_reindexing_ops)
    params.dimensions = dimensions
    params.best_monoclinic_beta = best_monoclinic_beta
    if use_known_space_group:
        params.space_group = expected_space_group.info()
    if use_known_lattice_group:
        params.lattice_group = expected_space_group.info()

    params.normalisation = None
    cosym = CosymAnalysis(datasets, params)
    cosym.run()
    d = cosym.as_dict()
    if not use_known_space_group:
        assert d["subgroup_scores"][0]["likelihood"] > 0.89
        assert (
            sgtbx.space_group(d["subgroup_scores"][0]["patterson_group"])
            == sgtbx.space_group_info(space_group)
            .group()
            .build_derived_patterson_group()
        )

    space_groups = {}
    reindexing_ops = {}
    for dataset_id in cosym.reindexing_ops.keys():
        if 0 in cosym.reindexing_ops[dataset_id]:
            cb_op = cosym.reindexing_ops[dataset_id][0]
            reindexing_ops.setdefault(cb_op, set())
            reindexing_ops[cb_op].add(dataset_id)
        if dataset_id in cosym.space_groups:
            space_groups.setdefault(cosym.space_groups[dataset_id], set())
            space_groups[cosym.space_groups[dataset_id]].add(dataset_id)

    assert len(reindexing_ops) == len(expected_reindexing_ops)
    assert len(space_groups) == 1
    space_group_info = list(space_groups)[0].info()

    if use_known_space_group:
        expected_sg = sgtbx.space_group_info(space_group).group()
    else:
        expected_sg = (
            sgtbx.space_group_info(space_group).group().build_derived_patterson_group()
        )
    assert cosym.best_subgroup["best_subsym"].space_group() == expected_sg

    for cb_op, ridx_set in reindexing_ops.items():
        for expected_set in expected_reindexing_ops.values():
            assert (len(ridx_set.symmetric_difference(expected_set)) == 0) or (
                len(ridx_set.intersection(expected_set)) == 0
            )
        for d_id in ridx_set:
            reindexed = (
                datasets[d_id]
                .change_basis(sgtbx.change_of_basis_op(cb_op))
                .customized_copy(
                    space_group_info=space_group_info.change_basis(
                        cosym.cb_op_inp_min.inverse()
                    )
                )
            )
            assert reindexed.is_compatible_unit_cell(), str(
                reindexed.crystal_symmetry()
            )
