from __future__ import absolute_import, division, print_function

import pytest

import libtbx
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
    ),
    [
        ("P2", None, None, 10, False, False),
        ("P3", None, None, 20, False, False),
        ("I23", None, libtbx.Auto, 10, False, False),
        ("I23", None, libtbx.Auto, 10, True, False),
        ("P422", (79, 79, 37, 90, 90, 90), None, 20, True, True),
    ],
)
def test_cosym(
    space_group,
    unit_cell,
    dimensions,
    sample_size,
    use_known_space_group,
    use_known_lattice_group,
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
    )
    expected_space_group = sgtbx.space_group_info(symbol=space_group).group()

    # Workaround fact that the minimum cell reduction can occassionally be unstable
    # The input *should* be already the minimum cell, but for some combinations of unit
    # cell parameters the change_of_basis_op_to_minimum_cell is never the identity.
    # Therefore apply this cb_op to the expected_reindexing_ops prior to the comparison.
    cb_op_inp_min = datasets[0].crystal_symmetry().change_of_basis_op_to_minimum_cell()
    expected_reindexing_ops = {
        (sgtbx.change_of_basis_op(cb_op) * cb_op_inp_min).as_xyz(): dataset_ids
        for cb_op, dataset_ids in expected_reindexing_ops.items()
    }

    params = phil_scope.extract()
    params.cluster.n_clusters = len(expected_reindexing_ops)
    params.dimensions = dimensions
    if use_known_space_group:
        params.space_group = expected_space_group.info()
    if use_known_lattice_group:
        params.lattice_group = expected_space_group.info()

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
    assert sorted(reindexing_ops.keys()) == sorted(expected_reindexing_ops.keys())

    for ridx_set in reindexing_ops.values():
        for expected_set in expected_reindexing_ops.values():
            assert (len(ridx_set.symmetric_difference(expected_set)) == 0) or (
                len(ridx_set.intersection(expected_set)) == 0
            )
