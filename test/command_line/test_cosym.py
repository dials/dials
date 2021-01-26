import os

import pytest

from cctbx import sgtbx, uctbx
from dxtbx.serialize import load

import dials.command_line.cosym as dials_cosym
from dials.algorithms.symmetry.cosym._generate_test_data import (
    generate_experiments_reflections,
)
from dials.array_family import flex
from dials.util import Sorry


@pytest.mark.parametrize("space_group", [None, "P 1", "P 4"])
def test_cosym(dials_data, run_in_tmpdir, space_group):
    mcp = dials_data("multi_crystal_proteinase_k")
    args = ["space_group=" + str(space_group)]
    for i in [1, 2, 3, 4, 5, 7, 8, 10]:
        args.append(mcp.join("experiments_%d.json" % i).strpath)
        args.append(mcp.join("reflections_%d.pickle" % i).strpath)
    dials_cosym.run(args=args)
    assert os.path.isfile("symmetrized.refl")
    assert os.path.isfile("symmetrized.expt")
    experiments = load.experiment_list("symmetrized.expt", check_format=False)
    if space_group is None:
        assert (
            experiments[0].crystal.get_space_group().type().lookup_symbol() == "P 4 2 2"
        )
    else:
        assert (
            experiments[0].crystal.get_space_group().type().lookup_symbol()
            == space_group
        )


def test_cosym_partial_dataset(dials_data, run_in_tmpdir):
    """Test how cosym handles partial/bad datasets."""
    mcp = dials_data("multi_crystal_proteinase_k")
    args = []
    for i in [1, 2]:
        args.append(mcp.join("experiments_%d.json" % i).strpath)
        args.append(mcp.join("reflections_%d.pickle" % i).strpath)
    # Make one dataset that will be removed in prefiltering
    r = flex.reflection_table.from_file(mcp.join("reflections_8.pickle").strpath)
    r["partiality"] = flex.double(r.size(), 0.1)
    r.as_file("renamed.refl")
    args.append("renamed.refl")
    args.append(mcp.join("experiments_8.json").strpath)
    # Add another good dataset at the end of the input list
    args.append(mcp.join("experiments_10.json").strpath)
    args.append(mcp.join("reflections_10.pickle").strpath)

    dials_cosym.run(args=args)
    assert os.path.exists("symmetrized.refl")
    assert os.path.exists("symmetrized.expt")
    experiments = load.experiment_list("symmetrized.expt", check_format=False)
    assert len(experiments) == 3


def test_cosym_partial_dataset_raises_sorry(dials_data, run_in_tmpdir, capsys):
    """Test how cosym handles partial/bad datasets."""
    mcp = dials_data("multi_crystal_proteinase_k")
    args = ["renamed.refl", mcp.join("experiments_8.json").strpath]
    r2 = flex.reflection_table.from_file(mcp.join("reflections_10.pickle").strpath)
    r2["partiality"] = flex.double(r2.size(), 0.1)
    r2.as_file("renamed2.refl")
    args.append("renamed2.refl")
    args.append(mcp.join("experiments_10.json").strpath)

    with pytest.raises(Sorry):
        dials_cosym.run(args=args)


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
        ("I23", None, 2, 10, False, False),
        ("P422", (79, 79, 37, 90, 90, 90), None, 10, True, False),
        ("P321", (59.39, 59.39, 28.35, 90, 90, 120), None, 5, False, False),
    ],
)
@pytest.mark.xfail("os.name == 'nt'", reason="UnicodeEncodeError in logging")
def test_synthetic(
    space_group,
    unit_cell,
    dimensions,
    sample_size,
    use_known_space_group,
    use_known_lattice_group,
    run_in_tmpdir,
):
    space_group = sgtbx.space_group_info(space_group).group()
    if unit_cell is not None:
        unit_cell = uctbx.unit_cell(unit_cell)
    experiments, reflections, _ = generate_experiments_reflections(
        space_group=space_group,
        unit_cell=unit_cell,
        unit_cell_volume=10000,
        sample_size=sample_size,
        map_to_p1=True,
        d_min=1.5,
    )

    experiments.as_json("tmp.expt")
    expt_file = "tmp.expt"
    joint_table = flex.reflection_table()
    for r in reflections:
        joint_table.extend(r)
    joint_table.as_file("tmp.refl")
    refl_file = "tmp.refl"

    args = [
        expt_file,
        refl_file,
        "output.experiments=symmetrized.expt",
        "output.reflections=symmetrized.refl",
        "output.html=cosym.html",
        "output.json=cosym.json",
    ]

    if use_known_space_group:
        args.append("space_group=%s" % space_group.info())
    if use_known_lattice_group:
        args.append("lattice_group=%s" % space_group.info())
    if dimensions is not None:
        args.append("dimensions=%s" % dimensions)

    dials_cosym.run(args=args)
    assert os.path.isfile("symmetrized.refl")
    assert os.path.isfile("symmetrized.expt")
    assert os.path.isfile("cosym.html")
    assert os.path.isfile("cosym.json")
    cosym_expts = load.experiment_list("symmetrized.expt", check_format=False)
    assert len(cosym_expts) == len(experiments)
    for expt in cosym_expts:
        if unit_cell is not None:
            assert expt.crystal.get_unit_cell().parameters() == pytest.approx(
                unit_cell.parameters()
            )
        if (
            str(expt.crystal.get_space_group().info()) == "P 6 2 2"
            and str(space_group.info()) == "P 3 2 1"
        ):
            # This is fine
            continue
        assert str(expt.crystal.get_space_group().info()) == str(space_group.info())
        assert expt.crystal.get_space_group() == space_group


def test_reindexing_identity(mocker):
    """
    Default to choosing the cluster that contains the most identity reindexing ops.

    This can be important if the lattice symmetry is only very approximately related by
    pseudosymmetry to the true symmetry. If all potential reindexing ops are genuine
    indexing ambiguities, then it doesn't matter which one is chosen, however if not,
    then choosing the wrong one will distort the true unit cell. In such cases it is
    likely that the input datasets were already indexed consistently, therefore default
    to choosing the cluster that contains the most identity reindexing ops.
    """
    # patch the cosym object, including the __init__
    mocker.patch.object(dials_cosym.cosym, "__init__", return_value=None, autospec=True)
    cosym = dials_cosym.cosym(mocker.ANY, mocker.ANY)
    cosym.observers = mocker.MagicMock()
    cosym.cosym_analysis = mocker.Mock()
    cosym._apply_reindexing_operators = mocker.Mock()
    cosym.params = dials_cosym.phil_scope.extract()

    # Mock cosym_analysis reindexing ops
    cosym.params.cluster.n_clusters = 6
    cosym.cosym_analysis.reindexing_ops = {
        0: {
            0: "-x,x+y,-z",
            1: "x,-x-y,-z",
            2: "-x,-y,z",
            3: "x,y,z",
            4: "-x-y,y,-z",
            5: "x+y,-y,-z",
        },
        1: {
            0: "-x,x+y,-z",
            1: "x,-x-y,-z",
            2: "-x,-y,z",
            3: "x,y,z",
            4: "-x-y,y,-z",
            5: "x+y,-y,-z",
        },
        2: {
            0: "-x,x+y,-z",
            1: "x,-x-y,-z",
            2: "-x,-y,z",
            3: "x,y,z",
            4: "-x-y,y,-z",
            5: "x+y,-y,-z",
        },
        3: {
            0: "-x,x+y,-z",
            1: "x,-x-y,-z",
            2: "-x,-y,z",
            3: "x,y,z",
            4: "-x-y,y,-z",
            5: "x+y,-y,-z",
        },
        4: {
            0: "-x,x+y,-z",
            1: "x,-x-y,-z",
            2: "-x,-y,z",
            3: "x,y,z",
            4: "-x-y,y,-z",
            5: "x+y,-y,-z",
        },
    }
    cosym.run()

    # Assert that chosen reindexing ops were the identity ops
    cosym._apply_reindexing_operators.assert_called_once_with(
        {"x,y,z": [0, 1, 2, 3, 4]},
        subgroup=mocker.ANY,
    )
