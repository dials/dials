from __future__ import absolute_import, division, print_function

import os
import pytest
import procrunner

from cctbx import sgtbx, uctbx
import scitbx.matrix
from dxtbx.serialize import load
from dxtbx.model import Crystal, Experiment, ExperimentList
from dials.array_family import flex
from dials.algorithms.symmetry.cosym._generate_test_data import (
    generate_experiments_reflections,
)
from dials.command_line.cosym import cosym


@pytest.mark.parametrize("space_group", [None, "P 1", "P 4"])
def test_cosym(dials_data, tmpdir, space_group):
    mcp = dials_data("multi_crystal_proteinase_k")
    command = ["dials.cosym", "space_group=" + str(space_group)]
    for i in [1, 2, 3, 4, 5, 7, 8, 10]:
        command.append(mcp.join("experiments_%d.json" % i).strpath)
        command.append(mcp.join("reflections_%d.pickle" % i).strpath)

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check(file=1)
    assert tmpdir.join("symmetrized.expt").check(file=1)
    experiments = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    if space_group is None:
        assert (
            experiments[0].crystal.get_space_group().type().lookup_symbol() == "P 4 2 2"
        )
    else:
        assert (
            experiments[0].crystal.get_space_group().type().lookup_symbol()
            == space_group
        )


def test_cosym_partial_dataset(dials_data, tmpdir):
    """Test how cosym handles partial/bad datasets."""
    mcp = dials_data("multi_crystal_proteinase_k")
    command = ["dials.cosym"]
    for i in [1, 2]:
        command.append(mcp.join("experiments_%d.json" % i).strpath)
        command.append(mcp.join("reflections_%d.pickle" % i).strpath)
    # Make one dataset that will be removed in prefiltering
    r = flex.reflection_table.from_file(mcp.join("reflections_8.pickle").strpath)
    r["partiality"] = flex.double(r.size(), 0.1)
    r.as_file(tmpdir.join("renamed.refl").strpath)
    command.append(tmpdir.join("renamed.refl").strpath)
    command.append(mcp.join("experiments_8.json").strpath)
    # Add another good dataset at the end of the input list
    command.append(mcp.join("experiments_10.json").strpath)
    command.append(mcp.join("reflections_10.pickle").strpath)

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check(file=1)
    assert tmpdir.join("symmetrized.expt").check(file=1)
    experiments = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    assert len(experiments) == 3


def test_cosym_partial_dataset_raises_sorry(dials_data, tmpdir):
    """Test how cosym handles partial/bad datasets."""
    mcp = dials_data("multi_crystal_proteinase_k")
    command = ["dials.cosym"]
    command.append(tmpdir.join("renamed.refl").strpath)
    command.append(mcp.join("experiments_8.json").strpath)
    r2 = flex.reflection_table.from_file(mcp.join("reflections_10.pickle").strpath)
    r2["partiality"] = flex.double(r2.size(), 0.1)
    r2.as_file(tmpdir.join("renamed2.refl").strpath)
    command.append(tmpdir.join("renamed2.refl").strpath)
    command.append(mcp.join("experiments_10.json").strpath)

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    # Sorry exceptions are only raised as text at the system level
    assert result.stderr.startswith(b"Sorry")


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
def test_synthetic(
    space_group,
    unit_cell,
    dimensions,
    sample_size,
    use_known_space_group,
    use_known_lattice_group,
    tmpdir,
):
    os.chdir(tmpdir.strpath)
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
    expt_file = tmpdir.join("tmp.expt").strpath
    joint_table = flex.reflection_table()
    for r in reflections:
        joint_table.extend(r)
    joint_table.as_pickle("tmp.refl")
    refl_file = tmpdir.join("tmp.refl").strpath

    command = [
        "dials.cosym",
        expt_file,
        refl_file,
        "output.experiments=symmetrized.expt",
        "output.reflections=symmetrized.refl",
        "output.html=cosym.html",
        "output.json=cosym.json",
    ]

    if use_known_space_group:
        command.append("space_group=%s" % space_group.info())
    if use_known_lattice_group:
        command.append("lattice_group=%s" % space_group.info())
    if dimensions is not None:
        command.append("dimensions=%s" % dimensions)

    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check(file=1)
    assert tmpdir.join("symmetrized.expt").check(file=1)
    assert tmpdir.join("cosym.html").check(file=1)
    assert tmpdir.join("cosym.json").check(file=1)
    """experiments, reflections = assign_unique_identifiers(experiments, reflections)
    cosym_instance = cosym(experiments, reflections, params=params)
    register_default_cosym_observers(cosym_instance)
    cosym_instance.run()
    cosym_instance.export()



    assert os.path.exists(params.output.experiments)
    assert os.path.exists(params.output.reflections)
    assert os.path.exists(params.output.html)
    assert os.path.exists(params.output.json)"""
    cosym_expts = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    assert len(cosym_expts) == len(experiments)
    for expt in cosym_expts:
        if unit_cell is not None:
            assert expt.crystal.get_unit_cell().parameters() == pytest.approx(
                unit_cell.parameters()
            )
        assert str(expt.crystal.get_space_group().info()) == str(space_group.info())
        assert expt.crystal.get_space_group() == space_group


def test_eliminate_sys_absent():
    refl = flex.reflection_table()
    refl["miller_index"] = flex.miller_index(
        [(-31, -5, -3), (-25, -3, -3), (0, 1, 0), (-42, -8, -2)]
    )
    sgi = sgtbx.space_group_info("C121")
    uc = sgi.any_compatible_unit_cell(volume=1000)
    B = scitbx.matrix.sqr(uc.fractionalization_matrix()).transpose()
    expt = Experiment(crystal=Crystal(B, space_group=sgi.group(), reciprocal=True))
    reflections = cosym._eliminate_sys_absent([expt], [refl])
    assert list(reflections[0]["miller_index"]) == [
        (-31, -5, -3),
        (-25, -3, -3),
        (-42, -8, -2),
    ]


def test_map_to_minimum_cell():
    # Input and expected output
    input_ucs = [
        (39.7413, 183.767, 140.649, 90, 90, 90),
        (40.16, 142.899, 92.4167, 90, 102.48, 90),
        (180.613, 40.1558, 142.737, 90, 90.0174, 90),
    ]
    input_sgs = ["C 2 2 21", "P 1 2 1", "C 1 2 1"]
    input_hkl = [
        [(1, -75, -71), (1, -73, -70), (1, -71, -69)],
        [(14, -37, -36), (-2, -35, -46), (-3, -34, -47)],
        [(-31, -5, -3), (-25, -3, -3), (-42, -8, -2)],
    ]
    expected_ucs = [
        (39.7413, 94.00755450320204, 140.649, 90.0, 90.0, 77.79717980856927),
        (40.16, 92.46399390642911, 142.899, 90.0, 90.0, 77.3882749092846),
        (
            40.1558,
            92.51154528306184,
            142.73699999999997,
            89.9830147351441,
            90.0,
            77.46527404307477,
        ),
    ]
    expected_output_hkl = [
        [(-1, 37, -71), (-1, 36, -70), (-1, 35, -69)],
        [(-14, 22, 37), (2, 48, 35), (3, 50, 34)],
        [(-5, 13, -3), (-3, 11, -3), (-8, 17, -2)],
    ]

    # Setup the input experiments and reflection tables
    expts = ExperimentList()
    reflections = []
    for uc, sg, hkl in zip(input_ucs, input_sgs, input_hkl):
        uc = uctbx.unit_cell(uc)
        sg = sgtbx.space_group_info(sg).group()
        B = scitbx.matrix.sqr(uc.fractionalization_matrix()).transpose()
        expts.append(Experiment(crystal=Crystal(B, space_group=sg, reciprocal=True)))
        refl = flex.reflection_table()
        refl["miller_index"] = flex.miller_index(hkl)
        reflections.append(refl)

    # Actually run the method we are testing
    expts_min, reflections_min = cosym._map_to_minimum_cell(
        expts, reflections, max_delta=5
    )

    # Verify that the unit cells have been transformed as expected
    for expt, uc in zip(expts, expected_ucs):
        assert expt.crystal.get_unit_cell().parameters() == pytest.approx(uc, abs=4e-2)

    # Space group should be set to P1
    assert [expt.crystal.get_space_group().type().number() for expt in expts_min] == [
        1,
        1,
        1,
    ]

    # Verify that the reflections have been reindexed as expected
    # Because the exact choice of minimum cell can be platform-dependent,
    # compare the magnitude, but not the sign of the output hkl values
    for refl, expected_hkl in zip(reflections, expected_output_hkl):
        for hkl, e_hkl in zip(refl["miller_index"], expected_hkl):
            assert [abs(h) for h in hkl] == [abs(eh) for eh in e_hkl]
