import procrunner
import pytest

from dials.array_family import flex


@pytest.fixture(scope="session")
def reflections(tmpdir_factory):
    # Make a dummy reflection table for the test setting some values and flags
    rt = flex.reflection_table.empty_standard(6)
    rt["iobs"] = flex.size_t_range(len(rt))
    rt["panel"] = flex.size_t_range(len(rt))
    rt["id"] = flex.int([0] * 5 + [1])
    rt["d"] = flex.double([50, 40, 3.0, 2.5, 2.0, 1.0])
    mask1 = flex.bool([True] * 3 + [False] * 3)
    mask2 = flex.bool([True, False] * 3)
    rt.set_flags(mask1, rt.flags.integrated)
    rt.set_flags(mask2, rt.flags.reference_spot)
    tmpdir = tmpdir_factory.mktemp("filter_reflections")
    rt_name = tmpdir.join("test_refs.refl")
    rt.as_file(rt_name.strpath)
    return rt_name


def test_filter_reflections_flag_expression(reflections, tmpdir):
    result = procrunner.run(
        [
            "dials.filter_reflections",
            reflections,
            "flag_expression='integrated & ~reference_spot'",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_file(tmpdir.join("filtered.refl").strpath)
    # The test selects only the 2nd reflection
    assert len(ref) == 1
    assert list(ref["iobs"]) == [1]


def test_filter_reflections_by_experiment_id(reflections, tmpdir):
    result = procrunner.run(
        ["dials.filter_reflections", reflections, "id=0"], working_directory=tmpdir
    )
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_file(tmpdir.join("filtered.refl").strpath)
    # The test selects only the first five reflections
    assert len(ref) == 5
    assert list(ref["iobs"]) == [0, 1, 2, 3, 4]


def test_filter_reflections_by_panel(reflections, tmpdir):
    result = procrunner.run(
        ["dials.filter_reflections", reflections, "panel=5"], working_directory=tmpdir
    )
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_file(tmpdir.join("filtered.refl").strpath)
    # The test selects only the last reflection
    assert len(ref) == 1
    assert list(ref["iobs"]) == [5]


def test_filter_reflections_by_resolution(reflections, tmpdir):
    result = procrunner.run(
        ["dials.filter_reflections", reflections, "d_max=3.0", "d_min=2.0"],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    ref = flex.reflection_table.from_file(tmpdir.join("filtered.refl").strpath)
    # The test selects only the 3rd, 4th and 5th reflections
    assert len(ref) == 3
    assert list(ref["iobs"]) == [2, 3, 4]


def test_filter_reflections_printing_analysis(reflections, tmpdir):
    result = procrunner.run(
        ["dials.filter_reflections", reflections], working_directory=tmpdir
    )
    assert not result.returncode and not result.stderr


@pytest.mark.parametrize(
    "dirname,expt_filename,refl_filename,args,expected",
    [
        (
            "centroid_test_data",
            "imported_experiments.json",
            "strong.pickle",
            ["ice_rings.filter=True"],
            615,
        ),
        (
            "centroid_test_data",
            "imported_experiments.json",
            "strong.pickle",
            ["d_min=3.1", "d_max=20"],
            153,
        ),
        (
            "centroid_test_data",
            "experiments.json",
            "integrated.refl",
            ["d_min=2", "d_max=20"],
            371,
        ),
    ],
)
def test_filter_reflections(
    dirname, expt_filename, refl_filename, args, expected, dials_data, tmpdir
):
    directory = dials_data(dirname)
    result = procrunner.run(
        [
            "dials.filter_reflections",
            directory.join(expt_filename).strpath,
            directory.join(refl_filename).strpath,
        ]
        + args,
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    filtered_refl = flex.reflection_table.from_file(tmpdir.join("filtered.refl"))
    assert len(filtered_refl) == expected
