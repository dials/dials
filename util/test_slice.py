import copy
import glob
import os
import pytest

from dxtbx.model import Experiment, ExperimentList
from dxtbx.model import Scan
from dxtbx.model.experiment_list import ExperimentListFactory

from dials.array_family import flex
from dials.util.slice import slice_experiments, slice_reflections


def test_slice_experiments():
    image_range = (0, 1000)
    oscillation = (0, 0.1)
    scan = Scan(image_range, oscillation)
    experiments = ExperimentList([Experiment(scan=scan)])
    sliced_image_range = [(1, 5)]
    sliced_experiments = slice_experiments(experiments, sliced_image_range)
    assert sliced_experiments[0].scan.get_image_range() == sliced_image_range[0]
    copy.deepcopy(sliced_experiments)


@pytest.mark.xfail
def test_slice_experiments_centroid_test_data(dials_regression):
    filenames = sorted(
        glob.glob(os.path.join(dials_regression, "centroid_test_data", "*.cbf"))
    )
    experiments = ExperimentListFactory.from_filenames(filenames)
    sliced_image_range = [(1, 3)]
    sliced_experiments = slice_experiments(experiments, sliced_image_range)
    assert sliced_experiments[0].scan.get_image_range() == sliced_image_range[0]
    # for some reason the sliced_experiments is not copyable
    copy.deepcopy(sliced_experiments)


def test_slice_reflections():
    r = flex.reflection_table()
    r["id"] = flex.int([0, 0, 0, 1, 1, 1, 2, 2, 2])
    image_number = [0, 1, 2, 0, 1, 2, 0, 1, 2]
    print((0 * len(r), 0 * len(r), image_number))
    print(zip([0] * len(r), [0] * len(r), image_number))

    r["xyzobs.px.value"] = flex.vec3_double(
        zip([0] * len(r), [0] * len(r), image_number)
    )
    sliced_r = slice_reflections(r, [(1, 2), (1, 1), (2, 3)])
    print(list(sliced_r["id"]))
    print(list(sliced_r["xyzobs.px.value"].parts()[2]))
