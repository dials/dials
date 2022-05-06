from __future__ import annotations

import copy

from dxtbx.model import Experiment, ExperimentList, Scan
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


def test_slice_experiments_centroid_test_data(dials_data):
    files = sorted(dials_data("centroid_test_data", pathlib=True).glob("*.cbf"))
    experiments = ExperimentListFactory.from_filenames(files)
    sliced_image_range = [(1, 3)]
    sliced_experiments = slice_experiments(experiments, sliced_image_range)
    assert sliced_experiments[0].scan.get_image_range() == sliced_image_range[0]
    # for some reason the sliced_experiments is not copyable
    assert copy.deepcopy(sliced_experiments)


def test_slice_experiments_centroid_test_data_starting_from_2(dials_data):
    files = sorted(dials_data("centroid_test_data", pathlib=True).glob("*.cbf"))[1:]
    experiments = ExperimentListFactory.from_filenames(files)
    sliced_image_range = [(2, 4)]
    sliced_experiments = slice_experiments(experiments, sliced_image_range)
    assert sliced_experiments[0].scan.get_image_range() == sliced_image_range[0]


def test_slice_reflections():
    r = flex.reflection_table()
    r["id"] = flex.int([0, 0, 0, 1, 1, 1, 2, 2, 2])
    image_number = [0, 1, 2, 0, 1, 2, 0, 1, 2]
    r["xyzobs.px.value"] = flex.vec3_double(
        zip([0] * len(r), [0] * len(r), image_number)
    )
    sliced_r = slice_reflections(r, [(1, 2), (1, 1), (2, 3)])
    assert list(sliced_r["id"]) == [0, 0, 1, 2, 2]
    assert list(sliced_r["xyzobs.px.value"].parts()[2]) == [0.0, 1.0, 0.0, 1.0, 2.0]
