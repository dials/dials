from __future__ import annotations

from os.path import join

import numpy as np

from dxtbx import flumpy
from dxtbx.model.experiment_list import ExperimentListFactory

from dials.array_family import flex
from dials_algorithms_tof_integration_ext import (
    tof_calculate_ellipse_shoebox_mask,
    tof_calculate_seed_skewness_shoebox_mask,
)
from dials_tof_scaling_ext import (
    tof_extract_shoeboxes_to_reflection_table,
)


def test_tof_calculate_shoebox_mask(dials_data):
    def pad_bboxes(bboxes, max_size, padding):
        for i in range(len(bboxes)):
            b = list(bboxes[i])
            b[0] -= padding[0]
            b[1] += padding[0]
            b[2] -= padding[1]
            b[3] += padding[1]
            b[4] -= padding[2]
            b[5] += padding[2]

            if b[0] < 0:
                b[0] = 0
            if b[1] > max_size[0]:
                b[1] = max_size[0]

            if b[2] < 0:
                b[2] = 0
            if b[3] > max_size[1]:
                b[3] = max_size[1]

            if b[4] < 0:
                b[4] = 0
            if b[5] > max_size[2]:
                b[5] = max_size[2]

            bboxes[i] = tuple(b)

        return bboxes

    image_file = join(dials_data("isis_sxd_example_data"), "sxd_nacl_run.nxs")
    experiments = ExperimentListFactory.from_filenames([image_file])
    reflections = flex.reflection_table.from_msgpack_file(
        join(dials_data("isis_sxd_nacl_processed"), "strong.refl")
    )
    _, _, pz = reflections["xyzobs.px.value"].parts()
    reflections = reflections.select(pz < 100)
    experiments[0].imageset = experiments[0].imageset[:105]
    experiments[0].scan = experiments[0].scan[:105]

    expt_data = experiments[0].imageset
    bboxes = pad_bboxes(reflections["bbox"], max_size=(64, 64, 105), padding=(1, 1, 1))

    reflections["shoebox"] = flex.shoebox(
        reflections["panel"],
        bboxes,
        allocate=False,
        flatten=False,
    )

    tof_extract_shoeboxes_to_reflection_table(
        reflections, experiments[0], expt_data, False
    )

    tof_calculate_ellipse_shoebox_mask(reflections, experiments[0])

    mask_arr = flumpy.to_numpy(reflections["shoebox"][0].mask)
    expected_ellipse_mask = np.array(
        [
            [[3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 3, 3]],
            [[3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 3, 5], [3, 3, 3, 3], [3, 3, 3, 3]],
            [[3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 3, 5], [3, 3, 3, 3], [3, 3, 3, 3]],
            [[3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 5, 5], [3, 3, 3, 3], [3, 3, 3, 3]],
        ],
        dtype=np.int32,
    )

    assert np.array_equal(mask_arr, expected_ellipse_mask)

    tof_calculate_seed_skewness_shoebox_mask(reflections, experiments[0], 0.01, 10)

    mask_arr = flumpy.to_numpy(reflections["shoebox"][0].mask)

    expected_ss_mask = np.array(
        [
            [[3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 5, 5], [3, 3, 3, 5], [3, 3, 3, 3]],
            [[3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 5, 5], [3, 3, 5, 3], [3, 3, 3, 3]],
            [[3, 3, 3, 3], [3, 3, 3, 3], [3, 3, 5, 5], [3, 3, 3, 5], [3, 3, 3, 3]],
            [[3, 3, 3, 3], [3, 3, 5, 3], [3, 3, 3, 5], [3, 3, 3, 3], [3, 3, 3, 3]],
        ],
        dtype=np.int32,
    )

    assert np.array_equal(mask_arr, expected_ss_mask)
