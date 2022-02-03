import pytest

from dxtbx.serialize import load

from dials.algorithms.profile_model.ellipsoid.algorithm import initial_integrator
from dials.algorithms.profile_model.ellipsoid.indexer import reindex
from dials.algorithms.profile_model.ellipsoid.refiner import RefinerData
from dials.array_family import flex


def test_initial_integrator(dials_data):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    _ = dials_data("cunir_serial", pathlib=True)
    refls = flex.reflection_table.from_file(str(ssx / "indexed.refl"))
    expts = load.experiment_list(str(ssx / "indexed.expt"))[0:1]
    refls = refls.select_on_experiment_identifiers([expts[0].identifier])

    refls = reindex(refls, expts[0])
    refls, sigma_d = initial_integrator(expts, refls)

    assert sigma_d == pytest.approx(0.00062, abs=1e-5)
    assert list(refls["intensity.sum.value"]) == pytest.approx(
        [20.51, 183.67, 132.11, 103.50, 150.80, 418.075, 412.45, 524.93, 138.86]
        + [161.51, 77.45, 218.50, 109.02, 375.47, 64.86, 81.46, 136.36, 391.07]
        + [107.60, 408.10, 348.73, 143.41, 115.20, 72.93],
        abs=1e-2,
    )

    data = RefinerData.from_reflections(expts[0], refls)

    assert list(data.s0) == pytest.approx([0.0, 0.0, -0.72669], abs=1e-5)
    assert list(data.h_list) == (
        [(-11, 16, 30), (-19, 6, 18), (-18, 6, 17), (-26, -11, 6), (-19, -15, -2)]
        + [(-13, -8, -1), (-7, -15, -8), (3, 3, 2), (-4, -17, -10), (4, -5, -5)]
        + [(-3, -21, -12), (8, -7, -7), (5, -19, -13), (20, 7, 5), (15, -1, -3)]
        + [(12, -8, -8), (16, -1, -3), (21, 6, 4), (25, 8, 7), (18, -8, -8)]
        + [(27, 3, 2), (21, -10, -9), (33, 5, 6), (30, 0, 0)]
    )
    assert list(data.ctot_list) == pytest.approx(
        [91.912, 211.155, 189.053, 164.854, 179.438, 466.371, 450.249, 584.141]
        + [172.747, 208.006, 126.703, 255.102, 157.813, 426.994, 124.453, 130.475]
        + [167.967, 455.962, 155.182, 463.591, 387.441, 196.885, 153.05, 138.08],
        abs=1e-3,
    )
    expected_mobs = [
        (1e-08, -2.22e-07),
        (-8e-09, -2.78e-07),
        (1.5e-08, -9.4e-08),
        (-1.4e-08, -2.01e-07),
        (5.3e-08, -1.78e-07),
        (-1e-09, -2.5e-08),
        (4e-09, -3.7e-08),
        (-1e-09, -7e-09),
        (-4e-09, -7.2e-08),
        (8e-09, -6.5e-08),
        (0.0, -8e-08),
        (1.4e-08, -3e-08),
        (2e-08, -1.48e-07),
        (8e-09, -7.2e-08),
        (2.2e-08, -1.14e-07),
        (3.3e-08, -8.6e-08),
        (-1.8e-08, -5.9e-08),
        (-1.5e-08, -4.1e-08),
        (-4e-09, -1.67e-07),
        (1e-09, -2.4e-08),
        (-4.2e-08, -1.29e-07),
        (-2e-09, -7.5e-08),
        (3.3e-08, -2.55e-07),
        (-5.5e-08, -2e-07),
    ]
    for i in range(len(expected_mobs)):
        assert data.mobs_list[i] == pytest.approx(expected_mobs[i], abs=1e-9)
    # print(list((round(i[0], 9), round(i[1], 9)) for i in data.mobs_list))
    assert list(data.sobs_list) == pytest.approx(
        [
            1.85e-07,
            -1.5e-08,
            -1.5e-08,
            2.86e-07,
            4.9e-07,
            1.6e-08,
            1.6e-08,
            5.02e-07,
            1.68e-07,
            -2.7e-08,
            -2.7e-08,
            1.81e-07,
            2.21e-07,
            2.3e-08,
            2.3e-08,
            3.36e-07,
            3.48e-07,
            -1.09e-07,
            -1.09e-07,
            3.6e-07,
            1.29e-07,
            3e-09,
            3e-09,
            8e-08,
            9.3e-08,
            -1.2e-08,
            -1.2e-08,
            1.02e-07,
            7.7e-08,
            1.6e-08,
            1.6e-08,
            7.5e-08,
            1.54e-07,
            1.2e-08,
            1.2e-08,
            1.75e-07,
            3.28e-07,
            -4.5e-08,
            -4.5e-08,
            3.99e-07,
            2.12e-07,
            1e-09,
            1e-09,
            1.57e-07,
            2.04e-07,
            -5.1e-08,
            -5.1e-08,
            1.2e-07,
            3.42e-07,
            -4.3e-08,
            -4.3e-08,
            3.06e-07,
            1.39e-07,
            -1.6e-08,
            -1.6e-08,
            1.65e-07,
            2.96e-07,
            -7e-08,
            -7e-08,
            3.71e-07,
            4.45e-07,
            -1.05e-07,
            -1.05e-07,
            2.55e-07,
            3.76e-07,
            5.5e-08,
            5.5e-08,
            1.78e-07,
            1.24e-07,
            3.4e-08,
            3.4e-08,
            9.1e-08,
            3.82e-07,
            8e-09,
            8e-09,
            2.93e-07,
            1.95e-07,
            -1e-09,
            -1e-09,
            5.5e-08,
            2.39e-07,
            7.5e-08,
            7.5e-08,
            2.28e-07,
            2.98e-07,
            4e-09,
            4e-09,
            1.47e-07,
            4.6e-07,
            -4.3e-08,
            -4.3e-08,
            3.45e-07,
            2.75e-07,
            8.7e-08,
            8.7e-08,
            3.14e-07,
        ],
        abs=1e-9,
    )
