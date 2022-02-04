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
