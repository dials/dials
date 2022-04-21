from __future__ import annotations


def test_run(dials_data):
    experiments = dials_data("centroid_test_data", pathlib=True) / "experiments.json"

    from dxtbx.model.experiment_list import ExperimentListFactory

    experiments = ExperimentListFactory.from_json_file(str(experiments))

    from dials.algorithms.profile_model.modeller import EwaldSphereSampler

    beam = experiments[0].beam
    detector = experiments[0].detector
    goniometer = experiments[0].goniometer
    scan = experiments[0].scan

    sampler = EwaldSphereSampler(beam, detector, goniometer, scan, 1)

    assert sorted(sampler.nearest_n(0)) == sorted([0, 1, 2, 3, 4, 5, 6, 7, 8])

    assert sorted(sampler.nearest_n(1)) == sorted([0, 1, 2, 8, 9, 10])
    assert sorted(sampler.nearest_n(2)) == sorted([0, 2, 3, 1, 11, 12])
    assert sorted(sampler.nearest_n(3)) == sorted([0, 3, 4, 2, 13, 14])
    assert sorted(sampler.nearest_n(4)) == sorted([0, 4, 5, 3, 15, 16])
    assert sorted(sampler.nearest_n(5)) == sorted([0, 5, 6, 4, 17, 18])
    assert sorted(sampler.nearest_n(6)) == sorted([0, 6, 7, 5, 19, 20])
    assert sorted(sampler.nearest_n(7)) == sorted([0, 7, 8, 6, 21, 22])
    assert sorted(sampler.nearest_n(8)) == sorted([0, 8, 1, 7, 23, 24])

    assert sorted(sampler.nearest_n(9)) == sorted([1, 9, 10, 24, 25, 26])
    assert sorted(sampler.nearest_n(10)) == sorted([1, 10, 11, 9, 27, 28])
    assert sorted(sampler.nearest_n(24)) == sorted([8, 24, 9, 23, 55, 56])

    assert sorted(sampler.nearest_n(25)) == sorted([9, 25, 26, 56])
    assert sorted(sampler.nearest_n(26)) == sorted([9, 26, 27, 25])
    assert sorted(sampler.nearest_n(56)) == sorted([24, 56, 25, 55])
