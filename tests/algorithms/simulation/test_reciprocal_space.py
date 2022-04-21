from __future__ import annotations

import math

from dxtbx.model.experiment_list import ExperimentListFactory

from dials.algorithms.simulation.reciprocal_space import Simulator


def test_simulation(dials_data):
    experiments = ExperimentListFactory.from_json_file(
        dials_data("centroid_test_data", pathlib=True) / "experiments.json",
        check_format=False,
    )

    sigma_b = 0.058 * math.pi / 180
    sigma_m = 0.157 * math.pi / 180
    n_sigma = 3

    N = 100
    In = 1000
    B = 10
    simulate = Simulator(experiments[0], sigma_b, sigma_m, n_sigma)
    simulate.with_random_intensity(N, In, B, 0, 0, 0)
