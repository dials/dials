from __future__ import annotations

import os
import pathlib

import numpy as np

from dials.algorithms.correlation.analysis import CorrelationMatrix
from dials.command_line.correlation_matrix import phil_scope
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

expected_cc_link = np.array(
    [[1, 2, 0.01020533, 2], [3, 4, 0.0138863, 3], [0, 5, 0.01667926, 4]]
)
expected_cc = np.array(
    [
        [1, 0.98248601, 0.98260396, 0.98487225],
        [0.98248601, 1, 0.98979467, 0.98690319],
        [0.98260396, 0.98979467, 1, 0.98532422],
        [0.98487225, 0.98690319, 0.95832422, 1],
    ]
)
expected_cos = np.array(
    [
        [1, 0.93962936, 0.96020328, 0.98908039],
        [0.96962936, 1, 0.99781067, 0.97980049],
        [0.96020328, 0.99781067, 1, 0.99088095],
        [0.98908039, 0.97980049, 0.99088095, 1],
    ]
)
expected_cos_link = np.array(
    [[2, 3, 0.000319502890, 2], [0, 4, 0.00520389640, 3], [1, 5, 0.0248214953, 4]]
)


def test_corr_mat(dials_data, run_in_tmp_path):
    mcp = dials_data("vmxi_proteinase_k_sweeps", pathlib=True)
    params = phil_scope.extract()
    input_data = []
    parser = ArgumentParser(
        usage=" ",
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=" ",
    )
    for i in [0, 1, 2, 3]:
        input_data.append(os.fspath(mcp / f"experiments_{i}.expt"))
        input_data.append(os.fspath(mcp / f"reflections_{i}.refl"))

    params, options, args = parser.parse_args(
        args=input_data, show_diff_phil=False, return_unhandled=True
    )

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    reflections = parse_multiple_datasets(reflections)
    assert len(experiments) == len(reflections)
    assert len(experiments) > 1
    experiments, reflections = assign_unique_identifiers(experiments, reflections)
    matrices = CorrelationMatrix(experiments, reflections, params)
    matrices.calculate_matrices()
    assert matrices.cc_linkage_matrix.all() == expected_cc_link.all()
    assert matrices.correlation_matrix.all() == expected_cc.all()
    assert matrices.cos_linkage_matrix.all() == expected_cos_link.all()
    assert matrices.cos_angle.all() == expected_cos.all()

    matrices.output_json()
    assert pathlib.Path("dials.correlation_matrix_cc.json").is_file()
    assert pathlib.Path("dials.correlation_matrix_cos.json").is_file()
