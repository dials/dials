from __future__ import annotations

import os
import pathlib

from dials.algorithms.correlation.analysis import CorrelationMatrix
from dials.command_line.correlation_matrix import phil_scope
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files


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

    params.output.json = "dials.correlation_matrix.json"

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    reflections = parse_multiple_datasets(reflections)
    assert len(experiments) == len(reflections)
    assert len(experiments) > 1
    experiments, reflections = assign_unique_identifiers(experiments, reflections)
    matrices = CorrelationMatrix(experiments, reflections, params)
    matrices.calculate_matrices()
    matrices.convert_to_json()
    matrices.output_json()
    assert pathlib.Path("dials.correlation_matrix.json").is_file()
