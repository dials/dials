from __future__ import annotations

import json
import pathlib

import pytest

from dials.algorithms.correlation.analysis import CorrelationMatrix
from dials.command_line.correlation_matrix import phil_scope
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files


@pytest.fixture()
def proteinase_k(dials_data):
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
        input_data.append(str(mcp / f"experiments_{i}.expt"))
        input_data.append(str(mcp / f"reflections_{i}.refl"))

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
    yield experiments, reflections, params


def test_corr_mat(proteinase_k, run_in_tmp_path):
    experiments, reflections, params = proteinase_k
    matrices = CorrelationMatrix(experiments, reflections, params)
    matrices.calculate_matrices()
    matrices.output_json()
    assert pathlib.Path("dials.correlation_matrix.json").is_file()


def test_filtered_corr_mat(proteinase_k, run_in_tmp_path):
    experiments, reflections, params = proteinase_k
    ids_to_identifiers_map = {}
    for table in reflections:
        ids_to_identifiers_map.update(table.experiment_identifiers())

    # Simulate filtered dataset by multiplex
    id_to_remove = [ids_to_identifiers_map[2]]
    ids_to_identifiers_map.pop(2)
    reflections.pop(2)
    experiments.remove_on_experiment_identifiers(id_to_remove)

    matrices = CorrelationMatrix(
        experiments, reflections, params, ids_to_identifiers_map
    )
    matrices.calculate_matrices()
    matrices.output_json()
    assert pathlib.Path("dials.correlation_matrix.json").is_file()

    expected_ids = [[0, 3], [0, 1, 3]]

    # Check main algorithm correct with filtering
    for i, j in zip(matrices.correlation_clusters, expected_ids):
        assert i.labels == j

    # Check json output also correct
    with open(pathlib.Path("dials.correlation_matrix.json")) as f:
        data = json.load(f)

    assert len(data["correlation_matrix_clustering"]) == len(expected_ids)
    for i, j in zip(data["correlation_matrix_clustering"], expected_ids):
        assert len(data["correlation_matrix_clustering"][i]["datasets"]) == len(j)
        for a, e in zip(data["correlation_matrix_clustering"][i]["datasets"], j):
            assert a == ids_to_identifiers_map[e]
