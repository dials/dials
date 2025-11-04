from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

from dials.algorithms.correlation.analysis import CorrelationMatrix
from dials.command_line.correlation_matrix import phil_scope
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

parent_path = pathlib.Path(__file__).parent.resolve()

# cows + people from Thompson, A. J. et al. (2025) Acta Cryst. D81, 278-290 - two clear clusters
data_1 = np.loadtxt(parent_path / "test_coords/data_1.txt")
expected_1 = np.loadtxt(parent_path / "test_coords/labels_1.txt")

# cryo cows + pigs + people from Thompson, A. J. et al. (2025) Acta Cryst. D81, 278-290 - three clear clusters
data_2 = np.loadtxt(parent_path / "test_coords/data_2.txt")
expected_2 = np.loadtxt(parent_path / "test_coords/labels_2.txt")

# 4 x CPVs from VMXm - four clear clusters
data_3 = np.loadtxt(parent_path / "test_coords/data_3.txt")
expected_3 = np.loadtxt(parent_path / "test_coords/labels_3.txt")

# room temp cows + pigs + people from Thompson, A. J. et al. (2025) Acta Cryst. D81, 278-290 - three clusters + noise
data_4 = np.loadtxt(parent_path / "test_coords/data_4.txt")
expected_4 = np.loadtxt(parent_path / "test_coords/labels_4.txt")

# example made up coordinates - one clear cluster + noise
data_5 = np.loadtxt(parent_path / "test_coords/data_5.txt")
expected_5 = np.loadtxt(parent_path / "test_coords/labels_5.txt")

# single cluster + noise, high dimension
data_6 = np.loadtxt(parent_path / "test_coords/data_6.txt")
expected_6 = np.loadtxt(parent_path / "test_coords/labels_6.txt")

# One large cluster + one small cluster
data_7 = np.loadtxt(parent_path / "test_coords/data_7.txt")
expected_7 = np.loadtxt(parent_path / "test_coords/labels_7.txt")

# Two tight clusters with one obvious outlier
data_8 = np.loadtxt(parent_path / "test_coords/data_8.txt")
expected_8 = np.loadtxt(parent_path / "test_coords/labels_8.txt")


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


# Test for very clearcut cases
# initial guesses taken from heuristic (Thompson, A.J. et al 2025) for known datasets
@pytest.mark.parametrize(
    "coordinates,expected_labels,initial_min_samples",
    [
        (data_1, expected_1, 5),
        (data_2, expected_2, 6),
        (data_3, expected_3, 8),
        (data_5, expected_5, 10),
        (data_8, expected_8, 5),
    ],
)
def test_optics_classification_definitive(
    coordinates, expected_labels, initial_min_samples
):
    _, _, _, actual_labels, _ = CorrelationMatrix.optimise_clustering(
        coordinates, initial_min_samples=initial_min_samples
    )

    assert np.array_equal(actual_labels, expected_labels)


# Test for more ambiguous datasets - set up to add more in future if needed
@pytest.mark.parametrize(
    "coordinates,expected_labels,initial_min_samples",
    [
        (data_4, expected_4, 27),
        (data_6, expected_6, 5),
        (data_7, expected_7, 40),
    ],
)
def test_optics_classification_variable(
    coordinates, expected_labels, initial_min_samples
):
    _, _, _, actual_labels, _ = CorrelationMatrix.optimise_clustering(
        coordinates, initial_min_samples=initial_min_samples
    )

    differences = actual_labels != expected_labels

    difference_indices = np.argwhere(differences)

    values = [
        (actual_labels[tuple(idx)], expected_labels[tuple(idx)])
        for idx in difference_indices
    ]

    noise_change_count = 0
    changed_classification = []

    # Expected as algorithms change that a small number of datasets may switch between noise/cluster
    for i in values:
        if i[0] == -1 or i[1] == -1:
            noise_change_count += 1
        else:
            changed_classification.append(i)

    assert noise_change_count < 0.03 * len(coordinates)  # 3% error in noise permitted
    assert (
        len(changed_classification) == 0
    )  # do not want any clusters to change identity
