from __future__ import annotations

import pytest

from cctbx import sgtbx, uctbx
from dxtbx.imageset import ImageSet
from dxtbx.serialize import load

from dials.algorithms.indexing import lattice_search, stills_indexer
from dials.array_family import flex
from dials.command_line.index import phil_scope
from dials.command_line.slice_sequence import slice_experiments, slice_reflections


@pytest.fixture
def i04_weak_data(dials_data):
    data_dir = dials_data("i04_weak_data")
    reflections_path = data_dir / "full.pickle"
    experiments_path = data_dir / "experiments_import.json"

    reflections = flex.reflection_table.from_file(reflections_path)
    experiments = load.experiment_list(experiments_path, check_format=False)

    return {
        "reflections": slice_reflections(reflections, [(1, 20)]),
        "experiments": slice_experiments(experiments, [(1, 20)]),
    }


@pytest.mark.parametrize(
    "indexing_method,space_group,unit_cell",
    (
        ("fft1d", None, None),
        ("fft3d", None, None),
        ("real_space_grid_search", "P422", (57.8, 57.8, 150.0, 90, 90, 90)),
    ),
)
def test_BasisVectorSearch_i04_weak_data(
    i04_weak_data, indexing_method, space_group, unit_cell
):
    reflections = i04_weak_data["reflections"]
    experiments = i04_weak_data["experiments"]
    params = phil_scope.fetch().extract()
    params.indexing.refinement_protocol.n_macro_cycles = 2
    params.indexing.basis_vector_combinations.max_refine = 5
    params.indexing.method = indexing_method
    if unit_cell is not None:
        params.indexing.known_symmetry.unit_cell = uctbx.unit_cell(unit_cell)
    if space_group is not None:
        params.indexing.known_symmetry.space_group = sgtbx.space_group_info(space_group)
    idxr = lattice_search.BasisVectorSearch(reflections, experiments, params)
    idxr.index()

    indexed_experiments = idxr.refined_experiments
    assert len(indexed_experiments) == 1
    assert indexed_experiments[0].crystal.get_unit_cell().parameters() == pytest.approx(
        (57.752, 57.776, 150.013, 90.0101, 89.976, 90.008), rel=1e-3
    )


@pytest.mark.parametrize(
    "indexing_method,space_group,unit_cell",
    (
        ("fft3d", "P422", (57.8, 57.8, 150.0, 90, 90, 90)),
        ("fft1d", "P422", (57.8, 57.8, 150.0, 90, 90, 90)),
        ("real_space_grid_search", "P422", (57.8, 57.8, 150.0, 90, 90, 90)),
        ("low_res_spot_match", "P422", (57.8, 57.8, 150.0, 90, 90, 90)),
    ),
)
def test_stills_indexer_methods_i04_weak_data(
    i04_weak_data, indexing_method, space_group, unit_cell
):
    reflections = slice_reflections(i04_weak_data["reflections"], [(1, 2)])
    experiments = slice_experiments(i04_weak_data["experiments"], [(1, 2)])
    for experiment in experiments:
        experiment.imageset = ImageSet(
            experiment.imageset.data(), experiment.imageset.indices()
        )
        experiment.imageset.set_scan(None)
        experiment.imageset.set_goniometer(None)
        experiment.scan = None
        experiment.goniometer = None
    params = phil_scope.fetch().extract()
    params.indexing.method = indexing_method
    params.indexing.basis_vector_combinations.max_refine = 5
    if unit_cell is not None:
        params.indexing.known_symmetry.unit_cell = uctbx.unit_cell(unit_cell)
    if space_group is not None:
        params.indexing.known_symmetry.space_group = sgtbx.space_group_info(space_group)
    try:
        idxr = stills_indexer.StillsIndexerBasisVectorSearch(
            reflections, experiments, params
        )
    except RuntimeError:
        idxr = stills_indexer.StillsIndexerLatticeSearch(
            reflections, experiments, params
        )
    idxr.index()

    indexed_experiments = idxr.refined_experiments
    assert len(indexed_experiments) == 1
    assert indexed_experiments[0].crystal.get_unit_cell().parameters() == pytest.approx(
        (57.752, 57.776, 150.013, 90.0101, 89.976, 90.008), rel=1e-2
    )
