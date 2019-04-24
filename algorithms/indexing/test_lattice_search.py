from __future__ import absolute_import, division, print_function

import pytest
import py.path

from cctbx import sgtbx, uctbx
from dxtbx.serialize import load
from dials.algorithms.indexing import lattice_search
from dials.command_line.index import phil_scope
from dials.command_line.slice_sweep import slice_experiments, slice_reflections
from dials.array_family import flex


@pytest.fixture
def i04_weak_data(dials_regression):
    data_dir = py.path.local(dials_regression).join(
        "indexing_test_data", "i04_weak_data"
    )
    reflections_path = data_dir.join("full.pickle").strpath
    experiments_path = data_dir.join("experiments_import.json").strpath

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
    indexed_reflections = idxr.refined_reflections
    assert len(indexed_experiments) == 1
    assert indexed_experiments[0].crystal.get_unit_cell().parameters() == pytest.approx(
        (57.752, 57.776, 150.013, 90.0101, 89.976, 90.008), rel=1e-3
    )
