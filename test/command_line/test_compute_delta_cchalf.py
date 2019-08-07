"""Tests for dials.compute_delta_cchalf."""
from __future__ import absolute_import, division, print_function

import copy
import os

import procrunner
from dials.algorithms.statistics.delta_cchalf import PerImageCChalfStatistics
from iotbx.reflection_file_reader import any_reflection_file


def test_compute_delta_cchalf_scaled_data(dials_data, tmpdir):
    """Test dials.compute_delta_cchalf on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled")
    refls = location.join("scaled_20_25.refl").strpath
    expts = location.join("scaled_20_25.expt").strpath

    # set cutoff to 0.0 to force one to be 'rejected'
    command = [
        "dials.compute_delta_cchalf",
        refls,
        expts,
        "stdcutoff=0.0",
        "output.reflections=filtered.refl",
        "output.experiments=filtered.expt",
    ]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("filtered.expt").check()
    assert tmpdir.join("filtered.refl").check()


def test_compute_delta_cchalf(dials_regression):
    """Test compute delta cchalf on an integrated mtz."""
    filename = os.path.join(
        dials_regression, "delta_cchalf_test_data", "test.XDS_ASCII.mtz"
    )

    # Read the mtz file
    reader = any_reflection_file(filename)

    # Get the columns as miller arrays
    miller_arrays = reader.as_miller_arrays(merge_equivalents=False)

    # Select the desired columns
    intensities = None
    batches = None
    for array in miller_arrays:
        if array.info().labels == ["I", "SIGI"]:
            intensities = array
        if array.info().labels == ["BATCH"]:
            batches = array
    assert intensities is not None
    assert batches is not None
    assert len(batches.data()) == len(intensities.data())

    # Get the unit cell and space group
    unit_cell = intensities.unit_cell()
    space_group = intensities.crystal_symmetry().space_group()

    # The reflection data
    miller_index = intensities.indices()
    batch = batches.data()
    intensity = intensities.data()
    variance = intensities.sigmas() ** 2

    # Create unit cell list
    min_batch = min(batch)
    dataset = batch - min_batch
    num_datasets = max(dataset) + 1
    unit_cell_list = [unit_cell for _ in range(num_datasets)]
    identifiers = list(set(dataset))
    # Add in dummy images for now
    images = copy.deepcopy(batch)

    # Compute the CC 1/2 Stats
    statistics = PerImageCChalfStatistics(
        miller_index,
        identifiers,
        dataset,
        images,
        intensity,
        variance,
        unit_cell_list,
        space_group,
        nbins=1,
    )

    mean_cchalf = statistics.mean_cchalf()
    cchalf_i = statistics.cchalf_i()

    assert abs(100 * mean_cchalf - 94.582) < 1e-3
    assert abs(100 * cchalf_i[0] - 79.587) < 1e-3
    assert abs(100 * cchalf_i[1] - 94.238) < 1e-3
