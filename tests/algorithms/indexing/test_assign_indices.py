from __future__ import annotations

import copy
import os
import random

import pytest

from cctbx import crystal, sgtbx
from cctbx.sgtbx import bravais_types
from dxtbx.format import Format
from dxtbx.imageset import ImageSequence, ImageSetData
from dxtbx.model import (
    BeamFactory,
    Crystal,
    DetectorFactory,
    GoniometerFactory,
    ScanFactory,
)
from dxtbx.model.experiment_list import Experiment, ExperimentList
from dxtbx.serialize import load
from scitbx import matrix
from scitbx.math import euler_angles_as_matrix

from dials.algorithms.indexing.assign_indices import (
    AssignIndicesGlobal,
    AssignIndicesLocal,
)
from dials.array_family import flex


def random_rotation(angle_min=0, angle_max=360):
    angles = [random.uniform(angle_min, angle_max) for i in range(3)]
    print("Rotation: ", angles)
    return euler_angles_as_matrix(angles, deg=True)


@pytest.fixture
def crystal_factory():
    def _crystal_factory(space_group_symbol):
        # set random seeds so tests more reliable
        seed = 54321
        random.seed(seed)
        flex.set_random_seed(seed)

        space_group_info = sgtbx.space_group_info(symbol=space_group_symbol)
        space_group = space_group_info.group()
        unit_cell = space_group_info.any_compatible_unit_cell(
            volume=random.uniform(1e4, 1e6)
        )

        crystal_symmetry = crystal.symmetry(
            unit_cell=unit_cell, space_group=space_group
        )
        crystal_symmetry.show_summary()

        # the reciprocal matrix
        B = matrix.sqr(unit_cell.fractionalization_matrix()).transpose()
        U = random_rotation()
        A = U * B

        direct_matrix = A.inverse()
        return Crystal(
            direct_matrix[0:3],
            direct_matrix[3:6],
            direct_matrix[6:9],
            space_group=space_group,
        )

    return _crystal_factory


@pytest.fixture
def experiment():
    beam = BeamFactory.make_beam(wavelength=0.97625, sample_to_source=(0, 0, 1))

    detector = DetectorFactory.simple(
        sensor="PAD",
        distance=265.27,
        beam_centre=(210.7602, 205.27684),
        fast_direction="+x",
        slow_direction="-y",
        pixel_size=(0.172, 0.172),
        image_size=(2463, 2527),
        trusted_range=(-1, 1e8),
    )

    goniometer = GoniometerFactory.single_axis()

    scan = ScanFactory.make_scan(
        image_range=(1, 20),
        exposure_times=0.067,
        oscillation=(82, 0.15),
        epochs=[0] * 20,
    )

    isetdata = ImageSetData(
        reader=Format.Reader(None, ["path"] * len(scan)), masker=None
    )
    iset = ImageSequence(
        isetdata, beam=beam, detector=detector, goniometer=goniometer, scan=scan
    )

    return Experiment(
        imageset=iset, beam=beam, detector=detector, goniometer=goniometer, scan=scan
    )


@pytest.mark.parametrize("space_group_symbol", bravais_types.acentric)
def test_assign_indices(space_group_symbol, experiment, crystal_factory):
    cryst_model = crystal_factory(space_group_symbol)
    experiment.crystal = cryst_model

    predicted_reflections = flex.reflection_table.from_predictions(experiment)
    use_fraction = 0.3
    use_sel = flex.random_selection(
        len(predicted_reflections), int(use_fraction * len(predicted_reflections))
    )
    predicted_reflections = predicted_reflections.select(use_sel)
    miller_indices = predicted_reflections["miller_index"]
    predicted_reflections["xyzobs.mm.value"] = predicted_reflections["xyzcal.mm"]
    predicted_reflections["id"] = flex.int(len(predicted_reflections), 0)
    predicted_reflections.map_centroids_to_reciprocal_space(
        ExperimentList([experiment])
    )

    # check that local and global indexing worked equally well in absence of errors
    result = CompareGlobalLocal(experiment, predicted_reflections, miller_indices)
    assert result.misindexed_local == 0
    assert result.misindexed_global == 0

    a, b, c = map(matrix.col, experiment.crystal.get_real_space_vectors())
    relative_error = 0.02
    a *= 1 + relative_error
    b *= 1 + relative_error
    c *= 1 + relative_error

    cryst_model2 = Crystal(a, b, c, space_group=cryst_model.get_space_group())
    experiment.crystal = cryst_model2

    result = CompareGlobalLocal(experiment, predicted_reflections, miller_indices)

    # check that the local indexing did a better job given the errors in the basis vectors
    # assert result.misindexed_local < result.misindexed_global
    assert result.misindexed_local == 0
    assert result.correct_local > result.correct_global
    # usually the number misindexed is much smaller than this
    assert result.misindexed_local < (0.001 * len(result.reflections_local))

    # the reciprocal matrix
    A = matrix.sqr(cryst_model.get_A())
    A = random_rotation(angle_max=0.5) * A

    direct_matrix = A.inverse()
    cryst_model2 = Crystal(
        direct_matrix[0:3],
        direct_matrix[3:6],
        direct_matrix[6:9],
        space_group=cryst_model.get_space_group(),
    )
    experiment.crystal = cryst_model2

    result = CompareGlobalLocal(experiment, predicted_reflections, miller_indices)

    # check that the local indexing did a better job given the errors in the basis vectors
    assert result.misindexed_local <= result.misindexed_global, (
        result.misindexed_local,
        result.misindexed_global,
    )
    assert result.misindexed_local < 0.01 * result.correct_local
    assert result.correct_local >= result.correct_global
    # usually the number misindexed is much smaller than this
    assert result.misindexed_local < (0.001 * len(result.reflections_local))


class CompareGlobalLocal:
    def __init__(self, experiment, reflections, expected_miller_indices):

        index_reflections_global = AssignIndicesGlobal()
        index_reflections_local = AssignIndicesLocal()

        # index reflections using simple "global" method
        self.reflections_global = copy.deepcopy(reflections)
        self.reflections_global["id"] = flex.int(len(self.reflections_global), -1)
        self.reflections_global["imageset_id"] = flex.int(
            len(self.reflections_global), 0
        )
        index_reflections_global(self.reflections_global, ExperimentList([experiment]))
        non_zero_sel = self.reflections_global["miller_index"] != (0, 0, 0)
        assert self.reflections_global["id"].select(~non_zero_sel).all_eq(-1)
        self.misindexed_global = (
            (expected_miller_indices == self.reflections_global["miller_index"])
            .select(non_zero_sel)
            .count(False)
        )
        self.correct_global = (
            expected_miller_indices == self.reflections_global["miller_index"]
        ).count(True)

        # index reflections using xds-style "local" method
        self.reflections_local = copy.deepcopy(reflections)
        self.reflections_local["id"] = flex.int(len(self.reflections_local), -1)
        index_reflections_local(self.reflections_local, ExperimentList([experiment]))
        non_zero_sel = self.reflections_local["miller_index"] != (0, 0, 0)
        assert self.reflections_local["id"].select(~non_zero_sel).all_eq(-1)
        self.misindexed_local = (
            (expected_miller_indices == self.reflections_local["miller_index"])
            .select(non_zero_sel)
            .count(False)
        )
        self.correct_local = (
            expected_miller_indices == self.reflections_local["miller_index"]
        ).count(True)

        print(
            "Global misindexed: %d, correct: %d, total: %d"
            % (
                self.misindexed_global,
                self.correct_global,
                len(self.reflections_global),
            )
        )
        print(
            " Local misindexed: %d, correct: %d, total: %d"
            % (self.misindexed_local, self.correct_local, len(self.reflections_local))
        )


def test_index_reflections(dials_regression):
    experiments_json = os.path.join(
        dials_regression, "indexing_test_data", "i04_weak_data", "experiments.json"
    )
    experiments = load.experiment_list(experiments_json, check_format=False)
    reflections = flex.reflection_table.from_file(
        os.path.join(
            dials_regression, "indexing_test_data", "i04_weak_data", "full.pickle"
        )
    )
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    reflections["imageset_id"] = flex.int(len(reflections), 0)
    reflections["id"] = flex.int(len(reflections), -1)
    AssignIndicesGlobal(tolerance=0.3)(reflections, experiments)
    assert "miller_index" in reflections
    counts = reflections["id"].counts()
    assert dict(counts) == {-1: 1390, 0: 114692}


def test_local_multiple_rotations(dials_data):
    """Test the fix for https://github.com/dials/dials/issues/1458"""

    experiments = load.experiment_list(
        dials_data("insulin_processed", pathlib=True) / "indexed.expt",
        check_format=False,
    )
    # Override the scan range to ensure we have 4 full rotations
    experiments[0].scan.set_image_range((1, 1440))

    # Generate some predicted reflections
    reflections = flex.reflection_table.from_predictions(experiments[0], dmin=4)
    reflections["imageset_id"] = flex.int(len(reflections), 0)
    reflections["id"] = flex.int(len(reflections), -1)
    reflections["xyzobs.px.value"] = reflections["xyzcal.px"]
    reflections["xyzobs.mm.value"] = reflections["xyzcal.mm"]
    predicted_miller_indices = reflections["miller_index"]

    # Prepare reflections for indexing
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)

    # Reset miller indices to ensure they are reindexed
    reflections["miller_index"] = flex.miller_index(len(reflections), (0, 0, 0))

    # Assign indices with the correct scan oscillation
    AssignIndicesLocal()(reflections, experiments)

    # Assert we have correctly indexed all reflections
    assert (reflections["miller_index"] == (0, 0, 0)).count(True) == 0
    assert (reflections["miller_index"] == predicted_miller_indices).count(False) == 0

    # Modify the scan oscillation such that we are out by 1 degree per rotation
    experiments[0].scan.set_oscillation((0, 1 - 1 / 360), deg=True)

    # Reset miller indices and re-map to reciprocal space
    reflections["miller_index"] = flex.miller_index(len(reflections), (0, 0, 0))
    reflections["id"] = flex.int(len(reflections), -1)
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)

    # Assign indices, this time with the incorrect scan oscillation
    AssignIndicesLocal()(reflections, experiments)

    # Assert that most reflections have been indexed
    indexed_sel = reflections["miller_index"] == (0, 0, 0)
    assert indexed_sel.count(True) < 10

    # Assert that all indexed miller indices are correct
    assert (
        (reflections["miller_index"] != predicted_miller_indices) & ~indexed_sel
    ).count(True) == 0
