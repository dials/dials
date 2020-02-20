from __future__ import absolute_import, division, print_function

import dials.pychef
import pytest
from cctbx import sgtbx
from cctbx.array_family import flex
from iotbx.reflection_file_reader import any_reflection_file
from dxtbx.model import Experiment, ExperimentList, Scan


def test_observations():
    miller_indices = flex.miller_index(
        ((1, 2, 3), (-1, 2, 3), (1, -2, 3), (4, 5, 6), (4, 5, -6))
    )
    sg = sgtbx.space_group_info(symbol="I222").group()
    anomalous_flag = True
    observations = dials.pychef.Observations(miller_indices, sg, anomalous_flag)
    groups = observations.observation_groups()
    assert list(groups[(1, 2, 3)].iplus()) == [0]
    assert list(groups[(1, 2, 3)].iminus()) == [1, 2]


def test_accumulators(dials_data):
    f = dials_data("pychef").join("insulin_dials_scaled_unmerged.mtz").strpath
    reader = any_reflection_file(f)
    assert reader.file_type() == "ccp4_mtz"
    arrays = reader.as_miller_arrays(merge_equivalents=False)
    for ma in arrays:
        if ma.info().labels == ["BATCH"]:
            batches = ma
        elif ma.info().labels == ["I", "SIGI"]:
            intensities = ma
        elif ma.info().labels == ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]:
            intensities = ma
    assert intensities is not None
    assert batches is not None

    anomalous_flag = True
    if anomalous_flag:
        intensities = intensities.as_anomalous_array()

    stats = dials.pychef.Statistics(intensities, batches.data())

    # test completeness
    assert stats.iplus_comp_overall.size() == 46
    assert stats.iplus_comp_overall[45] == pytest.approx(0.9428352196431997)
    assert stats.iminus_comp_overall[44] == pytest.approx(0.07769038941108766)
    assert stats.ieither_comp_overall[1] == pytest.approx(0.03721465566852101)
    assert stats.iboth_comp_overall[-1] == pytest.approx(0.0779781315940917)

    # test rcp,scp
    print(list(stats.rcp))
    print(list(stats.scp))
    assert stats.rcp[45] == pytest.approx(0.04844584637191411)
    assert stats.scp[45] == pytest.approx(0.9201295853457298)

    # test Rd
    print(list(stats.rd))
    assert stats.rd[0] == pytest.approx(0.05234416616316846)


def test_interpret_images_to_doses_options():
    """Test handling of command line options for experiments input."""
    params = dials.pychef.phil_scope.extract()
    experiments = ExperimentList()
    experiments.append(Experiment(scan=Scan(image_range=(1, 10), oscillation=(0, 1.0))))
    experiments.append(Experiment(scan=Scan(image_range=(1, 20), oscillation=(0, 1.0))))
    experiments.append(Experiment(scan=Scan(image_range=(1, 10), oscillation=(0, 1.0))))

    # Default
    starting_doses, dpi = dials.pychef.interpret_images_to_doses_options(
        params, experiments
    )
    assert starting_doses == [0, 0, 0]
    assert dpi == [1.0, 1.0, 1.0]

    # Multi-sweep measurements on same crystal
    params.dose.experiments.shared_crystal = True
    starting_doses, dpi = dials.pychef.interpret_images_to_doses_options(
        params, experiments
    )
    assert starting_doses == [0, 10, 30]
    assert dpi == [1.0, 1.0, 1.0]

    # Specify starting doses
    params.dose.experiments.shared_crystal = False
    params.dose.experiments.starting_doses = [0, 20, 0]
    starting_doses, dpi = dials.pychef.interpret_images_to_doses_options(
        params, experiments
    )
    assert starting_doses == [0, 20, 0]
    assert dpi == [1.0, 1.0, 1.0]

    # Specify doses per image and shared crystal
    params.dose.experiments.starting_doses = None
    params.dose.experiments.shared_crystal = True
    params.dose.experiments.dose_per_image = [1.0, 2.0, 1.0]
    starting_doses, dpi = dials.pychef.interpret_images_to_doses_options(
        params, experiments
    )
    assert starting_doses == [0, 10, 50]
    assert dpi == [1.0, 2.0, 1.0]

    # Test error is raised if bad input values for starting doses or dose per image.
    params.dose.experiments.shared_crystal = False
    params.dose.experiments.starting_doses = [0, 1]
    with pytest.raises(ValueError):
        _, __ = dials.pychef.interpret_images_to_doses_options(params, experiments)
    params.dose.experiments.starting_doses = None
    params.dose.experiments.dose_per_image = [1.0, 2.0]
    with pytest.raises(ValueError):
        _, __ = dials.pychef.interpret_images_to_doses_options(params, experiments)
