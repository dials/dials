"""
Test combination of multiple experiments and reflections files.
"""

from __future__ import absolute_import, division, print_function

import copy
import os
import procrunner

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from dials.util import Sorry
from dials.array_family import flex
import dials.command_line.combine_experiments as combine_experiments


def test(dials_regression, run_in_tmpdir):
    data_dir = os.path.join(
        dials_regression, "refinement_test_data", "multi_narrow_wedges"
    )

    input_range = list(range(2, 49))
    for i in (8, 10, 15, 16, 34, 39, 45):
        input_range.remove(i)

    phil_input = "\n".join(
        (
            "  input.experiments={0}/data/sweep_%03d/experiments.json\n"
            + "  input.reflections={0}/data/sweep_%03d/reflections.pickle"
        )
        % (i, i)
        for i in input_range
    )
    # assert phil_input == "\n" + phil_input2 + "\n "

    input_phil = (
        phil_input.format(data_dir)
        + """
 reference_from_experiment.beam=0
 reference_from_experiment.scan=0
 reference_from_experiment.goniometer=0
 reference_from_experiment.detector=0
 """
    )

    with open("input.phil", "w") as phil_file:
        phil_file.writelines(input_phil)

    result = procrunner.run(["dials.combine_experiments", "input.phil"])
    assert not result.returncode and not result.stderr

    # load results
    exp = ExperimentListFactory.from_json_file("combined.expt", check_format=False)
    ref = flex.reflection_table.from_pickle("combined.refl")

    # test the experiments
    assert len(exp) == 103
    assert len(exp.crystals()) == 103
    assert len(exp.beams()) == 1
    assert len(exp.scans()) == 1
    assert len(exp.detectors()) == 1
    assert len(exp.goniometers()) == 1
    for e in exp:
        assert e.imageset is not None

    # test the reflections
    assert len(ref) == 11689

    result = procrunner.run(
        ["dials.split_experiments", "combined.expt", "combined.refl"]
    )
    assert not result.returncode and not result.stderr

    for i, e in enumerate(exp):
        assert os.path.exists("split_%03d.expt" % i)
        assert os.path.exists("split_%03d.refl" % i)

        exp_single = ExperimentListFactory.from_json_file(
            "split_%03d.expt" % i, check_format=False
        )
        ref_single = flex.reflection_table.from_pickle("split_%03d.refl" % i)

        assert len(exp_single) == 1
        assert exp_single[0].crystal == e.crystal
        assert exp_single[0].beam == e.beam
        assert exp_single[0].detector == e.detector
        assert exp_single[0].scan == e.scan
        assert exp_single[0].goniometer == e.goniometer
        assert exp_single[0].imageset == e.imageset
        assert len(ref_single) == len(ref.select(ref["id"] == i))
        assert ref_single["id"].all_eq(0)

    result = procrunner.run(
        ["dials.split_experiments", "combined.expt", "output.experiments_prefix=test"]
    )
    assert not result.returncode and not result.stderr

    for i in range(len(exp)):
        assert os.path.exists("test_%03d.expt" % i)

    # Modify a copy of the detector
    detector = copy.deepcopy(exp.detectors()[0])
    panel = detector[0]
    x, y, z = panel.get_origin()
    panel.set_frame(panel.get_fast_axis(), panel.get_slow_axis(), (x, y, z + 10))
    # Set half of the experiments to the new detector
    for i in range(len(exp) // 2):
        exp[i].detector = detector
    exp.as_json("modded.expt")

    result = procrunner.run(
        [
            "dials.split_experiments",
            "modded.expt",
            "combined.refl",
            "output.experiments_prefix=test_by_detector",
            "output.reflections_prefix=test_by_detector",
            "by_detector=True",
        ]
    )
    assert not result.returncode and not result.stderr

    for i in range(2):
        assert os.path.exists("test_by_detector_%03d.expt" % i)
        assert os.path.exists("test_by_detector_%03d.refl" % i)
    assert not os.path.exists("test_by_detector_%03d.expt" % 2)
    assert not os.path.exists("test_by_detector_%03d.refl" % 2)

    # Now do test when input has identifiers set
    reflections = flex.reflection_table().from_pickle("combined.refl")
    explist = ExperimentListFactory.from_json_file("combined.expt", check_format=False)
    # set string identifiers as nonconsecutive 0,2,4,6....
    for i, exp in enumerate(explist):
        assert i in reflections["id"]
        reflections.experiment_identifiers()[i] = str(i * 2)
        exp.identifier = str(i * 2)
    reflections.as_pickle("assigned.refl")
    explist.as_json("assigned.expt")

    result = procrunner.run(
        ["dials.split_experiments", "assigned.expt", "assigned.refl"]
    )
    assert not result.returncode and not result.stderr

    for i in range(len(explist)):
        assert os.path.exists("split_%03d.expt" % i)
        assert os.path.exists("split_%03d.refl" % i)

        exp_single = ExperimentListFactory.from_json_file(
            "split_%03d.expt" % i, check_format=False
        )
        ref_single = flex.reflection_table.from_pickle("split_%03d.refl" % i)

        assert len(exp_single) == 1
        # resets all ids to 0, but keeps mapping to unique identifier.
        # doesn't have to be set to 0 but doing this to keep more consistent with
        # other dials programs
        assert ref_single["id"].all_eq(0)
        assert ref_single.experiment_identifiers()[0] == str(i * 2)

    # update modded experiments to have same identifiers as assigned_experiments
    moddedlist = ExperimentListFactory.from_json_file("modded.expt", check_format=False)
    for i, exp in enumerate(moddedlist):
        exp.identifier = str(i * 2)
    moddedlist.as_json("modded.expt")

    result = procrunner.run(
        [
            "dials.split_experiments",
            "modded.expt",
            "assigned.refl",
            "output.experiments_prefix=test_by_detector",
            "output.reflections_prefix=test_by_detector",
            "by_detector=True",
        ]
    )
    assert not result.returncode and not result.stderr

    # Expect each datasets to have ids from 0..50 with experiment identifiers
    # all kept from before 0,2,4,6,...
    current_exp_id = 0
    for i in range(2):
        assert os.path.exists("test_by_detector_%03d.expt" % i)
        assert os.path.exists("test_by_detector_%03d.refl" % i)
        explist = ExperimentListFactory.from_json_file(
            "test_by_detector_%03d.expt" % i, check_format=False
        )
        refl = flex.reflection_table.from_pickle("test_by_detector_%03d.refl" % i)

        for k in range(len(explist)):
            assert refl.experiment_identifiers()[k] == str(current_exp_id)
            current_exp_id += 2

    assert not os.path.exists("test_by_detector_%03d.expt" % 2)
    assert not os.path.exists("test_by_detector_%03d.refl" % 2)


def test_failed_tolerance_error(dials_regression, monkeypatch):
    """Test that we get a sensible error message on tolerance failures"""
    # Select some experiments to use for combining
    jsons = os.path.join(
        dials_regression,
        "refinement_test_data",
        "multi_narrow_wedges",
        "data",
        "sweep_{:03d}",
        "{}",
    )
    files = [
        jsons.format(2, "experiments.json"),
        jsons.format(2, "reflections.pickle"),
        jsons.format(3, "experiments.json"),
        jsons.format(3, "reflections.pickle"),
    ]

    # Use the combine script
    script = combine_experiments.Script()
    # Disable writing output
    monkeypatch.setattr(script, "_save_output", lambda *args: None)
    # Parse arguments and configure
    params, options = script.parser.parse_args(files)
    params.reference_from_experiment.beam = 0

    # Validate that these pass
    script.run_with_preparsed(params, options)

    # Now, alter the beam to check it doesn't pass
    exp_2 = params.input.experiments[1].data[0]
    exp_2.beam.set_wavelength(exp_2.beam.get_wavelength() * 2)

    with pytest.raises(Sorry) as exc:
        script.run_with_preparsed(params, options)
    assert "Beam" in str(exc.value)
    print("Got (expected) error message:", exc.value)
