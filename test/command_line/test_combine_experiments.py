"""
Test combination of multiple experiments and reflections files.
"""

from __future__ import absolute_import, division, print_function

import copy
import os
import procrunner

def test(dials_regression, tmpdir):
  from dials.array_family import flex
  from dxtbx.model.experiment_list import ExperimentListFactory

  tmpdir.chdir()

  data_dir = os.path.join(dials_regression, "refinement_test_data",
                          "multi_narrow_wedges")

  input_range = range(2, 49)
  for i in (8, 10, 15, 16, 34, 39, 45):
    input_range.remove(i)

  phil_input = "\n".join(
    ("  input.experiments={0}/data/sweep_%03d/experiments.json\n" +
     "  input.reflections={0}/data/sweep_%03d/reflections.pickle")
    % (i, i) for i in input_range)
# assert phil_input == "\n" + phil_input2 + "\n "

  input_phil = phil_input.format(data_dir) + """
 reference_from_experiment.beam=0
 reference_from_experiment.scan=0
 reference_from_experiment.goniometer=0
 reference_from_experiment.detector=0
 """

  with open("input.phil","w") as phil_file:
    phil_file.writelines(input_phil)

  result = procrunner.run_process(["dials.combine_experiments", "input.phil"])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  # load results
  exp = ExperimentListFactory.from_json_file("combined_experiments.json",
              check_format=False)
  ref = flex.reflection_table.from_pickle("combined_reflections.pickle")

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

  result = procrunner.run_process([
      "dials.split_experiments",
      "combined_experiments.json",
      "combined_reflections.pickle",
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  for i, e in enumerate(exp):
    assert os.path.exists("experiments_%03d.json" % i)
    assert os.path.exists("reflections_%03d.pickle" % i)

    exp_single = ExperimentListFactory.from_json_file(
      "experiments_%03d.json" % i, check_format=False)
    ref_single = flex.reflection_table.from_pickle("reflections_%03d.pickle" % i)

    assert len(exp_single) == 1
    assert exp_single[0].crystal == e.crystal
    assert exp_single[0].beam == e.beam
    assert exp_single[0].detector == e.detector
    assert exp_single[0].scan == e.scan
    assert exp_single[0].goniometer == e.goniometer
    assert exp_single[0].imageset == e.imageset
    assert len(ref_single) == len(ref.select(ref['id'] == i))
    assert ref_single['id'].all_eq(0)

  result = procrunner.run_process([
      "dials.split_experiments",
      "combined_experiments.json",
      "output.experiments_prefix=test",
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  for i in range(len(exp)):
    assert os.path.exists("test_%03d.json" %i)

  # Modify a copy of the detector
  detector = copy.deepcopy(exp.detectors()[0])
  panel = detector[0]
  x, y, z = panel.get_origin()
  panel.set_frame(panel.get_fast_axis(),
                  panel.get_slow_axis(),
                  (x, y, z+10))
  # Set half of the experiments to the new detector
  for i in xrange(len(exp)//2):
    exp[i].detector = detector
  from dxtbx.serialize import dump
  dump.experiment_list(exp, "modded_experiments.json")

  result = procrunner.run_process([
      "dials.split_experiments",
      "modded_experiments.json",
      "combined_reflections.pickle",
      "output.experiments_prefix=test_by_detector",
      "output.reflections_prefix=test_by_detector",
      "by_detector=True",
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  for i in range(2):
    assert os.path.exists("test_by_detector_%03d.json" % i)
    assert os.path.exists("test_by_detector_%03d.pickle" % i)
  assert not os.path.exists("test_by_detector_%03d.json" % 2)
  assert not os.path.exists("test_by_detector_%03d.pickle" % 2)
