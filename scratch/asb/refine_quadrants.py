#!/usr/bin/env dials.python
from __future__ import division

import sys, os

from libtbx.phil import command_line, parse

if len(sys.argv) != 2: exit("please pass the path to a phil file")
#with(open(sys.argv[1])) as f:
#  phil = f.read()

#print phil
phil = sys.argv[1]

master_phil = parse("""
  input
    .multiple = True
  {
    experiments = None
      .type = path
    reflections = None
      .type = path
  }
  """)

cmd_line = command_line.argument_interpreter(master_params=master_phil)
working_phil = cmd_line.process_and_fetch(args=(phil,))
working_params = working_phil.extract()

for input in working_params.input:
  print input.experiments, input.reflections

from dials.model.serialize import load
def load_input(exp_path, ref_path):
  #from dials.model.experiment.experiment_list import ExperimentListFactory
  #e = ExperimentListFactory().from_args(working_params.input[0].experiments)
  refs = load.reflections(ref_path)
  exp = load.experiment_list(exp_path , check_format=False)[0]
  return refs, exp

from dials.model.experiment.experiment_list import ExperimentList, Experiment
class ExperimentFromCrystal(object):

  def __init__(self, reference_beam, reference_detector):

    self.reference_beam = reference_beam
    self.reference_detector = reference_detector
    return

  def __call__(self, crystal):

    return Experiment(beam=self.reference_beam,
                      detector=self.reference_detector,
                      crystal=crystal)

assert len(working_params.input)
e = enumerate(working_params.input)
i, line = e.next()
reflections, exp = load_input(line.experiments, line.reflections)
assert reflections['id'].all_eq(0)
experiment_from_crystal=ExperimentFromCrystal(exp.beam, exp.detector)

from dials.model.experiment.experiment_list import ExperimentList

experiments=ExperimentList()
experiments.append(experiment_from_crystal(exp.crystal))

from scitbx.array_family import flex
for i, line in e:
  print i

  refs, exp = load_input(line.experiments, line.reflections)
  refs['id'] = flex.int(len(refs),i)
  reflections.extend(refs)
  experiments.append(experiment_from_crystal(exp.crystal))



