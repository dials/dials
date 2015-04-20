#!/usr/bin/env dials.python
from __future__ import division
import sys
from libtbx.phil import command_line, parse

if len(sys.argv) != 2: exit("please pass the path to a phil file")
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

from dials.model.serialize import load as load_dials
from dxtbx.serialize import load as load_dxtbx
def load_input(exp_path, ref_path):
  refs = load_dials.reflections(ref_path)
  exp = load_dxtbx.experiment_list(exp_path , check_format=False)[0]
  return refs, exp

from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment
class ExperimentFromCrystal(object):

  def __init__(self, reference_beam, reference_detector):

    self.reference_beam = reference_beam
    self.reference_detector = reference_detector
    return

  def __call__(self, crystal):

    return Experiment(beam=self.reference_beam,
                      detector=self.reference_detector,
                      crystal=crystal)

assert len(working_params.input) > 1
print len(working_params.input), "datasets specified as input"

e = enumerate(working_params.input)
i, line = e.next()
reflections, exp = load_input(line.experiments, line.reflections)
assert reflections['id'].all_eq(0)
experiment_from_crystal=ExperimentFromCrystal(exp.beam, exp.detector)

from dxtbx.model.experiment.experiment_list import ExperimentList
experiments=ExperimentList()
experiments.append(experiment_from_crystal(exp.crystal))

from scitbx.array_family import flex
for i, line in e:
  refs, exp = load_input(line.experiments, line.reflections)
  refs['id'] = flex.int(len(refs),i)
  reflections.extend(refs)
  experiments.append(experiment_from_crystal(exp.crystal))

# analysis of panel sampling
#TODO

# refinement - limit Jacobian calculation to 100'000 reflections at a time
from libtbx.phil import parse
user_phil=parse("""
refinement{
  parameterisation {
    beam.fix=all
    detector.hierarchy_level=1
    sparse=True
  }
  target.gradient_calculation_blocksize=100000
  reflections.outlier.algorithm=tukey
  reflections.minimum_number_of_reflections=1
}""")
from dials.data.refinement import phil_scope as master_phil
working_phil = master_phil.fetch(
  sources=[user_phil])
working_phil.show()
params = working_phil.extract()
from dials.algorithms.refinement import RefinerFactory
refiner = RefinerFactory.from_parameters_data_experiments(
    params, reflections, experiments,
    verbosity=2)
refiner.run()

# save the refined experiments
from dxtbx.model.experiment.experiment_list import ExperimentListDumper
experiments = refiner.get_experiments()
dump = ExperimentListDumper(experiments)
experiments_filename = "refined_experiments.json"
dump.as_json(experiments_filename)
print "refined geometry written to {0}".format(experiments_filename)

# save reflections used in refinement
matches = refiner.get_matches()
reflections_filename = "refined_reflections.json"
matches.as_pickle(reflections_filename)
print "reflections used in refinement written to {0}".format(reflections_filename)

