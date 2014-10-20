#!/usr/bin/env dials.python
#
#
#  Copyright (C) 2014 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK, Lawrence Berkeley National Laboratory, USA
#
#  Author: David Waterman and Aaron Brewster
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""This script is intended for panel group refinement using still shot data
collected on a CSPAD detector. This version of the script uses a 'hybrid
minimiser'. Rather than a single joint refinement job of all crystals and the
detector, only the detector parameters are refined at first (using all data)
then each crystal is refined individually. This forms one macrocycle."""

from __future__ import division
from math import sqrt

from libtbx.phil import parse
from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment
from dials.algorithms.indexing.indexer import indexer_base

from dials.array_family import flex
from dials.algorithms.refinement import RefinerFactory

from libtbx.utils import Sorry
from libtbx import easy_mp

class ExperimentFromCrystal(object):

  def __init__(self, reference_beam, reference_detector):

    self.reference_beam = reference_beam
    self.reference_detector = reference_detector
    return

  def __call__(self, crystal):

    return Experiment(beam=self.reference_beam,
                      detector=self.reference_detector,
                      crystal=crystal)

def check_experiment(experiment, reflections):

  # predict reflections in place
  from dials.algorithms.spot_prediction import StillsReflectionPredictor
  sp = StillsReflectionPredictor(experiment)
  UB = experiment.crystal.get_U() * experiment.crystal.get_B()
  try:
    sp.for_reflection_table(reflections, UB)
  except RuntimeError:
    return False

  # calculate unweighted RMSDs
  x_obs, y_obs, _ = reflections['xyzobs.px.value'].parts()
  delpsi = reflections['delpsical.rad']
  x_calc, y_calc, _ = reflections['xyzcal.px'].parts()

  # calculate residuals and assign columns
  x_resid = x_calc - x_obs
  x_resid2 = x_resid**2
  y_resid = y_calc - y_obs
  y_resid2 = y_resid**2
  delpsical2 = delpsi**2
  r_x = flex.sum(x_resid2)
  r_y = flex.sum(y_resid2)
  r_z = flex.sum(delpsical2)

  # rmsd calculation
  n = len(reflections)
  rmsds = (sqrt(r_x / n),
           sqrt(r_y / n),
           sqrt(r_z / n))

  # check positional RMSDs are within 5 pixels
  if rmsds[0] > 5: return False
  if rmsds[1] > 5: return False

  return True

def detector_refiner(params, experiments, reflections):

  print "Refining detector at hierarchy_level=" + \
    str(params.refinement.parameterisation.detector.hierarchy_level)
  refiner = RefinerFactory.from_parameters_data_experiments(
        params, reflections, experiments)
  refiner.run()
  return refiner.get_experiments()

def crystals_refiner(params, experiments, reflections):

  def do_work(item):
    iexp, exp = item

    print "Refining crystal", iexp
    # reflection subset for a single experiment
    refs = reflections.select(reflections['id'] == iexp)
    refs['id'] = flex.size_t(len(refs),0)
    # experiment list for a single experiment
    exps=ExperimentList()
    exps.append(exp)
    refiner = RefinerFactory.from_parameters_data_experiments(
      params, refs, exps)
    # do refinement
    refiner.run()
    refined_exps = refiner.get_experiments()
    # replace this experiment with the refined one
    experiments[iexp] = refined_exps[0]

  print "Beginning crystal refinement with %d processor(s)"%params.mp.nproc
  easy_mp.parallel_map(
    func = do_work,
    iterable = enumerate(experiments),
    processes = params.mp.nproc,
    method = params.mp.method,
    asynchronous=True,
    preserve_exception_message=True)

  return experiments


class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''
      output {
        experiments_filename = refined_experiments.json
          .type = str
          .help = "The filename for refined experimental models"

        reflections_filename = None
          .type = str
          .help = "The filename for output of refined reflections"
      }

      n_macrocycles = 1
        .type = int(value_min=1)

      detector_phase {
        include scope dials.algorithms.refinement.refiner.phil_scope
      }

      crystals_phase {
        include scope dials.algorithms.refinement.refiner.phil_scope
        include scope dials.data.multiprocessing.phil_scope
      }

    ''', process_includes=True)

    # Set new defaults for detector and crystals refinement phases
    default_phil = parse('''
    crystals_phase.refinement {
        parameterisation {
          beam.fix=all
          detector.fix=all
        }
      reflections.do_outlier_rejection=True
      refinery.engine=LBFGScurvs
      verbosity=1
    }
    detector_phase.refinement {
        parameterisation {
          beam.fix=all
          crystal.fix=all
          detector.hierarchy_level=1
          sparse=True
        }
      target.jacobian_max_nref=100000
      reflections{
        do_outlier_rejection=True
        weighting_strategy.override=stills
        weighting_strategy.delpsi_constant=1000000
      }
      refinery.engine=LBFGScurvs
      verbosity=2
    }
    ''')

    # combine these
    working_phil = phil_scope.fetch(source=default_phil)

    # The script usage
    usage  = "usage: %s [options] [param.phil] " \
             "experiments.json reflections.pickle" \
               % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_reflections=True,
      read_experiments=True,
      check_format=False)

  def run(self):

    print "Parsing input"
    params, options = self.parser.parse_args(show_diff_phil=True)

    # Try to obtain the models and data
    if len(params.input.experiments) == 0:
      raise Sorry("No Experiments found in the input")
    if len(params.input.reflections) == 0:
      raise Sorry("No reflection data found in the input")
    try:
      assert len(params.input.reflections) == len(params.input.experiments)
    except AssertionError:
      raise Sorry("The number of input reflections files does not match the "
        "number of input experiments")

    # set up global experiments and reflections lists
    from dials.array_family import flex
    reflections = flex.reflection_table()
    global_id = 0
    from dxtbx.model.experiment.experiment_list import ExperimentList
    experiments=ExperimentList()

    # set the experiment factory that combines a crystal with the reference beam
    # and detector from the first experiment of the first experiment list
    ref_exp = params.input.experiments[0].data[0]
    experiment_from_crystal=ExperimentFromCrystal(ref_exp.beam, ref_exp.detector)

    # keep track of the number of refl per accepted experiment for a table
    nrefs_per_exp = []

    # loop through the input, building up the global lists
    for ref_wrapper, exp_wrapper in zip(params.input.reflections,
                                        params.input.experiments):
      refs = ref_wrapper.data
      exps = exp_wrapper.data

      # there might be multiple experiments already here. Loop through them
      for i, exp in enumerate(exps):

        # select the relevant reflections
        sel = refs['id'] == i
        sub_ref = refs.select(sel)

        if len(sub_ref) == 0:
          print "skipping experiment", i, "in", exp_wrapper.filename, "due to no reflections being found in", refs.filename
          continue

        # build an experiment with this crystal plus the reference models
        combined_exp = experiment_from_crystal(exp.crystal)

        # next experiment ID in series
        exp_id = len(experiments)

        # check this experiment
        if not check_experiment(combined_exp, sub_ref):
          print "skipping experiment", i, "in", exp_wrapper.filename, "due to poor RMSDs"
          continue

        # set reflections ID
        sub_ref['id'] = flex.size_t(len(sub_ref), exp_id)

        # keep number of reflections for the table
        nrefs_per_exp.append(len(sub_ref))

        # obtain mm positions on the reference detector
        sub_ref = indexer_base.map_spots_pixel_to_mm_rad(sub_ref,
          combined_exp.detector, combined_exp.scan)

        # extend refl and experiments lists
        reflections.extend(sub_ref)
        experiments.append(combined_exp)

    # print number of reflections per accepted experiment
    from libtbx.table_utils import simple_table
    header = ["Experiment", "Nref"]
    rows = [(str(i), str(n)) for (i, n) in enumerate(nrefs_per_exp)]
    st = simple_table(rows, header)
    print st.format()

    for cycle in range(params.n_macrocycles):

      print "MACROCYCLE %02d" % (cycle + 1)
      print "=============\n"
      # first run: multi experiment joint refinement of detector with fixed beam and
      # crystals
      print "PHASE 1"
      experiments = detector_refiner(params.detector_phase, experiments, reflections)

      # second run
      print "PHASE 2"
      experiments = crystals_refiner(params.crystals_phase, experiments, reflections)

    # Save the refined experiments to file
    output_experiments_filename = params.output.experiments_filename
    print 'Saving refined experiments to {0}'.format(output_experiments_filename)
    from dxtbx.model.experiment.experiment_list import ExperimentListDumper
    dump = ExperimentListDumper(experiments)
    dump.as_json(output_experiments_filename)

    # Write out refined reflections, if requested
    if params.output.reflections_filename:
      print 'Saving refined reflections to {0}'.format(
        params.output.reflections_filename)
      reflections.as_pickle(params.output.reflections_filename)

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
