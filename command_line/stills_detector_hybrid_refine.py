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

# LIBTBX_SET_DISPATCHER_NAME dev.dials.stills_detector_hybrid_refine

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
from dials.algorithms.refinement.stills_detector_metrology import \
  StillsDetectorRefinerFactory
from dials.algorithms.refinement.refinement_helpers import \
  get_panel_groups_at_depth, get_panel_ids_at_root
from dials.util import log

from libtbx.utils import Sorry
from libtbx import easy_mp
import copy

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

  # Here use the specialised faster refiner
  refiner = StillsDetectorRefinerFactory.from_parameters_data_experiments(
        params, reflections, experiments)
  refiner.run()
  return refiner.get_experiments()

def detector_parallel_refiners(params, experiments, reflections):

  print "Refining detector at hierarchy_level=" + \
    str(params.refinement.parameterisation.detector.hierarchy_level), "\n"
  orig_detector = experiments.detectors()[0]
  try:
    h = orig_detector.hierarchy()
  except AttributeError:
    print "This detector does not have a hierarchy"
    raise

  # get the panel groups at the chosen level
  level = params.refinement.parameterisation.detector.hierarchy_level
  try:
    groups = get_panel_groups_at_depth(h, level)
  except AttributeError:
    print "Cannot access the hierarchy at the depth level={0}".format(level)
    raise

  # collect the panel ids for each Panel within the groups
  panels = [p for p in orig_detector]
  panel_ids_by_group = [get_panel_ids_at_root(panels, g) for g in groups]

  print "The detector will be divided into", len(panel_ids_by_group), \
    "groups consisting of the following panels:"
  for i, g in enumerate(panel_ids_by_group):
    print "Group%02d:" % (i+1), g
  print

  # now construct sub-detectors
  def recursive_add_child(d, parent, child):
    """ Creates either a panel group or a panel on the parent,
        and sets it up to match the child """
    if child.is_group():
      newchild = parent.add_group()
    else:
      newchild = parent.add_panel()
      newchild.set_image_size(child.get_image_size())
      newchild.set_trusted_range(child.get_trusted_range())
      newchild.set_pixel_size(child.get_pixel_size())
      newchild.set_px_mm_strategy(child.get_px_mm_strategy())

    m = child.get_local_d_matrix()
    newchild.set_local_frame(m[0::3],m[1::3],m[2::3])
    newchild.set_name(child.get_name())
    if child.is_group():
      for c in child.children():
        recursive_add_child(d, newchild, c)

  from dxtbx.model import Detector
  sub_detectors = [Detector() for e in groups]
  for d, g in zip(sub_detectors, groups):
    d.hierarchy().set_name(g.get_name())
    d.hierarchy().set_frame(g.get_fast_axis(),
                            g.get_slow_axis(),
                            g.get_origin())
    if g.is_group():
      for c in g.children():
        recursive_add_child(d, d.hierarchy(), c)
    else: # at the bottom of the hierarchy. Note the new panel's frame will be the identity matrix.
      p = d.hierarchy().add_panel()
      p.set_image_size(g.get_image_size())
      p.set_trusted_range(g.get_trusted_range())
      p.set_pixel_size(g.get_pixel_size())
      p.set_px_mm_strategy(g.get_px_mm_strategy())
      p.set_name(g.get_name())

  # set experiment lists for each sub-detector
  sub_det_expts = [copy.deepcopy(experiments) for e in groups]
  for d, exp in zip(sub_detectors, sub_det_expts):
    exp.replace(exp.detectors()[0], d)

  # divide the reflections by sub-detector
  sub_reflections = []
  for pnls in panel_ids_by_group:
    isels = [(reflections['panel'] == pnl).iselection() for pnl in pnls]
    isel = flex.size_t()
    for s in isels: isel.extend(s)
    gp_refs = reflections.select(isel)
    # reset panel number to match the sub-detector
    for new_id, old_id in enumerate(pnls):
      sel = gp_refs['panel'] == old_id
      gp_refs['panel'].set_selected(sel, new_id)
    sub_reflections.append(gp_refs)

  # We wish to refine each whole sub-detector as a single group. Therefore
  # we must use hierarchy_level=0 for these jobs
  tmplevel = params.refinement.parameterisation.detector.hierarchy_level
  params.refinement.parameterisation.detector.hierarchy_level=0

  # do refinements and collect the refined experiments
  def do_work(item):
    refs, exps = item

    if len(refs) < 20:
      print "Cannot refine detector", exps[0].detector.hierarchy().get_name(), "due to too few reflections (", len(refs), ")"
      return exps # do not refine this detector element

    # Here use the specialised faster refiner
    refiner = StillsDetectorRefinerFactory.from_parameters_data_experiments(
        params, refs, exps)
    refiner.run()
    return refiner.get_experiments()

  refined_exps = easy_mp.parallel_map(
    func = do_work,
    iterable = zip(sub_reflections, sub_det_expts),
    processes = params.mp.nproc,
    method = params.mp.method,
    asynchronous=True,
    preserve_exception_message=True)

  # update the full detector
  for group, refined_exp in zip(groups, refined_exps):
    refined_det = refined_exp.detectors()[0]
    local_root = refined_det[0]
    f = local_root.get_fast_axis()
    s = local_root.get_slow_axis()
    o = local_root.get_origin()
    group.set_frame(f, s, o) # propagates local frame changes

  # refine the full detector to get RMSDs per panel
  print
  print "Refining full recombined detector"
  print "---------------------------------"
  experiments = detector_refiner(params, experiments, reflections)

  # reset hierarchy_level
  params.refinement.parameterisation.detector.hierarchy_level=tmplevel

  return experiments

def crystals_refiner(params, experiments, reflections):

  def do_work(item):
    iexp, exp = item

    print "Refining crystal", iexp
    # reflection subset for a single experiment
    refs = reflections.select(reflections['id'] == iexp)
    refs['id'] = flex.int(len(refs),0)

    # DGW commented out as reflections.minimum_number_of_reflections no longer exists
    #if len(refs) < params.refinement.reflections.minimum_number_of_reflections:
    #  print "Not enough reflections to refine experiment"
    #  return

    # experiment list for a single experiment
    exps=ExperimentList()
    exps.append(exp)
    try:
      refiner = RefinerFactory.from_parameters_data_experiments(
        params, refs, exps)
      # do refinement
      refiner.run()
    except Exception, e:
      print "Error,", str(e)
      return

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
          .help = The filename for refined experimental models

        reflections_filename = None
          .type = str
          .help = The filename for output of refined reflections
      }

      n_macrocycles = 1
        .type = int(value_min=1)

      detector_phase {
        include scope dials.algorithms.refinement.refiner.phil_scope
        include scope dials.data.multiprocessing.phil_scope
      }

      crystals_phase {
        include scope dials.algorithms.refinement.refiner.phil_scope
        include scope dials.data.multiprocessing.phil_scope
      }

      reference_detector = *first average
        .type = choice
        .help = First: use the first detector found in the experiment \
                Average: create an average detector from all experiments

    ''', process_includes=True)

    # Set new defaults for detector and crystals refinement phases
    default_phil = parse('''
    crystals_phase.refinement {
        parameterisation {
          beam.fix=all
          detector.fix=all
        }
      reflections.outlier.algorithm=null
      refinery.engine=LevMar
      verbosity=1
    }
    detector_phase.refinement {
        parameterisation {
          beam.fix=all
          crystal.fix=all
          detector.hierarchy_level=1
          sparse=True
        }
      target.gradient_calculation_blocksize=100000
      reflections{
        outlier.algorithm=tukey
        outlier.separate_experiments=False
        weighting_strategy.override=stills
        weighting_strategy.delpsi_constant=1000000
      }
      refinery.engine=LevMar
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
      phil=working_phil,
      read_reflections=True,
      read_experiments=True,
      check_format=False)

  def run(self):

    print "Parsing input"
    params, options = self.parser.parse_args(show_diff_phil=True)

    #Configure the logging
    log.config(params.detector_phase.refinement.verbosity,
      info='dials.refine.log', debug='dials.refine.debug.log')

    # Try to obtain the models and data
    if not params.input.experiments:
      raise Sorry("No Experiments found in the input")
    if not params.input.reflections:
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

    if params.reference_detector == "first":
      # Use the first experiment of the first experiment list as the reference detector
      ref_exp = params.input.experiments[0].data[0]
    else:
      # Average all the detectors to generate a reference detector
      assert params.detector_phase.refinement.parameterisation.detector.hierarchy_level == 0
      from scitbx.matrix import col
      panel_fasts = []
      panel_slows = []
      panel_oris = []
      for exp_wrapper in params.input.experiments:
        exp = exp_wrapper.data[0]
        if panel_oris:
          for i, panel in enumerate(exp.detector):
            panel_fasts[i] += col(panel.get_fast_axis())
            panel_slows[i] += col(panel.get_slow_axis())
            panel_oris[i] += col(panel.get_origin())
        else:
          for i, panel in enumerate(exp.detector):
            panel_fasts.append(col(panel.get_fast_axis()))
            panel_slows.append(col(panel.get_slow_axis()))
            panel_oris.append(col(panel.get_origin()))

      ref_exp = copy.deepcopy(params.input.experiments[0].data[0])
      for i, panel in enumerate(ref_exp.detector):
        # Averaging the fast and slow axes can make them be non-orthagonal. Fix by finding
        # the vector that goes exactly between them and rotate
        # around their cross product 45 degrees from that vector in either direction
        vf = panel_fasts[i]/len(params.input.experiments)
        vs = panel_slows[i]/len(params.input.experiments)
        c = vf.cross(vs)
        angle = vf.angle(vs, deg=True)
        v45 = vf.rotate(c, angle/2, deg=True)
        vf = v45.rotate(c, -45, deg=True)
        vs = v45.rotate(c, 45, deg=True)
        panel.set_frame(vf, vs,
                        panel_oris[i]/len(params.input.experiments))

      print "Reference detector (averaged):", str(ref_exp.detector)

    # set the experiment factory that combines a crystal with the reference beam
    # and the reference detector
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

        ## DGW commented out as reflections.minimum_number_of_reflections no longer exists
        #if len(sub_ref) < params.crystals_phase.refinement.reflections.minimum_number_of_reflections:
        #  print "skipping experiment", i, "in", exp_wrapper.filename, "due to insufficient strong reflections in", ref_wrapper.filename
        #  continue

        # build an experiment with this crystal plus the reference models
        combined_exp = experiment_from_crystal(exp.crystal)

        # next experiment ID in series
        exp_id = len(experiments)

        # check this experiment
        if not check_experiment(combined_exp, sub_ref):
          print "skipping experiment", i, "in", exp_wrapper.filename, "due to poor RMSDs"
          continue

        # set reflections ID
        sub_ref['id'] = flex.int(len(sub_ref), exp_id)

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
    print "Number of reflections per experiment"
    print st.format()

    for cycle in range(params.n_macrocycles):

      print "MACROCYCLE %02d" % (cycle + 1)
      print "=============\n"
      # first run: multi experiment joint refinement of detector with fixed beam and
      # crystals
      print "PHASE 1"

      # SET THIS TEST TO FALSE TO REFINE WHOLE DETECTOR AS SINGLE JOB
      if params.detector_phase.refinement.parameterisation.detector.hierarchy_level > 0:
        experiments = detector_parallel_refiners(params.detector_phase, experiments, reflections)
      else:
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
