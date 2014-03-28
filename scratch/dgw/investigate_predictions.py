#!/usr/bin/env python
#
# dials.refine.py
#
#  Copyright (C) 2013 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK.
#
#  Author: James Parkhurst and David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner
from dials.util.command_line import Importer
from dials.model.experiment.experiment_list import \
  ExperimentList, ExperimentListDumper
from math import pi
from scitbx import matrix

from dials.algorithms.spot_prediction import ray_intersection

TWO_PI = 2.0*pi

class Script(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = "usage: %prog [options] [param.phil] " \
             "experiments.json reflections.pickle"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Output experiments filename option
    self.config().add_option(
        '--output-experiments-filename',
        dest = 'output_experiments_filename',
        type = 'string', default = 'refined_experiments.json',
        help = 'Set the filename for refined experimental models.')

    self.config().add_option(
        '--output-centroids-filename',
        dest = 'output_centroids_filename',
        type = 'string',
        help = 'Set the filename for the table of centroids at the end of refinement.')

    self.config().add_option(
        '--output-parameters-filename',
        dest = 'output_parameters_filename',
        type = 'string',
        help = 'Set the filename for the table of scan varying parameter values'
               ' at the end of refinement.')

    # Add a verbosity option
    self.config().add_option(
        "-v", "--verbosity",
        action="count", default=0,
        help="set verbosity level; -vv gives verbosity level 2.")


  def main(self, params, options, args):
    '''Execute the script.'''
    #from dials.algorithms.refinement import RefinerFactory
    import cPickle as pickle

    # Check the number of arguments
    if len(args) < 2:
      self.config().print_help()
      return

    importer = Importer(args, check_format=False, verbose=False)

    # Try to load the models and data
    experiments = importer.experiments
    if experiments is None:
      raise RuntimeError("No Experiments found in the input")
    reflections = importer.reflections
    if reflections is None:
      raise RuntimeError("No reflection data found in the input")
    if len(reflections) > 1:
      raise RuntimeError("Only one reflections list can be imported at present")
    reflections = importer.reflections[0]

    #from dials.util.command_line import interactive_console; interactive_console()
    # check that the beam vectors are stored: if not, compute them
    nrefs_wo_s1 = (reflections['s1'].norms() < 1.e-6).count(True)
    if nrefs_wo_s1 > 0:
      print "Setting scattering vectors for", nrefs_wo_s1, "reflections"
    for i in xrange(reflections.nrows()):
      ref = reflections[i]
      if ref['s1'] != (0.0, 0.0, 0.0):
        continue
      beam = experiments[ref['id']].beam
      detector = experiments[ref['id']].detector
      panel = detector[ref['panel']]
      impact = ref['xyzobs.mm.value'][0:2]
      x, y = panel.millimeter_to_pixel(impact)
      s1 = matrix.col(panel.get_pixel_lab_coord(
          (x, y))).normalize() / beam.get_wavelength()
      row = { 's1' : s1 }
      reflections[i] = row

    print "how many reflections?", len(reflections)
    from copy import deepcopy
    old_reflections = deepcopy(reflections)

    # make ReflectionManagers
    from dials.algorithms.refinement.reflection_manager import ReflectionManager
    from dials.algorithms.refinement.target import ReflectionManager as OldReflectionManager
    print "Build new ReflectionManager"
    new_refman = ReflectionManager(reflections, experiments, verbosity=2, iqr_multiplier=None)
    print "Build old ReflectionManager"
    old_refman = OldReflectionManager(reflections, experiments, verbosity=2, iqr_multiplier=None)

    # make reflection predictors
    from dials.algorithms.refinement.prediction import ExperimentsPredictor
    from dials.algorithms.refinement.prediction import ScansRayPredictor
    new_ref_predictor = ExperimentsPredictor(experiments)
    sweep_range = experiments[0].scan.get_oscillation_range(deg=False)
    #old_ref_predictor = ScansRayPredictor(experiments, sweep_range)
    old_ref_predictor = ScansRayPredictor(experiments)

    # make a new target object
    # pass None as a prediction_parameterisation - this doesn't matter because
    # we don't use it to predict reflections
    from dials.algorithms.refinement.target_new import LeastSquaresPositionalResidualWithRmsdCutoff
    new_target = LeastSquaresPositionalResidualWithRmsdCutoff(experiments, new_ref_predictor, new_refman,
                   prediction_parameterisation=None)

    print
    print "Compare sizes: input, accepted and sample"
    print "New refman:", (new_refman.get_input_size(),
                          new_refman.get_accepted_refs_size(),
                          new_refman.get_sample_size())
    print "Old refman:", (old_refman.get_input_size(),
                          old_refman.get_accepted_refs_size(),
                          old_refman.get_sample_size())

    # compare predictions
    print "OLD TARGET PREDICT"
    self.old_target_predict(old_ref_predictor, old_refman, experiments)
    print "NEW TARGET PREDICT"
    new_target.predict()

    print "number of matches"
    print "new = ", len(new_refman.get_matches())
    print "old = ", len(old_refman.get_matches())

    return

  @staticmethod
  def old_target_predict(reflection_predictor, reflection_manager, experiments):
    """standalone method copied out of old Target class to do prediction
    as loop over hkl and update the reflection manager"""

    # update the reflection_predictor and the prediction parameterisation
    # with the scan-independent part of the current geometry
    reflection_predictor.update()
    #self._prediction_parameterisation.prepare()

    # reset the 'use' flag for all observations
    reflection_manager.reset_accepted_reflections()

    # loop over all reflections in the manager
    for obs in reflection_manager.get_obs():

      # get data from the observation
      h = obs.miller_index
      frame_id = obs.frame_obs
      panel_id = obs.panel

      # FIXME current reflection container does not track experiment id. Here
      # we use the crystal_id as its proxy. User code must ensure that the order
      # of experiments matches the index given by the crystal_id. Unless
      # 'crystal_id' is abused, this effectively means that only multiple
      # crystal experiments are supported, not a multiplicity of other models.
      experiment_id = obs.crystal_id


      # predict for this hkl
      predictions = reflection_predictor.predict(
                             h, experiment_id=experiment_id)

      # obtain the impact positions
      impacts = ray_intersection(experiments[experiment_id].detector,
                                 predictions,
                                 panel=panel_id)

      print h, len(impacts)
      # find the prediction with the right 'entering' flag
      try:
        i = [x.entering == obs.entering \
             for x in impacts].index(True)
      except ValueError:
        # we don't have a prediction for this obs
        print "NO PREDICTION FOR", h
        continue

      ref = impacts[i]
      x_calc, y_calc = ref.image_coord_mm

      # do not wrap around multiples of 2*pi; keep the full rotation
      # from zero to differentiate repeat observations.
      resid = ref.rotation_angle - (obs.phi_obs % TWO_PI)

      # ensure this is the smaller of two possibilities
      resid = (resid + pi) % TWO_PI - pi

      phi_calc = obs.phi_obs + resid
      s_calc = matrix.col(ref.beam_vector)

      # calculate gradients for this reflection
      #grads = self._prediction_parameterisation.get_gradients(
      #                            h, s_calc, phi_calc, panel_id, frame_id,
      #                            experiment_id=experiment_id)

      # store all this information in the matched obs-pred pair
      obs.update_prediction(x_calc, y_calc, phi_calc, s_calc, gradients=None)

    if reflection_manager.first_update:

      # print summary before outlier rejection
      reflection_manager.print_stats_on_matches()

      # flag potential outliers
      rejection_occurred = reflection_manager.reject_outliers()

      # delete all obs-pred pairs from the manager that do not
      # have a prediction or were flagged as outliers
      reflection_manager.strip_unmatched_observations()

      # print summary after outlier rejection
      if rejection_occurred: reflection_manager.print_stats_on_matches()

      reflection_manager.first_update = False

    return

if __name__ == '__main__':
  script = Script()
  script.run()
