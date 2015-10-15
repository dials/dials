#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
A test script to run refinement against a synthetic rotation series from some
input experiments.json

"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi
from libtbx.phil import parse
from libtbx.test_utils import not_approx_equal
from scitbx import matrix
from scitbx.array_family import flex
from logging import info, debug
from copy import deepcopy
import random

from dials.algorithms.refinement import RefinerFactory

# Model parameterisations
from dials.algorithms.refinement.parameterisation.detector_parameters import \
    DetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.beam_parameters import \
    BeamParameterisation
from dials.algorithms.refinement.parameterisation.crystal_parameters import \
    CrystalOrientationParameterisation, CrystalUnitCellParameterisation

# Symmetry constrained parameterisation for the unit cell
from cctbx.uctbx import unit_cell
from rstbx.symmetry.constraints.parameter_reduction import \
    symmetrize_reduce_enlarge

# Reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dials.algorithms.refinement.prediction import ScansRayPredictor, \
  ExperimentsPredictor
from dials.algorithms.spot_prediction import ray_intersection
from cctbx.sgtbx import space_group, space_group_symbols

class ExperimentsPerturber(object):
  '''Perturb the models in an experiment list. For simplicity create a complete
  Refiner for this, which will just be used as a carrier for the parameterised
  models'''

  def __init__(self, experiments):

    self.dummy_reflections = generate_reflections(experiments)

    # use default options for Refiners built internally
    phil_scope = parse('''
    include scope dials.algorithms.refinement.refiner.phil_scope
    ''', process_includes=True)
    self.params = phil_scope.extract()

    # changes for building internal Refiners
    self.params.refinement.reflections.outlier.algorithm="null"
    self.params.refinement.reflections.random_seed=None
    self.params.refinement.verbosity=0

    self.original_experiments = experiments
    return

  def random_perturbation(self, fraction=0.05):
    '''randomly perturb each model parameter value by an amount drawn from a
    normal distribution with sigma equal to fraction*value'''

    refiner = RefinerFactory.from_parameters_data_experiments(
      self.params, self.dummy_reflections, self.original_experiments)

    pp = refiner._pred_param
    vals = pp.get_param_vals()
    sigmas = [fraction * val for val in vals]
    new_vals = [random.gauss(mu, sig) for (mu, sig) in zip(vals, sigmas)]
    pp.set_param_vals(new_vals)

    return refiner.get_experiments()

  def known_perturbation(self, fraction=0.05):
    '''Perturb each model parameter value by a known relative amount given by
    fraction*value'''

    refiner = RefinerFactory.from_parameters_data_experiments(
      self.params, self.dummy_reflections, self.original_experiments)

    pp = refiner._pred_param
    vals = pp.get_param_vals()
    shifts = [fraction * val for val in vals]
    new_vals = [val + shift for (val, shift) in zip(vals, shifts)]
    pp.set_param_vals(new_vals)

    return refiner.get_experiments()

def generate_reflections(experiments):

  refs = []
  for iexp, exp in enumerate(experiments):

    info("Generating reflections for experiment {0}".format(iexp))

    # All indices in a 1.5 Angstrom sphere
    resolution = 1.5
    index_generator = IndexGenerator(exp.crystal.get_unit_cell(),
                    space_group(space_group_symbols(1).hall()).type(), resolution)
    indices = index_generator.to_array()

    # Predict rays within the sweep range
    ray_predictor = ScansRayPredictor(experiments,
      exp.scan.get_oscillation_range(deg=False))
    obs_refs = ray_predictor.predict(indices, experiment_id=iexp)

    info("Total number of reflections excited: {0}".format(len(obs_refs)))

    # Take only those rays that intersect the detector
    intersects = ray_intersection(exp.detector, obs_refs)
    obs_refs = obs_refs.select(intersects)
    obs_refs['id'] = flex.size_t(len(obs_refs), iexp)
    refs.append(obs_refs)

    info("Total number of impacts: {0}".format(len(obs_refs)))

  # Concatenate reflections
  obs_refs = reduce(lambda x, y: x.extend(y), refs)

  # Make a reflection predictor and re-predict for all these reflections. The
  # result is the same, but we gain also the flags and xyzcal.px columns
  ref_predictor = ExperimentsPredictor(experiments)
  obs_refs = ref_predictor.predict(obs_refs)

  # Set 'observed' centroids from the predicted ones
  obs_refs['xyzobs.mm.value'] = obs_refs['xyzcal.mm']

  # Invent some variances for the centroid positions of the simulated data
  im_width = exp.scan.get_oscillation()[1] * pi / 180.
  px_size = exp.detector[0].get_pixel_size()
  var_x = flex.double(len(obs_refs), (px_size[0] / 2.)**2)
  var_y = flex.double(len(obs_refs), (px_size[1] / 2.)**2)
  var_phi = flex.double(len(obs_refs), (im_width / 2.)**2)
  obs_refs['xyzobs.mm.variance'] = flex.vec3_double(var_x, var_y, var_phi)

  info("Total number of observations made: {0}".format(len(obs_refs)))

  return obs_refs


if __name__ == "__main__":

  from dials.util.options import OptionParser
  # The phil scope
  phil_scope = parse('''
  include scope dials.algorithms.refinement.refiner.phil_scope
  ''', process_includes=True)

  # Initialise the option parser
  usage = "usage: %prog experiments.json [params]"
  parser = OptionParser(
    phil=phil_scope,
    usage=usage,
    read_experiments=True,
    check_format=False)

  # Parse the command line arguments
  params, options = parser.parse_args(show_diff_phil=True)
  #print params.input.experiments
  from dials.util.options import flatten_experiments
  experiments = flatten_experiments(params.input.experiments)

  from dials.util import log
  log.config(params.refinement.verbosity,
      info='synthetic_refine.log', debug='synthetic_refine.debug.log')

  info("Perturbing original experiments")
  exp_perturber = ExperimentsPerturber(experiments)
  perturbed_experiments = exp_perturber.random_perturbation()

  info("Generating 'observations' to refine against")
  reflections = generate_reflections(perturbed_experiments)

  info("Running refinement starting from original experiments")
  refiner = RefinerFactory.from_parameters_data_experiments(
    params, reflections, experiments)
  refiner.run()
  refined_experiments = refiner.get_experiments()

  # quick check on refined detector geometry using panel 0
  old_detector = experiments[0].detector
  new_detector = refined_experiments[0].detector
  old_origin = matrix.col(old_detector[0].get_origin())
  new_origin = matrix.col(new_detector[0].get_origin())
  dorigin = new_origin - old_origin
  print "origin offset is", dorigin.length(), "mm"

  old_fast = matrix.col(old_detector[0].get_fast_axis())
  old_slow = matrix.col(old_detector[0].get_slow_axis())
  new_fast = matrix.col(new_detector[0].get_fast_axis())
  new_slow = matrix.col(new_detector[0].get_slow_axis())

  print "offset angle between fast axes is", old_fast.accute_angle(
    new_fast, deg=True), "degrees"
  print "offset angle between slow axes is", old_slow.accute_angle(
    new_slow, deg=True), "degrees"
