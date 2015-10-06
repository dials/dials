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

def generate_reflections(experiment):

  # copy the experiment to avoid modifying the original models
  from copy import deepcopy
  experiment = deepcopy(experiment)

  detector = experiment.detector
  beam = experiment.beam
  goniometer = experiment.goniometer
  crystal = experiment.crystal
  scan = experiment.scan
  assert [goniometer, scan].count(None) == 0

  from dxtbx.model.experiment.experiment_list import \
    ExperimentList
  experiments = ExperimentList()
  experiments.append(experiment)

  ###########################
  # Parameterise the models #
  ###########################

  det_param = DetectorParameterisationSinglePanel(detector)
  s0_param = BeamParameterisation(beam, goniometer)
  xlo_param = CrystalOrientationParameterisation(crystal)
  xluc_param = CrystalUnitCellParameterisation(crystal)

  # Fix beam to the X-Z plane (imgCIF geometry), fix wavelength
  s0_param.set_fixed([True, False, True])

  ################################
  # Apply known parameter shifts #
  ################################

  # shift detector by 1.0 mm each translation and 2 mrad each rotation
  det_p_vals = det_param.get_param_vals()
  p_vals = [a + b for a, b in zip(det_p_vals,
                                  [1.0, 1.0, 1.0, 2., 2., 2.])]
  det_param.set_param_vals(p_vals)

  # shift beam by 2 mrad in free axis
  s0_p_vals = s0_param.get_param_vals()
  p_vals = list(s0_p_vals)

  p_vals[0] += 2.
  s0_param.set_param_vals(p_vals)

  # rotate crystal a bit (=10 mrad each rotation)
  xlo_p_vals = xlo_param.get_param_vals()
  p_vals = [a + b for a, b in zip(xlo_p_vals, [10., 10., 10.])]
  xlo_param.set_param_vals(p_vals)

  # work out what unit cell parameters we can vary by apparent lattice symmetry
  xluc_p_vals = xluc_param.get_param_vals()
  (a,b,c,aa,bb,cc) = crystal.get_unit_cell().parameters()
  a += 0.1
  b += 0.1
  c += 0.1
  if not_approx_equal(aa, 90) or not_approx_equal(aa, 120): aa += 0.1
  if not_approx_equal(bb, 90) or not_approx_equal(bb, 120): bb += 0.1
  if not_approx_equal(cc, 90) or not_approx_equal(cc, 120): cc += 0.1
  new_uc = unit_cell((a,b,c,aa,bb,cc))
  newB = matrix.sqr(new_uc.fractionalization_matrix()).transpose()
  S = symmetrize_reduce_enlarge(crystal.get_space_group())
  S.set_orientation(orientation=newB)
  X = tuple([e * 1.e5 for e in S.forward_independent_parameters()])
  xluc_param.set_param_vals(X)

  #############################
  # Generate some reflections #
  #############################

  # All indices in a 1.0 Angstrom sphere
  resolution = 1.5
  index_generator = IndexGenerator(crystal.get_unit_cell(),
                  space_group(space_group_symbols(1).hall()).type(), resolution)
  indices = index_generator.to_array()

  # Predict rays within the sweep range
  ray_predictor = ScansRayPredictor(experiments,
    scan.get_oscillation_range(deg=False))
  obs_refs = ray_predictor.predict(indices)

  print "Total number of reflections excited", len(obs_refs)

  # Take only those rays that intersect the detector
  intersects = ray_intersection(detector, obs_refs)
  obs_refs = obs_refs.select(intersects)

  # Make a reflection predictor and re-predict for all these reflections. The
  # result is the same, but we gain also the flags and xyzcal.px columns
  ref_predictor = ExperimentsPredictor(experiments)
  obs_refs['id'] = flex.size_t(len(obs_refs), 0)
  obs_refs = ref_predictor.predict(obs_refs)

  # Set 'observed' centroids from the predicted ones
  obs_refs['xyzobs.mm.value'] = obs_refs['xyzcal.mm']

  # Invent some variances for the centroid positions of the simulated data
  im_width = scan.get_oscillation()[1] * pi / 180.
  px_size = detector[0].get_pixel_size()
  var_x = flex.double(len(obs_refs), (px_size[0] / 2.)**2)
  var_y = flex.double(len(obs_refs), (px_size[1] / 2.)**2)
  var_phi = flex.double(len(obs_refs), (im_width / 2.)**2)
  obs_refs['xyzobs.mm.variance'] = flex.vec3_double(var_x, var_y, var_phi)

  print "Total number of observations made", len(obs_refs)

  return obs_refs, experiment


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
  print params.input.experiments
  from dials.util.options import flatten_experiments
  experiments = flatten_experiments(params.input.experiments)

  try:
    # only allow single experiment at the moment
    assert len(experiments) == 1
  except AssertionError:
    print usage
    raise

  from dials.util import log
  log.config(params.refinement.verbosity,
      info='synthetic_refine.log', debug='synthetic_refine.debug.log')

  reflections, perturbed_exp = generate_reflections(experiments[0])

  from dials.algorithms.refinement import RefinerFactory
  refiner = RefinerFactory.from_parameters_data_experiments(
    params, reflections, experiments)
  refiner.run()

  refined_experiments = refiner.get_experiments()

  # quick check on refined detector geometry using panel 0
  old_detector = perturbed_exp.detector
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
