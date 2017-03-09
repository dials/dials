#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Test hierarchical detector refinement.

"""

# python imports
from __future__ import absolute_import, division
import os
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx.test_utils import approx_equal
from math import pi
#from libtbx.test_utils import open_tmp_directory

def generate_reflections(experiments):

  from dials.algorithms.spot_prediction import IndexGenerator
  from dials.algorithms.refinement.prediction import \
    ScansRayPredictor, ExperimentsPredictor
  from dials.algorithms.spot_prediction import ray_intersection
  from cctbx.sgtbx import space_group, space_group_symbols
  from scitbx.array_family import flex

  detector = experiments[0].detector
  crystal = experiments[0].crystal

  # All indices in a 2.0 Angstrom sphere
  resolution = 2.0
  index_generator = IndexGenerator(crystal.get_unit_cell(),
                  space_group(space_group_symbols(1).hall()).type(), resolution)
  indices = index_generator.to_array()

  # Predict rays within the sweep range
  scan = experiments[0].scan
  sweep_range = scan.get_oscillation_range(deg=False)
  ray_predictor = ScansRayPredictor(experiments, sweep_range)
  obs_refs = ray_predictor(indices)

  # Take only those rays that intersect the detector
  intersects = ray_intersection(detector, obs_refs)
  obs_refs = obs_refs.select(intersects)

  # Make a reflection predictor and re-predict for all these reflections. The
  # result is the same, but we gain also the flags and xyzcal.px columns
  ref_predictor = ExperimentsPredictor(experiments)
  obs_refs['id'] = flex.int(len(obs_refs), 0)
  obs_refs = ref_predictor(obs_refs)

  # Set 'observed' centroids from the predicted ones
  obs_refs['xyzobs.mm.value'] = obs_refs['xyzcal.mm']

  # Invent some variances for the centroid positions of the simulated data
  im_width = 0.1 * pi / 180.
  px_size = detector[0].get_pixel_size()
  var_x = flex.double(len(obs_refs), (px_size[0] / 2.)**2)
  var_y = flex.double(len(obs_refs), (px_size[1] / 2.)**2)
  var_phi = flex.double(len(obs_refs), (im_width / 2.)**2)
  obs_refs['xyzobs.mm.variance'] = flex.vec3_double(var_x, var_y, var_phi)

  return obs_refs, ref_predictor

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # use a datablock that contains a CS-PAD detector description
  data_dir = os.path.join(dials_regression, "refinement_test_data",
                          "hierarchy_test")
  datablock_path = os.path.join(data_dir, "datablock.json")
  assert os.path.exists(datablock_path)

  # load models
  from dxtbx.datablock import DataBlockFactory
  datablock = DataBlockFactory.from_serialized_format(datablock_path, check_format=False)
  im_set = datablock[0].extract_imagesets()[0]
  from copy import deepcopy
  detector = deepcopy(im_set.get_detector())
  beam = im_set.get_beam()

  # we'll invent a crystal, goniometer and scan for this test
  from dxtbx.model import Crystal
  crystal = Crystal((40.,0.,0.) ,(0.,40.,0.), (0.,0.,40.),
                          space_group_symbol = "P1")

  from dxtbx.model import goniometer_factory
  goniometer = goniometer_factory.known_axis((1., 0., 0.))

  # Build a mock scan for a 180 degree sweep
  from dxtbx.model.scan import scan_factory
  sf = scan_factory()
  scan = sf.make_scan(image_range = (1,1800),
                      exposure_times = 0.1,
                      oscillation = (0, 0.1),
                      epochs = range(1800),
                      deg = True)
  sweep_range = scan.get_oscillation_range(deg=False)
  im_width = scan.get_oscillation(deg=False)[1]
  assert sweep_range == (0., pi)
  assert approx_equal(im_width, 0.1 * pi / 180.)

  from dxtbx.model.experiment_list import ExperimentList, Experiment

  # Build an experiment list
  experiments = ExperimentList()
  experiments.append(Experiment(
        beam=beam, detector=detector, goniometer=goniometer,
        scan=scan, crystal=crystal, imageset=None))

  # simulate some reflections
  refs, ref_predictor = generate_reflections(experiments)

  # move the detector quadrants apart by 2mm both horizontally and vertically
  from dials.algorithms.refinement.parameterisation \
    import DetectorParameterisationHierarchical
  det_param = DetectorParameterisationHierarchical(detector, level=1)
  det_p_vals = det_param.get_param_vals()
  p_vals = list(det_p_vals)
  p_vals[1] += 2
  p_vals[2] -= 2
  p_vals[7] += 2
  p_vals[8] += 2
  p_vals[13] -= 2
  p_vals[14] += 2
  p_vals[19] -= 2
  p_vals[20] -= 2
  det_param.set_param_vals(p_vals)

  # reparameterise the detector at the new perturbed geometry
  det_param = DetectorParameterisationHierarchical(detector, level=1)

  # parameterise other models
  from dials.algorithms.refinement.parameterisation.beam_parameters import \
      BeamParameterisation
  from dials.algorithms.refinement.parameterisation.crystal_parameters import \
      CrystalOrientationParameterisation, CrystalUnitCellParameterisation
  beam_param = BeamParameterisation(beam, goniometer)
  xlo_param = CrystalOrientationParameterisation(crystal)
  xluc_param = CrystalUnitCellParameterisation(crystal)

  # fix beam
  beam_param.set_fixed([True]*3)

  # fix crystal
  xluc_param.set_fixed([True]*6)
  xlo_param.set_fixed([True]*3)

  # parameterisation of the prediction equation
  from dials.algorithms.refinement.parameterisation.prediction_parameters import \
      XYPhiPredictionParameterisation
  from dials.algorithms.refinement.parameterisation.parameter_report import \
      ParameterReporter
  pred_param = XYPhiPredictionParameterisation(experiments,
    [det_param], [beam_param], [xlo_param], [xluc_param])
  param_reporter = ParameterReporter([det_param], [beam_param],
                                     [xlo_param], [xluc_param])

  # reflection manager and target function
  from dials.algorithms.refinement.target import \
    LeastSquaresPositionalResidualWithRmsdCutoff
  from dials.algorithms.refinement.reflection_manager import ReflectionManager
  refman = ReflectionManager(refs, experiments, nref_per_degree=20)

  # set a very tight rmsd target of 1/10000 of a pixel
  target = LeastSquaresPositionalResidualWithRmsdCutoff(experiments,
      ref_predictor, refman, pred_param, restraints_parameterisation=None,
      frac_binsize_cutoff=0.0001)

  # minimisation engine
  from dials.algorithms.refinement.engine \
    import LevenbergMarquardtIterations as Refinery
  refinery = Refinery(target = target,
                      prediction_parameterisation = pred_param,
                      log = None,
                      verbosity = 0,
                      track_step = False,
                      track_gradient = False,
                      track_parameter_correlation = False,
                      max_iterations = 20)

  # Refiner
  from dials.algorithms.refinement.refiner import Refiner
  refiner = Refiner(reflections=refs,
                    experiments=experiments,
                    pred_param=pred_param,
                    param_reporter=param_reporter,
                    refman=refman,
                    target=target,
                    refinery=refinery,
                    verbosity=0)

  history = refiner.run()
  assert history.reason_for_termination == "RMSD target achieved"

  #compare detector with original detector
  orig_det = im_set.get_detector()
  refined_det = refiner.get_experiments()[0].detector

  from scitbx import matrix
  import math
  for op, rp in zip(orig_det, refined_det):
    # compare the origin vectors by...
    o1 = matrix.col(op.get_origin())
    o2 = matrix.col(rp.get_origin())
    # ...their relative lengths
    assert approx_equal(
      math.fabs(o1.length() - o2.length()) / o1.length(), 0, eps=1e-5)
    # ...the angle between them
    assert approx_equal(o1.accute_angle(o2), 0, eps=1e-5)

  print "OK"
  return

def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return

  test1()

if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()
