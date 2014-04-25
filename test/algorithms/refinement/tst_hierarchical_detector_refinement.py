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
from __future__ import division
import os
import shutil
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from math import pi
#from libtbx.test_utils import open_tmp_directory

def generate_reflections(experiments):

  from dials.algorithms.spot_prediction import IndexGenerator
  from dials.algorithms.refinement.prediction import \
    ScansRayPredictor, ExperimentsPredictor
  from dials.algorithms.spot_prediction import ray_intersection
  from cctbx.sgtbx import space_group, space_group_symbols

  #beam = experiments[0].beam
  detector = experiments[0].detector
  crystal = experiments[0].crystal

  # All indices in a 2.0 Angstrom sphere
  resolution = 2.0
  index_generator = IndexGenerator(crystal.get_unit_cell(),
                  space_group(space_group_symbols(1).hall()).type(), resolution)
  indices = index_generator.to_array()

  # Build a reflection predictor
  scan = experiments[0].scan
  sweep_range = scan.get_oscillation_range(deg=False)
  ref_predictor = ScansRayPredictor(experiments, sweep_range)

  obs_refs = ref_predictor.predict(indices)

  # Invent some variances for the centroid positions of the simulated data
  im_width = 0.1 * pi / 180.
  px_size = detector[0].get_pixel_size()
  var_x = (px_size[0] / 2.)**2
  var_y = (px_size[1] / 2.)**2
  var_phi = (im_width / 2.)**2

  obs_refs = ray_intersection(detector, obs_refs)
  for ref in obs_refs:

    # set the 'observed' centroids
    ref.centroid_position = ref.image_coord_mm + (ref.rotation_angle, )

    # set the centroid variance
    ref.centroid_variance = (var_x, var_y ,var_phi)

    # set the frame number, calculated from rotation angle
    ref.frame_number = scan.get_image_index_from_angle(
        ref.rotation_angle, deg=False)

  # redefine the reflection predictor to use the ExperimentsPredictor class
  ref_predictor = ExperimentsPredictor(experiments)

  return obs_refs.to_table(centroid_is_mm=True), ref_predictor

def test1():

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # use the i04_weak_data for this test
  data_dir = os.path.join(dials_regression, "refinement_test_data", "hierarchy_test")
  datablock_path = os.path.join(data_dir, "datablock.json")

  assert os.path.exists(datablock_path)

  # load models
  from dxtbx.datablock import DataBlockFactory
  datablock = DataBlockFactory.from_serialized_format(datablock_path, check_format=False)
  im_set = datablock[0].extract_imagesets()[0]
  from copy import deepcopy
  detector = deepcopy(im_set.get_detector())
  beam = im_set.get_beam()

  # we'll make a crystal, goniometer and scan for this test
  from cctbx.crystal.crystal_model import crystal_model
  crystal = crystal_model((40.,0.,0.) ,(0.,40.,0.), (0.,0.,40.),
                          space_group_symbol = "P1")

  from dials.model.experiment import goniometer_factory
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
  temp = scan.get_oscillation(deg=False)
  im_width = temp[1] - temp[0]
  assert sweep_range == (0., pi)
  assert approx_equal(im_width, 0.1 * pi / 180.)

  from dials.model.experiment.experiment_list import ExperimentList, Experiment

  # Build an experiment list
  experiments = ExperimentList()
  experiments.append(Experiment(
        beam=beam, detector=detector, goniometer=goniometer,
        scan=scan, crystal=crystal, imageset=None))

  # simulate some reflections
  refs, ref_predictor = generate_reflections(experiments)

  # move the detector quadrants apart by 2mm both horizontally and vertically
  from dials.algorithms.refinement.parameterisation import DetectorParameterisationHierarchical
  from dials.algorithms.refinement.parameterisation.detector_parameters \
   import DetectorParameterisationHierarchical2
  #det_param = DetectorParameterisationHierarchical(detector, beam, level=1)
  det_param = DetectorParameterisationHierarchical2(detector, level=1)
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
  det_param = DetectorParameterisationHierarchical(detector, beam, level=1)

  # Parameterise other models
  from dials.algorithms.refinement.parameterisation.beam_parameters import \
      BeamParameterisationOrientation
  from dials.algorithms.refinement.parameterisation.crystal_parameters import \
      CrystalOrientationParameterisation, CrystalUnitCellParameterisation
  beam_param = BeamParameterisationOrientation(beam, goniometer)
  xlo_param = CrystalOrientationParameterisation(crystal)
  xluc_param = CrystalUnitCellParameterisation(crystal)

  # Fix beam
  beam_param.set_fixed([True]*2)

  # Fix crystal
  xluc_param.set_fixed([True]*6)
  xlo_param.set_fixed([True]*3)

  # Parameterisation of the prediction equation
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
  # very tight rmsd target of 1/1000 of a pixel
  target = LeastSquaresPositionalResidualWithRmsdCutoff(experiments,
      ref_predictor, refman, pred_param,frac_binsize_cutoff=0.001)

  # minimisation engine
  from dials.algorithms.refinement.engine import LevenbergMarquardtIterations as Refinery
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
                    crystal_ids=[0],
                    pred_param=pred_param,
                    param_reporter=param_reporter,
                    refman=refman,
                    target=target,
                    refinery=refinery,
                    verbosity=0)

  history = refiner.run()
  assert history.reason_for_termination == "RMSD target achieved"

  #from dials.util.command_line import interactive_console; interactive_console()
  #f=open("temp2.dat","w")
  #for ref in refs:
  #  p = ref['panel']
  #  impact = ref['xyzobs.mm.value'][0:2]
  #  coord = detector[p].get_lab_coord(impact)
  #  #f.write("{0}\n".format(ref['panel']))
  #  f.write("{0} {1} {2}\n".format(*coord))
  #f.close()

  # compare detector with original detector
  #orig_det = im_set.get_detector()
  #refined_det = refiner.get_experiments()[0].detector
  #
  #from dials.util.command_line import interactive_console; interactive_console()
  #from scitbx import matrix
  #for op, rp in zip(orig_det, refined_det):
  #  # difference between origin vectors
  #  o1 = matrix.col(op.get_origin())
  #  o2 = matrix.col(rp.get_origin())
  #  test = (o1 - o2).length()
  #  # scale by length of origin vector
  #  test /= o1.length()
  #  approx_equal(test, 0., eps=1e-5)

  print "OK"
  return

def run():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tests in " + __file__ + " as dials_regression not present"
    return

  test1()

if __name__ == '__main__':
  import sys
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()
