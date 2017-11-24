#!/usr/bin/env python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division
import sys
from cctbx.sgtbx import space_group, space_group_symbols
from libtbx.test_utils import approx_equal
from libtbx.phil import parse
from math import pi
from scitbx.array_family import flex
from dials.test.algorithms.refinement.setup_geometry import Extract
from dxtbx.model import ScanFactory
from dials.algorithms.spot_prediction import IndexGenerator, ray_intersection
from dxtbx.model.experiment_list import ExperimentList, Experiment
from dials.algorithms.refinement.prediction import ScansRayPredictor, \
  ExperimentsPredictor
from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters import \
    ScanVaryingPredictionParameterisation
from dials.algorithms.refinement.parameterisation.scan_varying_crystal_parameters \
    import ScanVaryingCrystalOrientationParameterisation, \
           ScanVaryingCrystalUnitCellParameterisation
from dials.algorithms.refinement.parameterisation.scan_varying_beam_parameters \
    import ScanVaryingBeamParameterisation
from dials.algorithms.refinement.parameterisation.scan_varying_detector_parameters \
    import ScanVaryingDetectorParameterisationSinglePanel
from dials.algorithms.refinement.parameterisation.scan_varying_goniometer_parameters \
    import ScanVaryingGoniometerParameterisation

class Test(object):

  def create_models(self, cmdline_overrides=None):

    if cmdline_overrides is None:
      cmdline_overrides = []
    overrides = """geometry.parameters.crystal.a.length.range = 10 50
geometry.parameters.crystal.b.length.range = 10 50
geometry.parameters.crystal.c.length.range = 10 50"""

    master_phil = parse("""
    include scope dials.test.algorithms.refinement.geometry_phil
    """, process_includes=True)

    # Extract models
    models = Extract(master_phil, overrides, cmdline_args = cmdline_overrides)
    self.detector = models.detector
    self.goniometer = models.goniometer
    self.crystal = models.crystal
    self.beam = models.beam

    # Make a scan of 1-360 * 0.5 deg images
    sf = ScanFactory()
    self.scan = sf.make_scan((1,360), 0.5, (0, 0.5), range(360))

    # Generate an ExperimentList
    self.experiments = ExperimentList()
    self.experiments.append(Experiment(
          beam=self.beam, detector=self.detector, goniometer=self.goniometer,
          scan=self.scan, crystal=self.crystal, imageset=None))

    # Create a reflection predictor for the experiments
    self.ref_predictor = ExperimentsPredictor(self.experiments)

    # Create scan-varying parameterisations of these models, with 5 samples
    self.det_param = ScanVaryingDetectorParameterisationSinglePanel(
            self.detector, self.scan.get_array_range(), 5)
    self.s0_param = ScanVaryingBeamParameterisation(
            self.beam, self.scan.get_array_range(), 5, self.goniometer)
    self.xlo_param = ScanVaryingCrystalOrientationParameterisation(
            self.crystal, self.scan.get_array_range(), 5)
    self.xluc_param = ScanVaryingCrystalUnitCellParameterisation(
            self.crystal, self.scan.get_array_range(), 5)
    self.gon_param = ScanVaryingGoniometerParameterisation(
            self.goniometer, self.scan.get_array_range(), 5, self.beam)

    return

  def generate_reflections(self):
    sweep_range = self.scan.get_oscillation_range(deg=False)
    resolution = 2.0
    index_generator = IndexGenerator(self.crystal.get_unit_cell(),
                          space_group(space_group_symbols(1).hall()).type(),
                          resolution)
    indices = index_generator.to_array()

    # Predict rays within the sweep range
    ray_predictor = ScansRayPredictor(self.experiments, sweep_range)
    obs_refs = ray_predictor(indices)

    # Take only those rays that intersect the detector
    intersects = ray_intersection(self.detector, obs_refs)
    obs_refs = obs_refs.select(intersects)

    # Re-predict using the Experiments predictor for all these reflections. The
    # result is the same, but we gain also the flags and xyzcal.px columns
    obs_refs['id'] = flex.int(len(obs_refs), 0)
    obs_refs = self.ref_predictor(obs_refs)

    # Set 'observed' centroids from the predicted ones
    obs_refs['xyzobs.mm.value'] = obs_refs['xyzcal.mm']

    # Invent some variances for the centroid positions of the simulated data
    im_width = 0.1 * pi / 180.
    px_size = self.detector[0].get_pixel_size()
    var_x = flex.double(len(obs_refs), (px_size[0] / 2.)**2)
    var_y = flex.double(len(obs_refs), (px_size[1] / 2.)**2)
    var_phi = flex.double(len(obs_refs), (im_width / 2.)**2)
    obs_refs['xyzobs.mm.variance'] = flex.vec3_double(var_x, var_y, var_phi)

    # set the flex random seed to an 'uninteresting' number
    flex.set_random_seed(12407)

    # take 5 random reflections for speed
    reflections = obs_refs.select(flex.random_selection(len(obs_refs), 5))

    # use a BlockCalculator to calculate the blocks per image
    from dials.algorithms.refinement.reflection_manager import BlockCalculator
    block_calculator = BlockCalculator(self.experiments, reflections)
    reflections = block_calculator.per_image()

    return reflections

  def __call__(self, cmdline_overrides):

    self.create_models(cmdline_overrides)
    reflections = self.generate_reflections()

    # use a ReflectionManager to exclude reflections too close to the spindle,
    # plus set the frame numbers
    from dials.algorithms.refinement.reflection_manager import ReflectionManager
    refman = ReflectionManager(reflections, self.experiments,
      outlier_detector=None)

    # create prediction parameterisation of the requested type

    # FIXME At the moment the test of scan-varying goniometer parameters will
    # fail, because ScanVaryingReflectionPredictor.for_reflection_table, used
    # in managed_predictors.py, cannot take a varying rotation axis. This has
    # to be changed before this test can pass. Until then, remove the
    # goniometer parameterisation from the test
    #pred_param = ScanVaryingPredictionParameterisation(self.experiments,
    #    [self.det_param], [self.s0_param], [self.xlo_param], [self.xluc_param],
    #    [self.gon_param])
    pred_param = ScanVaryingPredictionParameterisation(self.experiments,
        [self.det_param], [self.s0_param], [self.xlo_param], [self.xluc_param])

    # make a target to ensure reflections are predicted and refman is finalised
    from dials.algorithms.refinement.target import \
      LeastSquaresPositionalResidualWithRmsdCutoff
    target = LeastSquaresPositionalResidualWithRmsdCutoff(self.experiments,
        self.ref_predictor, refman, pred_param, restraints_parameterisation=None)

    # keep only those reflections that pass inclusion criteria and have predictions
    reflections = refman.get_matches()

    # get analytical gradients
    pred_param.compose(reflections)
    an_grads = pred_param.get_gradients(reflections)

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    p_names = pred_param.get_param_names()
    deltas = [1.e-7] * len(p_vals)

    for i in range(len(deltas)):

      val = p_vals[i]

      p_vals[i] -= deltas[i] / 2.
      pred_param.set_param_vals(p_vals)
      pred_param.compose(reflections)

      self.ref_predictor(reflections)

      rev_state = reflections['xyzcal.mm'].deep_copy()

      p_vals[i] += deltas[i]
      pred_param.set_param_vals(p_vals)
      pred_param.compose(reflections)

      self.ref_predictor(reflections)

      fwd_state = reflections['xyzcal.mm'].deep_copy()
      p_vals[i] = val

      fd = (fwd_state - rev_state)
      x_grads, y_grads, phi_grads = fd.parts()
      x_grads /= deltas[i]
      y_grads /= deltas[i]
      phi_grads /= deltas[i]

      try:
        for n, (a,b) in enumerate(zip(x_grads, an_grads[i]["dX_dp"])):
          assert approx_equal(a, b, eps=1.e-5)
        for n, (a,b) in enumerate(zip(y_grads, an_grads[i]["dY_dp"])):
          assert approx_equal(a, b, eps=1.e-5)
        for n, (a,b) in enumerate(zip(phi_grads, an_grads[i]["dphi_dp"])):
          assert approx_equal(a, b, eps=1.e-5)
      except AssertionError:
        print "Failure for {0}".format(p_names[i])
        raise

    # return to the initial state
    pred_param.set_param_vals(p_vals)
    pred_param.compose(reflections)
    print "OK"

    return

if __name__ == "__main__":

  from time import time
  start_time = time()

  cmdline_overrides = sys.argv[1:]

  test1 = Test()
  test1(cmdline_overrides)

  finish_time = time()
  print "Time Taken: ",finish_time - start_time
