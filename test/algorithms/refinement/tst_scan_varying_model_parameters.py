#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# Python imports
from __future__ import absolute_import, division
from math import pi
import random

# CCTBX imports
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx.array_family import flex
from dxtbx.model.crystal import Crystal
from dxtbx.model.beam import Beam
from dxtbx.model.goniometer import Goniometer

# DIALS imports
from dials.algorithms.refinement.refinement_helpers \
    import get_fd_gradients, random_param_shift
from dials.algorithms.refinement.parameterisation.scan_varying_model_parameters \
    import ScanVaryingParameterSet, GaussianSmoother
from dials.algorithms.refinement.parameterisation.scan_varying_crystal_parameters \
    import ScanVaryingCrystalOrientationParameterisation, \
           ScanVaryingCrystalUnitCellParameterisation
from dials.algorithms.refinement.parameterisation.scan_varying_beam_parameters \
    import ScanVaryingBeamParameterisation
from dials.algorithms.refinement.parameterisation.scan_varying_detector_parameters \
    import ScanVaryingDetectorParameterisationSinglePanel

class SmootherTest(object):
  """Test a bare parameter set with the smoother"""
  def __init__(self, plots = False):

    # make scatterplots
    self.do_plots = plots

    # 7 values, all set to 1.0
    self.myparam = ScanVaryingParameterSet(1.0, 7)

    # Adjust a couple of the values
    self.myparam.value[3:4] = [2.0, 2.0]

    # Make a smoother with x_range as an 'image range', between 1 and 100.
    # This smoother needs 5 intervals (for 7 total values). The default
    # smoother uses an averaging window of 3 values
    self.smoother = GaussianSmoother((1, 100), 5)

  def run(self):

    assert self.smoother.num_values() == 7
    assert self.smoother.num_samples() == 5
    assert self.smoother.num_average() == 3

    # The smoother positions depend on the number of intervals but not on
    # the x_range
    assert self.smoother.positions() == [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]

    # By contrast the spacing is in units of the original, unnormalised
    # coordinate
    assert self.smoother.spacing() == 19.8

    # Use the single value smoother multiple times...
    smooth_at = [e for e in range(1, 101)]
    data = [self.smoother.value_weight(e, self.myparam) for e in smooth_at]
    vals, weights, sumweights = zip(*data)
    assert len(smooth_at) == len(vals)

    # ...and the multi value smoother once
    mvals, mweights, msumweights = self.smoother.multi_value_weight(smooth_at, self.myparam)

    # the results should be identical.
    assert (flex.double(vals) == mvals).all_eq(True)
    assert (flex.double(sumweights) == msumweights).all_eq(True)
    for v1, v2 in zip(weights, mweights.transpose().cols()):
      assert (v1.as_dense_vector() == v2.as_dense_vector()).all_eq(True)

    if self.do_plots:
      try:
        import matplotlib.pyplot as plt
        plt.ion()
        plt.scatter(smooth_at, vals)
        plt.draw()
      except ImportError as e:
        print "pyplot not available", e

    print "OK"

# Use classes to wrap up ScanVarying*Parameterisation classes with
# compose and get_state overridden, so they can be passed to the existing FD
# derivative code.
class TestOrientationModel(ScanVaryingCrystalOrientationParameterisation):

  def __init__(self, image_number, *args):

    ScanVaryingCrystalOrientationParameterisation.__init__(self, *args)

    # set overloads now, after construction of the base class
    self.compose = self._compose
    self.get_state = self._get_state

    self.set_time_point(image_number)

  def set_time_point(self, t):

    self.image_number = t
    self.compose()

  def _compose(self):
    """override for compose to pass in the requested t"""

    ScanVaryingCrystalOrientationParameterisation.compose(self,
            self.image_number)

  def _get_state(self):
    """override for get_state to do so only at the requested t"""

    # ensure the state is updated by re-composing
    self.compose()
    return ScanVaryingCrystalOrientationParameterisation.get_state(self)

class TestUnitCellModel(ScanVaryingCrystalUnitCellParameterisation):

  def __init__(self, image_number, *args):

    ScanVaryingCrystalUnitCellParameterisation.__init__(self, *args)

    # set overloads now, after construction of the base class
    self.compose = self._compose
    self.get_state = self._get_state

    self.set_time_point(image_number)

  def set_time_point(self, t):

    self.image_number = t
    self.compose()

  def _compose(self):
    """override for compose to pass in the requested t"""

    ScanVaryingCrystalUnitCellParameterisation.compose(self,
            self.image_number)

  def _get_state(self):
    """override for get_state to do so only at the requested t"""

    # ensure the state is updated by re-composing
    self.compose()
    return ScanVaryingCrystalUnitCellParameterisation.get_state(self)

class TestBeamModel(ScanVaryingBeamParameterisation):

  def __init__(self, image_number, *args):

    ScanVaryingBeamParameterisation.__init__(self, *args)

    # set overloads now, after construction of the base class
    self.compose = self._compose
    self.get_state = self._get_state

    self.set_time_point(image_number)

  def set_time_point(self, t):

    self.image_number = t
    self.compose()

  def _compose(self):
    """override for compose to pass in the requested t"""

    ScanVaryingBeamParameterisation.compose(self,
            self.image_number)

  def _get_state(self):
    """override for get_state to do so only at the requested t"""

    # ensure the state is updated by re-composing
    self.compose()
    return ScanVaryingBeamParameterisation.get_state(self)

class TestDetectorModel(ScanVaryingDetectorParameterisationSinglePanel):

  def __init__(self, image_number, *args):

    ScanVaryingDetectorParameterisationSinglePanel.__init__(self, *args)

    # set overloads now, after construction of the base class
    self.compose = self._compose
    self.get_state = self._get_state

    self.set_time_point(image_number)

  def set_time_point(self, t):

    self.image_number = t
    self.compose()

  def _compose(self):
    """override for compose to pass in the requested t"""

    ScanVaryingDetectorParameterisationSinglePanel.compose(self,
            self.image_number)

  def _get_state(self):
    """override for get_state to do so only at the requested t"""

    # ensure the state is updated by re-composing
    self.compose()
    return ScanVaryingDetectorParameterisationSinglePanel.get_state(self)

class TestScanVaryingModelParameterisation(object):

  def __init__(self, plots = False):

    # Do we make scatterplots?
    self.do_plots = plots

    # Let's say we have a scan of 100 images
    self.image_range = (1, 100)

    # Make a random P1 crystal
    a = random.uniform(10,50) * \
        self.random_direction_close_to(matrix.col((1, 0, 0)))
    b = random.uniform(10,50) * \
        self.random_direction_close_to(matrix.col((0, 1, 0)))
    c = random.uniform(10,50) * \
        self.random_direction_close_to(matrix.col((0, 0, 1)))
    self.xl = Crystal(a, b, c, space_group_symbol = "P 1")

    # Make a beam with wavelength in the range 0.8--1.2 and s0 direction close
    # to 0,0,1
    s0 = random.uniform(0.8, 1.2) * \
        self.random_direction_close_to(matrix.col((0, 0, 1)))
    self.beam = Beam(s0)

    # Make a standard goniometer model along X
    self.goniometer = Goniometer((1,0,0))

    # Make a simple single panel detector
    d1 = matrix.col((1, 0, 0))
    d2 = matrix.col((0, -1, 0))
    npx_fast = 1475
    npx_slow = 1679
    pix_size_f = pix_size_s = 0.172
    from dxtbx.model import detector_factory
    self.detector = detector_factory.make_detector("PAD", d1, d2,
      matrix.col((0, 0, -110)), (pix_size_f, pix_size_s),
      (npx_fast, npx_slow), (0, 2e20))

  def random_direction_close_to(self, vector):
    return vector.rotate_around_origin(matrix.col(
                (random.random(),
                 random.random(),
                 random.random())).normalize(),
                 random.gauss(0, 1.0),  deg = True)

class TestScanVaryingCrystalOrientationParameterisation(TestScanVaryingModelParameterisation):
  """Test a ScanVaryingCrystalOrientationParameterisation"""

  def test_num_intervals(self, nintervals):
    """Test a range of different numbers of intervals"""

    # Parameterise the crystal with the image range and five intervals. Init
    # TestOrientationModel to explore gradients at image 50, but actually
    # will try various time points in the test
    xl_op = TestOrientationModel(50, self.xl, self.image_range, nintervals)

    # How many parameters?
    num_param = xl_op.num_free()

    # shift the parameters away from zero
    p_vals = xl_op.get_param_vals()
    sigmas = [1.0] * len(p_vals)
    new_vals = random_param_shift(p_vals, sigmas)
    xl_op.set_param_vals(new_vals)

    # recalc state and gradients at image 50
    xl_op.compose()
    p_vals = xl_op.get_param_vals()
    #print "Shifted parameter vals", p_vals

    # compare analytical and finite difference derivatives at image 50
    an_ds_dp = xl_op.get_ds_dp()
    fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * num_param)
    param_names = xl_op.get_param_names()

    null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))
    for e, f in zip(an_ds_dp, fd_ds_dp):
      assert(approx_equal((e - f), null_mat, eps = 1.e-6))

    # Now test gradients at equally spaced time points across the whole
    # range
    num_points = 50
    smooth_at = []
    phi1_data = []
    phi2_data = []
    phi3_data = []
    step_size = (self.image_range[1] - self.image_range[0]) / num_points
    for t in [self.image_range[0] + e * step_size \
                for e in range(num_points + 1)]:

      # collect data for plot
      smooth_at.append(t)
      phi1_data.append(xl_op._smoother.value_weight(t, xl_op._param[0])[0])
      phi2_data.append(xl_op._smoother.value_weight(t, xl_op._param[1])[0])
      phi3_data.append(xl_op._smoother.value_weight(t, xl_op._param[2])[0])

      xl_op.set_time_point(t)
      an_ds_dp = xl_op.get_ds_dp()
      fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * num_param)
      #print t
      #print "Gradients:"
      #for s, a, f in zip(param_names, an_ds_dp, fd_ds_dp):
      #    print s
      #    print a
      #    print f
      #    print "diff:", a-f
      #    print
      #
      for e, f in zip(an_ds_dp, fd_ds_dp):
        assert(approx_equal((e - f), null_mat, eps = 1.e-6))

    if self.do_plots:
      try:
        import matplotlib.pyplot as plt
        plt.ion()
        plt.clf()
        plt.subplot(311)
        plt.cla()
        plt.scatter(smooth_at, phi1_data)
        plt.title("Phi1")
        plt.xlabel("image number")
        plt.ylabel("Phi1 (mrad)")
        plt.subplot(312)
        plt.cla()
        plt.scatter(smooth_at, phi2_data)
        plt.title("Phi2")
        plt.xlabel("image number")
        plt.ylabel("Phi2 (mrad)")
        plt.subplot(313)
        plt.cla()
        plt.scatter(smooth_at, phi3_data)
        plt.title("Phi3")
        plt.xlabel("image number")
        plt.ylabel("Phi3 (mrad)")
        plt.suptitle("Parameter smoothing with %d intervals" % nintervals)
        plt.draw()
      except ImportError as e:
        print "pyplot not available", e

    print "OK"

  def test_random(self):
    """Test random initial orientations, random parameter shifts and random
    times"""

    attempts = 100
    failures = 0
    null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))

    for i in range(attempts):

      # make a new P1 random crystal and parameterise it
      a = random.uniform(10,50) * \
              self.random_direction_close_to(matrix.col((1, 0, 0)))
      b = random.uniform(10,50) * \
              self.random_direction_close_to(matrix.col((0, 1, 0)))
      c = random.uniform(10,50) * \
              self.random_direction_close_to(matrix.col((0, 0, 1)))
      xl = Crystal(a, b, c, space_group_symbol="P 1")

      xl_op = TestOrientationModel(50, xl, self.image_range, 5)

      # How many parameters?
      num_param = xl_op.num_free()

      # apply random parameter shifts to the orientation (2.0 mrad each
      # checkpoint)
      p_vals = xl_op.get_param_vals()
      sigmas = [2.0] * len(p_vals)
      new_vals = random_param_shift(p_vals, sigmas)
      xl_op.set_param_vals(new_vals)

      # select random time point at which to make comparisons
      t = random.uniform(*self.image_range)
      xl_op.set_time_point(t)

      # compare analytical and finite difference derivatives
      xl_op_an_ds_dp = xl_op.get_ds_dp()
      xl_op_fd_ds_dp = get_fd_gradients(xl_op, [1.e-6 * pi/180] * \
                                        num_param)

      for j in range(num_param):
        try:
          assert(approx_equal((xl_op_fd_ds_dp[j] - xl_op_an_ds_dp[j]),
                              null_mat, eps = 1.e-6))
        except Exception:
          failures += 1
          print "for try", i
          print "failure for parameter number", j
          print "of the orientation parameterisation"
          print "with fd_ds_dp = "
          print xl_op_fd_ds_dp[j]
          print "and an_ds_dp = "
          print xl_op_an_ds_dp[j]
          print "so that difference fd_ds_dp - an_ds_dp ="
          print xl_op_fd_ds_dp[j] - xl_op_an_ds_dp[j]

    if failures == 0: print "OK"

  def run(self):

    for n in (1, 2, 3, 4, 5, 6, 7):
      self.test_num_intervals(n)

    self.test_random()

class TestScanVaryingCrystalUnitCellParameterisation(TestScanVaryingModelParameterisation):
  """Basic test of a ScanVaryingCrystalUnitCellParameterisation"""

  def test_num_intervals(self, nintervals):
    """Test a range of different numbers of intervals"""

    # Parameterise the crystal with the image range and five intervals. Init
    # TestOrientationModel to explore gradients at image 50
    xl_ucp = TestUnitCellModel(50, self.xl, self.image_range, nintervals)

    # How many parameters?
    num_param = xl_ucp.num_free()

    # apply a random parameter shift to the unit cell, on order of 2% of
    # the initial metrical matrix parameters
    p_vals = xl_ucp.get_param_vals()
    sigmas = [0.02  * p for p in p_vals]
    new_vals = random_param_shift(p_vals, sigmas)
    xl_ucp.set_param_vals(new_vals)

    # calculate state and gradients at image 50
    xl_ucp.compose()

    null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))

    # compare analytical and finite difference derivatives at image 50
    an_ds_dp = xl_ucp.get_ds_dp()
    fd_ds_dp = get_fd_gradients(xl_ucp, [1.e-7] * num_param)
    param_names = xl_ucp.get_param_names()

    for e, f in zip(an_ds_dp, fd_ds_dp):
      assert(approx_equal((e - f), null_mat, eps = 1.e-6))

    print "OK"
    return

  def run(self):

    for n in (1, 2, 3, 4, 5, 6, 7):
      self.test_num_intervals(n)

    return

class TestScanVaryingBeamParameterisation(TestScanVaryingModelParameterisation):
  """Basic test of a ScanVaryingBeamParameterisation"""

  def test_num_intervals(self, nintervals):
    """Test a range of different numbers of intervals"""

    # Parameterise the crystal with the image range and five intervals. Init
    # TestOrientationModel to explore gradients at image 50
    beam_p = TestBeamModel(50, self.beam, self.image_range, nintervals,
      self.goniometer)

    # How many parameters?
    num_param = beam_p.num_free()

    # apply a random parameter shift to the beam, on order of 2% of
    # the initial values
    p_vals = beam_p.get_param_vals()
    sigmas = [0.02  * p for p in p_vals]
    new_vals = random_param_shift(p_vals, sigmas)
    beam_p.set_param_vals(new_vals)

    # calculate state and gradients at image 50
    beam_p.compose()

    null_vec = matrix.col((0., 0., 0.))

    # compare analytical and finite difference derivatives at image 50
    an_ds_dp = beam_p.get_ds_dp()
    fd_ds_dp = get_fd_gradients(beam_p, [1.e-7] * num_param)
    param_names = beam_p.get_param_names()

    for e, f in zip(an_ds_dp, fd_ds_dp):
      assert(approx_equal((e - f), null_vec, eps = 1.e-6))

    print "OK"
    return

  def run(self):

    for n in (1, 2, 3, 4, 5, 6, 7):
      self.test_num_intervals(n)

    return

class TestScanVaryingDetectorParameterisation(TestScanVaryingModelParameterisation):
  """Basic test of a ScanVaryingDetectorParameterisationSinglePanel"""

  def test_num_intervals(self, nintervals):
    """Test a range of different numbers of intervals"""

    # Parameterise the detector with the image range and five intervals. Init
    # TestOrientationModel to explore gradients at image 50
    det_p = TestDetectorModel(50, self.detector, self.image_range, nintervals)

    # How many parameters?
    num_param = det_p.num_free()

    # apply a random parameter shift to the detector, on order of 2% of
    # the initial values
    p_vals = det_p.get_param_vals()
    sigmas = [0.02  * p for p in p_vals]
    new_vals = random_param_shift(p_vals, sigmas)
    det_p.set_param_vals(new_vals)

    # calculate state and gradients at image 50
    det_p.compose()

    null_mat = matrix.sqr((0., 0., 0., 0., 0., 0., 0., 0., 0.))

    # compare analytical and finite difference derivatives at image 50
    an_ds_dp = det_p.get_ds_dp()
    fd_ds_dp = get_fd_gradients(det_p, [1.e-7] * num_param)
    param_names = det_p.get_param_names()

    for e, f in zip(an_ds_dp, fd_ds_dp):
      assert(approx_equal((e - f), null_mat, eps = 1.e-6))

    print "OK"
    return

  def run(self):

    for n in (1, 2, 3, 4, 5, 6, 7):
      self.test_num_intervals(n)

    return


if __name__ == '__main__':

  print "Testing the GaussianSmoother and ScanVaryingParameterSet"
  test = SmootherTest(plots = False)
  test.run()
  print

  print "Testing the ScanVaryingCrystalOrientationParameterisation"
  test = TestScanVaryingCrystalOrientationParameterisation(plots = False)
  test.run()
  print

  print "Testing the ScanVaryingCrystalUnitCellParameterisation"
  test = TestScanVaryingCrystalUnitCellParameterisation(plots = False)
  test.run()

  print "Testing the ScanVaryingBeamParameterisation"
  test = TestScanVaryingBeamParameterisation(plots = False)
  test.run()

  print "Testing the ScanVaryingDetectorParameterisationSinglePanel"
  test = TestScanVaryingDetectorParameterisation(plots = False)
  test.run()
