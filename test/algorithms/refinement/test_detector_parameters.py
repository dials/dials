from __future__ import absolute_import, division, print_function
from math import sin, cos, pi, sqrt
import random

from scitbx import matrix
from libtbx.test_utils import approx_equal

from dxtbx.model import Panel, Detector
from dxtbx.model import DetectorFactory
from dxtbx.model import BeamFactory

from dials.algorithms.refinement.refinement_helpers \
    import get_fd_gradients, random_param_shift
from dials.algorithms.refinement.parameterisation.detector_parameters \
    import DetectorParameterisationSinglePanel, DetectorParameterisationMultiPanel


def random_panel(lim = (0, 50)):
  """For testing, return a square panel with a randomised position
  and orientation"""

  # start with a randomised origin vector
  o = matrix.col((random.uniform(-200, 200),
                  random.uniform(-200, 200),
                  random.uniform(-200, 200)))

  # two orthogonal unit vectors randomly oriented in the normal plane
  # of the origin vector
  u1 = o.ortho().normalize()
  u2 = o.cross(u1).normalize()
  #theta = random.uniform(0, 2. * pi)
  u1 = u1.rotate_around_origin(o, pi/12)
  u2 = u2.rotate_around_origin(o, pi/12)

  # offset the plane normal from the origin vector by random rotations
  # of up to 45 degrees for each direction
  u1 = u1.rotate_around_origin(u2, random.uniform(-pi/2., pi/2.))
  u2 = u2.rotate_around_origin(u1, random.uniform(-pi/2., pi/2.))

  return Panel("PAD", "Panel", u1, u2, o,
          (lim[1]/200, lim[1]/200), (200, 200), (0, 2e20), 0.0, "")

# local function required to make a 3x3 multi-panel detector
def make_multi_panel(single_panel_detector):
  """Create a 3x3 multi-panel detector filling the same space as
  a supplied single panel detector"""

  from dials.test.algorithms.refinement.test_multi_panel_detector_parameterisation \
      import make_panel_in_array
  from dials.test.algorithms.refinement.setup_geometry import \
      random_vector_close_to

  multi_panel_detector = Detector()
  for x in range(3):
    for y in range(3):
      new_panel = make_panel_in_array(
                      (x, y), single_panel_detector[0])
      multi_panel_detector.add_panel(new_panel)

  # apply small random shifts & rotations to each panel
  for p in multi_panel_detector:

    # perturb origin vector
    o_multiplier = random.gauss(1.0, 0.01)
    new_origin = random_vector_close_to(p.get_origin(), sd=0.1)
    new_origin *= o_multiplier

    # perturb fast direction vector
    new_dir1 = random_vector_close_to(p.get_fast_axis(), sd=0.5)

    # create vector in the plane of dir1-dir2
    dir1_dir2 = random_vector_close_to(p.get_slow_axis(), sd=0.5)

    # find normal to panel plane and thus new slow direction vector
    dn = new_dir1.cross(dir1_dir2)
    new_dir2 = dn.cross(new_dir1)

    # set panel frame
    p.set_frame(new_dir1, new_dir2, new_origin)

  return multi_panel_detector

def test():
  # set the random seed to make the test reproducible
  random.seed(1337)

  # set up a simple detector frame with directions aligned with
  # principal axes and sensor origin located on the z-axis at -110
  d1 = matrix.col((1, 0, 0))
  d2 = matrix.col((0, -1, 0))
  #lim = (0,50)
  npx_fast = 1475
  npx_slow = 1679
  pix_size_f = pix_size_s = 0.172
  detector = DetectorFactory.make_detector("PAD", d1, d2,
      matrix.col((0, 0, -110)), (pix_size_f, pix_size_s),
      (npx_fast, npx_slow), (0, 2e20))

  dp = DetectorParameterisationSinglePanel(detector)
  beam = BeamFactory().make_beam(
          sample_to_source=-1*(matrix.col((0, 0, -110)) + 10 * d1 + 10 * d2),
          wavelength=1.0)

  # Test change of parameters
  # =========================

  # 1. shift detector plane so that the z-axis intercepts its centre
  # at a distance of 100 along the initial normal direction. As the
  # initial normal is along -z, we expect the frame to intercept the
  # z-axis at -100.

  p_vals = dp.get_param_vals()
  p_vals[0:3] = [100., 0., 0.]
  dp.set_param_vals(p_vals)
  detector = dp._model
  assert(len(detector) == 1)
  panel = detector[0]
  v1 = matrix.col(panel.get_origin())
  v2 = matrix.col((0., 0., 1.))
  assert(approx_equal(v1.dot(v2), -100.))

  # 2. rotate frame around its initial normal by +90 degrees. Only d1
  # and d2 should change. As we rotate clockwise around the initial
  # normal (-z direction) then d1 should rotate onto the original
  # direction d2, and d2 should rotate to negative of the original
  # direction d1

  p_vals[3] = 1000. * pi/2 # set tau1 value
  dp.set_param_vals(p_vals)

  detector = dp._model
  assert(len(detector) == 1)
  panel = detector[0]
  assert(approx_equal(matrix.col(panel.get_fast_axis()).dot(dp._initial_state['d1']), 0.))
  assert(approx_equal(matrix.col(panel.get_slow_axis()).dot(dp._initial_state['d2']), 0.))
  assert(approx_equal(matrix.col(panel.get_normal()).dot(dp._initial_state['dn']), 1.))

  # 3. no rotation around initial normal, +10 degrees around initial
  # d1 direction and +10 degrees around initial d2. Check d1 and d2
  # match paper calculation

  p_vals[3] = 0.    # tau1
  p_vals[4] = 1000. * pi/18 # tau2
  p_vals[5] = 1000. * pi/18 # tau3
  dp.set_param_vals(p_vals)

  # paper calculation values
  v1 = matrix.col((cos(pi/18), 0, sin(pi/18)))
  v2 = matrix.col((sin(pi/18)**2,
                   -cos(pi/18),
                   sqrt((2*sin(pi/36)*sin(pi/18))**2 - sin(pi/18)**4) - sin(pi/18)))

  detector = dp._model
  assert(len(detector) == 1)
  panel = detector[0]
  assert(approx_equal(matrix.col(panel.get_fast_axis()).dot(v1), 1.))
  assert(approx_equal(matrix.col(panel.get_slow_axis()).dot(v2), 1.))

  # 4. Test fixing and unfixing of parameters
  p_vals = [100., 0., 0., 1000.*pi/18, 1000.*pi/18, 1000.*pi/18]
  dp.set_param_vals(p_vals)
  f = dp.get_fixed()
  f[0:3] = [True] * 3
  dp.set_fixed(f)
  p_vals2 = [0., 0., 0.]
  dp.set_param_vals(p_vals2)
  assert(dp.get_param_vals(only_free = False) == [100., 0., 0., 0., 0., 0.])

  an_ds_dp = dp.get_ds_dp()
  assert(len(an_ds_dp) == 3)

  f[0:3] = [False] * 3
  dp.set_fixed(f)
  p_vals = dp.get_param_vals()
  p_vals2 = [a + b for a, b in zip(p_vals, [-10., 1., 1., 0., 0., 0.])]
  dp.set_param_vals(p_vals2)
  assert(dp.get_param_vals() == [90., 1., 1., 0., 0., 0.])

  # 5. Tests of the calculation of derivatives

  # Now using parameterisation in mrad

  # random initial orientations with a random parameter shift at each
  attempts = 100
  for i in range(attempts):

    # create random initial position
    det = Detector(random_panel())
    dp = DetectorParameterisationSinglePanel(det)

    # apply a random parameter shift
    p_vals = dp.get_param_vals()
    p_vals = random_param_shift(p_vals, [10, 10, 10, 1000.*pi/18,
                                         1000.*pi/18, 1000.*pi/18])
    dp.set_param_vals(p_vals)

    # obtain current sensor state
    #state = dp.get_state()

    # compare analytical and finite difference derivatives.
    an_ds_dp = dp.get_ds_dp(multi_state_elt=0)
    fd_ds_dp = get_fd_gradients(dp,
                    [1.e-6] * 3 + [1.e-4 * pi/180] * 3)

    for j in range(6):
      assert(approx_equal((fd_ds_dp[j] - an_ds_dp[j]),
              matrix.sqr((0., 0., 0.,
                          0., 0., 0.,
                          0., 0., 0.)), eps = 1.e-6)), textwrap.dedent("""\
        Failure comparing analytical with finite difference derivatives.
        Failure in try {i}
        failure for parameter number {j}
        of the orientation parameterisation
        with fd_ds_dp =
        {fd}
        and an_ds_dp =
        {an}
        so that difference fd_ds_dp - an_ds_dp =
        {diff}
        """).format(i=i, j=j, fd=fd_ds_dp[j], an=an_ds_dp[j], diff=fd_ds_dp[j] - an_ds_dp[j])

  # 5. Test a multi-panel detector with non-coplanar panels.

  # place a beam at the centre of the single panel detector (need a
  # beam to initialise the multi-panel detector parameterisation)
  lim = det[0].get_image_size_mm()
  shift1 = lim[0] / 2.
  shift2 = lim[1] / 2.
  beam_centre = matrix.col(det[0].get_origin()) + \
                shift1 * matrix.col(det[0].get_fast_axis()) + \
                shift2 * matrix.col(det[0].get_slow_axis())
  beam = BeamFactory().make_beam(sample_to_source=-1.*beam_centre,
                                  wavelength=1.0)

  multi_panel_detector = make_multi_panel(det)

  # parameterise this detector
  dp = DetectorParameterisationMultiPanel(multi_panel_detector,
                                          beam)

  # ensure the beam still intersects the central panel
  intersection = multi_panel_detector.get_ray_intersection(
                                                      beam.get_s0())
  assert intersection[0] == 4

  # record the offsets and dir1s, dir2s
  offsets_before_shift = dp._offsets
  dir1s_before_shift = dp._dir1s
  dir2s_before_shift = dp._dir2s

  # apply a random parameter shift (~10 mm distances, ~50 mrad angles)
  p_vals = dp.get_param_vals()
  p_vals = random_param_shift(p_vals, [10, 10, 10, 50,
                                       50, 50])

  # reparameterise the detector
  dp = DetectorParameterisationMultiPanel(multi_panel_detector,
                                          beam)

  # record the offsets and dir1s, dir2s
  offsets_after_shift = dp._offsets
  dir1s_after_shift = dp._dir1s
  dir2s_after_shift = dp._dir2s

  # ensure the offsets, dir1s and dir2s are the same. This means that
  # each panel in the detector moved with the others as a rigid body
  for a, b in zip(offsets_before_shift, offsets_after_shift):
    assert(approx_equal(a, b, eps=1.e-10))

  for a, b in zip(dir1s_before_shift, dir1s_after_shift):
    assert(approx_equal(a, b, eps=1.e-10))

  for a, b in zip(dir2s_before_shift, dir2s_after_shift):
    assert(approx_equal(a, b, eps=1.e-10))

  attempts = 5
  for i in range(attempts):

    multi_panel_detector = make_multi_panel(det)

    # parameterise this detector
    dp = DetectorParameterisationMultiPanel(multi_panel_detector,
                                            beam)
    p_vals = dp.get_param_vals()

    # apply a random parameter shift
    p_vals = random_param_shift(p_vals, [10, 10, 10, 1000.*pi/18,
                                         1000.*pi/18, 1000.*pi/18])
    dp.set_param_vals(p_vals)

    # obtain current state of the 1st panel
    state = dp.get_state()

    # compare analytical and finite difference derivatives
    # get_fd_gradients will implicitly only get gradients for the
    # 1st panel in the detector, so explicitly get the same for the
    # analytical gradients

    for j in range(9):

      an_ds_dp = dp.get_ds_dp(multi_state_elt=j)
      fd_ds_dp = get_fd_gradients(dp, [1.e-7] * dp.num_free(),
                                  multi_state_elt=j)

      for k in range(6):
        assert approx_equal((fd_ds_dp[k] - matrix.sqr(an_ds_dp[k])),
                matrix.sqr((0., 0., 0.,
                            0., 0., 0.,
                            0., 0., 0.)),
                            eps = 1.e-5,
                            out = None), textwrap.dedent("""\
        Failure comparing analytical with finite difference derivatives.
        Failure in try {i}
        for panel number {j]
        failure for parameter number {k}
        of the orientation parameterisation
        with fd_ds_dp =
        {fd}
        and an_ds_dp =
        {an}
        so that difference fd_ds_dp - an_ds_dp =
        {diff}
        """).format(i=i, j=j, k=k, fd=fd_ds_dp[k], an=an_ds_dp[k], diff=fd_ds_dp[k] - matrix.sqr(an_ds_dp[k]))
