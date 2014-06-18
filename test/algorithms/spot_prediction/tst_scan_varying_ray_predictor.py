class Test:

  def __init__(self):
    from scitbx import matrix
    # cubic, 50A cell, 1A radiation, 1 deg osciillation, everything ideal
    a = 50.0
    self.ub = matrix.sqr((1.0 / a, 0.0, 0.0,
                         0.0, 1.0 / a, 0.0,
                         0.0, 0.0, 1.0 / a))

    self.axis = matrix.col((0, 1, 0))
    self.s0 = matrix.col((0, 0, 1))

    self.dmin = 1.5
    self.margin = 1
    self.image_range = (1, 10)
    self.dphi = (0.0, 1.0)

  def run(self):
    from dials.algorithms.spot_prediction import ReekeIndexGenerator
    i0, i1 = self.image_range

    # Check indices are the same for each frame
    for frame in range(i0-1,i1):

      # Get the setting matrices
      ub_beg, ub_end = self.get_ub(frame)

      # Get the miller indices
      r = ReekeIndexGenerator(ub_beg, ub_end, self.axis,
                              self.s0, self.dmin, self.margin)

      hkl = r.to_array()
      py_rays = self.predict_py(hkl, frame, ub_beg, ub_end)
      cc_rays = self.predict_cc(hkl, frame, ub_beg, ub_end)

      assert(len(py_rays) == len(cc_rays))

      eps = 1e-7
      for r1, r2 in zip(py_rays, cc_rays):
        count = [r1,r2].count(None)
        if count == 0:
          a1 = r1.rotation_angle
          a2 = r2.angle
          e1 = r1.entering
          e2 = r2.entering
          s11 = r1.beam_vector
          s12 = r2.s1
          assert(e1 == e2)
          assert(abs(a1 - a2) < eps)
          assert(all(abs(a - b) < eps for a, b in zip(s11, s12)))
        else:
          assert(count == 2)

      print 'OK'

  def get_ub(self, frame):
    from scitbx import matrix
    import scitbx.math

    angle_beg = self.dphi[0] + frame * self.dphi[1]
    angle_end = self.dphi[0] + (frame+1) * self.dphi[1]

    r_osc_beg = matrix.sqr(
      scitbx.math.r3_rotation_axis_and_angle_as_matrix(
      axis = self.axis, angle = angle_beg, deg=True))

    r_osc_end = matrix.sqr(
      scitbx.math.r3_rotation_axis_and_angle_as_matrix(
      axis = self.axis, angle = angle_end, deg=True))

    ub_beg = r_osc_beg * self.ub
    ub_end = r_osc_end * self.ub
    return ub_beg, ub_end

  def predict_py(self, hkl, frame, A1, A2):
    from dials.algorithms.refinement.prediction.predictors import \
      ScanVaryingReflectionPredictor
    from dxtbx.model import Beam, Goniometer, Scan
    from dxtbx.model.crystal import crystal_model
    beam = Beam(self.s0)
    gonio = Goniometer(self.axis)
    scan = Scan(self.image_range, self.dphi)
    crystal = crystal_model((1, 0, 0), (0, 1, 0), (0, 0, 1), 1)
    crystal.set_A_at_scan_points([self.ub] * (scan.get_num_images() + 1))
    predict = ScanVaryingReflectionPredictor(crystal, beam, gonio, scan, self.dmin)
    predict.prepare(frame, 1)
    eps = 1e-7
    assert(all(abs(a1-a2) < eps for a1, a2 in zip(A1, predict.get_A1())))
    assert(all(abs(a1-a2) < eps for a1, a2 in zip(A2, predict.get_A2())))
    result = [predict.predict(h) for h in hkl]
    return result

  def predict_cc(self, hkl, frame, A1, A2):
    from dials.algorithms.spot_prediction import ScanVaryingRayPredictor
    from math import pi
    dphi = (self.dphi[0] * pi / 180.0, self.dphi[1] * pi / 180.0)
    predict = ScanVaryingRayPredictor(self.s0, self.axis, dphi, self.dmin)
    result = [predict(h, A1, A2, frame, 1) for h in hkl]
    return result

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
