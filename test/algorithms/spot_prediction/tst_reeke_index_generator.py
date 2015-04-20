from __future__ import division


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

  def run(self):

    num = 20

    # Check indices are the same for each frame. Exceptionally there may be
    # small differences that we can accept (e.g. 32 bit Linux C++ generate_cpp
    # has been seen to generate 3 fewer indices for frame=0).
    for frame in range(num):

      py_index = self.generate_python(frame)
      cp_index = self.generate_cpp(frame)

      py_set = set(py_index)
      cp_set = set(cp_index)

      # Sanity check that we're not generating zero length lists, for example,
      # that would pass other tests.
      assert len(py_index) > 3800

      # Check indices are unique in list
      assert(len(py_set) == len(py_index))
      assert(len(cp_set) == len(cp_index))

      # Check there are 6 or fewer differences between the two lists. A
      # difference of 6 indices may be caused by just 2 additional Ewald sphere
      # intersections in one list compared to the other. Each of these may
      # be increased to 3 extra indices by the addition of a margin +/1 around
      # that index.
      diff = py_set ^ cp_set
      assert len(diff) <= 6

      print 'OK'

  def get_ub(self, frame):
    from scitbx import matrix
    import scitbx.math

    angle_beg = frame * 1
    angle_end = (frame+1) * 1

    r_osc_beg = matrix.sqr(
      scitbx.math.r3_rotation_axis_and_angle_as_matrix(
      axis = self.axis, angle = angle_beg, deg=True))

    r_osc_end = matrix.sqr(
      scitbx.math.r3_rotation_axis_and_angle_as_matrix(
      axis = self.axis, angle = angle_end, deg=True))

    ub_beg = r_osc_beg * self.ub
    ub_end = r_osc_end * self.ub
    return ub_beg, ub_end

  def generate_python(self, frame):
    from dials.algorithms.spot_prediction.reeke import reeke_model
    ub_beg, ub_end = self.get_ub(frame)
    r = reeke_model(ub_beg, ub_end, self.axis, self.s0, self.dmin, self.margin)
    hkl = r.generate_indices()
    hkl = [h for h in hkl if h != (0, 0, 0)]
    return sorted(hkl)

  def generate_cpp(self, frame):
    from dials.algorithms.spot_prediction import ReekeIndexGenerator
    from cctbx.sgtbx import space_group_info
    space_group_type = space_group_info("P 1").group().type()
    ub_beg, ub_end = self.get_ub(frame)
    r = ReekeIndexGenerator(ub_beg, ub_end, space_group_type, self.axis, self.s0, self.dmin, self.margin)
    hkl = r.to_array()
    return sorted(list(hkl))

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
