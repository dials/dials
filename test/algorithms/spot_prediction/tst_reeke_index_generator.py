

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

    # Check indices are the same for each frame
    for frame in range(num):

      py_index = self.generate_python(frame)
      cp_index = self.generate_cpp(frame)

      # Check indices are unique in list
      assert(len(set(cp_index)) == len(cp_index))

      ## Check indices are in order
      #for idx1, idx2 in zip(cp_index[0:-1], cp_index[1:]):
        #print idx1, idx2
        #assert(idx2[0] >= idx1[0])
        #if idx2[0] == idx1[0]:
          #assert(idx2[1] >= idx1[1])
        #if idx2[1] == idx2[0]:
          #assert(idx2[2] >= idx1[2])

      # Check all the indices are the same
      assert(len(py_index) == len(cp_index))
      for idx1, idx2 in zip(py_index, cp_index):
        assert(idx1 == idx2)

      print 'OK'

  def get_ub(self, frame):
    import scitbx
    from scitbx import matrix

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
    from dials.algorithms.refinement.prediction.reeke import reeke_model
    ub_beg, ub_end = self.get_ub(frame)
    r = reeke_model(ub_beg, ub_end, self.axis, self.s0, self.dmin, self.margin)
    hkl = r.generate_indices()
    hkl = [h for h in hkl if h != (0, 0, 0)]
    return sorted(hkl)

  def generate_cpp(self, frame):
    from dials.algorithms.spot_prediction import ReekeIndexGenerator
    ub_beg, ub_end = self.get_ub(frame)
    r = ReekeIndexGenerator(ub_beg, ub_end, self.axis, self.s0, self.dmin, self.margin)
    hkl = r.to_array()
    return sorted(list(hkl))

if __name__ == '__main__':
  test = Test()
  test.run()
